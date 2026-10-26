#!/usr/bin/env python3
"""
Build the revised cytochrome P450 dataset and download structures.
Standard library only (no pandas / requests needed). Python >= 3.6.

1. UniProtKB query: Swiss-Prot entries with the P450 domain (Pfam PF00067)
   and at least one experimental structure. All result pages are retrieved.
2. PDB structures of each UniProt entry  -> pdb_new/
   (legacy .pdb; entries without a .pdb file are downloaded as .cif)
3. AlphaFill models of each UniProt entry -> alphafill_new/
   (model .cif + AlphaFill metadata .json with the transplanted ligands)

Outputs (in the working directory):
  cytochrome_p450_pdb_map_new.csv      uniprot, pdb, protein, organism
  download_manifest_pdb.csv            status of every PDB download
  download_manifest_alphafill.csv      status of every AlphaFill download (+ number of HEM transplants)
  QUERY.txt                            query text and retrieval date

Usage:
  python3 scripts/01_download_p450_structures.py --outdir work --limit 3    # quick test
  python3 scripts/01_download_p450_structures.py --outdir work              # everything

Already downloaded files are skipped, so the script can be re-run after an interruption.
"""

import argparse
import csv
import io
import json
import re
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from datetime import date
from pathlib import Path

UNIPROT_QUERY = "xref:pfam-PF00067 AND structure_3d:true AND reviewed:true"
UNIPROT_URL = "https://rest.uniprot.org/uniprotkb/search"
PDB_URL = "https://files.rcsb.org/download/{id}.{ext}"
AF_MODEL_URL = "https://alphafill.eu/v1/aff/{acc}"
AF_JSON_URL = "https://alphafill.eu/v1/aff/{acc}/json"

N_WORKERS = 6
TIMEOUT = 120


def http_get(url, retries=4):
    """Returns (status, body_bytes, headers). status None on repeated network failure."""
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "p450-tunnels/1.0"})
            with urllib.request.urlopen(req, timeout=TIMEOUT) as r:
                return r.status, r.read(), r.headers
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return 404, b"", e.headers
        except (urllib.error.URLError, TimeoutError, ConnectionError, OSError):
            pass
        time.sleep(3 * (attempt + 1))
    return None, b"", {}


def next_link(headers):
    link = headers.get("Link") or headers.get("link") or ""
    m = re.search(r'<([^>]+)>;\s*rel="next"', link)
    return m.group(1) if m else None


def fetch_uniprot():
    params = urllib.parse.urlencode({
        "query": UNIPROT_QUERY,
        "format": "tsv",
        "fields": "accession,protein_name,organism_name,xref_pdb",
        "size": 500,
    })
    url = f"{UNIPROT_URL}?{params}"
    rows, total = [], None
    while url:
        status, body, headers = http_get(url)
        if status != 200:
            raise SystemExit(f"UniProt query failed (status {status})")
        if total is None:
            total = headers.get("x-total-results")
        reader = csv.DictReader(io.StringIO(body.decode("utf-8")), delimiter="\t")
        rows.extend(reader)
        url = next_link(headers)
    print(f"UniProt query: {UNIPROT_QUERY}")
    print(f"  hits reported: {total}, retrieved: {len(rows)}")
    if total is not None and int(total) != len(rows):
        raise SystemExit("ERROR: pagination incomplete!")
    return rows


def build_map(rows):
    out = []
    for row in rows:
        pdbs = (row.get("PDB") or "").strip()
        for pdb in pdbs.split(";"):
            pdb = pdb.strip().upper()
            if pdb:
                out.append({"uniprot": row["Entry"], "pdb": pdb,
                            "protein": row.get("Protein names", ""), "organism": row.get("Organism", "")})
    return out


def write_csv(path, rows, fields):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)


def nonempty(p):
    return p.exists() and p.stat().st_size > 0


def download_pdb(pdb, outdir):
    for ext in ("pdb", "cif"):
        out = outdir / f"{pdb.lower()}.{ext}"
        if nonempty(out):
            return {"pdb": pdb, "file": out.name, "status": "exists"}
    status, body, _ = http_get(PDB_URL.format(id=pdb, ext="pdb"))
    if status == 200:
        (outdir / f"{pdb.lower()}.pdb").write_bytes(body)
        return {"pdb": pdb, "file": f"{pdb.lower()}.pdb", "status": "downloaded"}
    status, body, _ = http_get(PDB_URL.format(id=pdb, ext="cif"))   # large entries: no legacy .pdb
    if status == 200:
        (outdir / f"{pdb.lower()}.cif").write_bytes(body)
        return {"pdb": pdb, "file": f"{pdb.lower()}.cif", "status": "downloaded_cif"}
    return {"pdb": pdb, "file": "", "status": "failed" if status is None else f"http_{status}"}


def download_alphafill(acc, outdir):
    model = outdir / f"AF-{acc}-F1.cif"
    meta = outdir / "json" / f"AF-{acc}-F1.json"
    status = "exists"
    if not nonempty(model):
        st, body, _ = http_get(AF_MODEL_URL.format(acc=acc))
        if st != 200 or not body.startswith(b"data_"):
            return {"uniprot": acc, "file": "", "status": "not_in_alphafill" if st is not None else "failed",
                    "n_hem_transplants": 0}
        model.write_bytes(body)
        status = "downloaded"
    if not nonempty(meta):
        st, body, _ = http_get(AF_JSON_URL.format(acc=acc))
        if st == 200:
            meta.write_bytes(body)
    n_hem = 0
    if nonempty(meta):
        try:
            d = json.loads(meta.read_text())
            n_hem = sum(t.get("analogue_id") == "HEM" or t.get("compound_id") == "HEM"
                        for h in (d.get("hits") or []) for t in (h.get("transplants") or []))
        except json.JSONDecodeError:
            pass
    return {"uniprot": acc, "file": model.name, "status": status, "n_hem_transplants": n_hem}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=".", help="base folder (default: current directory)")
    ap.add_argument("--limit", type=int, default=None, help="only first N UniProt entries (testing)")
    args = ap.parse_args()

    base = Path(args.outdir)
    pdb_dir, af_dir = base / "pdb_new", base / "alphafill_new"
    pdb_dir.mkdir(parents=True, exist_ok=True)
    (af_dir / "json").mkdir(parents=True, exist_ok=True)

    mapping = build_map(fetch_uniprot())
    if args.limit:
        keep = []
        for r in mapping:
            if r["uniprot"] not in keep:
                keep.append(r["uniprot"])
        keep = set(keep[:args.limit])
        mapping = [r for r in mapping if r["uniprot"] in keep]
    write_csv(base / "cytochrome_p450_pdb_map_new.csv", mapping, ["uniprot", "pdb", "protein", "organism"])
    (base / "QUERY.txt").write_text(f"UniProtKB query: {UNIPROT_QUERY}\nretrieved: {date.today()}\n")

    accs = sorted({r["uniprot"] for r in mapping})
    pdbs = sorted({r["pdb"] for r in mapping})
    print(f"  UniProt entries with PDB: {len(accs)}, unique PDB entries: {len(pdbs)}", flush=True)

    print(f"\nDownloading {len(pdbs)} PDB entries -> {pdb_dir}", flush=True)
    res_pdb = []
    with ThreadPoolExecutor(N_WORKERS) as ex:
        for i, res in enumerate(ex.map(lambda p: download_pdb(p, pdb_dir), pdbs), 1):
            res_pdb.append(res)
            if i % 100 == 0 or i == len(pdbs):
                print(f"  {i}/{len(pdbs)}", flush=True)
    pdb2acc = {}
    for r in mapping:
        pdb2acc.setdefault(r["pdb"], []).append(r["uniprot"])
    for r in res_pdb:
        r["uniprot"] = ";".join(pdb2acc.get(r["pdb"], []))
    write_csv(base / "download_manifest_pdb.csv", res_pdb, ["pdb", "file", "status", "uniprot"])

    print(f"\nDownloading {len(accs)} AlphaFill models -> {af_dir}", flush=True)
    with ThreadPoolExecutor(N_WORKERS) as ex:
        res_af = list(ex.map(lambda a: download_alphafill(a, af_dir), accs))
    write_csv(base / "download_manifest_alphafill.csv", res_af, ["uniprot", "file", "status", "n_hem_transplants"])

    print("\n=== SUMMARY ===")
    print("PDB:", dict(Counter(r["status"] for r in res_pdb)))
    print("AlphaFill:", dict(Counter(r["status"] for r in res_af)))
    print("AlphaFill models without transplanted HEM:",
          sum(1 for r in res_af if r["file"] and r["n_hem_transplants"] == 0))


if __name__ == "__main__":
    main()
