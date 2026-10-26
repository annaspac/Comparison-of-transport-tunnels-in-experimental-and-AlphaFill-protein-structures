#!/usr/bin/env python3
"""
Workflow for the revision with the RECALCULATED cytochrome P450 tunnels.

Step 1  collect : reads the new MOLE outputs (<structure>/json/data.json, same parsing as the
                  original scripts) and writes the channel table + UniProt->PDB mapping in the
                  format of the original scripts.
Step 2  compare : runs tunnel_pipeline.py (SIFTS residue mapping, best matches in all directions,
                  1:1 control, statistical tests, figures) for the new CYP data together with the
                  ChannelsDB-annotated ("reviewed") dataset.

Figures: Fig. 4 and Fig. 5 have the median as a dashed line; Fig. 5 contains both datasets.
Tests  : asymmetry_tests.csv (group-level Wilcoxon + Cliff's delta + bootstrap CI, 1:1 control,
         pathway-level Mann-Whitney U and Kolmogorov-Smirnov).

Usage (from the repository root):

  python3 scripts/03_cyp_workflow.py --pdb-out MOLE_OUT/pdb --af-out MOLE_OUT/alphafill

Original usage notes:

  # tunnels calculated with mole_from_heme.py (uses its manifest)
  python3 cyp_workflow.py --manifest cytochrome_p450/mole_heme/manifest.csv

  # tunnels calculated with your own scripts: give the output folders
  python3 cyp_workflow.py --pdb-out cytochrome_p450/out_pdb_new --af-out cytochrome_p450/out_alphafill_new

Optional: --map (default cytochrome_p450/cytochrome_p450_pdb_map_new.csv),
          --annotated-csv / --annotated-map (reviewed dataset, defaults as in the original layout),
          --out (default results_revision_new).
Requires pandas, numpy, scipy, matplotlib, requests and internet for the first run (SIFTS, UniProt).
"""

import argparse
import csv
import json
import re
from pathlib import Path

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
import tunnel_pipeline as tp

SECTIONS = ["Tunnels", "Paths", "Pores", "MergedPores"]


def extract_residue_flow(item):
    layers = item.get("Layers")
    if isinstance(layers, dict) and isinstance(layers.get("ResidueFlow"), list):
        return layers["ResidueFlow"]
    rf = item.get("ResidueFlow")
    return rf if isinstance(rf, list) else None


def channel_rows(data_json, structure_id, source, chain_rename=None):
    data = json.loads(Path(data_json).read_text())
    rows = []
    for section in SECTIONS:
        for item in (data.get("Channels", {}) or {}).get(section, []) or []:
            rf = extract_residue_flow(item)
            if not rf:
                continue
            if chain_rename:   # chain renamed in cleaned input (multi-character chain IDs) -> original
                rf = [re.sub(r"^(\S+\s+-?\d+\S*\s+)(\S+)", lambda m: m.group(1) + chain_rename.get(m.group(2), m.group(2)), t)
                      for t in rf]
            rows.append({
                "structure_id": structure_id, "structure_source": source,
                "channel_type": section, "channel_id": str(item.get("Id")),
                "cavity": item.get("Cavity"), "auto": item.get("Auto"),
                "residue_flow": repr(rf), "n_residues": len(rf),
            })
    return rows


def af_accession(name):
    m = re.search(r"([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})", name.upper())
    return m.group(1) if m else None


def collect(args, outdir):
    u2p = {}
    if args.map.endswith(".json"):          # {uniprot: [pdb, ...]}
        for k, v in json.loads(Path(args.map).read_text()).items():
            u2p[k.strip().upper()] = {p.strip().lower() for p in v}
    else:                                   # csv with columns uniprot, pdb
        for r in csv.DictReader(open(args.map)):
            u2p.setdefault(r["uniprot"].strip().upper(), set()).add(r["pdb"].strip().lower())
    rows, missing = [], []

    if args.manifest:
        base = Path(args.manifest).parent
        for m in csv.DictReader(open(args.manifest)):
            if m["status"] != "ok":
                continue
            js = Path(m["xml"]).parent / "out" / "json" / "data.json"
            if not js.exists():   # manifest moved: resolve relative to it
                js = base / ("pdb" if m["source"] == "PDB" else "alphafill") / m["id"] / "out" / "json" / "data.json"
            if not js.exists():
                missing.append(m["id"])
                continue
            if m["source"] == "PDB":
                ren = {m["chain_in_file"]: m["chain_original"]} if m["chain_in_file"] != m["chain_original"] else None
                rows += channel_rows(js, m["structure"].lower(), "pdb", ren)
            else:
                rows += channel_rows(js, m["uniprot"].upper(), "uniprot")
    else:
        for d in sorted(Path(args.pdb_out).iterdir()):
            js = d / "json" / "data.json"
            if js.exists():
                rows += channel_rows(js, d.name.lower(), "pdb")
            else:
                missing.append(d.name)
        for d in sorted(Path(args.af_out).iterdir()):
            js = d / "json" / "data.json"
            acc = af_accession(d.name)
            if js.exists() and acc:
                rows += channel_rows(js, acc, "uniprot")
            else:
                missing.append(d.name)

    pdb_with = {r["structure_id"] for r in rows if r["structure_source"] == "pdb"}
    af_with = {r["structure_id"] for r in rows if r["structure_source"] == "uniprot"}
    print(f"Collected {len(rows)} pathways: {len(pdb_with)} PDB structures, {len(af_with)} AlphaFill models; "
          f"{len(missing)} structures without MOLE output")
    outdir.mkdir(parents=True, exist_ok=True)
    ch_csv = outdir / "cyp_channels_new.csv"
    with open(ch_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["structure_id", "structure_source", "channel_type", "channel_id",
                                          "cavity", "auto", "residue_flow", "n_residues"])
        w.writeheader()
        w.writerows(rows)
    map_json = outdir / "cyp_uniprot_to_pdb_new.json"
    map_json.write_text(json.dumps({k: sorted(v) for k, v in u2p.items()}, indent=1))
    (outdir / "structures_without_mole_output.txt").write_text("\n".join(missing) + "\n")
    return ch_csv, map_json


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", help="manifest.csv written by mole_from_heme.py")
    ap.add_argument("--pdb-out", help="folder with MOLE outputs of PDB structures (<id>/json/data.json)")
    ap.add_argument("--af-out", help="folder with MOLE outputs of AlphaFill models")
    ap.add_argument("--map", default="data/cyp_uniprot_to_pdb.json", help="UniProt->PDB mapping (.json or .csv)")
    ap.add_argument("--annotated-csv", default="data/reviewed/reviewed_pdb_vs_uni_channels_uniprot_and_pdb.csv")
    ap.add_argument("--annotated-map", default="data/reviewed/reviewed_pdb_vs_uni_uniprot_to_pdb.json")
    ap.add_argument("--out", default="results_new")
    ap.add_argument("--reviewed-label", default="Reviewed",
                    help='name of the ChannelsDB dataset in figures and tables (default "Reviewed")')
    args = ap.parse_args()
    if not args.manifest and not (args.pdb_out and args.af_out):
        ap.error("give --manifest, or both --pdb-out and --af-out")

    out = Path(args.out)
    ch_csv, map_json = collect(args, out / "input_tables")

    tp.DATASETS = {   # order = order in figures and tables: Reviewed first, then CYP
        args.reviewed_label: (args.annotated_csv, args.annotated_map, None),
        "CYP": (str(ch_csv), str(map_json), "PF00067"),
    }
    tp.OUT = out
    tp.RUN_RAW = False
    tp.ORDER = [(ds, d) for ds in tp.DATASETS for d in ["PDB-PDB", "PDB-AlphaFill", "AlphaFill-PDB"]]
    tp.main()


if __name__ == "__main__":
    main()
