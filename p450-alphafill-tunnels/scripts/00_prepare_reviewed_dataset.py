#!/usr/bin/env python3
"""
Step 0: build the ChannelsDB reviewed dataset table (standard library only).

Reads the ChannelsDB 2.0 data archives and writes the pathway table used by 03_cyp_workflow.py:
  - PDB entries:    <pdb-base>/*/<pdb id>/data.zip      -> authors.json   (literature-curated pathways)
  - UniProt models: <uniprot-base>/*/<UniProt>/data.zip -> alphafill.json (pathways in AlphaFill models)
and a UniProt -> PDB mapping (JSON {uniprot: [pdb, ...]}).

Outputs (in --out):
  reviewed_pdb_vs_uni_channels_uniprot_and_pdb.csv
  reviewed_pdb_vs_uni_uniprot_to_pdb.json

Usage:
  python3 scripts/00_prepare_reviewed_dataset.py --map uniprot_to_pdb_all.json \
      --pdb-base ChannelsDB/pdb --uniprot-base ChannelsDB/uniprot --out data/reviewed
The processed result used in the article is already provided in data/reviewed/.
"""

import argparse
import csv
import json
import zipfile
from pathlib import Path

SECTIONS = ["Tunnels", "Paths", "Pores", "MergedPores"]
FIELDS = ["structure_id", "structure_source", "pdb_id", "uniprot_id", "mapped_pdb_ids", "channel_type",
          "channel_id", "channel_name", "cavity", "auto", "residue_flow", "n_residues", "residue_flow_str",
          "has_profile"]


def load_json_from_zip(zip_path, name):
    with zipfile.ZipFile(zip_path) as z:
        if name not in z.namelist():
            return None
        with z.open(name) as f:
            return json.load(f)


def residue_flow(item):
    layers = item.get("Layers")
    if isinstance(layers, dict) and isinstance(layers.get("ResidueFlow"), list):
        return layers["ResidueFlow"]
    rf = item.get("ResidueFlow")
    return rf if isinstance(rf, list) else None


def channel_rows(sid, data, source, mapped=None):
    rows = []
    channels = data.get("Channels", {}) if isinstance(data, dict) else {}
    for section in SECTIONS:
        for item in channels.get(section, []) or []:
            rf = residue_flow(item)
            name = (item.get("mapped_name") or item.get("Name") or item.get("name") or item.get("Label")
                    or item.get("label") or item.get("Type") or item.get("type"))
            rows.append({
                "structure_id": sid, "structure_source": source,
                "pdb_id": sid if source == "pdb" else "", "uniprot_id": sid if source == "uniprot" else "",
                "mapped_pdb_ids": ";".join(mapped) if mapped else "",
                "channel_type": section, "channel_id": "" if item.get("Id") is None else str(item.get("Id")),
                "channel_name": name, "cavity": item.get("Cavity"), "auto": item.get("Auto"),
                "residue_flow": repr(rf) if rf is not None else "", "n_residues": len(rf) if rf else "",
                "residue_flow_str": "; ".join(rf) if rf else "", "has_profile": isinstance(item.get("Profile"), list),
            })
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--map", required=True, help="JSON {UniProt: [PDB IDs]} of the reviewed dataset")
    ap.add_argument("--pdb-base", required=True, help="ChannelsDB PDB folder (*/<pdb>/data.zip)")
    ap.add_argument("--uniprot-base", required=True, help="ChannelsDB UniProt folder (*/<UniProt>/data.zip)")
    ap.add_argument("--out", default="data/reviewed")
    a = ap.parse_args()

    u2p = {str(k).upper(): sorted({str(x).lower() for x in v}) for k, v in json.load(open(a.map)).items()}
    pdb_zip = {p.parent.name.lower(): p for p in Path(a.pdb_base).glob("*/*/data.zip")}
    uni_zip = {p.parent.name.upper(): p for p in Path(a.uniprot_base).glob("*/*/data.zip")}

    rows, missing_pdb, missing_af = [], [], []
    for pdb in sorted({p for v in u2p.values() for p in v}):
        data = load_json_from_zip(pdb_zip[pdb], "authors.json") if pdb in pdb_zip else None
        if data is None:
            missing_pdb.append(pdb)
            continue
        rows += channel_rows(pdb, data, "pdb", [pdb])
    for uni, pdbs in sorted(u2p.items()):
        data = load_json_from_zip(uni_zip[uni], "alphafill.json") if uni in uni_zip else None
        if data is None:
            missing_af.append(uni)
            continue
        rows += channel_rows(uni, data, "uniprot", pdbs)

    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    with open(out / "reviewed_pdb_vs_uni_channels_uniprot_and_pdb.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS)
        w.writeheader()
        w.writerows(rows)
    (out / "reviewed_pdb_vs_uni_uniprot_to_pdb.json").write_text(json.dumps(u2p, indent=2))
    print(f"UniProt groups: {len(u2p)}; pathways: {len(rows)}; "
          f"PDB without authors.json: {len(missing_pdb)}; UniProt without alphafill.json: {len(missing_af)}")


if __name__ == "__main__":
    main()
