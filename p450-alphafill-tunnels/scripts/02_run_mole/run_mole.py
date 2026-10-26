#!/usr/bin/env python3
"""
Step 2: run MOLE 2.5 on all structures in a folder (standard library only).

For every structure (*.pdb, *.cif) in --input, a copy of the template XML is written to
<output>/<structure id>/run.xml with <Input> and <WorkingDirectory> replaced, and
    mono <mole2.exe> run.xml
is executed. Results end up in <output>/<structure id>/json/data.json, the layout read by
03_cyp_workflow.py. Structures that already have data.json are skipped, so the script can be
re-run after an interruption.

Usage:
  python3 scripts/02_run_mole/run_mole.py --mole /path/to/mole2.exe \
      --input work/pdb_new --output mole_out/pdb --workers 8
  python3 scripts/02_run_mole/run_mole.py --mole /path/to/mole2.exe \
      --input work/alphafill_new --output mole_out/alphafill --workers 8

  --dry-run   only write the XML files, do not run MOLE

On MetaCentrum, load mono first (also inside PBS jobs):  module add mono
"""

import argparse
import subprocess
import sys
import xml.etree.ElementTree as ET
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent


def write_xml(template, path, structure, workdir):
    tree = ET.parse(template)
    root = tree.getroot()
    inp, wd = root.find("./Input"), root.find("./WorkingDirectory")
    if inp is None or wd is None:
        sys.exit("Template XML must contain <Input> and <WorkingDirectory>.")
    inp.text = str(structure)
    wd.text = str(workdir) + "/"
    tree.write(path, encoding="utf-8", xml_declaration=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mole", help="path to mole2.exe")
    ap.add_argument("--input", required=True, help="folder with *.pdb / *.cif structures")
    ap.add_argument("--output", required=True, help="output folder (one subfolder per structure)")
    ap.add_argument("--template", default=str(HERE / "mole_template.xml"))
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    if not a.dry_run and not a.mole:
        ap.error("--mole is required unless --dry-run is used")

    structures = sorted(p for p in Path(a.input).iterdir() if p.suffix.lower() in (".pdb", ".cif"))
    out = Path(a.output)
    out.mkdir(parents=True, exist_ok=True)
    todo = []
    for s in structures:
        wd = (out / s.stem).resolve()
        if (wd / "json" / "data.json").exists():
            continue
        wd.mkdir(parents=True, exist_ok=True)
        write_xml(a.template, wd / "run.xml", s.resolve(), wd)
        todo.append(wd)
    print(f"{len(structures)} structures, {len(structures) - len(todo)} already done, {len(todo)} to run", flush=True)
    if a.dry_run:
        return

    def run(wd):
        with open(wd / "mole.log", "w") as log:
            rc = subprocess.run(["mono", a.mole, str(wd / "run.xml")], stdout=log, stderr=subprocess.STDOUT).returncode
        return wd.name, rc, (wd / "json" / "data.json").exists()

    failed = []
    with ThreadPoolExecutor(a.workers) as ex:
        for i, (name, rc, ok) in enumerate(ex.map(run, todo), 1):
            if not ok:
                failed.append(name)
                print(f"  FAILED {name} (return code {rc}, see {name}/mole.log)", flush=True)
            if i % 50 == 0 or i == len(todo):
                print(f"  {i}/{len(todo)}", flush=True)
    (out / "failed.txt").write_text("\n".join(failed) + "\n")
    print(f"Done. Failed: {len(failed)} (list in {out / 'failed.txt'})")


if __name__ == "__main__":
    main()
