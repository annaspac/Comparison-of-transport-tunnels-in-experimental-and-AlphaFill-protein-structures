# Step 2: tunnel calculation with MOLE 2.5

MOLE 2.5 (command-line version, `mole2.exe`) is run on every downloaded PDB structure
(`work/pdb_new/`) and AlphaFill model (`work/alphafill_new/`) with the settings in
`mole_template.xml`:

- origins: heme residues (`Residues("HEM")`), automatic origin detection disabled (`Auto="0"`);
- probe radius 5 Å, interior threshold 1.4 Å, origin radius 5 Å, surface cover radius 10 Å,
  bottleneck radius 1.25 Å, maximum tunnel similarity 0.7;
- heteroatoms ignored during cavity detection (`IgnoreHETAtoms="1"`).

`run_mole.py` writes a copy of the template for each structure with `<Input>` and `<WorkingDirectory>`
replaced and calls `mono mole2.exe run.xml`. The output of each structure ends up in
`<output>/<structure id>/json/data.json` (PDB: folder named by PDB ID; AlphaFill: e.g. `AF-P08684-F1`),
which is what `03_cyp_workflow.py` reads. Finished structures are skipped on re-run; failures are listed
in `<output>/failed.txt` with the MOLE log in `<structure id>/mole.log`.

```bash
python3 run_mole.py --mole /path/to/mole2.exe --input ../../work/pdb_new --output ../../mole_out/pdb --workers 8
python3 run_mole.py --input ../../work/pdb_new --output ../../mole_out/pdb --dry-run   # only write XML files
```

On MetaCentrum, load mono inside the job (`module add mono`); `run_mole.pbs` is an example job script.
