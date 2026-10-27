# Transport tunnels in experimental and AlphaFill protein structures

Code and processed data for the article

> Špačková A., Bazgier V., Berka K. *Comparison of transport tunnels in experimental and AlphaFill protein structures.* Frontiers in Chemistry (2026).

Transport tunnels (MOLE 2.5) in experimentally determined PDB structures are compared with tunnels in
ligand-completed AlphaFill models of the same protein. Tunnels are represented by their lining residues
(ResidueFlow), mapped to UniProt sequence positions with SIFTS, and compared with the Jaccard index,
overlap coefficient and Dice coefficient in three directions: PDB→PDB, PDB→AlphaFill and AlphaFill→PDB.

## Datasets

| Dataset | Source | UniProt groups | PDB structures | AlphaFill models |
|---|---|---|---|---|
| CYP | UniProtKB query `xref:pfam-PF00067 AND structure_3d:true AND reviewed:true` (21 Oct 2026) | 105 | 1,079 | 98 |
| Reviewed | literature-curated pathways from ChannelsDB 2.0 | 133 | 304 | 17 groups with AlphaFill pathways |

## Repository layout

```
scripts/
  00_prepare_reviewed_dataset.py   pathway table of the ChannelsDB reviewed dataset (from ChannelsDB archives)
  01_download_p450_structures.py   UniProt query, download of PDB structures and AlphaFill models
  02_run_mole/                     MOLE 2.5 template XML, batch runner (run_mole.py) and PBS example
  03_cyp_workflow.py               collects MOLE outputs and runs the whole comparison
  tunnel_pipeline.py               SIFTS mapping, best matches, statistics, figures (used by 03 and 04)
  04_make_figures.py               redraws all figures from existing results
data/
  QUERY.txt                        UniProt query and retrieval date
  cyp_uniprot_to_pdb.json          UniProt -> PDB mapping of the CYP dataset
  reviewed/                        ChannelsDB reviewed dataset (pathways + mapping)
results/
  summary_statistics.csv           n, zero-overlap %, medians and quartiles per dataset and direction
  asymmetry_tests.csv              Wilcoxon, Cliff's delta, bootstrap CI, 1:1 control, Mann-Whitney U, KS
  table_S1_identifiers.csv         Supplementary Table S1 (all UniProt / PDB / AlphaFill identifiers)
  best_matches_all_directions.csv.gz   every best-matching pathway pair with all metrics
  pair_control_1to1.csv            1:1 comparison of each PDB structure with the AlphaFill model
  mapping_qc.csv                   residue mapping per PDB chain (SIFTS / sequence offset / excluded)
  excluded_groups.csv              groups removed by the Pfam filter
  cyp_channels.csv.gz              all calculated CYP pathways (ResidueFlow residue lists)
  figures/                         Figures 3-5 and supplementary figures (PNG 600 dpi + PDF)
```

## Requirements

- Python ≥ 3.6; step 1 needs only the standard library, steps 3–4 need `pip install -r requirements.txt`
- MOLE 2.5 command-line version and `mono` (step 2)
- internet access for steps 1 and 3 (UniProt, RCSB, AlphaFill, PDBe SIFTS API); API responses are cached in `cache/`

On MetaCentrum: `module add mono` (MOLE) and `module add py-pandas py-scipy py-matplotlib py-requests` (steps 3–4).

## Reproducing the analysis

```bash
# 0. (optional) rebuild the reviewed dataset table from ChannelsDB archives (result already in data/reviewed/)
python3 scripts/00_prepare_reviewed_dataset.py --map uniprot_to_pdb_all.json \
    --pdb-base ChannelsDB/pdb --uniprot-base ChannelsDB/uniprot --out data/reviewed

# 1. dataset and structures (-> work/pdb_new, work/alphafill_new)
python3 scripts/01_download_p450_structures.py --outdir work

# 2. MOLE 2.5 from the heme cofactor (-> mole_out/<structure id>/json/data.json)
module add mono      # MetaCentrum; or submit scripts/02_run_mole/run_mole.pbs with qsub
python3 scripts/02_run_mole/run_mole.py --mole /path/to/mole2.exe --input work/pdb_new --output mole_out/pdb --workers 8
python3 scripts/02_run_mole/run_mole.py --mole /path/to/mole2.exe --input work/alphafill_new --output mole_out/alphafill --workers 8

# 3. comparison, statistics and figures (-> results_new/)
python3 scripts/03_cyp_workflow.py --pdb-out mole_out/pdb --af-out mole_out/alphafill

# 4. only redraw the figures from existing results
python3 scripts/04_make_figures.py results
```

Step 3 runs in about 10–20 minutes. The PDB list depends on the date of the UniProt query; `data/` contains
the exact mapping used in the article.

## Method summary

1. **Residue mapping.** Lining residues of PDB pathways are mapped to UniProt positions with SIFTS (PDBe API)
   and validated against the UniProt sequence (≥ 90 % residue-name agreement; otherwise a sequence-based offset
   search, else the chain is excluded). AlphaFill models use UniProt numbering. Residues of other proteins and
   non-protein residues are discarded, backbone- and side-chain-lining entries of a residue are merged, and
   pathways with fewer than three residues are excluded.
2. **Best match.** Each query pathway is compared with all pathways of the same type in the target set of the
   same UniProt group; candidates are ranked by Jaccard, overlap, Dice and number of shared residues.
   PDB→PDB excludes pathways of the same PDB entry.
3. **Statistics.** Group-level paired Wilcoxon test on per-UniProt-group medians with Cliff's delta and bootstrap
   95 % CI (2,000 resamples); 1:1 control comparing each PDB structure with the AlphaFill model individually;
   descriptive pathway-level Mann–Whitney U and Kolmogorov–Smirnov tests.

## Data

Experimental structures are available from the Protein Data Bank and AlphaFill models from https://alphafill.eu; the exact identifiers and the UniProt query with its retrieval date are provided in data/. Processed results (all calculated pathways, best-match tables and statistics) are included in results/. Supplementary figures are deposited in Zenodo: https://doi.org/10.5281/zenodo.23666013.

## License

Code: MIT (see `LICENSE`). Structures are subject to the terms of the PDB and AlphaFill.

## Contact

Karel Berka, Department of Physical Chemistry, Faculty of Science, Palacký University Olomouc – karel.berka@upol.cz
