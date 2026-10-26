#!/usr/bin/env python3
"""Redraw all figures at print size from existing results (no recalculation).
Usage (from the repository root): python3 scripts/04_make_figures.py [results_dir] [output_dir]
"""
import sys
from pathlib import Path
import pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parent))
import tunnel_pipeline as tp

R = Path(sys.argv[1] if len(sys.argv) > 1 else "results")
O = Path(sys.argv[2] if len(sys.argv) > 2 else R / "figures")
O.mkdir(parents=True, exist_ok=True)
bm = R / "best_matches_all_directions.csv"
b = pd.read_csv(bm if bm.exists() else bm.with_suffix(".csv.gz"))
p = pd.read_csv(R / "pair_control_1to1.csv")
order = [ds for ds in b.dataset.unique() if ds != "CYP"] + ["CYP"]   # Reviewed first, then CYP
tp.DATASETS = {ds: None for ds in order}
tp.ORDER = [(ds, d) for ds in tp.DATASETS for d in ["PDB-PDB", "PDB-AlphaFill", "AlphaFill-PDB"]]
tp.fig3(b, O / "Fig3_overlap_and_dice.png")
tp.fig4(b, O / "Fig4_PDB_PDB.png")
tp.fig5(b, O / "Fig5_AlphaFill_vs_PDB.png", only=["CYP"])
for ds in tp.DATASETS:
    if ds != "CYP":
        tp.fig5(b, O / f"FigS_{ds}_AlphaFill_vs_PDB.png", only=[ds])
for m in ["jaccard", "overlap"]:
    tp.fig_si_metric_box(b, m, O / f"FigS_{m}_boxplot.png")
    tp.fig4(b, O / f"FigS_{m}_PDB_PDB.png", metric=m)
    tp.fig5(b, O / f"FigS_{m}_AlphaFill_vs_PDB.png", metric=m)
tp.fig_pair_control(p, O / "FigS_1to1_control.png")
print("Figures (PNG 600 dpi + PDF) in", O)
