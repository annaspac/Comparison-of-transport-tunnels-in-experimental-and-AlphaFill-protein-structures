#!/usr/bin/env python3
"""
Tunnel comparison pipeline: experimental (PDB) vs AlphaFill transport pathways.

Revised version addressing reviewer comments. Main changes vs. the original scripts:
  * residues are mapped to UniProt positions via SIFTS (PDBe API), so that
    chain IDs (A/B/C..., symmetry copies "C-2") and author numbering offsets
    no longer prevent matching of identical residues;
  * every mapping is validated against the UniProt sequence; chains with poor
    agreement fall back to a sequence-based offset search or are dropped;
  * backbone- and side-chain-lining contributions of a residue are merged
    (one residue = one element), configurable;
  * the CYP dataset is restricted to UniProt entries with the P450 Pfam domain
    (PF00067); this removes NADPH-cytochrome P450 reductases and the
    accidental "UNIPROT" header group;
  * a symmetric 1:1 structure-pair control and group-level statistics
    (paired Wilcoxon, Cliff's delta, bootstrap CI) test the PDB/AlphaFill asymmetry;
  * figures include median lines, the annotated dataset in Fig. 5,
    and SI figures for Jaccard and overlap coefficients.

Inputs are the channel tables already produced by the original scripts
(they contain the MOLE ResidueFlow of every pathway), so MOLE does not need to be re-run.

Run:  python tunnel_pipeline.py            (SIFTS mode + raw mode for comparison)
"""

import ast
import json
import re
import sys
import time
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
import requests
from scipy.stats import wilcoxon, mannwhitneyu, ks_2samp

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ============================================================
# CONFIG
# ============================================================
DATASETS = {
    # name: (channels csv with PDB + AlphaFill ResidueFlows, uniprot->pdb json, require Pfam)
    "CYP": (
        "cyt_uni_pdb/cytochrome_p450_channels_uniprot_and_pdb.csv",
        "cyt_uni_pdb/cytochrome_p450_uniprot_to_pdb.json",
        "PF00067",
    ),
    "Annotated": (
        "reviewed_pdb_vs_uni/reviewed_pdb_vs_uni_channels_uniprot_and_pdb.csv",
        "reviewed_pdb_vs_uni/reviewed_pdb_vs_uni_uniprot_to_pdb.json",
        None,
    ),
}

OUT = Path("results_revision")
CACHE = Path("cache")
RUN_RAW = True               # also run the original residue handling (raw strings) for comparison
MERGE_BACKBONE = True          # one residue = one element regardless of backbone/side-chain lining
MIN_RESIDUES_IN_GROUP = 3      # pathways with fewer residues mapped to the group protein are dropped
NAME_AGREEMENT_OK = 0.90       # minimum residue-name agreement with UniProt sequence
OFFSET_SEARCH = range(-150, 151)
N_WORKERS = 8
N_BOOT = 2000
RNG = np.random.RandomState(42)

COLORS = {"PDB-PDB": "#5E8B7E", "PDB-AlphaFill": "#B85C5C", "AlphaFill-PDB": "#7DA0B6"}

AA3 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
    "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
    "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "MSE": "M", "SEC": "U", "PYL": "O", "HYP": "P", "CSO": "C", "CME": "C", "OCS": "C",
    "SEP": "S", "TPO": "T", "PTR": "Y", "KCX": "K", "MLY": "K", "CSD": "C",
}


def log(*a):
    print(*a, flush=True)


# ============================================================
# EXTERNAL DATA (cached)
# ============================================================
def cached_json(url, path):
    path = Path(path)
    if path.exists():
        return json.loads(path.read_text())
    path.parent.mkdir(parents=True, exist_ok=True)
    for attempt in range(4):
        try:
            r = requests.get(url, timeout=60)
            if r.status_code == 404:
                data = {}
                break
            r.raise_for_status()
            data = r.json()
            break
        except Exception as e:  # noqa
            if attempt == 3:
                log(f"  ! failed {url}: {e}")
                return None
            time.sleep(2 * (attempt + 1))
    path.write_text(json.dumps(data))
    return data


def get_sifts(pdb):
    d = cached_json(f"https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{pdb.lower()}",
                    CACHE / "sifts" / f"{pdb.lower()}.json")
    return (d or {}).get(pdb.lower(), {}).get("UniProt", {})


def get_uniprot(acc):
    d = cached_json(f"https://rest.uniprot.org/uniprotkb/{acc}.json",
                    CACHE / "uniprot" / f"{acc}.json")
    if not d or "sequence" not in d:
        return None
    desc = d.get("proteinDescription", {})
    name = (desc.get("recommendedName") or (desc.get("submissionNames") or [{}])[0]).get("fullName", {}).get("value", "")
    genes = d.get("genes") or [{}]
    return {
        "accession": acc,
        "name": name,
        "gene": genes[0].get("geneName", {}).get("value", ""),
        "organism": d.get("organism", {}).get("scientificName", ""),
        "sequence": d["sequence"]["value"],
        "pfam": sorted(x["id"] for x in d.get("uniProtKBCrossReferences", []) if x.get("database") == "Pfam"),
    }


# ============================================================
# RESIDUE PARSING AND MAPPING
# ============================================================
TOKEN_RE = re.compile(r"^(\S+)\s+(-?\d+)([A-Za-z]?)\s+(\S+)(\s+Backbone)?$")


def parse_token(tok):
    """'GLN 245 A Backbone' -> (resname, number, ins_code, base_chain, is_backbone)."""
    m = TOKEN_RE.match(tok.strip())
    if not m:
        return None
    resname, num, ins, chain, bb = m.groups()
    base_chain = chain.split("-")[0]          # MOLE symmetry copies: C-2 -> C
    return resname.upper(), int(num), ins, base_chain, bool(bb)


def name_agreement(pairs, seq, offset=0):
    """pairs: list of (resname, uniprot_pos). Fraction matching UniProt sequence."""
    n = ok = 0
    for rn, pos in pairs:
        aa = AA3.get(rn)
        p = pos + offset
        if aa is None or p < 1 or p > len(seq):
            continue
        n += 1
        ok += seq[p - 1] == aa
    return (ok / n if n else np.nan), n


def _nz(x):
    return 0.0 if (x is None or (isinstance(x, float) and np.isnan(x))) else x


def build_pdb_chain_maps(pdb, residues_by_chain, uniprot_info, group_accs=()):
    """
    Returns {chain: (accession, segments, method, agreement)} where segments is a list of
    (author_lo, author_hi, offset) and uniprot_pos = author_number + offset of the segment
    containing the residue. Mapping is validated against the UniProt sequence.
    residues_by_chain: {chain: set((resname, author_num))} observed in pathways.
    """
    sifts = get_sifts(pdb)
    segs = defaultdict(lambda: defaultdict(list))   # chain -> acc -> [(lo, hi, off)]
    for acc, v in sifts.items():
        for m in v.get("mappings", []):
            s, e = m.get("start", {}), m.get("end", {})
            sa, ea = s.get("author_residue_number"), e.get("author_residue_number")
            sr, er = s.get("residue_number"), e.get("residue_number")
            if sa is not None:
                off = m["unp_start"] - sa
            elif ea is not None:
                off = m["unp_end"] - ea
            elif sr is not None:
                off = m["unp_start"] - sr
            else:
                continue
            lo = sa if sa is not None else -10**9
            hi = ea if ea is not None else 10**9
            segs[m["chain_id"]][acc].append((lo, hi, off))
    out = {}
    for chain, res in residues_by_chain.items():
        if chain not in segs:
            # chain absent from SIFTS (e.g. very recent entries): sequence-based offset vs group proteins
            cands = []
            for acc in group_accs:
                if acc in uniprot_info and len(res) >= 5:
                    seq = uniprot_info[acc]["sequence"]
                    ag, o = max(((_nz(name_agreement(list(res), seq, o)[0]), o) for o in OFFSET_SEARCH),
                                key=lambda x: (x[0], -abs(x[1])))
                    cands.append((ag, acc, o))
            if cands and max(cands)[0] >= NAME_AGREEMENT_OK:
                ag, acc, o = max(cands)
                out[chain] = (acc, [(-10**9, 10**9, o)], "sequence_offset_no_sifts", ag)
            else:
                out[chain] = (None, [], "no_sifts_chain", np.nan)
            continue
        best = None
        for acc, seglist in segs[chain].items():
            info = uniprot_info.get(acc)
            if info is None:
                continue
            pairs = []
            for rn, num in res:
                for lo, hi, off in seglist:
                    if lo <= num <= hi:
                        pairs.append((rn, num + off))
                        break
            ag, n = name_agreement(pairs, info["sequence"])
            score = (n, _nz(ag))          # chain is assigned to the accession covering most residues
            if best is None or score > best[0]:
                best = (score, acc, seglist, ag)
        if best is None:
            out[chain] = (None, [], "uniprot_unavailable", np.nan)
            continue
        (n, _), acc, seglist, ag = best
        method = "sifts"
        if n == 0 or (n >= 5 and _nz(ag) < NAME_AGREEMENT_OK):
            seq = uniprot_info[acc]["sequence"]
            ref = seglist[0][2]
            fb = max(((_nz(name_agreement(list(res), seq, o)[0]), o) for o in OFFSET_SEARCH),
                     key=lambda x: (x[0], -abs(x[1] - ref)))
            if fb[0] >= NAME_AGREEMENT_OK:
                seglist, ag, method = [(-10**9, 10**9, fb[1])], fb[0], "sequence_offset"
            else:
                method = "dropped_low_agreement"
        out[chain] = (acc, seglist, method, ag)
    return out


def map_residue(seglist, num):
    for lo, hi, off in seglist:
        if lo <= num <= hi:
            return num + off
    return None


# ============================================================
# LOADING
# ============================================================
def load_dataset(name, channels_csv, map_json, require_pfam):
    log(f"\n=== Loading {name} ===")
    df = pd.read_csv(channels_csv, dtype={"structure_id": str, "channel_id": str})
    df["tokens"] = df["residue_flow"].apply(lambda x: ast.literal_eval(x) if isinstance(x, str) else [])
    df = df[df.tokens.str.len() > 0].copy()
    df["source"] = np.where(df.structure_source.str.lower() == "pdb", "PDB", "AlphaFill")
    df["sid"] = np.where(df.source == "PDB", df.structure_id.str.lower(), df.structure_id.str.upper())

    u2p_raw = json.loads(Path(map_json).read_text())
    u2p = {}
    for k, v in u2p_raw.items():
        k = str(k).strip().upper()
        if not re.fullmatch(r"[A-Z0-9]{6,10}", k) or k == "UNIPROT":
            log(f"  dropping invalid group key: {k}")
            continue
        u2p[k] = sorted({str(x).lower() for x in v})

    # UniProt info for groups + every accession appearing in SIFTS of used PDBs
    pdb_ids = sorted(df.loc[df.source == "PDB", "sid"].unique())
    log(f"  fetching SIFTS for {len(pdb_ids)} PDB entries ...")
    with ThreadPoolExecutor(N_WORKERS) as ex:
        list(ex.map(get_sifts, pdb_ids))
    accs = set(u2p)
    for p in pdb_ids:
        accs |= set(get_sifts(p).keys())
    log(f"  fetching UniProt for {len(accs)} accessions ...")
    with ThreadPoolExecutor(N_WORKERS) as ex:
        infos = dict(zip(sorted(accs), ex.map(get_uniprot, sorted(accs))))
    infos = {k: v for k, v in infos.items() if v}

    excluded_groups = []
    if require_pfam:
        for g in list(u2p):
            if g not in infos or require_pfam not in infos[g]["pfam"]:
                excluded_groups.append((g, infos.get(g, {}).get("name", "?"), f"no {require_pfam}"))
                del u2p[g]
    for g, nm, why in excluded_groups:
        log(f"  excluded group {g} ({nm}): {why}")
    return df, u2p, infos, excluded_groups


def map_dataset(df, u2p, infos, mode):
    """
    Returns pathway records: list of dicts with uid, group, source, sid, ctype, set, n.
    mode 'sifts' : residues -> UniProt positions of the group protein
    mode 'raw'   : original behaviour (raw ResidueFlow strings), for reproducing manuscript numbers
    """
    qc_rows = []
    p2groups = defaultdict(list)
    for g, pl in u2p.items():
        for p in pl:
            p2groups[p].append(g)

    chain_maps = {}
    if mode == "sifts":
        for sid, sub in df[df.source == "PDB"].groupby("sid"):
            rbc = defaultdict(set)
            for toks in sub.tokens:
                for t in toks:
                    pt = parse_token(t)
                    if pt and pt[0] in AA3:          # protein residues only (heme etc. are not mapped)
                        rbc[pt[3]].add((pt[0], pt[1]))
            chain_maps[sid] = build_pdb_chain_maps(sid, rbc, infos, p2groups.get(sid, []))
            for ch, (acc, segl, meth, ag) in chain_maps[sid].items():
                offs = ";".join(str(s[2]) for s in segl)
                qc_rows.append({"pdb": sid, "chain": ch, "uniprot": acc, "offset_author_to_uniprot": offs,
                                "method": meth, "name_agreement": ag, "groups": ";".join(p2groups.get(sid, []))})

    recs = []
    dropped = defaultdict(int)
    for r in df.itertuples():
        if r.source == "PDB":
            groups = p2groups.get(r.sid, [])
        else:
            groups = [r.sid] if r.sid in u2p else []
        for g in groups:
            elems = set()
            for t in r.tokens:
                pt = parse_token(t)
                if pt is None:
                    continue
                rn, num, ins, ch, bb = pt
                if mode == "sifts" and rn not in AA3:   # heme and other non-polymer residues
                    continue
                if mode == "raw":
                    elems.add(t)
                    continue
                if r.source == "PDB":
                    acc, segl, meth, _ = chain_maps[r.sid].get(ch, (None, [], "", 0))
                    if acc != g or meth.startswith("dropped") or ins:
                        continue
                    pos = map_residue(segl, num)
                    if pos is None:
                        continue
                else:
                    pos = num  # AlphaFill models use UniProt numbering
                elems.add((pos, "BB") if (bb and not MERGE_BACKBONE) else pos)
            if mode == "sifts" and len(elems) < MIN_RESIDUES_IN_GROUP:
                dropped[(g, r.source)] += 1
                continue
            recs.append({
                "uid": f"{r.source}:{r.sid}__{r.channel_type}__{r.channel_id}",
                "group": g, "source": r.source, "sid": r.sid, "ctype": r.channel_type,
                "set": frozenset(elems), "n": len(elems),
            })
    # AlphaFill sequence check (QC only)
    if mode == "sifts":
        for sid, sub in df[df.source == "AlphaFill"].groupby("sid"):
            if sid not in infos:
                continue
            pairs = set()
            for toks in sub.tokens:
                for t in toks:
                    pt = parse_token(t)
                    if pt:
                        pairs.add((pt[0], pt[1]))
            ag, n = name_agreement(list(pairs), infos[sid]["sequence"])
            qc_rows.append({"pdb": f"AlphaFill:{sid}", "chain": "A", "uniprot": sid, "offset_author_to_uniprot": 0,
                            "method": "alphafill_uniprot_numbering", "name_agreement": ag, "groups": sid})
    return recs, pd.DataFrame(qc_rows), dropped


# ============================================================
# COMPARISONS
# ============================================================
def metrics(a, b):
    i = len(a & b)
    u = len(a) + len(b) - i
    return (i / u if u else 0.0,
            i / min(len(a), len(b)) if min(len(a), len(b)) else 0.0,
            2 * i / (len(a) + len(b)) if (len(a) + len(b)) else 0.0,
            i)


def best_match(query, targets):
    """Hierarchical ranking: Jaccard, overlap, Dice, shared residues."""
    qs, qn = query["set"], query["n"]
    best, best_t = (-1, -1, -1, -1), None
    for t in targets:
        i = len(qs & t["set"])
        tn = t["n"]
        m = (i / (qn + tn - i), i / min(qn, tn), 2 * i / (qn + tn), i)
        if m > best:
            best, best_t = m, t
    return best, best_t


def run_comparisons(recs, dataset):
    by_group = defaultdict(lambda: {"PDB": defaultdict(list), "AlphaFill": defaultdict(list)})
    for r in recs:
        by_group[r["group"]][r["source"]][r["sid"]].append(r)

    rows = []
    pair_rows = []  # 1:1 structure-pair control
    t0 = time.time()
    for gi, (g, d) in enumerate(by_group.items(), 1):
        pdb, af = d["PDB"], d["AlphaFill"]
        log(f"    [{dataset}] group {gi}/{len(by_group)} {g}: {len(pdb)} PDB, "
            f"{sum(len(v) for v in pdb.values())} pathways ({time.time() - t0:.0f}s)")
        af_all = [t for v in af.values() for t in v]
        pdb_ids = sorted(pdb)

        def add(direction, q, best, t):
            j, o, dc, i = best
            rows.append({"dataset": dataset, "direction": direction, "group": g,
                         "uid1": q["uid"], "structure1": q["sid"], "ctype": q["ctype"], "n1": q["n"],
                         "uid2": t["uid"] if t else None, "structure2": t["sid"] if t else None,
                         "n2": t["n"] if t else None, "intersection": i,
                         "jaccard": j, "overlap": o, "dice": dc})

        # PDB -> PDB (other structures of the same group)
        if len(pdb_ids) >= 2:
            for s in pdb_ids:
                others = defaultdict(list)
                for s2 in pdb_ids:
                    if s2 != s:
                        for t in pdb[s2]:
                            others[t["ctype"]].append(t)
                for q in pdb[s]:
                    if others[q["ctype"]]:
                        add("PDB-PDB", q, *best_match(q, others[q["ctype"]]))
        if not af_all:
            continue
        af_by_type = defaultdict(list)
        for t in af_all:
            af_by_type[t["ctype"]].append(t)
        # PDB -> AlphaFill
        for s in pdb_ids:
            for q in pdb[s]:
                if af_by_type[q["ctype"]]:
                    add("PDB-AlphaFill", q, *best_match(q, af_by_type[q["ctype"]]))
        # AlphaFill -> PDB (pooled over all experimental structures, as in manuscript)
        pdb_by_type = defaultdict(list)
        for s in pdb_ids:
            for t in pdb[s]:
                pdb_by_type[t["ctype"]].append(t)
        for q in af_all:
            if pdb_by_type[q["ctype"]]:
                add("AlphaFill-PDB", q, *best_match(q, pdb_by_type[q["ctype"]]))
        # 1:1 control: each PDB structure vs the AlphaFill model, both directions
        for s in pdb_ids:
            P = [q for q in pdb[s] if af_by_type[q["ctype"]]]
            A = [q for q in af_all if any(t["ctype"] == q["ctype"] for t in pdb[s])]
            if not P or not A:
                continue
            pd_ = [best_match(q, [t for t in A if t["ctype"] == q["ctype"]])[0][2] for q in P]
            ad_ = [best_match(q, [t for t in P if t["ctype"] == q["ctype"]])[0][2] for q in A]
            pair_rows.append({"dataset": dataset, "group": g, "pdb": s, "n_pdb": len(P), "n_af": len(A),
                              "pdb_matched": np.mean(np.array(pd_) > 0), "af_matched": np.mean(np.array(ad_) > 0),
                              "pdb_median_dice": np.median(pd_), "af_median_dice": np.median(ad_)})
    return pd.DataFrame(rows), pd.DataFrame(pair_rows)


# ============================================================
# STATISTICS
# ============================================================
def cliffs_delta(x, y):
    x, y = np.asarray(x), np.asarray(y)
    gt = sum((xi > y).sum() for xi in x)
    lt = sum((xi < y).sum() for xi in x)
    return (gt - lt) / (len(x) * len(y))


def boot_ci_median_diff(diff):
    diff = np.asarray(diff)
    bs = [np.median(RNG.choice(diff, len(diff), replace=True)) for _ in range(N_BOOT)]
    return np.percentile(bs, [2.5, 97.5])


def asymmetry_tests(best, pairs):
    out = []
    for ds in best.dataset.unique():
        b = best[best.dataset == ds]
        # (1) pooled best matches, aggregated per group
        g1 = b[b.direction == "PDB-AlphaFill"].groupby("group").dice.median()
        g2 = b[b.direction == "AlphaFill-PDB"].groupby("group").dice.median()
        gg = pd.concat([g1.rename("pdb_af"), g2.rename("af_pdb")], axis=1).dropna()
        # (0) pathway-level distributions (descriptive; pathways within a group are not independent)
        d_pa = b[b.direction == "PDB-AlphaFill"].dice.values
        d_ap = b[b.direction == "AlphaFill-PDB"].dice.values
        if len(d_pa) and len(d_ap):
            out.append({"dataset": ds, "test": "pathway-level Dice distributions (Mann-Whitney U)",
                        "n_groups": b[b.direction != "PDB-PDB"].group.nunique(), "n_AlphaFill_to_PDB": len(d_ap), "n_PDB_to_AlphaFill": len(d_pa),
                        "median_AlphaFill_to_PDB": np.median(d_ap), "median_PDB_to_AlphaFill": np.median(d_pa),
                        "median_difference": np.median(d_ap) - np.median(d_pa),
                        "mannwhitney_p": mannwhitneyu(d_ap, d_pa, alternative="two-sided").pvalue,
                        "ks_statistic": ks_2samp(d_ap, d_pa).statistic, "ks_p": ks_2samp(d_ap, d_pa).pvalue,
                        "cliffs_delta": cliffs_delta(d_ap, d_pa) if len(d_ap) * len(d_pa) < 5e7 else np.nan,
                        "note": "pathways not independent; see group-level tests"})
        # (2) 1:1 control aggregated per group
        p = pairs[pairs.dataset == ds].groupby("group")[["pdb_matched", "af_matched",
                                                         "pdb_median_dice", "af_median_dice"]].median()
        tests = [
            ("pooled best-match Dice (group medians)", gg.af_pdb, gg.pdb_af),
            ("1:1 fraction of pathways with counterpart", p.af_matched, p.pdb_matched),
            ("1:1 best-match Dice (group medians)", p.af_median_dice, p.pdb_median_dice),
        ]
        for label, a, b_ in tests:
            n = len(a)
            if n < 5:
                out.append({"dataset": ds, "test": label, "n_groups": n, "note": "too few groups"})
                continue
            diff = np.asarray(a) - np.asarray(b_)
            try:
                pval = wilcoxon(a, b_).pvalue
            except ValueError:
                pval = np.nan
            lo, hi = boot_ci_median_diff(diff)
            out.append({"dataset": ds, "test": label, "n_groups": n,
                        "median_AlphaFill_to_PDB": np.median(a), "median_PDB_to_AlphaFill": np.median(b_),
                        "median_difference": np.median(diff), "ci95_low": lo, "ci95_high": hi,
                        "wilcoxon_p": pval, "cliffs_delta": cliffs_delta(a, b_),
                        "groups_AF_higher": int((diff > 0).sum())})
    return pd.DataFrame(out)


def summary_table(best):
    rows = []
    for (ds, d), sub in best.groupby(["dataset", "direction"]):
        rows.append({"dataset": ds, "direction": d, "n_pathways": len(sub), "n_groups": sub.group.nunique(),
                     "zero_overlap_pct": 100 * (sub.intersection == 0).mean(),
                     **{f"median_{m}": sub[m].median() for m in ["dice", "jaccard", "overlap", "intersection"]},
                     **{f"q1_{m}": sub[m].quantile(.25) for m in ["dice"]},
                     **{f"q3_{m}": sub[m].quantile(.75) for m in ["dice"]}})
    return pd.DataFrame(rows)


# ============================================================
# FIGURES  (drawn at print size: Frontiers double column = 180 mm = 7.09 in)
# ============================================================
ORDER = [(ds, d) for ds in DATASETS for d in ["PDB-PDB", "PDB-AlphaFill", "AlphaFill-PDB"]]
FIG_W = 7.09
FIG_DPI = 600
ARROW = {"PDB-PDB": "PDB→PDB", "PDB-AlphaFill": "PDB→AlphaFill", "AlphaFill-PDB": "AlphaFill→PDB"}
plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 10, "axes.titlesize": 10.5, "axes.labelsize": 10,
    "xtick.labelsize": 9, "ytick.labelsize": 9, "legend.fontsize": 9,
    "axes.linewidth": 0.8, "lines.linewidth": 1.2, "pdf.fonttype": 42, "ps.fonttype": 42,
})


TICK = {"PDB-PDB": "PDB→PDB", "PDB-AlphaFill": "PDB→\nAlphaFill", "AlphaFill-PDB": "AlphaFill→\nPDB"}


def label(ds, d):
    return TICK[d]


def group_labels(ax, present, y=-0.30):
    """Dataset name centred under its group of x positions (1-based positions)."""
    groups = {}
    for i, (ds, d) in enumerate(present, 1):
        groups.setdefault(ds, []).append(i)
    tr = ax.get_xaxis_transform()
    for ds, xs in groups.items():
        lo, hi = min(xs) - 0.4, max(xs) + 0.4
        ax.plot([lo, hi], [y + 0.07, y + 0.07], color="black", lw=0.8, transform=tr, clip_on=False)
        ax.text((lo + hi) / 2, y, ds, ha="center", va="top", fontsize=10, fontweight="bold", transform=tr)


def clean(ax):
    for k in ("top", "right"):
        ax.spines[k].set_visible(False)


def panel_letter(ax, letter, x=-0.12, y=1.06):
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=13, fontweight="bold", va="bottom", ha="left")


def save(fig, path):
    path = Path(path)
    fig.savefig(path, dpi=FIG_DPI, bbox_inches="tight")
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def present_pairs(best, order=None):
    return [(ds, d) for ds, d in (order or ORDER) if ((best.dataset == ds) & (best.direction == d)).any()]


def fig3(best, path):
    present = present_pairs(best)
    fig, ax = plt.subplots(2, 1, figsize=(FIG_W, 7.0), gridspec_kw={"hspace": 0.75})
    x = np.arange(1, len(present) + 1)
    zero = [100 * (best[(best.dataset == ds) & (best.direction == d)].intersection == 0).mean() for ds, d in present]
    ax[0].bar(x, zero, color="#B85C5C", width=0.7, label="Zero overlap")
    ax[0].bar(x, [100 - z for z in zero], bottom=zero, color="#5E8B7E", width=0.7, label="Non-zero overlap")
    for i, z in enumerate(zero, 1):
        ax[0].text(i, z + (100 - z) / 2, f"{100 - z:.1f}%", ha="center", va="center", color="white", fontsize=9)
        ax[0].text(i, z + 2, f"{z:.1f}%", ha="center", va="bottom", fontsize=9,
                   color="white" if z > 8 else "black")
    ax[0].set_xticks(x)
    ax[0].set_xticklabels([label(*p) for p in present])
    ax[0].set_ylabel("Best matches (%)")
    ax[0].set_ylim(0, 100)
    ax[0].legend(frameon=False, ncol=2, loc="lower right", bbox_to_anchor=(1.0, 1.0))
    vals = [best[(best.dataset == ds) & (best.direction == d)].dice.values for ds, d in present]
    bp = ax[1].boxplot(vals, patch_artist=True, widths=0.6,
                       flierprops={"marker": "o", "markersize": 2, "alpha": 0.4, "markeredgewidth": 0.5})
    for p_, (ds, d) in zip(bp["boxes"], present):
        p_.set_facecolor(COLORS[d]); p_.set_linewidth(0.8)
    for m in bp["medians"]:
        m.set_color("black"); m.set_linewidth(1.6)
    for i, v in enumerate(vals, 1):
        ax[1].text(i, 1.13, f"{np.median(v):.2f}", ha="center", va="bottom", fontsize=9.5, fontweight="bold")
        ax[1].text(i, 1.04, f"n = {len(v):,}", ha="center", va="bottom", fontsize=8, color="#444444")
    ax[1].set_xticks(x)
    ax[1].set_xticklabels([label(*p) for p in present])
    ax[1].set_ylabel("Dice similarity")
    ax[1].set_ylim(-0.03, 1.22)
    ax[1].set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    for a_, l_ in zip(ax, "AB"):
        clean(a_); panel_letter(a_, l_, x=-0.09); group_labels(a_, present)
    save(fig, path)


def hist_scatter_panel(ax_h, ax_s, sub, metric, color, title, nonzero=False):
    v = sub[sub.intersection > 0] if nonzero else sub
    counts, _, _ = ax_h.hist(v[metric], bins=30, color=color, edgecolor="none")
    ax_h.set_ylim(0, (counts.max() if len(counts) else 1) * 1.35)   # headroom for the median label
    med = v[metric].median()
    ax_h.axvline(med, ls="--", color="black", lw=1.2)
    ax_h.text(0.03, 0.95, f"median = {med:.2f}\nn = {len(v):,}", transform=ax_h.transAxes,
              ha="left", va="top", fontsize=9,
              bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.85, "pad": 1.5})
    ax_h.set_title(title)
    ax_h.set_xlabel(f"{metric.capitalize()} coefficient"); ax_h.set_ylabel("Count")
    ax_h.set_xlim(-0.02, 1.02)
    ax_s.scatter(v.intersection, v[metric], alpha=0.35, s=5, color=color, linewidths=0, rasterized=True)
    ax_s.set_xlabel("Shared residues"); ax_s.set_ylabel(f"{metric.capitalize()} coefficient")
    ax_s.set_ylim(-0.03, 1.05)
    clean(ax_h); clean(ax_s)


def fig4(best, path, metric="dice"):
    dsets = [ds for ds in DATASETS if ((best.dataset == ds) & (best.direction == "PDB-PDB")).any()]
    fig, ax = plt.subplots(2, len(dsets), figsize=(FIG_W, 5.6), squeeze=False,
                           gridspec_kw={"hspace": 0.45, "wspace": 0.35})
    k = 0
    for j, ds in enumerate(dsets):
        sub = best[(best.dataset == ds) & (best.direction == "PDB-PDB")]
        hist_scatter_panel(ax[0, j], ax[1, j], sub, metric, COLORS["PDB-PDB"], f"{ds}: PDB→PDB")
    for i in range(2):
        for j in range(len(dsets)):
            panel_letter(ax[i, j], "ABCDEFGH"[k], x=-0.2); k += 1
    save(fig, path)


def fig5(best, path, metric="dice", only=None):
    dsets = [ds for ds in DATASETS if ((best.dataset == ds) & (best.direction != "PDB-PDB")).any()]
    if only is not None:
        dsets = [ds for ds in dsets if ds in only]
    dirs = ["AlphaFill-PDB", "PDB-AlphaFill"]
    fig, ax = plt.subplots(2 * len(dsets), 2, figsize=(FIG_W, 4.4 * len(dsets)), squeeze=False,
                           gridspec_kw={"hspace": 0.6, "wspace": 0.32})
    k = 0
    for i, ds in enumerate(dsets):
        for j, d in enumerate(dirs):
            sub = best[(best.dataset == ds) & (best.direction == d)]
            if len(sub):
                hist_scatter_panel(ax[2 * i, j], ax[2 * i + 1, j], sub, metric, COLORS[d],
                                   f"{ds}: {ARROW[d]}", nonzero=True)
            else:
                ax[2 * i, j].axis("off"); ax[2 * i + 1, j].axis("off")
    for i in range(2 * len(dsets)):
        for j in range(2):
            if ax[i, j].axison:
                panel_letter(ax[i, j], "ABCDEFGH"[k], x=-0.2); k += 1
    save(fig, path)


def fig_si_metric_box(best, metric, path):
    present = present_pairs(best)
    vals = [best[(best.dataset == ds) & (best.direction == d)][metric].values for ds, d in present]
    fig, ax = plt.subplots(figsize=(FIG_W, 3.4))
    bp = ax.boxplot(vals, patch_artist=True, widths=0.6,
                    flierprops={"marker": "o", "markersize": 2, "alpha": 0.4, "markeredgewidth": 0.5})
    for p_, (ds, d) in zip(bp["boxes"], present):
        p_.set_facecolor(COLORS[d]); p_.set_linewidth(0.8)
    for m in bp["medians"]:
        m.set_color("black"); m.set_linewidth(1.6)
    for i, v in enumerate(vals, 1):
        ax.text(i, 1.04, f"{np.median(v):.2f}", ha="center", va="bottom", fontsize=9.5, fontweight="bold")
    ax.set_xticks(range(1, len(present) + 1)); ax.set_xticklabels([label(*p) for p in present])
    ax.set_ylabel(f"{metric.capitalize()} coefficient")
    ax.set_ylim(-0.03, 1.12)
    ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    clean(ax); group_labels(ax, present)
    save(fig, path)


def fig_pair_control(pairs, path):
    dsets = [ds for ds in DATASETS if (pairs.dataset == ds).any()]
    fig, ax = plt.subplots(len(dsets), 2, figsize=(FIG_W, 3.5 * len(dsets)), squeeze=False,
                           gridspec_kw={"hspace": 0.55, "wspace": 0.35})
    k = 0
    for i, ds in enumerate(dsets):
        g = pairs[pairs.dataset == ds].groupby("group")[["pdb_matched", "af_matched",
                                                         "pdb_median_dice", "af_median_dice"]].median()
        for j, (a, b, ttl) in enumerate([("pdb_matched", "af_matched", "Fraction with counterpart"),
                                         ("pdb_median_dice", "af_median_dice", "Median best-match Dice")]):
            ax_ = ax[i, j]
            ax_.plot([0, 1], [0, 1], ls="--", color="grey", lw=0.8)
            ax_.scatter(g[a], g[b], color="#555555", alpha=0.7, s=14, linewidths=0)
            ax_.set_xlabel("PDB→AlphaFill"); ax_.set_ylabel("AlphaFill→PDB")
            ax_.set_title(f"{ds}: {ttl}")
            ax_.set_xlim(-0.02, 1.02); ax_.set_ylim(-0.02, 1.02); ax_.set_aspect("equal")
            clean(ax_); panel_letter(ax_, "ABCDEFGH"[k], x=-0.25); k += 1
    save(fig, path)


# ============================================================
# MAIN
# ============================================================
def main():
    OUT.mkdir(exist_ok=True)
    all_best, all_pairs, all_qc, raw_best, s1_rows, excl = [], [], [], [], [], []
    for ds, (csv, mj, pf) in DATASETS.items():
        df, u2p, infos, excluded = load_dataset(ds, csv, mj, pf)
        excl += [{"dataset": ds, "group": g, "name": n, "reason": w} for g, n, w in excluded]

        # raw mode = original manuscript behaviour (for the response letter)
        if RUN_RAW:
            log("  [raw mode = original behaviour]")
            recs_raw, _, _ = map_dataset(df, u2p, infos, "raw")
            b_raw, _ = run_comparisons(recs_raw, ds)
            raw_best.append(b_raw)

        log("  [SIFTS mode = corrected]")
        recs, qc, dropped = map_dataset(df, u2p, infos, "sifts")
        qc["dataset"] = ds
        all_qc.append(qc)
        log(f"  pathways kept: {len(recs)}; dropped (<{MIN_RESIDUES_IN_GROUP} residues in group protein): "
            f"{sum(dropped.values())}")
        b, pairs = run_comparisons(recs, ds)
        all_best.append(b); all_pairs.append(pairs)

        # Table S1
        rec_df = pd.DataFrame(recs)
        for g, pl in sorted(u2p.items()):
            sub = rec_df[rec_df.group == g] if len(rec_df) else rec_df
            info = infos.get(g, {})
            pdb_used = sorted(sub[sub.source == "PDB"].sid.unique()) if len(sub) else []
            s1_rows.append({"dataset": ds, "uniprot": g, "protein": info.get("name", ""), "gene": info.get("gene", ""),
                            "organism": info.get("organism", ""),
                            "alphafill_model": f"AF-{g}-F1" if (len(sub) and (sub.source == "AlphaFill").any()) else "",
                            "n_pdb_in_group": len(pl), "n_pdb_with_pathways": len(pdb_used),
                            "n_pdb_pathways": int((sub.source == "PDB").sum()) if len(sub) else 0,
                            "n_alphafill_pathways": int((sub.source == "AlphaFill").sum()) if len(sub) else 0,
                            "pdb_ids": " ".join(pdb_used)})

    best = pd.concat(all_best, ignore_index=True)
    pairs = pd.concat(all_pairs, ignore_index=True)
    qc = pd.concat(all_qc, ignore_index=True)
    rawb = pd.concat(raw_best, ignore_index=True) if raw_best else None

    best.to_csv(OUT / "best_matches_all_directions.csv", index=False)
    pairs.to_csv(OUT / "pair_control_1to1.csv", index=False)
    qc.to_csv(OUT / "mapping_qc.csv", index=False)
    pd.DataFrame(s1_rows).to_csv(OUT / "table_S1_identifiers.csv", index=False)
    pd.DataFrame(excl).to_csv(OUT / "excluded_groups.csv", index=False)

    summ = summary_table(best)
    summ.to_csv(OUT / "summary_statistics.csv", index=False)
    if rawb is not None:
        summ_raw = summary_table(rawb)
        comp = summ_raw[["dataset", "direction", "n_pathways", "zero_overlap_pct", "median_dice"]].merge(
            summ[["dataset", "direction", "n_pathways", "zero_overlap_pct", "median_dice"]],
            on=["dataset", "direction"], suffixes=("_raw_strings", "_sifts_corrected"))
        comp.to_csv(OUT / "original_vs_corrected.csv", index=False)
    else:
        comp = summ[["dataset", "direction", "n_pathways", "zero_overlap_pct", "median_dice"]]
    tests = asymmetry_tests(best, pairs)
    tests.to_csv(OUT / "asymmetry_tests.csv", index=False)

    fig3(best, OUT / "Fig3_overlap_and_dice.png")
    fig4(best, OUT / "Fig4_PDB_PDB.png")
    fig5(best, OUT / "Fig5_AlphaFill_vs_PDB.png", only=["CYP"])
    for ds in DATASETS:
        if ds != "CYP":
            fig5(best, OUT / f"FigS_{ds}_AlphaFill_vs_PDB.png", only=[ds])
    for m in ["jaccard", "overlap"]:
        fig_si_metric_box(best, m, OUT / f"FigS_{m}_boxplot.png")
        fig4(best, OUT / f"FigS_{m}_PDB_PDB.png", metric=m)
        fig5(best, OUT / f"FigS_{m}_AlphaFill_vs_PDB.png", metric=m)
    fig_pair_control(pairs, OUT / "FigS_1to1_control.png")

    pd.set_option("display.width", 200)
    log("\n=== SUMMARY ===\n" + comp.round(3).to_string(index=False))
    log("\n=== ASYMMETRY TESTS ===\n" + tests.round(4).to_string(index=False))
    log("\n=== MAPPING QC (PDB chains) ===\n" + qc[~qc.pdb.str.startswith("AlphaFill")]
        .groupby(["dataset", "method"]).size().to_string())
    log(f"\nAll outputs in {OUT.resolve()}")


if __name__ == "__main__":
    main()
