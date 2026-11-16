#!/usr/bin/env python3
"""Preranked GSEA of feather promoter coverage (down, flight) and follicle FPKM.

Genes are ranked three ways: by fragment midpoints in each gene's promoter window
(TSS -1000/+500 bp) summed over the down or the flight libraries, and by mean FPKM in
the Zheng et al. feather-follicle RNA-seq. Redundant GO BP terms are merged into
higher-order sets, then each ranking is run through gseapy prerank.

Outputs (in --out-dir):
  ranking_{down,flight,fpkm}.tsv     gene rankings (percentile, highest first)
  gsea_{profile}.tsv                 all tested terms
  gsea_{profile}_significant.tsv     FDR < 0.05 terms (Table S4 = down)
  gsea_dotplot.{png,pdf}             Figure 3
  supp_nes_correlation.{png,pdf}     Figure S3
  supp_term_venn.{png,pdf}           Figure S4

usage: gsea.py --bams feather_bams.tsv --genes Galgal6_refseq.csv --gmt GOBP_chicken.gmt \
               --fpkm GSE146956_expression_all_genes.fpkm.txt
  feather_bams.tsv: bam_path <TAB> down|flight;  Galgal6_refseq.csv: UCSC refGene-format table
"""
import argparse
import re
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import gseapy as gp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pysam
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform
from scipy.stats import spearmanr

UP, DOWN = 1000, 500                                  # promoter window around the TSS
SKIP_CONTIGS = {"Z", "W", "M", "MT"}                  # sex chromosomes, mitochondrion
MIN_SIZE, MAX_SIZE, PERMS, SEED, MERGE_AT, FDR = 15, 1000, 1000, 1, 0.7, 0.05
TERM_CHARS = 38                                       # dot plot label width
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9,
                     "axes.spines.top": False, "axes.spines.right": False,
                     "axes.linewidth": 0.8, "xtick.major.width": 0.8,
                     "ytick.major.width": 0.8, "savefig.dpi": 300})


# ------------------------------------------------------------------ promoter coverage

def load_promoters(path):
    """One promoter window per gene: its locus with the most transcripts, TSS at the
    outermost transcript start, strand-aware."""
    df = pd.read_csv(path, sep="\t").rename(columns={"name2": "gene"})
    df = df[~df["chrom"].str.removeprefix("chr").str.upper().isin(SKIP_CONTIGS)]
    loci = (df.groupby(["gene", "chrom", "strand"], as_index=False)
              .agg(n=("txStart", "size"), lo=("txStart", "min"), hi=("txEnd", "max"))
              .sort_values(["gene", "n", "hi"], ascending=[True, False, False])
              .drop_duplicates("gene"))
    plus = (loci["strand"] == "+").to_numpy()
    start = np.where(plus, loci["lo"] - UP, loci["hi"] - 1 - DOWN)
    win = pd.DataFrame({"gene": loci["gene"].to_numpy(), "chrom": loci["chrom"].to_numpy(),
                        "start": start, "end": start + UP + DOWN + 1})
    return win[win["start"] >= 0].sort_values(["chrom", "start"]).reset_index(drop=True)


def midpoint(read, merged_groups):
    """Fragment midpoint: a merged read's centre, or for an unmerged pair the centre of
    the insert, counted once from the leftmost mate."""
    if read.get_tag("RG") in merged_groups:
        return (read.reference_start + read.reference_end) // 2
    return read.reference_start + read.template_length // 2 if read.template_length > 0 else None


def count_bam(bam_path, win):
    """Fragment midpoints per promoter window for one BAM."""
    counts = np.zeros(len(win), dtype=np.int64)
    with pysam.AlignmentFile(bam_path, threads=2) as bam:
        merged = {rg["ID"] for rg in bam.header.to_dict()["RG"] if rg["ID"].endswith("_merged")}
        for chrom, w in win.groupby("chrom"):
            if chrom not in bam.references:
                continue
            starts, ends, idx = w["start"].to_numpy(), w["end"].to_numpy(), w.index.to_numpy()
            for read in bam.fetch(chrom):
                if (read.is_secondary or read.is_supplementary or read.is_duplicate
                        or read.is_qcfail or read.reference_end is None):
                    continue
                mid = midpoint(read, merged)
                if mid is not None:
                    # fixed-width, start-sorted windows: both edges are a bisect
                    counts[idx[np.searchsorted(ends, mid, "right"):
                               np.searchsorted(starts, mid, "right")]] += 1
    return counts


def rank_genes(values, genes):
    """Percentile within the profile (ties averaged), highest first."""
    return (pd.DataFrame({"gene": genes, "score": pd.Series(values).rank(pct=True) * 100})
              .groupby("gene", as_index=False)["score"].mean()
              .sort_values("score", ascending=False).reset_index(drop=True))


def mean_fpkm(path, universe):
    """Mean FPKM per gene symbol (top locus per symbol), genes with a promoter window."""
    df = pd.read_csv(path, sep="\t", low_memory=False)
    df["gene"] = df["gene_short_name"].astype(str).str.strip()
    df["fpkm"] = df.filter(like="_FPKM").mean(axis=1)
    df = (df[~df["gene"].isin(["-", "nan"])].sort_values("fpkm", ascending=False)
            .drop_duplicates("gene").dropna(subset=["fpkm"]))
    return df[df["gene"].isin(universe)]


# ------------------------------------------------------------------ gene sets and GSEA

def read_gmt(path):
    return {f[0]: [g for g in f[2:] if g]
            for f in (line.rstrip("\n").split("\t") for line in open(path)) if len(f) > 2}


def merge_sets(sets, universe):
    """Merge redundant GO terms by gene membership. Similarity is the mean of the
    Jaccard and overlap coefficients; average-linkage clusters cut at MERGE_AT take the
    union of their members' genes and the name of their medoid term."""
    sets = {n: g for n, g in sets.items() if MIN_SIZE <= len(set(g) & universe) <= MAX_SIZE}
    names = list(sets)
    col = {g: i for i, g in enumerate(sorted({g for v in sets.values() for g in v} & universe))}
    M = np.zeros((len(names), len(col)), dtype=np.float32)
    for r, n in enumerate(names):
        M[r, [col[g] for g in sets[n] if g in col]] = 1.0

    inter = M @ M.T
    size = np.diag(inter).copy()
    union = size[:, None] + size[None, :] - inter
    smaller = np.minimum(size[:, None], size[None, :])
    with np.errstate(invalid="ignore", divide="ignore"):
        S = (np.where(union > 0, inter / union, 0) + np.where(smaller > 0, inter / smaller, 0)) / 2
    S = np.clip(S, 0, 1)
    S = (S + S.T) / 2
    np.fill_diagonal(S, 1.0)

    labels = fcluster(linkage(squareform(1 - S, checks=False), method="average"),
                      t=1 - MERGE_AT, criterion="distance")
    merged = {}
    for cl in np.unique(labels):
        members = np.flatnonzero(labels == cl)
        medoid = names[members[np.argmax(S[np.ix_(members, members)].mean(axis=1))]]
        name = medoid if members.size == 1 else f"{medoid} (+{members.size - 1} terms)"
        merged[name] = sorted({g for i in members for g in sets[names[i]]})
    return merged


def run_gsea(rank, sets):
    res = gp.prerank(rnk=rank, gene_sets=sets, min_size=MIN_SIZE, max_size=MAX_SIZE,
                     permutation_num=PERMS, seed=SEED, outdir=None, no_plot=True,
                     verbose=False).res2d
    tbl = res.rename(columns={"FDR q-val": "fdr", "NOM p-val": "pval",
                              "Lead_genes": "leading_edge"})
    tbl[["go_id", "term"]] = tbl["Term"].str.split("|", n=1, expand=True)
    for c in ("NES", "ES", "fdr", "pval"):
        tbl[c] = pd.to_numeric(tbl[c])
    lead, size = tbl["Tag %"].str.split("/", expand=True).astype(float).T.to_numpy()
    tbl["pct_in_set"] = 100 * lead / size             # share of the set in the leading edge
    return (tbl[["go_id", "term", "ES", "NES", "pval", "fdr", "pct_in_set", "leading_edge"]]
              .sort_values("fdr").reset_index(drop=True))


def clean_term(term):
    return re.sub(r"\s*\(\+\d+ terms\)$", "", str(term))


def write_significant(tbl, path):
    sig = tbl[tbl["fdr"] < FDR].copy()
    sig["term"] = sig["term"].map(clean_term)
    sig.insert(2, "direction", np.where(sig["NES"] > 0, "up", "down"))
    sig.sort_values("NES", ascending=False).to_csv(path, sep="\t", index=False,
                                                    float_format="%.4g")


# ------------------------------------------------------------------ figures

def save(fig, stem):
    fig.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def plot_dots(tables, out, top_n=10):
    """Figure 3: flight's top upregulated terms (all also significant in down), down on
    the same rows; terms also upregulated in the follicle FPKM are bold and starred."""
    def scored(tbl):
        s = tbl.set_index("go_id")
        s["logq"] = np.log10(1 / s["fdr"].clip(lower=1 / (PERMS + 1)))
        s["combined"] = s["NES"] * s["logq"]
        return s

    lead = scored(tables["flight"])
    top = lead[(lead["NES"] > 0) & (lead["fdr"] < FDR)].nlargest(top_n, "combined")
    ids = list(top.sort_values("combined").index)                 # strongest at the top
    fpkm = scored(tables["fpkm"])
    starred = {g for g in ids if g in fpkm.index
               and fpkm.at[g, "fdr"] < FDR and fpkm.at[g, "NES"] > 0}

    panels = {"Flight Feathers": scored(tables["flight"]).reindex(ids),
              "Down Feathers": scored(tables["down"]).reindex(ids)}
    allv = pd.concat(panels.values())
    vmin, vmax, smax = allv["logq"].min(), allv["logq"].max(), allv["pct_in_set"].max()

    fig = plt.figure(figsize=(6.6 * 2 + 2.0, 0.46 * len(ids) + 3.4), layout="constrained")
    gs = fig.add_gridspec(2, 3, width_ratios=[1, 1, 0.34], height_ratios=[1, 1])
    axes = [fig.add_subplot(gs[:, i]) for i in range(2)]
    axes[1].sharey(axes[0])
    y = np.arange(len(ids))
    for ax, (title, s) in zip(axes, panels.items()):
        dots = ax.scatter(s["combined"], y, c=s["logq"], s=30 + 320 * s["pct_in_set"] / smax,
                          cmap="autumn_r", vmin=vmin, vmax=vmax, edgecolors="none")
        ax.set_title(title, fontsize=18)
        ax.tick_params(labelsize=12)
        ax.grid(axis="y", lw=0.4, color="#e8e8e8")
        ax.set_axisbelow(True)
        lo, hi = min(0.0, s["combined"].min()), s["combined"].max()
        pad = 0.12 * (hi - lo)
        ax.set_xlim(lo - (pad if lo < 0 else 0), hi + pad)
    axes[1].tick_params(labelleft=False)

    labels = [clean_term(top.at[g, "term"]) for g in ids]
    labels = [t if len(t) <= TERM_CHARS else t[:TERM_CHARS - 1].rstrip() + "…" for t in labels]
    axes[0].set_yticks(y, [t + (" *" if g in starred else "") for t, g in zip(labels, ids)],
                       fontsize=15)
    for lab, g in zip(axes[0].get_yticklabels(), ids):
        if g in starred:
            lab.set_fontweight("bold")

    cell = fig.add_subplot(gs[0, 2]); cell.axis("off")
    bar = fig.colorbar(dots, cax=cell.inset_axes([0.10, 0.06, 0.20, 0.88]))
    bar.ax.tick_params(labelsize=11)
    bar.set_label("log$_{10}$ 1/FDR", size=13)
    keys = np.unique(np.round(np.linspace(allv["pct_in_set"].min(), smax, 3)))
    key_ax = fig.add_subplot(gs[1, 2]); key_ax.axis("off")
    key_ax.legend([plt.scatter([], [], s=30 + 320 * k / smax, color="#999999") for k in keys],
                  [f"{k:.0f}" for k in keys], title="% genes\nin set", loc="upper left",
                  frameon=False, fontsize=12, title_fontsize=12, labelspacing=1.6, borderpad=0.6)

    # one x label centred on the two panels (not the figure, which includes the keys)
    fig.canvas.draw()
    fig.set_layout_engine("none")
    boxes = [a.get_position() for a in axes]
    bottom = min(a.get_tightbbox(fig.canvas.get_renderer()).y0 for a in axes)
    fig.text((boxes[0].x0 + boxes[1].x1) / 2, bottom / fig.bbox.height - 0.012,
             "Combined score  (NES x log$_{10}$ 1/FDR)", ha="center", va="top", fontsize=15)
    save(fig, out)


def plot_nes_correlation(tables, out):
    """Figure S3: NES of every shared term, one panel per pair of profiles."""
    pairs = [("down", "flight"), ("down", "fpkm"), ("flight", "fpkm")]
    nes = {n: t.set_index("go_id")["NES"] for n, t in tables.items()}
    fig, axes = plt.subplots(1, 3, figsize=(13.8, 4.6), layout="constrained")
    for ax, (a, b) in zip(axes, pairs):
        shared = nes[a].index.intersection(nes[b].index).sort_values()
        x, y = nes[a][shared].to_numpy(), nes[b][shared].to_numpy()
        lim = np.abs(np.concatenate([x, y])).max() * 1.08
        ax.axhline(0, color="#d0d0d0", lw=0.6)
        ax.axvline(0, color="#d0d0d0", lw=0.6)
        ax.plot([-lim, lim], [-lim, lim], color="#c94040", lw=0.8, ls="--")
        ax.scatter(x, y, s=9, alpha=0.45, color="#4a4a4a", lw=0, rasterized=True)
        ax.set(xlim=(-lim, lim), ylim=(-lim, lim), aspect="equal")
        ax.set_xlabel(f"NES, {a}", fontsize=11)
        ax.set_ylabel(f"NES, {b}", fontsize=11)
        ax.set_title(f"{a.capitalize()} vs {b.capitalize()}\nSpearman rho = "
                     f"{spearmanr(x, y).statistic:+.3f}  ({len(shared):,} terms)", fontsize=11)
    fig.suptitle("NES agreement between profiles (dashed line = equality)", fontsize=12)
    save(fig, out)


def plot_venn(tables, out):
    """Figure S4: three-way overlap of significant terms, up- and downregulated."""
    centres = ((0.36, 0.42), (0.64, 0.42), (0.50, 0.66))
    colours = ("#4472C4", "#ED7D31", "#4E9A51")
    fig, axes = plt.subplots(1, 2, figsize=(11, 5.6), layout="constrained")
    for ax, (label, sign) in zip(axes, (("Upregulated", 1), ("Downregulated", -1))):
        a, b, c = ({t.go_id for t in tbl.itertuples() if t.fdr < FDR and sign * t.NES > 0}
                   for tbl in tables.values())
        for xy, col in zip(centres, colours):
            ax.add_patch(plt.Circle(xy, 0.26, alpha=0.35, color=col, lw=0))
        regions = {(0.22, 0.34): a - b - c, (0.78, 0.34): b - a - c, (0.50, 0.82): c - a - b,
                   (0.50, 0.30): (a & b) - c, (0.31, 0.60): (a & c) - b,
                   (0.69, 0.60): (b & c) - a, (0.50, 0.50): a & b & c}
        for (x, y), members in regions.items():
            ax.text(x, y, f"{len(members):,}", ha="center", va="center",
                    fontsize=12, fontweight="bold")
        for (x, y), name, col in zip(((0.16, 0.18), (0.84, 0.18), (0.50, 0.95)),
                                     tables, colours):
            ax.text(x, y, name.capitalize(), ha="center", va="center", fontsize=12,
                    fontweight="bold", color=col)
        ax.set(xlim=(0, 1), ylim=(0.05, 1.02), aspect="equal", title=f"{label}  (FDR < {FDR})")
        ax.axis("off")
    save(fig, out)


# ------------------------------------------------------------------ driver

def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--bams", default="feather_bams.tsv", help="TSV: bam_path <TAB> down|flight")
    p.add_argument("--genes", default="Galgal6_refseq.csv", help="UCSC refGene-format table")
    p.add_argument("--gmt", default="GOBP_chicken.gmt")
    p.add_argument("--fpkm", default="GSE146956_expression_all_genes.fpkm.txt")
    p.add_argument("--out-dir", default="gsea_out")
    p.add_argument("--threads", type=int, default=6)
    args = p.parse_args()
    out = Path(args.out_dir)
    out.mkdir(exist_ok=True)

    bams = pd.read_csv(args.bams, sep="\t", header=None, names=["bam", "category"])
    win = load_promoters(args.genes)
    genes = win["gene"].to_numpy()
    print(f"[genes] {len(win):,} promoter windows")
    sets = merge_sets(read_gmt(args.gmt), set(genes))

    with ProcessPoolExecutor(args.threads) as ex:
        counts = list(ex.map(count_bam, bams["bam"], [win] * len(bams)))
    rankings = {cat: rank_genes(np.sum([c for c, k in zip(counts, bams["category"]) if k == cat],
                                       axis=0), genes)
                for cat in ("down", "flight")}
    fpkm = mean_fpkm(args.fpkm, set(genes))
    rankings["fpkm"] = rank_genes(fpkm["fpkm"].to_numpy(), fpkm["gene"].to_numpy())

    tables = {}
    for name, rank in rankings.items():
        rank.to_csv(out / f"ranking_{name}.tsv", sep="\t", index=False, float_format="%.4f")
        tables[name] = run_gsea(rank, sets)
        tables[name].to_csv(out / f"gsea_{name}.tsv", sep="\t", index=False, float_format="%.4g")
        write_significant(tables[name], out / f"gsea_{name}_significant.tsv")
        print(f"[gsea] {name}: {len(rank):,} genes, {len(tables[name]):,} sets, "
              f"{(tables[name]['fdr'] < FDR).sum()} at FDR < {FDR}")

    plot_dots(tables, out / "gsea_dotplot")
    plot_nes_correlation(tables, out / "supp_nes_correlation")
    plot_venn(tables, out / "supp_term_venn")


if __name__ == "__main__":
    main()
