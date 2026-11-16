#!/usr/bin/env python3
"""Figure 2 and Table S2: feather cfDNA promoter CG/NCG and TSS coverage by expression septile.

Genes (Ensembl GRCg6a GTF, chrZ/chrW excluded) with mean FPKM > 0 across the feather
follicle samples (Zheng et al., GSE146956) are split into seven equal-size septiles.
  A: per sample and septile, CG and NCG termini are pooled over the septile's promoters
     (TSS -1500..+500 bp, strand-aware) and their ratio plotted.
  B: mean per-base depth in a 2 kb window centred on each TSS (strand-oriented),
     summed over the feather BAMs and averaged over the septile's genes.

Inputs: the GTF, the FPKM table, one whole-genome fragkon file per feather library
(columns CpGID + termini counts, from 02_cpg_termini) and the feather .rmdup.bam files.

usage: septiles.py --gtf Gallus_gallus.GRCg6a.104.gtf \
           --fpkm GSE146956_expression_all_genes.fpkm.txt \
           --fragkons SCS_*.CpG_allCpG.fragkon --bams SCS_*.rmdup.bam --out-dir .
"""
import argparse
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import polars as pl
import pysam

N_BINS = 7
HALF = 1000                         # TSS window is HALF bp either side
SEX_CHROMS = {"Z", "W"}
# fragkon termini columns generating 5'-CG / 5'-NCG ends (5' and 3', both strands)
CG_COLS = ["forStartC", "forEnd5", "revStart8", "revEndG"]
NCG_COLS = ["forStart5", "forEnd4", "revStart9", "revEnd8"]
COLOURS = [matplotlib.colors.to_hex(plt.get_cmap("coolwarm")(i / (N_BINS - 1)))
           for i in range(N_BINS)]
LABELS = [f"S{i}" + (" (lowest)" if i == 1 else " (highest)" if i == N_BINS else "")
          for i in range(1, N_BINS + 1)]


def load_genes(gtf):
    """One promoter window per gene (outermost transcript span), and the chrZ/chrW genes."""
    df = pd.read_csv(gtf, sep="\t", header=None, comment="#", usecols=[0, 2, 3, 4, 6, 8],
                     names=["chrom", "feature", "start", "end", "strand", "attr"],
                     dtype={"chrom": str})
    df = df[df["feature"].isin(["transcript", "exon"])]
    df["gene"] = df["attr"].str.extract(r'gene_id "([^"]+)"', expand=False)
    df["chrom"] = df["chrom"].str.replace(r"(?i)^chr", "", regex=True)
    sex = set(df.loc[df["chrom"].str.upper().isin(SEX_CHROMS), "gene"])

    g = (df[df["feature"] == "transcript"]
         .assign(start=lambda d: d["start"] - 1)                       # to 0-based
         .groupby(["gene", "chrom", "strand"], as_index=False)
         .agg(start=("start", "min"), end=("end", "max")))
    plus = g["strand"] == "+"
    g["tss"] = np.where(plus, g["start"], g["end"] - 1)
    g["prom_lo"] = np.where(plus, g["tss"] - 1500, g["tss"] - 500).clip(min=0)
    g["prom_hi"] = np.where(plus, g["tss"] + 500, g["tss"] + 1500)
    # keep genes extending more than 500 bp past the TSS (a non-empty gene body)
    body = np.where(plus, g["end"] - (g["tss"] + 500), (g["tss"] - 500) - g["start"])
    g = g[body > 0].drop_duplicates("gene")
    return g[~g["chrom"].str.upper().isin(SEX_CHROMS)].reset_index(drop=True), sex


def septiles(fpkm, exclude):
    """Expressed genes (mean FPKM > 0) in N_BINS equal-size bins, lowest first."""
    df = pd.read_csv(fpkm, sep="\t").dropna(subset=["tracking_id"])
    cols = [c for c in df.columns if c.endswith("_FPKM")]
    mean = df.groupby("tracking_id")[cols].mean().mean(axis=1)
    mean = mean[~mean.index.isin(exclude) & (mean > 0)].sort_values()
    return [mean.iloc[i] for i in np.array_split(np.arange(len(mean)), N_BINS)]


def promoter_ratios(fragkon, bins):
    """CG/NCG pooled over each septile's promoter windows for one library."""
    cpg = (pl.read_csv(fragkon, separator="\t", columns=["CpGID"] + CG_COLS + NCG_COLS,
                       infer_schema_length=500)
             .select(pl.col("CpGID").str.extract(r"^(?:chr)?(.+)_(\d+)$", 1).alias("chrom"),
                     pl.col("CpGID").str.extract(r"_(\d+)$", 1).cast(pl.Int64).alias("pos"),
                     pl.sum_horizontal(pl.col(CG_COLS).cast(pl.Float64).fill_null(0)).alias("cg"),
                     pl.sum_horizontal(pl.col(NCG_COLS).cast(pl.Float64).fill_null(0)).alias("ncg"))
             .sort("chrom", "pos").to_pandas())
    by_chrom = {c: (d["pos"].to_numpy(), np.r_[0, d["cg"].cumsum()], np.r_[0, d["ncg"].cumsum()])
                for c, d in cpg.groupby("chrom")}
    ratios = []
    for genes in bins:
        cg = ncg = 0.0
        for chrom, w in genes.groupby("chrom"):
            if chrom not in by_chrom:
                continue
            pos, cg_cum, ncg_cum = by_chrom[chrom]
            lo = np.searchsorted(pos, w["prom_lo"].to_numpy())
            hi = np.searchsorted(pos, w["prom_hi"].to_numpy())
            cg += (cg_cum[hi] - cg_cum[lo]).sum()
            ncg += (ncg_cum[hi] - ncg_cum[lo]).sum()
        ratios.append(cg / ncg)
    return ratios


def tss_depth(bam_path, genes):
    """Strand-oriented per-base depth around each gene's TSS (genes x 2*HALF)."""
    out = np.full((len(genes), 2 * HALF), np.nan)
    with pysam.AlignmentFile(bam_path) as bam:
        size = dict(zip(bam.references, bam.lengths))
        for i, (chrom, strand, tss) in enumerate(zip(genes["chrom"], genes["strand"], genes["tss"])):
            ref = chrom if chrom in size else "chr" + chrom
            if ref in size and tss - HALF >= 0 and tss + HALF <= size[ref]:
                d = np.sum(bam.count_coverage(ref, tss - HALF, tss + HALF), axis=0)
                out[i] = d[::-1] if strand == "-" else d
    return out


def plot(ratios, profiles, out):
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "axes.linewidth": 0.8, "xtick.major.width": 0.8,
                         "ytick.major.width": 0.8})
    fig, (a, b) = plt.subplots(1, 2, figsize=(10, 4.9), constrained_layout=True,
                               gridspec_kw={"width_ratios": [36, 64], "wspace": 0.08})
    pos = np.arange(1, N_BINS + 1)
    box = a.boxplot([ratios[c] for c in ratios.columns], positions=pos, patch_artist=True,
                    widths=0.45, showfliers=False, medianprops=dict(color="white", linewidth=1.8))
    for patch, col in zip(box["boxes"], COLOURS):
        patch.set_facecolor(col)
        patch.set_alpha(0.8)
    jitter = np.linspace(-0.14, 0.14, len(ratios))
    for p, c in zip(pos, ratios.columns):
        a.scatter(p + jitter, ratios[c], color="black", s=22, zorder=5, alpha=0.75, linewidths=0)
    a.set_xticks(pos, [l.split(" ")[0] for l in LABELS])
    a.set_ylabel("CG/NCG Ratio", fontsize=11)
    a.axhline(1.0, color="grey", linestyle=":", linewidth=0.9)

    x = np.arange(-HALF, HALF)
    for prof, lbl, col in zip(profiles, LABELS, COLOURS):
        b.plot(x, prof, color=col, linewidth=1.8, label=lbl)
    b.axvline(0, color="grey", linestyle="--", linewidth=0.8)
    b.set_xlabel("Position Relative to TSS", fontsize=11)
    b.set_ylabel("Mean Depth", fontsize=11)
    b.legend(fontsize=9, frameon=False, title="Expression", title_fontsize=9, loc="upper right")
    for ax, letter in ((a, "A."), (b, "B.")):
        ax.tick_params(axis="both", labelsize=11)
        ax.annotate(letter, xy=(0, 1), xycoords="axes fraction", xytext=(-42, 6),
                    textcoords="offset points", fontsize=16, ha="left", va="bottom")
    for ext in (".png", ".svg"):
        fig.savefig(out.with_suffix(ext), dpi=300)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--gtf", required=True)
    p.add_argument("--fpkm", required=True)
    p.add_argument("--fragkons", nargs="+", required=True)
    p.add_argument("--bams", nargs="+", required=True)
    p.add_argument("--out-dir", default=".")
    p.add_argument("--threads", type=int, default=6)
    args = p.parse_args()
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    genes, sex = load_genes(args.gtf)
    expr = septiles(args.fpkm, sex)
    pd.DataFrame({"septile": [f"S{i}" for i in range(1, N_BINS + 1)],
                  "n_genes": [len(s) for s in expr],
                  "fpkm_min": [s.min() for s in expr], "fpkm_max": [s.max() for s in expr],
                  "fpkm_mean": [s.mean() for s in expr]}
                 ).to_csv(out / "table_s2_septiles.tsv", sep="\t", index=False, float_format="%.4g")
    bins = [genes[genes["gene"].isin(s.index)] for s in expr]

    fragkons = sorted(args.fragkons)
    with ProcessPoolExecutor(args.threads) as ex:
        ratios = list(ex.map(promoter_ratios, fragkons, [bins] * len(fragkons)))
        depth = sum(ex.map(tss_depth, args.bams, [genes] * len(args.bams)))
    ratios = pd.DataFrame(ratios, columns=[f"S{i}" for i in range(1, N_BINS + 1)],
                          index=[Path(f).name.split(".")[0] for f in fragkons])
    ratios.to_csv(out / "promoter_cg_ncg_by_septile.tsv", sep="\t", float_format="%.5f")
    gene_bin = np.full(len(genes), -1)
    for i, b in enumerate(bins):
        gene_bin[b.index] = i
    profiles = [np.nanmean(depth[(gene_bin == i) & ~np.isnan(depth[:, 0])], axis=0)
                for i in range(N_BINS)]
    pd.DataFrame(np.array(profiles).T, index=pd.Index(np.arange(-HALF, HALF), name="pos"),
                 columns=ratios.columns).to_csv(out / "tss_depth_by_septile.tsv", sep="\t",
                                                float_format="%.5f")
    plot(ratios, profiles, out / "figure_2")


if __name__ == "__main__":
    main()
