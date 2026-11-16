#!/usr/bin/env python3
"""CG/NCG ratio by genomic compartment: Figure 1C and Figure S2.

Per CpG, CG counts termini whose fragmentation leaves a 5'-CG end (forward starts at
the C, forward ends just before it, reverse starts just after the G, reverse ends at
the G) and NCG those one base further out (5'-NCG). Each sample's ratio pools all CpGs
of a compartment; its 95% interval comes from resampling CpGs with replacement (2,000
replicates). The two snakeskin libraries are pooled into one sample.

Fig 1C: per-sample ratios by compartment, down / flight / snakeskin, with two-sided
Wilcoxon rank-sum p-values between compartments (down and flight pooled, n = 6).
Fig S2: one row per sample, one marker per compartment, bootstrap intervals.

Input: TSV with columns sample, group (down|flight|snakeskin), compartment
(whole_genome|CR1|CpG_island), fragkon (path from parse_termini.py).

usage: cg_ncg_compartments.py samples.tsv --out-dir results
"""
import argparse
import zlib
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from scipy.stats import ranksums

CG = ["forStartC", "forEnd5", "revStart8", "revEndG"]
NCG = ["forStart5", "forEnd4", "revStart9", "revEnd8"]
REGIONS = {"whole_genome": "Whole Genome", "CR1": "CR1 Elements", "CpG_island": "CpG Islands"}
COLOUR = {"down": "#ED7D31", "flight": "#4472C4", "snakeskin": "black"}
LABEL = {"down": "Down", "flight": "Flight", "snakeskin": "Snakeskin"}
MARKER = {"whole_genome": "o", "CR1": "s", "CpG_island": "^"}
N_BOOT = 2000


def cpg_counts(paths):
    """Per-CpG (cg, ncg) pairs over the given fragkon files, as {pair: number of CpGs}."""
    parts = []
    for path in paths:
        for chunk in pd.read_csv(path, sep="\t", usecols=CG + NCG, chunksize=2_000_000):
            parts.append(pd.DataFrame({"cg": chunk[CG].sum(axis=1),
                                       "ncg": chunk[NCG].sum(axis=1)}).value_counts())
    return pd.concat(parts).groupby(level=["cg", "ncg"]).sum()


def ratio_ci(counts, seed):
    """CG/NCG and its 95% CpG-bootstrap interval. Drawing N CpGs with replacement is a
    multinomial over the distinct (cg, ncg) pairs, so the bootstrap is exact."""
    pairs = np.array(counts.index.tolist(), dtype=float)
    mult = counts.to_numpy(dtype=np.int64)
    n = mult.sum()
    draws = np.random.default_rng(seed).multinomial(n, mult / n, size=N_BOOT)
    boot = (draws @ pairs[:, 0]) / (draws @ pairs[:, 1])
    return (pairs[:, 0] @ mult) / (pairs[:, 1] @ mult), *np.percentile(boot, [2.5, 97.5])


def sample_ratios(sheet):
    rows = []
    for (group, sample, region), files in sheet.groupby(["group", "sample", "compartment"]):
        counts = cpg_counts(files["fragkon"])
        ratio, lo, hi = ratio_ci(counts, zlib.crc32(sample.encode()))
        rows.append(dict(group=group, sample=sample, region=region, ratio=ratio, lo=lo, hi=hi,
                         n_cpg=int(counts.sum())))
        print(f"{sample} {region}: {ratio:.4f} [{lo:.4f}, {hi:.4f}]", flush=True)
    return pd.DataFrame(rows)


def compartment_tests(res):
    feathers = res[res.group.isin(["down", "flight"])]
    regions = list(REGIONS)
    return [(a, b, ranksums(feathers.ratio[feathers.region == a],
                            feathers.ratio[feathers.region == b]).pvalue)
            for i, a in enumerate(regions) for b in regions[i + 1:]]


def figure_1c(res, out):
    """Down, flight and snakeskin slots per compartment; samples staggered in order."""
    fig, ax = plt.subplots(figsize=(7.2, 8.6))
    width, gap, group_gap = 0.35, 0.10, 0.45
    slot = {}
    x = 1.0
    for region in REGIONS:
        for k, group in enumerate(("down", "flight", "snakeskin")):
            slot[region, group] = x + k * (width + gap)
        x = slot[region, "snakeskin"] + width + group_gap
    for (region, group), sub in res.sort_values("sample").groupby(["region", "group"]):
        offsets = np.linspace(-0.10, 0.10, len(sub)) if len(sub) > 1 else [0.0]
        for off, r in zip(offsets, sub.itertuples()):
            xp = slot[region, group] + off
            ax.errorbar([xp], [r.ratio], yerr=[[r.ratio - r.lo], [r.hi - r.ratio]], fmt="none",
                        ecolor=COLOUR[group], elinewidth=1.2, capsize=3, capthick=1.2, zorder=6)
            ax.scatter([xp], [r.ratio], color=COLOUR[group], s=95, alpha=0.9,
                       edgecolor="black", linewidth=0.5, zorder=5)

    ax.set_xticks([slot[r, "flight"] for r in REGIONS], list(REGIONS.values()), fontsize=14)
    ax.set_ylabel("CG/NCG Ratio", fontsize=14)
    ax.tick_params(axis="y", labelsize=13)
    ax.set_xlim(slot["whole_genome", "down"] - 0.5, slot["CpG_island", "snakeskin"] + 0.5)
    ax.spines[["top", "right"]].set_visible(False)
    ax.axhline(1.0, color="gray", linestyle="--", linewidth=1, alpha=0.7)
    ax.legend(handles=[Line2D([], [], linestyle="none", marker="o", markersize=14,
                              markerfacecolor=COLOUR[g], markeredgecolor="black",
                              markeredgewidth=0.5, alpha=0.9, label=LABEL[g]) for g in COLOUR],
              loc="center right", frameon=False, fontsize=17)

    # brackets between compartments, shortest lowest, labelled with the p-value
    ymin, ymax = res[["ratio", "lo", "hi"]].min().min(), res[["ratio", "lo", "hi"]].max().max()
    span = ymax - ymin
    order = list(REGIONS)
    tests = sorted(compartment_tests(res),
                   key=lambda t: (order.index(t[1]) - order.index(t[0]), order.index(t[0])))
    base, step, tick = ymax + 0.03 * span, 0.057 * span, 0.015 * span
    for tier, (a, b, p) in enumerate(tests):
        x1 = (slot[a, "down"] + slot[a, "flight"]) / 2
        x2 = (slot[b, "down"] + slot[b, "flight"]) / 2
        y = base + tier * step
        ax.plot([x1, x1, x2, x2], [y, y + tick, y + tick, y], lw=1.2, color="black")
        ax.text((x1 + x2) / 2, y + tick + 0.002 * span, f"p = {p:.2g}" + (" *" if p < 0.05 else ""),
                ha="center", va="bottom", fontsize=11)
    ax.set_ylim(ymin - 0.08 * span, base + (len(tests) - 1) * step + tick + 0.05 * span)
    fig.tight_layout()
    fig.savefig(out.with_suffix(".png"), dpi=200)
    fig.savefig(out.with_suffix(".svg"))
    plt.close(fig)


def figure_s2(res, out):
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "axes.linewidth": 0.8, "xtick.major.width": 0.8,
                         "ytick.major.width": 0.8, "savefig.dpi": 300})
    rank = {"down": 0, "flight": 1, "snakeskin": 2}
    rows = sorted({(rank[g], s, g) for g, s in zip(res.group, res["sample"])})
    y_of = {s: len(rows) - 1 - i for i, (_, s, _) in enumerate(rows)}
    fig, ax = plt.subplots(figsize=(8.0, 0.62 * len(rows) + 2.6), layout="constrained")
    for r in res.itertuples():
        y = y_of[r.sample]
        ax.plot([r.lo, r.hi], [y, y], color=COLOUR[r.group], lw=1.2, alpha=0.9, zorder=2,
                solid_capstyle="butt")
        ax.scatter([r.ratio], [y], marker=MARKER[r.region], s=70, color=COLOUR[r.group],
                   edgecolors="white", linewidths=0.6, zorder=3)
    for i in range(1, len(rows)):
        if rows[i][2] != rows[i - 1][2]:
            ax.axhline(len(rows) - 0.5 - i, color="#333333", lw=1.4, zorder=1)
    ax.set_yticks([y_of[s] for _, s, _ in rows], [s for _, s, _ in rows], fontsize=11)
    ax.set_xlabel("CG/NCG ratio", fontsize=13)
    ax.set_ylim(-0.6, len(rows) - 0.4)
    ax.tick_params(labelsize=11)
    ax.grid(axis="x", lw=0.4, color="#e8e8e8")
    ax.set_axisbelow(True)
    ax.set_title("CG/NCG by genomic compartment\nbars = 95% bootstrap interval over CpGs",
                 fontsize=13)
    handles = ([Line2D([], [], marker=MARKER[r], color="#555555", lw=0, markersize=8, label=n)
                for r, n in REGIONS.items()]
               + [Line2D([], [], marker="o", color=COLOUR[g], lw=0, markersize=8,
                         label=f"{LABEL[g]} feather" if g != "snakeskin" else LABEL[g])
                  for g in COLOUR])
    ax.legend(handles=handles, loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False,
              fontsize=10)
    for ext in (".png", ".svg"):
        fig.savefig(out.with_suffix(ext), bbox_inches="tight")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("samples")
    ap.add_argument("--out-dir", default="results")
    args = ap.parse_args()
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    res = sample_ratios(pd.read_csv(args.samples, sep="\t"))
    res.to_csv(out / "cg_ncg_by_compartment.tsv", sep="\t", index=False, float_format="%.6f")
    for a, b, p in compartment_tests(res):
        print(f"Wilcoxon rank-sum {a} vs {b} (down + flight, n = 6): p = {p:.4g}")
    figure_1c(res, out / "figure_1c_cg_ncg_by_compartment")
    figure_s2(res, out / "figure_s2_cg_ncg_per_sample")


if __name__ == "__main__":
    main()
