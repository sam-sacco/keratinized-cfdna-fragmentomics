#!/usr/bin/env python3
"""Figure 1A/1B: nuclear fragment length distributions of feather and snakeskin cfDNA.

Fragment length is the merged-read length, or the insert size of an unmerged pair
(counted once, from the forward mate). Unmapped, secondary, supplementary and
duplicate alignments are skipped, as are mitochondrial ones (chrM/MT). Prints the mean
length per group (30-599 bp) and the size under which 99% of fragments fall.

Input: TSV of bam<TAB>group, group in {down, flight, snakeskin} (see samples.tsv).
usage: fragment_lengths.py samples.tsv [--out figure_frag_panel]
"""
import argparse
from concurrent.futures import ProcessPoolExecutor

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pysam

LO, HI = 30, 600                  # histogram window (bp)
PLOT_LO, PLOT_HI = 30, 150        # plotted window (bp)
COLOUR = {"flight": "#4472C4", "down": "#ED7D31"}


def length_histogram(bam_path):
    """Counts per length in [LO, HI), plus the fragments below and at/above it."""
    hist = np.zeros(HI - LO, dtype=np.int64)
    below = above = 0
    with pysam.AlignmentFile(bam_path) as bam:
        mito = {"chrM", "MT"} & set(bam.references)
        for r in bam.fetch():
            if (r.is_unmapped or r.is_secondary or r.is_supplementary or r.is_duplicate
                    or r.reference_name in mito):
                continue
            if r.is_paired:
                if r.is_reverse or r.mate_is_unmapped:
                    continue
                n = abs(r.template_length)
            else:
                n = r.query_length or 0
            if n <= 0:
                continue
            if n < LO:
                below += 1
            elif n >= HI:
                above += 1
            else:
                hist[n - LO] += 1
    return hist, below, above


def cutoff(hist, below, above, q=0.99):
    """Smallest length under which more than a fraction q of all fragments fall."""
    cum = below + np.cumsum(hist)
    i = int(np.searchsorted(cum, q * (hist.sum() + below + above), side="right"))
    return LO + i + 1 if i < hist.size else None


def panel(ax, curves, combined, title, xlabel):
    x = np.arange(LO, HI)
    keep = (x >= PLOT_LO) & (x <= PLOT_HI)
    for colour, h in curves:
        ax.plot(x[keep], (h / h.sum())[keep], color=colour, linewidth=1.5, alpha=0.8)
    ax.plot(x[keep], (combined / combined.sum())[keep], color="black", linewidth=2.5, zorder=10)
    if xlabel:
        ax.set_xlabel("Fragment length (bp)", fontsize=17 * 1.15)
    ax.set_ylabel("Density", fontsize=17 * 1.15)
    ax.set_xlim(PLOT_LO, PLOT_HI)
    ax.set_title(title)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("samples")
    ap.add_argument("--out", default="figure_frag_panel")
    args = ap.parse_args()

    samples = pd.read_csv(args.samples, sep="\t", header=None, names=["bam", "group"])
    with ProcessPoolExecutor(max_workers=8) as ex:
        samples["hist"], samples["below"], samples["above"] = zip(*ex.map(length_histogram, samples["bam"]))

    x = np.arange(LO, HI)
    groups = {g: sub for g, sub in samples.groupby("group")}
    feather = samples[samples.group != "snakeskin"]
    for name, sub in [*groups.items(), ("feathers (pooled)", feather)]:
        h = np.sum(sub["hist"].tolist(), axis=0)
        print(f"{name}: n = {h.sum():,}, mean = {(x * h).sum() / h.sum():.1f} bp, "
              f"99% shorter than {cutoff(h, sub['below'].sum(), sub['above'].sum())} bp")

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 18, "axes.titlesize": 19,
                         "xtick.labelsize": 15, "ytick.labelsize": 15, "savefig.dpi": 300,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "axes.linewidth": 0.8, "xtick.major.width": 0.8, "ytick.major.width": 0.8})
    fig, (ax_f, ax_s) = plt.subplots(2, 1, figsize=(7, 10), constrained_layout=True)
    curves = [(COLOUR[g], h) for g in ("flight", "down") for h in groups[g]["hist"]]
    panel(ax_f, curves, np.sum(feather["hist"].tolist(), axis=0), "Feather", xlabel=False)
    panel(ax_s, [], np.sum(groups["snakeskin"]["hist"].tolist(), axis=0), "Snakeskin", xlabel=True)
    ax_f.legend(handles=[plt.Line2D([], [], color=COLOUR["flight"], linewidth=3.5, label="Flight"),
                         plt.Line2D([], [], color=COLOUR["down"], linewidth=3.5, label="Down"),
                         plt.Line2D([], [], color="black", linewidth=5.0, label="Combined")],
                loc="upper right", frameon=False, fontsize=20, handlelength=2.0)
    for ext in (".png", ".svg"):
        fig.savefig(args.out + ext)


if __name__ == "__main__":
    main()
