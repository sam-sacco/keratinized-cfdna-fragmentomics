#!/usr/bin/env python3
"""Figure S1: cfDNA coverage around CTCF motifs in feathers (A) and snakeskin (B).

Sites are FIMO hits for JASPAR MA0139.1 with p < 1e-6 (run_fimo.sh), centred on the
motif midpoint. For each library, per-base read depth (pysam count_coverage) is averaged
over sites within +/-1,000 bp, with minus-strand windows reversed to the motif
orientation, and divided by the mean depth of the outermost 250 bp at each end. Down and
flight are each the mean of their libraries (black: all six); snakeskin is the mean of
its libraries. Plotted lines are smoothed with a cubic 51-bp Savitzky-Golay filter.

Sample sheets: TSV of bam_path <TAB> group (down / flight for feathers; any label for
snakeskin). Writes <out>.png/.svg and per-library normalised profiles to <out>.npz.

usage: ctcf_profiles.py --feather-fimo fimo_out/fimo.tsv --feather-bams feather_bams.tsv \
                        --snake-fimo fimo_snake.tsv --snake-bams snake_bams.tsv \
                        [--out figure_ctcf_feather_snakeskin]
"""
import argparse
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pysam
from scipy.signal import savgol_filter

FLANK, BASELINE, PVAL_MAX, SAVGOL = 1000, 250, 1e-6, 51
COLOUR = {"down": "#ED7D31", "flight": "#4472C4"}


def load_sites(fimo_tsv):
    """FIMO hits (1-based) -> 0-based motif midpoints and strands."""
    df = pd.read_csv(fimo_tsv, sep="\t", comment="#")
    df = df[df["p-value"] < PVAL_MAX]
    mid = (df["start"] + df["stop"]) // 2 - 1
    sites = pd.DataFrame({"chrom": df["sequence_name"].astype(str), "mid": mid,
                          "strand": df["strand"]})
    return sites[sites["mid"] >= 0]


def profile(args):
    """Flank-normalised mean depth around the sites for one BAM."""
    bam_path, sites = args
    total, n = np.zeros(2 * FLANK + 1), 0
    with pysam.AlignmentFile(bam_path) as bam:
        refs = set(bam.references)
        for chrom, mid, strand in sites.itertuples(index=False):
            if chrom not in refs:
                continue
            try:
                depth = np.sum(bam.count_coverage(chrom, mid - FLANK, mid + FLANK + 1,
                                                  quality_threshold=0,
                                                  read_callback="nofilter"), axis=0)
            except ValueError:                    # window runs off the contig
                continue
            if depth.size != total.size:
                continue
            total += depth[::-1] if strand == "-" else depth
            n += 1
    total /= n
    return total / np.concatenate([total[:BASELINE], total[-BASELINE:]]).mean()


def profiles(fimo_tsv, sheet, workers):
    sites = load_sites(fimo_tsv)
    samples = pd.read_csv(sheet, sep="\t", header=None, names=["bam", "group"])
    print(f"[sites] {fimo_tsv}: {len(sites):,} CTCF sites (p < {PVAL_MAX:g})")
    with ProcessPoolExecutor(max_workers=workers) as ex:
        profs = list(ex.map(profile, [(b, sites) for b in samples["bam"]]))
    return {(Path(b).name.split(".")[0], g): p for b, g, p in zip(samples["bam"], samples["group"], profs)}


def smooth_mean(profs):
    return savgol_filter(np.mean(profs, axis=0), SAVGOL, 3)


def plot(feather, snake, out):
    plt.rcParams.update({"font.family": "DejaVu Sans", "axes.spines.top": False,
                         "axes.spines.right": False, "axes.titlesize": 19, "axes.labelsize": 17,
                         "xtick.labelsize": 15, "ytick.labelsize": 15})
    x = np.arange(-FLANK, FLANK + 1)
    fig, (ax_f, ax_s) = plt.subplots(1, 2, figsize=(15, 5.6), constrained_layout=True)
    for group, colour in COLOUR.items():
        ax_f.plot(x, smooth_mean([p for (_, g), p in feather.items() if g == group]),
                  color=colour, linewidth=1.8)
    ax_f.plot(x, smooth_mean(list(feather.values())), color="black", linewidth=2.5)
    ax_s.plot(x, smooth_mean(list(snake.values())), color="black", linewidth=2.5)

    ax_f.legend(handles=[plt.Line2D([], [], color=COLOUR["down"], linewidth=3.5, label="Down"),
                         plt.Line2D([], [], color=COLOUR["flight"], linewidth=3.5, label="Flight"),
                         plt.Line2D([], [], color="black", linewidth=5.0, label="Combined")],
                loc="upper right", frameon=False, fontsize=16, handlelength=2.0)
    for ax, title, letter in ((ax_f, "Feather", "A."), (ax_s, "Snakeskin", "B.")):
        ax.set_title(title)
        ax.set_xlim(-FLANK, FLANK)
        ax.axhline(1.0, color="grey", linestyle=":", linewidth=0.9)
        ax.axvline(0, color="grey", linestyle="--", linewidth=0.8)
        ax.set_xlabel("Position relative to CTCF motif (bp)")
        ax.annotate(letter, xy=(0, 1), xycoords="axes fraction", xytext=(-60, 10),
                    textcoords="offset points", fontsize=20, ha="left", va="bottom")
    ax_f.set_ylabel("Coverage / flank background")
    for ext in (".png", ".svg"):
        fig.savefig(f"{out}{ext}", dpi=300)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--feather-fimo", required=True)
    ap.add_argument("--feather-bams", required=True)
    ap.add_argument("--snake-fimo", required=True)
    ap.add_argument("--snake-bams", required=True)
    ap.add_argument("--out", default="figure_ctcf_feather_snakeskin")
    ap.add_argument("--workers", type=int, default=8)
    a = ap.parse_args()

    feather = profiles(a.feather_fimo, a.feather_bams, a.workers)
    snake = profiles(a.snake_fimo, a.snake_bams, a.workers)
    np.savez(f"{a.out}.npz", **{f"{n}__{g}": p for (n, g), p in {**feather, **snake}.items()})
    plot(feather, snake, a.out)


if __name__ == "__main__":
    main()
