#!/usr/bin/env python3
"""Figure S5: Spearman correlation of the feather promoter-coverage rankings (down,
flight) with embryonic chicken skin RNA-seq (Liang et al., GEO GSE136224).

Each RNA-seq replicate (GEO coverage bigWigs, galGal4) scores a gene as its coverage
averaged over the union of its RefSeq exons (UCSC galGal4 refGene). Replicates are turned
into percentile ranks and averaged per tissue; genes with no coverage in any replicate of
a tissue are dropped from it. Genes are matched to the feather rankings by symbol and
restricted to those present in both rankings and all five tissues.

usage: correlogram.py --refgene refGene.txt.gz --down ranking_down.tsv \
                      --flight ranking_flight.tsv [--bigwigs DIR] [--out-dir .]
"""
import argparse
import gzip
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import pyBigWig

HERE = Path(__file__).resolve().parent
TISSUE = {"E7 dorsal back skin epidermis": "E7 dorsal epidermis",       # plotting order
          "E9 dorsal back skin epidermis": "E9 dorsal epidermis",
          "E14 non-keratinized feather filament": "E14 feather filament",
          "E9 leg scale epidermis": "E9 leg scale epidermis",
          "E14 leg scale skin": "E14 leg scale skin"}


def merged_exons(path):
    """Gene symbol -> (chrom, merged exon intervals); a gene on several contigs keeps the
    contig carrying most of its transcripts."""
    per = defaultdict(lambda: defaultdict(list))
    with gzip.open(path, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            starts = map(int, f[9].rstrip(",").split(","))
            ends = map(int, f[10].rstrip(",").split(","))
            per[f[12]][f[2]].append(list(zip(starts, ends)))
    out = {}
    for gene, by_chrom in per.items():
        chrom = max(by_chrom, key=lambda c: len(by_chrom[c]))
        iv = sorted(x for tx in by_chrom[chrom] for x in tx)
        merged = [list(iv[0])]
        for s, e in iv[1:]:
            if s <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], e)
            else:
                merged.append([s, e])
        out[gene] = (chrom, merged)
    return out


def gene_expression(bw_path, genes):
    """Mean coverage over each gene's merged exons."""
    bw = pyBigWig.open(str(bw_path))
    chroms = bw.chroms()
    vals = {}
    for gene, (chrom, exons) in genes.items():
        if chrom not in chroms:
            continue
        tot = cov = 0.0
        for s, e in exons:
            e = min(e, chroms[chrom])
            if e > s:
                cov += (bw.stats(chrom, s, e, type="mean")[0] or 0.0) * (e - s)
                tot += e - s
        if tot:
            vals[gene] = cov / tot
    bw.close()
    return pd.Series(vals)


def tissue_ranks(samples, bigwigs, genes):
    """Gene x tissue: replicate percentile ranks averaged; NaN where no replicate has coverage."""
    reps = defaultdict(list)
    for r in samples.itertuples():
        reps[TISSUE[r.source]].append(gene_expression(bigwigs / r.url.rsplit("/", 1)[1], genes))
    rank = {t: pd.concat([e.rank(pct=True) for e in v], axis=1).mean(axis=1) for t, v in reps.items()}
    covered = {t: pd.concat([e > 0 for e in v], axis=1).fillna(False).any(axis=1) for t, v in reps.items()}
    return pd.DataFrame(rank).where(pd.DataFrame(covered))[list(TISSUE.values())]


def plot(rho, n, out):
    fig, ax = plt.subplots(figsize=(9, 3.3), constrained_layout=True)
    im = ax.imshow(rho, cmap="viridis", vmin=0, vmax=rho.values.max())
    for i, f in enumerate(rho.index):
        best = rho.loc[f].idxmax()
        for j, c in enumerate(rho.columns):
            ax.text(j, i, f"{rho.iat[i, j]:.2f}", ha="center", va="center", fontsize=12,
                    color="black" if rho.iat[i, j] > 0.6 * rho.values.max() else "white",
                    fontweight="bold" if c == best else "normal")
    ax.set_xticks(range(rho.shape[1]), [c.replace(" epidermis", "\nepidermis")
                                        .replace(" skin", "\nskin").replace(" filament", "\nfilament")
                                        for c in rho.columns], fontsize=10)
    ax.set_yticks(range(len(rho)), rho.index, fontsize=12)
    ax.set_ylabel("Feather promoter coverage", fontsize=12)
    ax.set_title(f"Spearman ρ, feather cfDNA vs reference RNA-seq  (n = {n:,} genes)")
    for s in ax.spines.values():
        s.set_visible(False)
    fig.colorbar(im, ax=ax, shrink=0.9, label="Spearman ρ")
    for ext in (".png", ".svg"):
        fig.savefig(out.with_suffix(ext), dpi=200)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--refgene", required=True, help="UCSC galGal4 refGene.txt.gz")
    p.add_argument("--down", required=True, help="down feather ranking (gene, score)")
    p.add_argument("--flight", required=True, help="flight feather ranking (gene, score)")
    p.add_argument("--samples", default=HERE / "samples.tsv", type=Path)
    p.add_argument("--bigwigs", default=Path("."), type=Path,
                   help="folder holding the bigWigs named in samples.tsv")
    p.add_argument("--out-dir", default=".", type=Path)
    a = p.parse_args()

    samples = pd.read_csv(a.samples, sep="\t", header=None, names=["gsm", "source", "url"])
    ref = tissue_ranks(samples, a.bigwigs, merged_exons(a.refgene)).dropna()
    feathers = {f: pd.read_csv(path, sep="\t").set_index("gene")["score"]
                for f, path in (("Down", a.down), ("Flight", a.flight))}
    genes = ref.index.intersection(feathers["Down"].index)
    rho = pd.DataFrame({c: {f: ref.loc[genes, c].corr(r.loc[genes], method="spearman")
                            for f, r in feathers.items()} for c in ref.columns})
    a.out_dir.mkdir(parents=True, exist_ok=True)
    rho.to_csv(a.out_dir / "correlogram_feathers_vs_rnaseq.tsv", sep="\t", float_format="%.4f")
    print(f"{len(genes):,} genes\n{rho.round(3).to_string()}")
    plot(rho, len(genes), a.out_dir / "figure_s5_correlogram")


if __name__ == "__main__":
    main()
