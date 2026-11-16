#!/bin/bash
# CpG windows for the three compartments of Fig 1C, for chicken (galGal6) and garter
# snake (thaSir1, NCBI Thamnophis_sirtalis-6.0). Needs bedtools and the inputs below.
#   galGal6.fa, thaSir1.fna                        genomes
#   galGal6.rmsk.txt.gz, galGal6.cpgIslandExt.txt.gz, thaSir1.cpgIslandExt.txt.gz
#                                                  UCSC tables (hgdownload .../database/)
#   thaSir.genome.fa.out                           RepeatMasker output of Gable et al. 2024
set -euo pipefail
here=$(dirname "$0")

# chicken: windows overlapping a CR1 / CpG island are kept whole
python "$here/make_cpg_windows.py" galGal6.fa > chicken.whole_genome.bed
zcat galGal6.rmsk.txt.gz | awk -v OFS='\t' '$11 ~ /CR1/ {print $6, $7, $8}' > chicken_cr1.bed
zcat galGal6.cpgIslandExt.txt.gz | cut -f2-4 > chicken_islands.bed
bedtools intersect -u -a chicken.whole_genome.bed -b chicken_cr1.bed > chicken.CR1.bed
bedtools intersect -u -a chicken.whole_genome.bed -b chicken_islands.bed > chicken.CpG_island.bed

# snake: windows are clipped to the CR1 / island boundaries (bedtools default). The
# island intersection was made with windows one base to the right and then shifted
# back, which is the same as intersecting with islands moved one base to the left.
python "$here/make_cpg_windows.py" thaSir1.fna > snake.whole_genome.bed
awk -v OFS='\t' '$11 == "LINE/CR1" {print $5, $6 - 1, $7}' thaSir.genome.fa.out > snake_cr1.bed
zcat thaSir1.cpgIslandExt.txt.gz |
    awk -v OFS='\t' '{sub(/v/, ".", $2); s = $3 - 1; print $2, (s < 0 ? 0 : s), $4 - 1}' \
    > snake_islands.bed
bedtools intersect -a snake.whole_genome.bed -b snake_cr1.bed > snake.CR1.bed
bedtools intersect -a snake.whole_genome.bed -b snake_islands.bed > snake.CpG_island.bed
