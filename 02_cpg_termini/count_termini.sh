#!/bin/bash
# Termini around CpGs for one library: collapse read pairs to fragments, pile up the
# CpG windows of each compartment (MAPQ >= 30, base quality >= 25, no BAQ) and write one
# <out_prefix>.<compartment>.fragkon per BED.
# usage: count_termini.sh <library.rmdup.bam> <genome.fa> <out_prefix> <compartment.bed>...
#   e.g. count_termini.sh SCS_DF1.rmdup.bam galGal6.fa fragkon/SCS_DF1 \
#            chicken.whole_genome.bed chicken.CR1.bed chicken.CpG_island.bed
set -euo pipefail
here=$(dirname "$0")
bam=$1; ref=$2; out=$3; shift 3

python "$here/make_fragment_bam.py" "$bam" "$out.frag.unsorted.bam"
samtools sort -o "$out.frag.bam" "$out.frag.unsorted.bam" && rm "$out.frag.unsorted.bam"
samtools index "$out.frag.bam"
for bed in "$@"; do
    compartment=$(basename "$bed" .bed); compartment=${compartment#*.}
    samtools mpileup -a -B -q 30 -Q 25 -l "$bed" -f "$ref" "$out.frag.bam" |
        python "$here/parse_termini.py" > "$out.$compartment.fragkon"
done
