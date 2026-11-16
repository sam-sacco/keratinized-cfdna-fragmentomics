#!/bin/bash
# CTCF motif sites (JASPAR MA0139.1, p < 1e-6) for ctcf_profiles.py, with MEME suite 5.5.9.
# usage: run_fimo.sh <galGal6.fa> <GCF_001077635.1_Thamnophis_sirtalis-6.0_genomic.fna>
set -euo pipefail
fimo --thresh 1e-6 MA0139.1.meme "$1"                                  # -> fimo_out/fimo.tsv
fimo --text --thresh 1e-6 MA0139.1.meme "$2" > fimo_snake.tsv
