#!/usr/bin/env python3
"""BED of 12-bp windows around every CpG in a genome: 5 bp before the C to 5 bp after
the G (0-based start = C - 5, end = C + 7). Windows that would start before the
beginning of a sequence are skipped (ones running past the end are kept; mpileup
clips them and the termini parser then drops the incomplete window).

usage: make_cpg_windows.py genome.fa > cpg_windows.bed
"""
import re
import sys

import pysam

with pysam.FastxFile(sys.argv[1]) as fa:
    for rec in fa:
        seq = rec.sequence.upper()
        for m in re.finditer("(?=CG)", seq):
            if m.start() >= 5:
                sys.stdout.write(f"{rec.name}\t{m.start() - 5}\t{m.start() + 7}\n")
