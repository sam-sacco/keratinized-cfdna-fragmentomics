#!/usr/bin/env python3
"""Count fragment termini around each CpG from a samtools mpileup (read from stdin) of
the 12-bp CpG windows, and write one row per CpG (a "fragkon" table) to stdout.

For each of the 12 window positions (1-5, C, G, 8-12) a row holds forward/reverse read
starts and ends, the same counts divided by depth, forward/reverse 'A' counts in the
pileup string and the depth. Reference-matching read starts/ends only ('^?.', '.$',
'^?,', ',$'). A CpG is a C followed by a G on the next pileup line; windows that run
off the pileup or are not 12 consecutive positions are skipped. The header names four
substitution blocks although rows carry two; it is kept as the downstream code expects.

usage: samtools mpileup ... | parse_termini.py > sample.fragkon
"""
import csv
import re
import sys
from collections import deque

POS = [1, 2, 3, 4, 5, "C", "G", 8, 9, 10, 11, 12]
BLOCKS = ["forStart", "forEnd", "revStart", "revEnd",
          "forStartN", "forEndN", "revStartN", "revEndN",
          "C->T_for", "C->T_rev", "G->A_for", "G->A_rev", "Coverage"]
TERMINI = [re.compile(p) for p in (r"\^.\.", r"\.\$", r"\^.,", r",\$")]


def parse(line):
    """(chrom, pos, base, depth, [for start, for end, rev start, rev end], A, a)."""
    chrom, pos, base, depth, pileup = line.split("\t")[:5]
    depth = int(depth)
    return (chrom, int(pos), base, depth, [len(p.findall(pileup)) for p in TERMINI],
            pileup.count("A"), pileup.count("a"))


def ratio(n, depth):
    return n / depth if depth else float("nan")


def fragkon_row(win):
    """Output row for one 12-position window."""
    chrom, pos = win[5][0], win[5][1]
    row = [f"{chrom}_{pos}", "".join(r[2] for r in win).upper()]
    for k in range(4):
        row += [r[4][k] for r in win]
    for k in range(4):
        row += [ratio(r[4][k], r[3]) for r in win]
    row += [r[5] for r in win] + [r[6] for r in win] + [r[3] for r in win]
    return row


def main():
    out = csv.writer(sys.stdout, delimiter="\t")
    out.writerow(["CpGID", "SeqContext"] + [f"{b}{p}" for b in BLOCKS for p in POS])
    win = deque(maxlen=12)              # the last 12 pileup lines
    for line in sys.stdin:
        win.append(parse(line))
        # the C of a candidate window sits 6 lines back, its G 5 lines back
        if (len(win) == 12 and win[5][2].upper() == "C" and win[6][2].upper() == "G"
                and win[11][1] - win[0][1] == 11):
            out.writerow(fragkon_row(win))


if __name__ == "__main__":
    main()
