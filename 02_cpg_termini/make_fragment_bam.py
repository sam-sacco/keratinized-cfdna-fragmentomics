#!/usr/bin/env python3
"""Collapse each proper read pair into one read spanning the whole fragment.

The termini parser counts every read start (^) and end ($) in a pileup as a fragment
terminus. That is right for merged reads, but for an unmerged pair each mate's inner end
is just where sequencing stopped. Here each proper pair becomes a single alignment from
the left mate's start to the right mate's end: both mates' bases are kept, the
unsequenced middle is a reference skip (N), and overlapping bases of the right mate are
dropped. The pileup then marks a start and an end only at the two outer coordinates.

Merged (unpaired) reads pass through unchanged. Improper pairs, mates on different
contigs and pairs clipped at an outer end (where the fragment end is uncertain) are
dropped; clipping on a mate's inner side is trimmed off. samtools mpileup already
skipped improper pairs.
The fragment takes read 1's strand (merged reads are in read-1 orientation) and the
lower of the two mapping qualities.

usage: make_fragment_bam.py <in.rmdup.bam> <out.bam>   (output is unsorted)
"""
import sys

import pysam

M, I, D, N, S, H, EQ, X = 0, 1, 2, 3, 4, 5, 7, 8
REF_OPS, QRY_OPS = {M, D, N, EQ, X}, {M, I, S, EQ, X}
SKIP = 0x4 | 0x100 | 0x200 | 0x400 | 0x800       # unmapped, secondary, qcfail, dup, supp


def trim_left(read, n_ref):
    """(cigar, seq, qual) of `read` with its first n_ref reference bases removed."""
    cigar, seq, qual = list(read.cigartuples), read.query_sequence, read.query_qualities
    q = 0
    while n_ref > 0 and cigar:
        op, ln = cigar[0]
        if op in REF_OPS:
            take = min(ln, n_ref)
            n_ref -= take
            if op in QRY_OPS:
                q += take
            cigar[0] = (op, ln - take)
            if cigar[0][1] == 0:
                cigar.pop(0)
        else:                                    # insertion before the cut: drop it
            q += ln if op in QRY_OPS else 0
            cigar.pop(0)
    while cigar and cigar[0][0] in (I, D, N):    # an alignment cannot open on these
        op, ln = cigar.pop(0)
        q += ln if op in QRY_OPS else 0
    return cigar, seq[q:], qual[q:]


def unclipped(read, side):
    """(cigar, seq, qual) with soft clips on `side` ("left"/"right") removed."""
    cigar, seq, qual = list(read.cigartuples), read.query_sequence, read.query_qualities
    if side == "left" and cigar[0][0] == S:
        n = cigar.pop(0)[1]
        seq, qual = seq[n:], qual[n:]
    if side == "right" and cigar[-1][0] == S:
        n = cigar.pop()[1]
        seq, qual = seq[:-n], qual[:-n]
    return cigar, seq, qual


def merge_pair(a, b, header):
    left, right = (a, b) if a.reference_start <= b.reference_start else (b, a)
    r1 = a if a.is_read1 else b
    if right.reference_end <= left.reference_end:   # right mate inside the left one
        return None
    lc, ls, lq = unclipped(left, "right")            # inner-side clips only
    if any(op in (S, H) for op, _ in lc):
        return None                                  # left mate's outer end is clipped
    overlap = left.reference_end - right.reference_start
    if overlap > 0:
        rc, rs, rq = trim_left(right, overlap)
        gap = []
    else:
        rc, rs, rq = unclipped(right, "left")
        gap = [(N, -overlap)] if overlap < 0 else []
    if not rc or any(op in (S, H) for op, _ in rc):
        return None                                  # right mate's outer end is clipped
    frag = pysam.AlignedSegment(header)
    frag.query_name = left.query_name
    frag.flag = 0x10 if r1.is_reverse else 0
    frag.reference_id = left.reference_id
    frag.reference_start = left.reference_start
    frag.mapping_quality = min(a.mapping_quality, b.mapping_quality)
    frag.cigartuples = lc + gap + rc
    frag.query_sequence = ls + rs
    frag.query_qualities = pysam.qualitystring_to_array(
        pysam.qualities_to_qualitystring(lq) + pysam.qualities_to_qualitystring(rq))
    return frag


def main(inp, out):
    stats = dict(merged=0, pairs=0, dropped_pairs=0, improper_reads=0)
    pending = {}
    with pysam.AlignmentFile(inp) as bam, \
            pysam.AlignmentFile(out, "wb", template=bam) as fo:
        for r in bam:
            if r.flag & SKIP:
                continue
            if not r.is_paired:
                fo.write(r)
                stats["merged"] += 1
                continue
            if (not r.is_proper_pair or r.mate_is_unmapped
                    or r.next_reference_id != r.reference_id):
                stats["improper_reads"] += 1
                continue
            mate = pending.pop(r.query_name, None)
            if mate is None:
                pending[r.query_name] = r
                continue
            frag = merge_pair(mate, r, bam.header)
            if frag is None:
                stats["dropped_pairs"] += 1
            else:
                fo.write(frag)
                stats["pairs"] += 1
    stats["unmatched_mates"] = len(pending)
    print("\t".join(f"{k}={v}" for k, v in stats.items()), file=sys.stderr)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
