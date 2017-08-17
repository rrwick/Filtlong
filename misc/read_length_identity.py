#!/usr/bin/env python3

"""
Copyright 2017 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/Filtlong

This script was used to get read length/identity stats from a minimap2 PAF alignment file. I used
it to make the scatterplot/histogram figures in the Filtlong README. It uses a particularly strict
definition of read identity: all bases are considered and unaligned bases are assigned an identity
of 0. So if a read had half of its bases align with an identity of 90% and the other half is
unaligned, then the read's final identity would be 45%.

Example commands:
  filtlong {OPTIONS} input.fastq.gz | gzip > filtered.fastq.gz
  minimap2 -x map10k -t 16 -c reference.fasta filtered.fastq.gz > alignments.paf
  read_length_identity.py alignments.paf > read_stats.tsv

This file is part of Filtlong. Filtlong is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Filtlong is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Filtlong. If
not, see <http://www.gnu.org/licenses/>.
"""

import sys
import collections
import statistics


def main():
    read_lengths = {}
    read_alignments = collections.defaultdict(list)

    paf_filename = sys.argv[1]
    with open(paf_filename, 'rt') as paf:
        for line in paf:
            paf_parts = line.strip().split('\t')
            if len(paf_parts) < 11:
                continue
            read_name = paf_parts[0]
            read_length = int(paf_parts[1])
            read_lengths[read_name] = read_length
            start = int(paf_parts[2])
            end = int(paf_parts[3])
            identity = 100.0 * int(paf_parts[9]) / int(paf_parts[10])
            read_alignments[read_name].append((start, end, identity))

    print('\t'.join(['Name', 'Length', 'Identity']))
    for read_name, read_length in read_lengths.items():
        alignments = read_alignments[read_name]
        identity_by_base = [0.0] * read_length
        for start, end, identity in alignments:
            for i in range(start, end):
                if identity > identity_by_base[i]:
                    identity_by_base[i] = identity
        whole_read_identity = statistics.mean(identity_by_base)
        print('\t'.join([read_name, str(read_length), str(whole_read_identity)]))


if __name__ == '__main__':
    main()
