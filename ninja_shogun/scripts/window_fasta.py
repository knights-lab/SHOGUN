#!/usr/bin/env python
from __future__ import print_function
import argparse
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def make_arg_parser():
    parser = argparse.ArgumentParser(description='Fragment a FASTA file into a database.')
    parser.add_argument('input', help='If nothing is given then stdin, else basename of the file for the'
                                              ' reference genome. Must be in FASTA format.')
    parser.add_argument('-o', '--output', help='If nothing is given, then stdout, else write to file')
    parser.add_argument('-w', '--window', help='The length of the sliding window.', default=100)
    parser.add_argument('-s', '--step', help='The step size for sliding the window.', default=1)
    return parser


def chunks(l, step=1, window=100):
    for i in range(0, len(l)-window, step):
        yield l[i:i+window], i, i+window


def main():
    parser = make_arg_parser()
    args = parser.parse_args()

    with open(args.input, 'r') if args.input else sys.stdin as inf:
        records = SeqIO.parse(inf, 'fasta')
        with open(args.output, 'w') if args.output else sys.stdout as outf:
            for record in records:
                for chunk, start, end in chunks(str(record.seq), step=args.step, window=args.window):
                    outf.write(SeqRecord(Seq(
                        chunk, record.seq.alphabet),
                        id=record.id + "_%d_%d" % (start, end), name=record.name,
                        description=record.description).format('fasta'))

if __name__ == '__main__':
    main()
