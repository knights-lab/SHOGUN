#!/usr/bin/env python
from __future__ import print_function
import argparse
import sys
import os


def make_arg_parser():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('input', help='If nothing is given, then stdin, else the input file')
    parser.add_argument('-o', '--output', help='If nothing is given, then stdout, else write to file')
    return parser


def read_fasta(f):
    title = None
    data = None
    for line in f:
        if line[0] == ">":
            if title:
                yield (title.strip(), data)
            title = line[1:]
            data = ''
        else:
            data += line.strip()
    if not title:
        yield None
    yield (title.strip(), data)


def main():
    parser = make_arg_parser()
    args = parser.parse_args()
    with open(args.input) if args.input else sys.stdin as inf:
        fasta_gen = read_fasta(inf)
        with open(args.output, 'w') if args.output else sys.stdout as outf:
            for title, data in fasta_gen:
                outf.write('>' + title + os.linesep)
                outf.write(data + os.linesep)

if __name__ == '__main__':
    main()
