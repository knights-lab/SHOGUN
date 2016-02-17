#!/usr/bin/env python
from __future__ import print_function
import argparse
import sys

# The arg parser for this wrapper
def make_arg_parser():
    parser = argparse.ArgumentParser(description='Simulate a metagenomic community from DWGSIM according to a powerlaw.')
    parser.add_argument('-i', '--input', type=str, help='The species name file. One species name per line.', required=True)
    parser.add_argument('-o', '--output', type=str, help='The output file.')
    return parser

def main(args):
    with sys.stdin if args.input == "-" else open(args.input, 'rb') as inf:
        with open(args.output, 'wb') if args.output else sys.stdout as outf:
            for line in inf:
                outf.write(line)

if __name__ == '__main__':
    parser = make_arg_parser()
    args = parser.parse_args()
    main(args)
