#!/usr/bin/env python
from __future__ import print_function
import argparse

# The arg parser for this wrapper
def make_arg_parser():
    parser = argparse.ArgumentParser(description='Simulate a metagenomic community from DWGSIM according to a powerlaw.')
    parser.add_argument('-i', '--input', type=str, help='The species name file. One species name per line.', required=True)
    return parser

def main(args):
    with sys.stdin if args.input == "-" else open(args.input, 'rb') as inf:
        for line in inf:
            print(line)


if __name__ == '__main__':
    parser = make_arg_parser()
    args = parser.parse_args()
    main(args)
