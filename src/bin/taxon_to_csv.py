#!/usr/bin/env python
from __future__ import print_function
import argparse
import sys
from sheq.parse_taxonomy import parse_taxonomy

def make_arg_parser():
    parser = argparse.ArgumentParser(description='Grab the Genus/Species from a taxon.tab.txt')
    parser.add_argument('input', type=str, help='The file to parse')
    parser.add_argument('-o', '--output', help='Write to the stdout if nothing is given, else write to file')
    return parser

def main():
    parser = make_arg_parser()
    args = parser.parse_args()
    df = parse_taxonomy(args.input)
    with open(args.output, 'w') if args.output else sys.stdout as outf:
        df.to_csv(outf, index=False)

if __name__ == '__main__':
    main()
