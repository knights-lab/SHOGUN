from __future__ import print_function
import argparse
import sys
import pandas as pd

def make_arg_parser():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('input', help='The input file')
    parser.add_argument('-o', '--output', help='If nothing is given, then stdout, else write to file')
    return parser

def index_stats(df, outf):
    outf.write('Domain hits:\t%d\n' % df.groupby('domain').sum().sum()[0])
    outf.write('Genus hits:\t%d\n' % df.groupby('genus').sum().sum()[0])
    outf.write('Species hits:\t%d\n' % df.groupby('species').sum().sum()[0])

def main():
    parser = make_arg_parser()
    args = parser.parse_args()
    df = pd.read_csv(args.input)
    with open(args.output, 'w') if args.output else sys.stdout as outf:
        index_stats(df, outf)

if __name__ == '__main__':
    main()
