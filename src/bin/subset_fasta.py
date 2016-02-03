from __future__ import print_function
import argparse
import random
import sys
import os


def make_arg_parser():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('input', help='If nothing is given, then stdin, else the input file')
    parser.add_argument('-o', '--output', help='If nothing is given, then stdout, else write to file')
    parser.add_argument('-n', '--number', help='The number of sequences', required=True, type=int)
    parser.add_argument('-k', '--keep', help='The number of sequences to keep', required=True, type=int)
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


def filter_fasta(fasta_gen, total, keep):
    lines_to_keep = set(random.sample(range(1, total), keep))
    for i, (title, data) in enumerate(fasta_gen):
        if i in lines_to_keep:
            yield title, data


def main():
    parser = make_arg_parser()
    args = parser.parse_args()
    with open(args.input) if args.input else sys.stdin as inf:
        fasta_gen = read_fasta(inf)
        filtered_fasta_gen = filter_fasta(fasta_gen, args.number, args.keep)
        with open(args.output, 'w') if args.output else sys.stdout as outf:
            for title, data in filtered_fasta_gen:
                outf.write('>' + title + os.linesep)
                outf.write(data + os.linesep)

if __name__ == '__main__':
    main()