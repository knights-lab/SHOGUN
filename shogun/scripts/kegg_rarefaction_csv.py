#!/usr/bin/env python
"""
Copyright 2015-2017 Knights Lab, Regents of the University of Minnesota.

This software is released under the GNU Affero General Public License (AGPL) v3.0 License.
"""

from __future__ import print_function
import argparse
import numpy as np
import csv
from multiprocessing import Pool, cpu_count


def make_arg_parser():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-i', '--input', help='The input file.', required=True)
    parser.add_argument('-o', '--output', help='If nothing is given, then stdout, else write to file')
    return parser


def count_lines(filename):
    lines = 0
    buffer = bytearray(2048)
    with open(filename, 'rb') as f:
        while f.readinto(buffer) > 0:
            lines += buffer.count(b'\n')
    return lines


def write_shuffles(args):
    size, num_lines, inf_path, outf_path = args
    sorted = np.random.shuffle(np.arange(num_lines))[:size].sort()
    with open(inf_path, 'r') as inf:
        csv_inf = csv.reader(inf)
        with open('%d.%s' % (size, args.output), 'a'):
            csv_outf = csv.writer(outf_path)
            for i, line in enumerate(csv_inf):
                ind = sorted.searchsorted(i)
                if sorted[ind] == i:
                    line[0] = '.'.join((line[0], size))
                    csv_outf.writerow(line)


def main():
    parser = make_arg_parser()
    args = parser.parse_args()

    num_lines = 100000000

    sizes = [np.power(10, i+1) for i in range(7)]

    pool = Pool(processes=cpu_count())
    pool.map(write_shuffles, zip(sizes, [num_lines]*7, [args.input]*7, [args.output]*7))


if __name__ == '__main__':
    main()
