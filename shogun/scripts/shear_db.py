#! /usr/bin/env python
# Shears fasta into substrings of length w step d
# usage
# *.py db.fasta 100 50
import sys

w = int(sys.argv[2])
d = int(sys.argv[3])
for line in open(sys.argv[1]):
    if line.startswith('>'):
        header = line.strip()
        header_words = header.split()
        header_id = header_words[0]
        header_comments = ''
        if len(header_words) > 1:
            header_comments = ' ' + ' '.join(header_words[1:])
    else:
        seq = line.strip()
        n = len(seq)
        start_pos = 0
        while start_pos < n:
            print(header_id + "_%09d" % (start_pos) + header_comments)
            if start_pos + w > n:
                start_pos = n - w
            print(seq[start_pos:(start_pos + w)])
            if start_pos == (n - w):
                start_pos = n
            else:
                start_pos += d
