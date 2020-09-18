#! /usr/bin/env python
# concatenates all DNA sequences from a single strain
# into one long genome with NNNNNNNNNNNNNNNNN between
# sequences
#
# uses taxonomy file to determine strain ID of each fasta
# sequence
#
# requires fasta in order by strain
#
# usage (NNNNNNNNN is what to separate strains with):
# *.py input.fna input.tax output.fna output.tax NNNNNNNNNN
import os
import sys

infna = sys.argv[1]
intax = sys.argv[2]
outfna = sys.argv[3]
outtax = sys.argv[4]
delim = sys.argv[5]

# load tax map
taxmap = {} # header:strain
taxa = set() # set of all strains
print('Reading/writing taxmap...')
count = 0
with open(intax,'r') as f:
    with open(outtax,'w') as g:
        for line in f:
            count += 1
            if count % 100000 == 0:
                sys.stdout.write(str(count) + ' ')
                sys.stdout.flush()
            words = line.strip().split('\t')
            header = words[0]
            tax = words[1]
            # write new taxmap line only if first instance of this strain
            if not tax in taxa:
                g.write(tax + '\t' + tax + '\n')
                taxa.add(tax)
            taxmap[header] = tax
sys.stdout.write('\n')

# read/write streaming fasta
# track strains found to see if a
# strain has reads that already showed up
# earlier in file (and fail in that case)
currentstrain = ''
currentdna = ''
count = 0
print('Reading/writing FASTA...')
with open(infna,'r') as f:
    with open(outfna, 'w') as g:
        for line in f:            
            if line.startswith('>'):
                count += 1
                if count % 100000 == 0:
                    sys.stdout.write(str(count) + ' ')
                    sys.stdout.flush()
                header = line.strip().split()[0][1:]
                strain = taxmap[header]
                if strain != currentstrain:
                    if currentstrain != '':
                        g.write('>' + currentstrain + '\n')
                        g.write(currentdna + '\n')
                        currentdna = ''
                    currentstrain = strain
            else:
                if currentdna != '':
                    currentdna += delim
                currentdna += line.strip()
        g.write('>' + currentstrain + '\n' + currentdna + '\n') # don't forget to write the last sequence
sys.stdout.write('\n')

