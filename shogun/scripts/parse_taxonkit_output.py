#! /usr/bin/env python

# Converts taxonkit output to greengenes taxonomy format
# usage :
# parse_taxonkit_output.py taxonkit_output.txt output.txt
#
# example input line:
# 372461    cellular organisms;Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacterales;Erwiniaceae;Buchnera;Buchnera aphidicola;Buchnera aphidicola (Cinara cedri);Buchnera aphidicola BCc 131567;2;1224;1236;91347;1903409;32199;9;261318;372461
import sys
import os

with open(sys.argv[1],'r') as f:
    with open(sys.argv[2],'w') as g:
        for line in f:
            cols = line.strip().split('\t')
            taxid = cols[0]
            taxraw = cols[1]
            taxraw = taxraw.split(';')
            tax = 'k__' + taxraw[1]
            tax += ';p__' + taxraw[2]
            tax += ';c__' + taxraw[3]
            tax += ';o__' + taxraw[4]
            tax += ';f__' + taxraw[5]
            tax += ';g__' + taxraw[6]
            tax += ';s__' + taxraw[7]
            if len(taxraw) == 10:
                if len(taxraw[9]) > 0:
                    tax += ';t__' + taxraw[9]
                else:
                    tax += ';t__' + taxraw[7]
            else:
                if len(taxraw[8]) > 0:
                    tax += ';t__' + taxraw[8]
                else:
                    tax += ';t__' + taxraw[7]
            g.write(taxid + '\t' + tax + '\n')
