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
            tax += ';p__'
            if len(taxraw) > 2:
                tax += taxraw[2]
            tax += ';c__'
            if len(taxraw) > 3:
                tax += taxraw[3]
            tax += ';o__'
            if len(taxraw) > 4:
                tax += taxraw[4]
            tax += ';f__'
            if len(taxraw) > 5:
                tax += taxraw[5]
            tax += ';g__'
            if len(taxraw) > 6:
                tax += taxraw[6]
            tax += ';s__'
            if len(taxraw) > 7:
                tax += taxraw[7]
            tax += ';t__'
            if len(taxraw) > 9:
                if len(taxraw[9]) > 0:
                    tax += taxraw[9]
                else:
                    tax += taxraw[7]
            elif len(taxraw) == 9:
                if len(taxraw[8]) > 0:
                    tax += taxraw[8]
                else:
                    tax += taxraw[7]
            else:
                tax += 'None'
            g.write(taxid + '\t' + tax + '\n')
