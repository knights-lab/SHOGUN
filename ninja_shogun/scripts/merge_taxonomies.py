#!/usr/bin/env python
# Parse taxonomy outputs and make one taxonomy table for qiime
#
# Note: delete headers add filename to each line and merge files:
# for i in *.txt; do sed '1d' $i > $i.bk; done
# for i in *.bk; do nawk '{print FILENAME","$0}' $i > $i.kk; mv $i.kk $i; done
# cat *.bk > taxonomies.txt
# rm *.bk
#
# usage:
# python merge_taxonomies.py -i taxonomies.txt -o taxonomy_table.txt

import sys
import os
import optparse


def get_opts():
    p = optparse.OptionParser()
    p.add_option("-i", "--input", type="string", \
                 default=None, help="Input usearch output file, where IMG was used as query [required].")
    p.add_option("-o", "--output", type="string", \
                 default=None, help="Output file [required].")
    p.add_option("--verbose", action="store_true", \
                 help="Print all output.")
    opts, args = p.parse_args(sys.argv)

    return opts, args


def check_opts(opts):
    if opts.input is None:
        raise ValueError('\n\nPlease include an input file.')
    if opts.output is None:
        raise ValueError('\n\nPlease include an output file.')


def create_taxonomy_dict(input_file):
    taxonomy_dict = {}
    # For each line in input file strip the lines
    # and split at comma
    # sampleID are the first column
    # Genus is third and species 4
    for line in input_file:
        words = line.strip().split(",")
        sampleID_trimmed = words[0].strip('"').split("_")[0]
        sampleID = sampleID_trimmed.split("-")[1]
        Genus = words[2].strip('"').replace(" ","")
        Species = words[3].strip('"').replace(" ","").replace("_","-")
        Taxonomy = Genus + "_" + Species
        Count = words[4].strip('"')

        # For each line in input_file check to see if
        # the bacteria is in the taxonomy_dict, if not, add it as a dict
        # then for that given taxonomy add the
        # sampleID on that line to the given to the dict
        if not taxonomy_dict.has_key(Taxonomy):
            taxonomy_dict[Taxonomy] = set()
            taxonomy_dict[Taxonomy] = {}
        if not taxonomy_dict[Taxonomy].has_key(sampleID):
            taxonomy_dict[Taxonomy][sampleID] = Count
    return taxonomy_dict


def create_sampleID_set(taxonomy_dict):
    sampleID_set = set()
    for taxon in taxonomy_dict:
        sampleID_set.update(taxonomy_dict[taxon])
    return sampleID_set


def write_header(sampleID_set, output_file):
    # In the output file write the each sampleID
    output_file.write('#OTU_ID')
    for IDs in sampleID_set:
        output_file.write('\t' + IDs)
    output_file.write('\n')


def write_rows(taxonomy_dict, sampleID_set, output_file):
    # For each taxon
    # print the taxonomy + tab
    # For sampleID in sampleID_set within the dict, print:
    # count of that taxon + tab
    # or 0 + tab
    for taxon in taxonomy_dict:
        output_file.write(taxon)
        for ID in sampleID_set:
            if ID in taxonomy_dict[taxon]:
                output_file.write('\t' + str(taxonomy_dict[taxon][ID]))
            else:
                output_file.write('\t0')
        output_file.write('\n')


def main(opts):
    # open input file
    input_file = open(opts.input, "r")

    # create an output file and open it
    output_file = open(opts.output, "w")

    # a dict mapping img IDs to a list of resfams_IDs
    taxonomy_dict = create_taxonomy_dict(input_file)

    # a set of all sample IDs (no duplicates)
    sampleID_set = create_sampleID_set(taxonomy_dict)

    # make sample_ID set and write header
    write_header(sampleID_set, output_file)

    # write the rows
    write_rows(taxonomy_dict, sampleID_set, output_file)

    # close the output file
    output_file.close

if __name__ == '__main__':
    opts, args = get_opts()
    check_opts(opts)
    main(opts)
