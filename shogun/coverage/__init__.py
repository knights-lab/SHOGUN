"""
Copyright 2015-2017 Knights Lab, Regents of the University of Minnesota.

This software is released under the GNU Affero General Public License (AGPL) v3.0 License.
"""

from collections import defaultdict
import numpy as np
import pandas as pd
import csv

from shogun import logger
from shogun.redistribute import summarize_bayes_at_level, parse_bayes

def zero_runs(a):
    # Stack Overflow:
    # https://stackoverflow.com/questions/24885092/finding-the-consecutive-zeros-in-a-numpy-array
    # Create an array that is 1 where a is 0, and pad each end with an extra 0.
    iszero = np.concatenate(([0], np.equal(a, 0).view(np.int8), [0]))
    absdiff = np.abs(np.diff(iszero))
    # Runs start and end where absdiff is 1.
    ranges = np.where(absdiff == 1)[0].reshape(-1, 2)
    return ranges

def get_coverage_of_microbes(infile, shear, level):
    #Load in the shear df at level
    shear_df = summarize_bayes_at_level(shear, level=level)

    samples_begin_map = dict()

    with open(infile) as utree_f:
        csv_embalm = csv.reader(utree_f, delimiter='\t')
        # qname, lca, confidence, support
        for line in csv_embalm:
            if line[-1] is not None:
                # TODO confidence/support filter
                begin = int(line[8])
                # samplename = '_'.join(line[0].split('_')[:-1])
                # if not samplename in samples_begin_map:
                #     samples_begin_map[samplename] = dict()
                taxaname = line[-1]
                taxa_level = taxaname.count(';') + 1
                if taxa_level >= level:
                    if taxa_level != level:
                        taxaname = ';'.join(taxaname.split(";")[:level])
                    if taxaname in shear_df.index:
                        indx = int(np.floor(begin/100.))
                        if not taxaname in samples_begin_map:
                            genome_length = shear_df['genome_length_median'][taxaname]
                            samples_begin_map[taxaname] = np.zeros(genome_length)
                        if indx == 0:
                            samples_begin_map[taxaname][0] += 1
                        elif indx >= shear_df['genome_length_median'][taxaname]:
                            samples_begin_map[taxaname][-1] += 1
                        else:
                            samples_begin_map[taxaname][indx] += 1
                            samples_begin_map[taxaname][indx+1] += 1
                    else:
                        logger.warning("The taxa %s not found." % taxaname)

    xx = np.zeros((len(samples_begin_map), 4))
    for i, taxaname in enumerate(sorted(samples_begin_map.keys())):
        hits = samples_begin_map[taxaname]
        coverages = zero_runs(hits)
        max_uncovered_region = np.max(coverages[:, 1] - coverages[:, 0])
        percent_max_unconvered = max_uncovered_region/shear_df['genome_length_median'][taxaname]
        percent_uncovered = np.sum(hits == 0)/shear_df['genome_length_median'][taxaname]
        xx[i] = np.array([max_uncovered_region, percent_max_unconvered, percent_uncovered, shear_df['genome_length_median'][taxaname]])

    df = pd.DataFrame(xx, columns=['max_unconvered_region', 'percent_max_uncovered', 'percent_uncovered', 'median_genome_size'], index=sorted(samples_begin_map.keys()))
    return df
