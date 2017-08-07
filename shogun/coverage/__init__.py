"""
Copyright 2015-2017 Knights Lab, Regents of the University of Minnesota.

This software is released under the GNU Affero General Public License (AGPL) v3.0 License.
"""

from collections import defaultdict
import numpy as np
import pandas as pd
import csv

from shogun import logger
from shogun.redistribute import summarize_bayes_at_level

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
    taxa_hits = defaultdict(int)

    logger.info("Started the coverage parsing.")
    with open(infile) as utree_f:
        csv_embalm = csv.reader(utree_f, delimiter='\t')
        # qname, lca, confidence, support
        for num, line in enumerate(csv_embalm):
            if num % 1000 == 0:
                logger.info("Parsed %d lines of b6." % num)
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
                        taxa_hits[taxaname] += 1
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

    xx = np.zeros((len(samples_begin_map), 8))
    for i, taxaname in enumerate(sorted(samples_begin_map.keys())):
        if i % 1000 == 0:
            logger.info("Calculated %d coverages." % i)
        unique_hits = taxa_hits[taxaname]
        hits = samples_begin_map[taxaname]
        coverages = zero_runs(hits)
        if coverages[0][0] == 0:
            if coverages[-1][-1] == hits.shape[0]:
                temp = coverages[:, 1] - coverages[:, 0]
                coverages = np.concatenate((coverages, np.atleast_2d(np.array([0, temp[0] + temp[-1]]))))
        max_uncovered_region = np.max(coverages[:, 1] - coverages[:, 0])
        percent_max_unconvered = max_uncovered_region/shear_df['genome_length_median'][taxaname]
        percent_covered = np.sum(hits > 0)/shear_df['genome_length_median'][taxaname]
        unique_counts = shear_df.iloc[:, level-1][taxaname]
        expected_c = expected_coverage(unique_counts, unique_hits)
        row = np.array([max_uncovered_region, percent_max_unconvered, percent_covered, shear_df['genome_length_median'][taxaname], unique_hits, unique_counts, expected_c, percent_covered/(expected_c)])
        row[np.isnan(row)] = 0
        xx[i] = row
    df = pd.DataFrame(xx, columns=['max_uncovered_region', 'percent_max_uncovered_region', 'percent_of_genome_covered', 'median_genome_size', 'hits_in_clade', 'unique_counts_of_clade', 'expected_coverage', 'ratio_covered_over_expected'], index=sorted(samples_begin_map.keys()))
    logger.info("Completed the coverage analysis.")
    return df

def expected_coverage(length_of_genome, number_of_trials):
    num = length_of_genome*(1-np.power((length_of_genome-1)/length_of_genome, number_of_trials))
    return num/length_of_genome
