"""
Copyright 2015-2017 Knights Lab, Regents of the University of Minnesota.

This software is released under the GNU Affero General Public License (AGPL) v3.0 License.
"""

from collections import defaultdict
import numpy as np
import pandas as pd
import csv

def get_coverage_of_microbes(infile, database):
    samples_begin_map = defaultdict(lambda: defaultdict(list))
    samples_end_map = defaultdict(lambda: defaultdict(list))
    with open(infile) as utree_f:
        csv_embalm = csv.reader(utree_f, delimiter='\t')
        # qname, lca, confidence, support
        for line in csv_embalm:
            if line[-1] is not None:
                # TODO confidence/support filter
                samples_begin_map['_'.join(line[0].split('_')[:-1])][line[-1]].append(int(line[8]))
                samples_end_map['_'.join(line[0].split('_')[:-1])][line[-1]].append(int(line[9]))

    df = pd.DataFrame(samples_end_map, dtype=np.int).fillna(0).astype(np.int)
    return df
