"""
Copyright 2015-2020 Knights Lab, Regents of the University of Minnesota.

This software is released under the GNU Affero General Public License (AGPL) v3.0 License.
"""

import os
from collections import defaultdict, Counter

import pandas as pd
from cytoolz import valfilter
import numpy as np

from shogun import logger
from shogun.parsers import yield_alignments_from_sam_inf
from ..utils.tree import Taxonomy
from shogun.utils.last_common_ancestor import build_lowest_common_ancestor_map
from shogun.wrappers import bowtie2_align
from shogun.utils.tree import build_tree_from_tax_file
from ._aligner import Aligner


class BowtieAligner(Aligner):
    _name = 'bowtie2'

    def __init__(self, database_dir, **kwargs):
        super().__init__(database_dir, **kwargs)

        # Setup the bowtie2 db
        self.prefix = os.path.join(database_dir, self.data_files[self._name])
        self.tree = build_tree_from_tax_file(self.tax)

    def align(self, infile, outdir, alignments_to_report=16):
        outfile = os.path.join(outdir, 'alignment.bowtie2.sam')

        #TODO: pie chart and coverage
        proc, out, err = bowtie2_align(infile, outfile, self.prefix,
                             num_threads=self.threads, alignments_to_report=alignments_to_report, shell=self.shell, percent_id=self.percent_id)
        if self.post_align:
            df = self._post_align(outfile)
            self.outfile = os.path.join(outdir, 'taxatable.bowtie2.txt')
            df.to_csv(self.outfile, sep='\t', float_format="%d", na_rep=0, index_label="#OTU ID")
        return proc, out, err

    def _post_align(self, sam_file: str, samples_iter=50) -> pd.DataFrame:
        logger.debug("Beginning post align with aligner %s" % self._name)
        align_gen = yield_alignments_from_sam_inf(sam_file)
        lca_map_gen = build_lowest_common_ancestor_map(align_gen, self.tree)
        sample_names_to_ix = dict()
        ix = 0
        mat_counts = np.zeros((self.tree.num_nodes, samples_iter), dtype=int)
        max_samples = samples_iter
        for rname, node_id in lca_map_gen:
            sample_name = rname.split('_')[0]
            if sample_name in sample_names_to_ix:
                c_ix = sample_names_to_ix[sample_name]
                mat_counts[node_id, c_ix] += 1
            else:
                if ix >= max_samples:
                    mat_counts = np.vstack((mat_counts, np.zeros(self.tree.num_nodes, dtype=int)))
                    b = np.zeros((max_samples, max_samples + samples_iter))
                    b[:, :-samples_iter] = mat_counts
                    mat_counts = b
                    max_samples += samples_iter
                sample_names_to_ix[sample_name] = ix
                mat_counts[node_id, ix] += 1
                ix += 1

        sample_names = [k for k, v in sorted(sample_names_to_ix.items(), key=lambda item: item[1])]

        df = pd.DataFrame(mat_counts[:, :ix], dtype=int, columns=sample_names)
        # drop all node ids of all zeros
        df = df.loc[~(df == 0).all(axis=1)].copy()
        df.index = [self.tree.node_id_to_taxa_name[node_id] for node_id in df.index]
        return df
