"""
Copyright 2015-2020 Knights Lab, Regents of the University of Minnesota.

This software is released under the GNU Affero General Public License (AGPL) v3.0 License.
"""

import os

import pandas as pd

from shogun import logger
from shogun.utils.lowest_common_ancestor import build_lca_df
from shogun.wrappers import bowtie2_align
from shogun.utils.tree import build_tree_from_tax_file
from shogun.aligners._aligner import Aligner


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

    def _post_align(self, sam_file: str, samples_iter: int = 50, confidence_threshold: float = 1.0, **kwargs) -> pd.DataFrame:
        logger.debug("Beginning post align with aligner %s" % self._name)
        df = build_lca_df(sam_file, self.tree, confidence_threshold=confidence_threshold, samples_iter=samples_iter)
        return df
