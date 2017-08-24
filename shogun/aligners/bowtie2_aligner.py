"""
Copyright 2015-2017 Knights Lab, Regents of the University of Minnesota.

This software is released under the GNU Affero General Public License (AGPL) v3.0 License.
"""

import os
from collections import defaultdict, Counter

import pandas as pd
from cytoolz import valfilter

from shogun import logger
from shogun.parsers import yield_alignments_from_sam_inf
from shogun.redistribute import Taxonomy
from shogun.utils import build_lca_map
from shogun.wrappers import bowtie2_align
from ._aligner import Aligner


class BowtieAligner(Aligner):
    _name = 'bowtie2'

    def __init__(self, database_dir, **kwargs):
        super().__init__(database_dir, **kwargs)

        # Setup the bowtie2 db
        self.prefix = os.path.join(database_dir, self.data_files[self._name])
        self.tree = Taxonomy(self.tax)

    def align(self, infile, outdir, alignments_to_report=16):
        outfile = os.path.join(outdir, 'bowtie2_results.sam')

        #TODO: pie chart and coverage
        proc, out, err = bowtie2_align(infile, outfile, self.prefix,
                             num_threads=self.threads, alignments_to_report=alignments_to_report, shell=self.shell)
        df = self._post_align(outfile)
        self.outfile = os.path.join(outdir, 'bowtie2_taxon_counts.txt')
        df.to_csv(self.outfile, sep='\t', float_format="%d",na_rep=0, index_label="#OTU ID")
        return proc, out, err

    def _post_align(self, sam_file: str) -> pd.DataFrame:
        logger.debug("Beginning post align with aligner %s" % self._name)
        align_gen = yield_alignments_from_sam_inf(sam_file)
        lca_map = build_lca_map(align_gen, self.tree)
        samples_lca_map = defaultdict(Counter)
        for key, value in valfilter(lambda x: x is not None, lca_map).items():
            samples_lca_map['_'.join(key.split('_')[:-1])].update([value])

        df = pd.DataFrame(samples_lca_map, dtype=int)
        return df
