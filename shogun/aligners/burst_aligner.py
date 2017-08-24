"""
Copyright 2015-2017 Knights Lab, Regents of the University of Minnesota.

This software is released under the GNU Affero General Public License (AGPL) v3.0 License.
"""

import csv
import os
from collections import defaultdict

import numpy as np
import pandas as pd

from shogun import logger
from shogun.aligners import Aligner
from shogun.redistribute import Taxonomy
from shogun.wrappers import embalmer_align


class EmbalmerAligner(Aligner):
    _name = 'embalmer'

    def __init__(self, database_dir, post_align='capitalist', **kwargs):
        super().__init__(database_dir,**kwargs)

        # Setup the embalmer database
        prefix = self.data_files[self._name]
        self.database = os.path.join(self.database_dir, prefix)

        if os.path.exists(self.database + '.acx'):
            self.accelerator = self.database + '.acx'
        else:
            self.accelerator = False
        self.tree = Taxonomy(self.tax)
        self.post_align = post_align


    def _post_align(self, outf):
        if self.post_align == 'capitalist':
            return self._post_align_capitalist(outf)
        else:
            return self._post_align_taxonomy(outf)


    def align(self, infile, outdir):
        if not os.path.exists(outdir):
            os.makedirs(outdir)

        self.outfile = os.path.join(outdir, 'embalmer_results.b6')

        #TODO: pie chart and coverage
        proc, out, err = embalmer_align(infile, self.outfile,
            self.database,tax=self.tax, accelerator=self.accelerator, shell=self.shell,
                                        taxa_ncbi=False, threads=self.threads)
        df = self._post_align(self.outfile)
        self.outfile = os.path.join(outdir, 'embalmer_taxatable.txt')
        df.to_csv(self.outfile, sep='\t', float_format="%d", na_rep=0, index_label="#OTU ID")
        return proc, out, err

    def _post_align_capitalist(self, outf):
        logger.debug("Beginning post align capitalist style with aligner %s" % self._name)
        # This alignment parsing assumes capitalist output
        samples_lca_map = defaultdict(lambda: defaultdict(int))
        with open(outf) as emb_inf:
            csv_embalm = csv.reader(emb_inf, delimiter='\t')
            # qname, lca, confidence, support
            for line  in csv_embalm:
                tax = self.tree(line[1])
                #TODO confidence/support filter
                samples_lca_map['_'.join(line[0].split('_')[:-1])][tax] += 1

        df = pd.DataFrame(samples_lca_map, dtype=np.int).fillna(0).astype(np.int)
        return df


    def _post_align_taxonomy(self, outf):
        logger.debug("Beginning post align taxonomy style with aligner %s" % self._name)
        samples_lca_map = defaultdict(lambda: defaultdict(int))
        with open(outf) as utree_f:
            csv_embalm = csv.reader(utree_f, delimiter='\t')
            # qname, lca, confidence, support
            for line  in csv_embalm:
                if line[-1] is not None:
                    #TODO confidence/support filter
                    samples_lca_map['_'.join(line[0].split('_')[:-1])][line[-1]] += 1

        df = pd.DataFrame(samples_lca_map, dtype=np.int).fillna(0).astype(np.int)
        return df
