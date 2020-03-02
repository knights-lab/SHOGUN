"""
Copyright 2015-2020 Knights Lab, Regents of the University of Minnesota.

This software is released under the GNU Affero General Public License (AGPL) v3.0 License.
"""

import csv
import os
import re
from collections import defaultdict, Counter

import pandas as pd

from shogun import logger
from shogun.wrappers import utree_search_gg
from ._aligner import Aligner


class UtreeAligner(Aligner):
    _name = 'utree'

    def __init__(self, database_dir, **kwargs):
        super().__init__(database_dir, **kwargs)

        # Setup the utree database
        prefix = self.data_files[self._name]
        self.compressed_tree = os.path.join(self.database_dir, prefix + ".ctr")

    def align(self, infile, outdir):
        if not os.path.exists(outdir):
            os.makedirs(outdir)

        outfile = os.path.join(outdir, 'alignment.utree.tsv')

        #TODO: pie chart and coverage
        proc, out, err = utree_search_gg(self.compressed_tree, infile, outfile, threads=self.threads, shell=self.shell)

        if self.post_align:
            df = self._post_align(outfile)
            self.outfile = os.path.join(outdir, 'taxatable.utree.txt')
            df.to_csv(self.outfile, sep='\t', float_format="%d",na_rep=0, index_label="#OTU ID")
        return proc, out, err

    def _post_align(self, utree_out: str) -> pd.DataFrame:
        logger.debug("Beginning post align with aligner %s" % self._name)
        samples_lca_map = defaultdict(Counter)
        with open(utree_out) as utree_f:
            csv_utree = csv.reader(utree_f, delimiter='\t')
            # qname, lca, confidence, support
            for line  in csv_utree:
                #TODO confidence/support filter
                taxonomy = split_utree_taxonomy(line[1])
                samples_lca_map['_'.join(line[0].split('_')[:-1])].update([taxonomy])
        df = pd.DataFrame(samples_lca_map, dtype=int)
        return df

def split_utree_taxonomy(tax):
    output = []
    for _ in itersplit(tax, sep=';'):
        if len(_) > 3:
            output.append(_)
        else:
            break
    return ";".join(output)


def itersplit(s, sep=None):
    # https://stackoverflow.com/questions/3862010/is-there-a-generator-version-of-string-split-in-python
    exp = re.compile(r'\s+' if sep is None else re.escape(sep))
    pos = 0
    while True:
        m = exp.search(s, pos)
        if not m:
            if pos < len(s) or sep is not None:
                yield s[pos:]
            break
        if pos < m.start() or sep is not None:
            yield s[pos:m.start()]
        pos = m.end()

