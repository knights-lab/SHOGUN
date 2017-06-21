"""
Copyright 2015-2017 Knights Lab, Regents of the University of Minnesota.

This software is released under the GNU Affero General Public License (AGPL) v3.0 License.
"""

import os
from yaml import load
from collections import defaultdict, Counter
import pandas as pd

from shogun.wrappers import embalmer_align, embalmulate, utree_search, bowtie2_align
from shogun.utils.last_common_ancestor import build_lca_map
from shogun.parsers import yield_alignments_from_sam_inf
from shogun.taxonomy import Taxonomy


class Aligner:
    _name = None

    def __init__(self, database_dir, threads=1, shell=False):
        self.threads = threads
        self.shell = shell
        check, msg = self.check_database(database_dir)

        with open(os.path.join(database_dir, 'metadata.yaml')) as stream:
            self.data_files = load(stream)

        if not check:
            raise Exception("Database %s is not formatted correctly: %s" % (database_dir, msg))

        self.database_dir = database_dir

        self.tax = os.path.join(database_dir, self.data_files['general']['taxonomy'])
        self.fasta = os.path.join(database_dir, self.data_files['general']['fasta'])

    @classmethod
    def check_database(cls, dir):
        with open(os.path.join(dir, 'metadata.yaml'), 'r') as stream:
            data_files = load(stream)

        SUFFICES = {
            'embalmer': ['.edx'],
            'utree': ['.ctr'],
            'bowtie2': ['.1.bt2']
        }
        for value in SUFFICES[cls._name]:
            if not os.path.exists(os.path.join(dir, data_files[cls._name] + value)):
                return False, '%s not found' % (os.path.join(data_files[cls._name] + value))
        return True, ''

    def align(self):
        raise NotImplementedError

    def _post_align(self):
        pass

class EmbalmerAligner(Aligner):
    _name = 'embalmer'

    def __init__(self, database_dir, **kwargs):
        super().__init__(database_dir,**kwargs)

        # Setup the embalmer database
        prefix = self.data_files[self._name]
        self.database = os.path.join(self.database_dir, prefix)

        if os.path.exists(self.database + '.acc'):
            self.accelerator = True
        else:
            self.accelerator = False

    def align(self, infile, outdir):
        if not os.path.exists(outdir):
            os.makedirs(outdir)

        outfile = os.path.join(outdir, 'results.b6')

        #TODO: pie chart and coverage
        proc, out, err = embalmer_align(infile, outfile,
            self.database,tax=self.tax, accelerator=self.accelerator, shell=self.shell,
                                        taxa_ncbi=False, threads=self.threads)
        proc2, out2, err2 = self._post_align(outfile, outdir)
        return (proc and proc2)

    def _post_align(self, outfile, outdir):
        return embalmulate(outfile, outdir)


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

        outfile = os.path.join(outdir, 'results.tsv')

        #TODO: pie chart and coverage
        return utree_search(self.compressed_tree, infile, outfile, shell=self.shell)


class BowtieAligner(Aligner):
    _name = 'bowtie2'

    def __init__(self, database_dir, **kwargs):
        super().__init__(database_dir, **kwargs)

        # Setup the bowtie2 db
        self.prefix = os.path.join(database_dir, self.data_files[self._name])
        self.tree = Taxonomy(self.tax)

    def align(self, infile, outdir, alignments_to_report=16):
        outfile = os.path.join(outdir, 'results.sam')


        #TODO: pie chart and coverage
        proc, out, err = bowtie2_align(infile, outfile, self.prefix,
                             num_threads=self.threads, alignments_to_report=alignments_to_report, shell=self.shell)
        df = self._post_align(outfile, outdir)
        print()

    def _post_align(self, sam_file: str, outdir: str) -> pd.DataFrame:
        align_gen = yield_alignments_from_sam_inf(sam_file)
        lca_map = build_lca_map(align_gen, self.tree)
        lca_map = filter(None, lca_map.values())
        samples_lca_map = defaultdict(Counter)
        for key, value in lca_map.items():
            samples_lca_map['_'.join(key.splt('_')[:-1])].update(value)

        return pd.DataFrame(samples_lca_map, index=samples_lca_map.keys())


__all__ = ["BowtieAligner", "UtreeAligner", "EmbalmerAligner"]
