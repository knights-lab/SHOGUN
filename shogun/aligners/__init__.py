"""
Copyright 2015-2017 Knights Lab, Regents of the University of Minnesota.

This software is released under the GNU Affero General Public License (AGPL) v3.0 License.
"""

import os
from yaml import load

from shogun.wrappers import embalmer_align, embalmulate, utree_search, bowtie2_align

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
        embalmer_align(infile, outfile,
            self.database, tax=self.tax, accelerator=self.accelerator, shell=self.shell, taxa_ncbi=False, threads=self.threads)
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

    def align(self, infile, outdir, alignments_to_report=16):
        outfile = os.path.join(outdir, 'results.sam')

        #TODO: pie chart and coverage
        return bowtie2_align(infile, outfile, self.prefix,
                             num_threads=self.threads, alignments_to_report=alignments_to_report, shell=self.shell)

__all__ = ["BowtieAligner", "UtreeAligner", "EmbalmerAligner"]
