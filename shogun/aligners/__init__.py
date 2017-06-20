"""
Copyright 2015-2017 Knights Lab, Regents of the University of Minnesota.

This software is released under the GNU Affero General Public License (AGPL) v3.0 License.
"""

import os
from yaml import load

from shogun.wrappers import embalmer_align, embalmulate

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
            'embalmer': ['.edb']
        }
        for value in SUFFICES[cls._name]:
            if not os.path.exists(os.path.join(dir, data_files[cls._name] + value)):
                return False, '%s not found' % (os.path.join(dir, value))
        return True, ''

    def align(self):
        raise NotImplementedError

class EmbalmerAligner(Aligner):
    _name = 'embalmer'

    def __init__(self, database_dir, threads=1, shell=False):
        super().__init__(database_dir, threads=threads, shell=shell)

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

        return embalmer_align(infile, outfile,
            self.database, tax=self.tax, accelerator=self.accelerator, shell=self.shell, taxa_ncbi=False, threads=self.threads)

        #TODO: Embalmulate

        #TODO: pie chart and coverage



class UtreeAligner(Aligner):
    def __init__(self, database_dir):
        super().__init__(database_dir)


class BowtieAligner(Aligner):
    def __init__(self, database_dir):
        super().__init__(database_dir)

__all__ = ["BowtieAligner", "UtreeAligner", "EmbalmerAligner"]
