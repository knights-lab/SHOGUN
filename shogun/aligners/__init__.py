"""
Copyright 2015-2017 Knights Lab, Regents of the University of Minnesota.

This software is released under the GNU Affero General Public License (AGPL) v3.0 License.
"""

import os
from yaml import load

from shogun.wrappers import embalmer_align

class Aligner:
    _name = None

    def __init__(self, database_dir, name):
        check, msg = self.check_database(database_dir, name)

        self.data_files = load(os.path.join(database_dir, 'metadata.yaml'))
        if not check:
            raise Exception("Database %s is not formatted correctly: %s" % (database_dir, msg))

        self.database_dir = database_dir

    @classmethod
    def check_database(cls, dir):
        stream = open(os.path.join(dir, 'metadata.yaml'))
        data_files = load(stream)
        for key, value in data_files[cls._name].items():
            if value:
                if not os.path.exists(os.path.join(dir, value)):
                    return False, '%s not found' % (os.path.join(dir, value))
        return True, ''

    def align(self):
        raise NotImplementedError

class EmbalmerAligner(Aligner):
    _name = 'embalmer'

    def __init__(self, database_dir):
        super().__init__(database_dir)

    def align(self, infile, outfile):
        embalmer_align(infile, outfile,
            self.database, tax=self.tax, accelerator=self.accelerator, shell=self.shell, taxa_ncbi=self.taxa_ncbi)




class UtreeAligner(Aligner):
    def __init__(self, database_dir):
        super().__init__(database_dir)


class BowtieAligner(Aligner):
    def __init__(self, database_dir):
        super().__init__(database_dir)

__all__ = ["BowtieAligner", "UtreeAligner", "EmbalmerAligner"]
