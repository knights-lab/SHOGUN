"""
Copyright 2015-2017 Knights Lab, Regents of the University of Minnesota.

This software is released under the GNU Affero General Public License (AGPL) v3.0 License.
"""

import os

from yaml import load

from shogun import logger


class Aligner:
    _name = None

    def __init__(self, database_dir, threads=1, post_align=True, shell=False, percent_id=.98):
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
        self.outfile = None
        logger.debug("Initiate Logger %s" % self._name)
        self.post_align = post_align
        self.percent_id = percent_id

    @classmethod
    def check_database(cls, dir):
        with open(os.path.join(dir, 'metadata.yaml'), 'r') as stream:
            data_files = load(stream)

        SUFFICES = {
            'burst': ['.edx'],
            'utree': ['.ctr'],
            'bowtie2': ['*'],
            'burst_best': ['.edx']
        }
        for value in SUFFICES[cls._name]:
            if cls._name == "bowtie2":
                from glob import glob
                files = glob(os.path.abspath(os.path.join(dir, data_files[cls._name] + value)))
                if len(files) == 0:
                    return False, "%s not found" % ("Prefix not found: %s"  % os.path.abspath(os.path.join(dir, data_files[cls._name])))
            elif not os.path.exists(os.path.abspath(os.path.join(dir, data_files[cls._name] + value))):
                return False, '%s not found' % (os.path.abspath(os.path.join(dir, data_files[cls._name] + value)))
        return True, ''

    def align(self, infile, outdir):
        raise NotImplementedError

    def _post_align(self, outf):
        raise NotImplementedError
