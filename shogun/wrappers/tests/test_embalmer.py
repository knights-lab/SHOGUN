"""
Copyright 2015-2017 Knights Lab, Regents of the University of Minnesota.

This software is released under the GNU Affero General Public License (AGPL) v3.0 License.
"""

import unittest
import shutil
import pkg_resources
import os
import tempfile

from shogun.utils import hash_file, read_checksums
from shogun.wrappers.embalmer import embalmer_align, embalmer_build


class TestEmbalmer(unittest.TestCase):
    def setUp(self):
        prefix = "shogun-test-temp-"

        self.checksums = read_checksums(pkg_resources.resource_filename(
           'shogun.wrappers.tests', os.path.join('data', 'embalmer', 'checksums.txt')))

        self.temp_dir = tempfile.TemporaryDirectory(prefix=prefix)

    def tearDown(self):
        self.temp_dir.cleanup()

    def test_embalmer_path(self):
        self.assertTrue(shutil.which("emb15") is not None)

    def test_embalmer_align(self):
        database = pkg_resources.resource_filename('shogun.wrappers.tests', os.path.join('data', 'embalmer', 'embalmer-test'))
        infile = pkg_resources.resource_filename('shogun.wrappers.tests', os.path.join('data', 'sims.fna'))
        outfile = os.path.join(self.temp_dir.name, 'embalmer-test-sims.b6')
        tax = pkg_resources.resource_filename('shogun.wrappers.tests', os.path.join('data', 'genomes.small.tax'))
        self.assertTrue(embalmer_align(infile, outfile, database, tax=tax)[0] == 0)

    def test_embalmer_build(self):
        fasta = pkg_resources.resource_filename('shogun.wrappers.tests', os.path.join('data', 'genomes.small.fna'))
        outfile = os.path.join(self.temp_dir.name, 'embalmer-test')
        embalmer_build(fasta, outfile, shell=False, cr=1050, s=500)

        for file in os.listdir(self.temp_dir.name):
           self.assertTrue(self.checksums[hash_file(os.path.join(self.temp_dir.name, file))] == file)


if __name__ == '__main__':
    unittest.main()
