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
from shogun.wrappers.burst import burst_align, burst_build, embalmulate


class TestBurst(unittest.TestCase):
    def setUp(self):
        prefix = "shogun-test-temp-"

        self.checksums = read_checksums(pkg_resources.resource_filename(
           'shogun.tests', os.path.join('data', 'burst', 'checksums.txt')))

        self.temp_dir = tempfile.TemporaryDirectory(prefix=prefix)

    def tearDown(self):
        self.temp_dir.cleanup()

    def test_burst_path(self):
        self.assertTrue(shutil.which("emb15") is not None)

    def test_burst_align(self):
        database = pkg_resources.resource_filename('shogun.tests', os.path.join('data', 'burst', 'genomes.small'))
        infile = pkg_resources.resource_filename('shogun.tests', os.path.join('data', 'combined_seqs.fna'))
        outfile = os.path.join(self.temp_dir.name, 'sims.b6')
        tax = pkg_resources.resource_filename('shogun.tests', os.path.join('data', 'genomes.small.tax'))
        self.assertTrue(burst_align(infile, outfile, database, tax=tax)[0] == 0)
        self.assertTrue(embalmulate(outfile, self.temp_dir.name)[0] == 0)

    def test_burst_build(self):
        fasta = pkg_resources.resource_filename('shogun.tests', os.path.join('data', 'genomes.small.fna'))
        outfile = os.path.join(self.temp_dir.name, 'genomes.small')
        burst_build(fasta, outfile, shell=False, clustradius=1050, shear=500)

        for file in os.listdir(self.temp_dir.name):
           # self.assertTrue(self.checksums[hash_file(os.path.join(self.temp_dir.name, file))] == file)
           continue

if __name__ == '__main__':
    unittest.main()
