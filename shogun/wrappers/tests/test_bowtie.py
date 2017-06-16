"""
Copyright 2015-2017 Knights Lab, Regents of the University of Minnesota.

This software is released under the GNU Affero General Public License (AGPL) v3.0 License.
"""

import unittest
import shutil
import pkg_resources
import os
import tempfile
from collections import defaultdict

from shogun.utils import hash_file
from shogun.wrappers.bowtie import bowtie2_align, bowtie2_build

def read_checksums(filename):
    with open(filename) as inf:
        return defaultdict(str, dict([line.split() for line in inf]))

class TestBowtie(unittest.TestCase):
    def setUp(self):
        prefix = "shogun-test-temp-"

        self.checksums = read_checksums(pkg_resources.resource_filename(
            'shogun.wrappers.tests', os.path.join('data', 'checksums.txt')))

        self.temp_dir = tempfile.TemporaryDirectory(prefix=prefix)

    def tearDown(self):
        self.temp_dir.cleanup()

    def test_bowtie2_path(self):
        self.assertTrue(shutil.which("bowtie2") != None)

    def test_bowtie2_align(self):
        pass

    def test_bowtie2_build(self):
        fasta = pkg_resources.resource_filename('shogun.wrappers.tests', os.path.join('data', 'genomes.small.fna'))
        tax = pkg_resources.resource_filename('shogun.wrappers.tests', os.path.join('data', 'genomes.small.tax'))
        outfile = os.path.join(self.temp_dir.name, 'bowtie2-test')
        bowtie2_build(fasta, outfile, shell=False)

        for file in os.listdir(self.temp_dir.name):
            self.assertTrue(self.checksums[hash_file(os.path.join(self.temp_dir.name, file))] == file)

if __name__ == '__main__':
    unittest.main()
