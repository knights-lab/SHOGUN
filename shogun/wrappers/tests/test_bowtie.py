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
from shogun.wrappers.bowtie2_wrapper import bowtie2_align, bowtie2_build


class TestBowtie(unittest.TestCase):
    def setUp(self):
        prefix = "shogun-test-temp-"

        self.checksums = read_checksums(pkg_resources.resource_filename(
            'shogun.tests', os.path.join('data', 'bowtie2', 'checksums.txt')))

        self.temp_dir = tempfile.TemporaryDirectory(prefix=prefix)

    def tearDown(self):
        self.temp_dir.cleanup()

    def test_bowtie2_path(self):
        self.assertTrue(shutil.which("bowtie2") is not None)

    def test_bowtie2_align(self):
        database = pkg_resources.resource_filename('shogun.tests', os.path.join('data', 'bowtie2', 'genomes.small'))
        infile = pkg_resources.resource_filename('shogun.tests', os.path.join('data', 'combined_seqs.fna'))
        outfile = os.path.join(self.temp_dir.name, 'sims.sam')
        self.assertTrue(bowtie2_align(infile, outfile, database)[0] == 0)
        self.assertTrue(os.path.isfile(outfile) and os.path.getsize(outfile) > 0)

    def test_bowtie2_build(self):
        fasta = pkg_resources.resource_filename('shogun.tests', os.path.join('data', 'genomes.small.fna'))
        outfile = os.path.join(self.temp_dir.name, 'genomes.small')
        # Test if system exit call is non-zero
        self.assertTrue(bowtie2_build(fasta, outfile, shell=False)[0] == 0)
        #TODO: Proper database file sniffing

if __name__ == '__main__':
    unittest.main()
