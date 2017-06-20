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
#TODO: Implement utree rank flexible testing
from shogun.wrappers.utree import utree_build, utree_build_gg, utree_compress, utree_search, utree_search_gg


class TestUtree(unittest.TestCase):
    def setUp(self):
        prefix = "shogun-test-temp-"

        self.checksums = read_checksums(pkg_resources.resource_filename(
            'shogun.tests', os.path.join('data', 'utree', 'checksums.txt')))

        self.temp_dir = tempfile.TemporaryDirectory(prefix=prefix)

    def tearDown(self):
        self.temp_dir.cleanup()

    def test_utree_path(self):
        self.assertTrue(shutil.which("utree-build") is not None)
        self.assertTrue(shutil.which("utree-compress") is not None)
        self.assertTrue(shutil.which("utree-search") is not None)

    def test_utree_build(self):
        fasta = pkg_resources.resource_filename('shogun.tests', os.path.join('data', 'genomes.small.fna'))
        tax = pkg_resources.resource_filename('shogun.tests', os.path.join('data', 'genomes.small.tax'))
        outfile_uncompressed = os.path.join(self.temp_dir.name, 'genomes.small.utr')
        outfile_compressed = os.path.join(self.temp_dir.name, 'genomes.small.ctr')
        utree_build(fasta, tax, outfile_uncompressed, shell=False)
        utree_compress(outfile_uncompressed, outfile_compressed, shell=False)
        self.assertTrue(self.checksums[hash_file(outfile_compressed)] == os.path.basename(outfile_compressed))

    def test_utree_build_gg(self):
        pass

    def test_utree_align(self):
        database = pkg_resources.resource_filename('shogun.tests', os.path.join('data', 'utree', 'genomes.small.ctr'))
        infile = pkg_resources.resource_filename('shogun.tests', os.path.join('data', 'combined_seqs.fna'))
        outfile = os.path.join(self.temp_dir.name, 'bowtie2-test-sims.txt')
        self.assertTrue(utree_search(database, infile, outfile)[0] == 0)

    def test_utree_align_gg(self):
        pass

if __name__ == '__main__':
    unittest.main()
