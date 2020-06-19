"""
Copyright 2015-2020 Knights Lab, Regents of the University of Minnesota.

This software is released under the GNU Affero General Public License (AGPL) v3.0 License.
"""

import unittest
import shutil
import pkg_resources
import os
import tempfile

from shogun.utils import hash_file, read_checksums
#TODO: Implement utree rank flexible testing
from shogun.wrappers.utree_wrapper import utree_build_gg, utree_compress, utree_search_gg


class TestUtree(unittest.TestCase):
    def setUp(self):
        prefix = "shogun-test-temp-"

        self.checksums = read_checksums(pkg_resources.resource_filename(
            'shogun.tests', os.path.join('data', 'utree', 'checksums.txt')))

        self.temp_dir = tempfile.TemporaryDirectory(prefix=prefix)

    def tearDown(self):
        self.temp_dir.cleanup()

    def test_utree_path(self):
        self.assertTrue(shutil.which("utree-build_gg") is not None)
        self.assertTrue(shutil.which("utree-compress") is not None)
        self.assertTrue(shutil.which("utree-search_gg") is not None)

    def test_utree_build_gg(self):
        fasta = pkg_resources.resource_filename('shogun.tests', os.path.join('data', 'genomes.small.fna'))
        tax = pkg_resources.resource_filename('shogun.tests', os.path.join('data', 'genomes.small.tax'))
        outfile_uncompressed = os.path.join(self.temp_dir.name, 'genomes.small.gg.utr')
        outfile_compressed = os.path.join(self.temp_dir.name, 'genomes.small.gg.ctr')
        utree_build_gg(fasta, tax, outfile_uncompressed, shell=False)
        self.assertTrue(os.path.isfile(outfile_uncompressed) and os.path.getsize(outfile_uncompressed) > 0)
        utree_compress(outfile_uncompressed, outfile_compressed, shell=False)
        self.assertTrue(os.path.isfile(outfile_compressed) and os.path.getsize(outfile_compressed) > 0)

    def test_utree_align_gg(self):
        database = pkg_resources.resource_filename('shogun.tests', os.path.join('data', 'utree', 'genomes.small.gg.ctr'))
        infile = pkg_resources.resource_filename('shogun.tests', os.path.join('data', 'combined_seqs.fna'))
        outfile = os.path.join(self.temp_dir.name, 'utree_gg-test-sims.txt')
        self.assertTrue(utree_search_gg(database, infile, outfile)[0] == 0)
        # TODO: Proper database outfile sniffing
        self.assertTrue(os.path.isfile(outfile) and os.path.getsize(outfile) > 0)

if __name__ == '__main__':
    unittest.main()
