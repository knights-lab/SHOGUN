"""
Copyright 2015-2017 Knights Lab, Regents of the University of Minnesota.

This software is released under the GNU Affero General Public License (AGPL) v3.0 License.
"""
import unittest
import pkg_resources
import os
import tempfile

from shogun.aligners import BurstAligner, UtreeAligner, BowtieAligner


class TestAligner(unittest.TestCase):
    def setUp(self):
        prefix = 'shogun-temp-dir-'
        self.temp_dir = tempfile.TemporaryDirectory(prefix=prefix)

    def tearDown(self):
        self.temp_dir.cleanup()

    def test_burst_db(self):
        self.assertTrue(BurstAligner.check_database(
            pkg_resources.resource_filename('shogun.tests', os.path.join('data')))[0])

    def test_burst_align(self):
        database = pkg_resources.resource_filename('shogun.tests', os.path.join('data'))
        aligner = BurstAligner(database, threads=1, shell=False)

        infile = pkg_resources.resource_filename('shogun.tests', os.path.join('data', 'combined_seqs.fna'))
        outdir = os.path.join(self.temp_dir.name)

        self.assertTrue(aligner.align(infile, outdir)[0] == 0)

    def test_burst_post_align(self):
        database = pkg_resources.resource_filename('shogun.tests', os.path.join('data'))
        aligner = BurstAligner(database, threads=1, shell=False)
        alignment_file = pkg_resources.resource_filename('shogun.tests', os.path.join('data', 'results', 'burst_results.b6'))
        df_capitalist = aligner._post_align(alignment_file)
        aligner.post_align = 'taxonomy'
        df_non_capitalist = aligner._post_align(alignment_file)
        self.assertTrue(df_non_capitalist.any().any())

    def test_utree_db(self):
        self.assertTrue(UtreeAligner.check_database(
            pkg_resources.resource_filename('shogun.tests', os.path.join('data')))[0])

    def test_utree_align(self):
        database = pkg_resources.resource_filename('shogun.tests', os.path.join('data'))
        aligner = UtreeAligner(database, shell=False)

        infile = pkg_resources.resource_filename('shogun.tests', os.path.join('data', 'combined_seqs.fna'))
        outdir = os.path.join(self.temp_dir.name)

        self.assertTrue(aligner.align(infile, outdir)[0] == 0)

    def test_bowtie2_db(self):
        self.assertTrue(BowtieAligner.check_database(
            pkg_resources.resource_filename('shogun.tests', os.path.join('data')))[0])

    def test_bowtie2_align(self):
        database = pkg_resources.resource_filename('shogun.tests', os.path.join('data'))
        aligner = BowtieAligner(database, shell=False)

        infile = pkg_resources.resource_filename('shogun.tests', os.path.join('data', 'combined_seqs.fna'))
        outdir = os.path.join(self.temp_dir.name)

        self.assertTrue(aligner.align(infile, outdir)[0] == 0)

    def test_bowtie2_post_align(self):
        database = pkg_resources.resource_filename('shogun.tests', os.path.join('data'))
        aligner = BowtieAligner(database, threads=1, shell=False)
        alignment_file = pkg_resources.resource_filename('shogun.tests', os.path.join('data', 'results', 'bowtie2_results.sam'))
        df = aligner._post_align(alignment_file)
        self.assertTrue(df.any().any())

if __name__ == '__main__':
    unittest.main()
