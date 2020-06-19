"""
Copyright 2015-2020 Knights Lab, Regents of the University of Minnesota.

This software is released under the GNU Affero General Public License (AGPL) v3.0 License.
"""
import unittest
import pkg_resources
from click.testing import CliRunner
import os
import tempfile
import pandas as pd

import glob

from shogun.__main__ import cli


class TestAligner(unittest.TestCase):
    def setUp(self):
        prefix = 'shogun-temp-dir-'
        self.temp_dir = tempfile.TemporaryDirectory(prefix=prefix)

    def tearDown(self):
        self.temp_dir.cleanup()

    def test_utree_pipeline(self):
        database = pkg_resources.resource_filename('shogun.tests', os.path.join('data'))
        infile = pkg_resources.resource_filename('shogun.tests', os.path.join('data', 'combined_seqs.fna'))

        outdir = os.path.join(self.temp_dir.name)

        runner = CliRunner()
        _log = runner.invoke(cli, ['--log', 'debug', 'pipeline', '--input', infile, '--database', database,
                                      '--output', outdir, '--aligner', 'utree', '--no-function'])

        outfile_ra = glob.glob(os.path.join(outdir, "*.ra.txt"))

        self.assertTrue(len(outfile_ra) == 1)

        outfile_ra = outfile_ra[0]

        df_infile = pd.read_csv(outfile_ra, sep="\t", index_col=0)

        # Assert the correct number of samples
        self.assertTrue(df_infile.shape[1] == 3)

        # Assert the type is float
        self.assertTrue(len(df_infile.select_dtypes(include=['float']).columns) == 3)

    def test_bowtie2_pipeline(self):
        database = pkg_resources.resource_filename('shogun.tests', os.path.join('data'))
        infile = pkg_resources.resource_filename('shogun.tests', os.path.join('data', 'combined_seqs.fna'))

        outdir = os.path.join(self.temp_dir.name)

        runner = CliRunner()
        _log = runner.invoke(cli, ['--log', 'debug', 'pipeline', '--input', infile, '--database', database,
                                      '--output', outdir, '--aligner', 'bowtie2', '--no-function'])

        outfile_ra = glob.glob(os.path.join(outdir, "*.ra.txt"))

        self.assertTrue(len(outfile_ra) == 1)

        outfile_ra = outfile_ra[0]

        df_infile = pd.read_csv(outfile_ra, sep="\t", index_col=0)

        # Assert the correct number of samples
        self.assertTrue(df_infile.shape[1] == 3)

        # Assert the type is float
        self.assertTrue(len(df_infile.select_dtypes(include=['float']).columns) == 3)

    def test_burst_pipeline(self):
        database = pkg_resources.resource_filename('shogun.tests', os.path.join('data'))
        infile = pkg_resources.resource_filename('shogun.tests', os.path.join('data', 'combined_seqs.fna'))

        outdir = os.path.join(self.temp_dir.name)

        runner = CliRunner()
        _log = runner.invoke(cli, ['--log', 'debug', 'pipeline', '--input', infile, '--database', database,
                                      '--output', outdir, '--aligner', 'burst', '--no-function'])

        outfile_ra = glob.glob(os.path.join(outdir, "*.ra.txt"))

        self.assertTrue(len(outfile_ra) == 1)

        outfile_ra = outfile_ra[0]

        df_infile = pd.read_csv(outfile_ra, sep="\t", index_col=0)

        # Assert the correct number of samples
        self.assertTrue(df_infile.shape[1] == 3)

        # Assert the type is float
        self.assertTrue(len(df_infile.select_dtypes(include=['float']).columns) == 3)


if __name__ == '__main__':
    unittest.main()
