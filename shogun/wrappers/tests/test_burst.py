"""
Copyright 2015-2017 Knights Lab, Regents of the University of Minnesota.

This software is released under the GNU Affero General Public License (AGPL) v3.0 License.
"""

import unittest
import shutil
import pkg_resources
import os
import tempfile

from shogun.wrappers.burst_wrapper import burst_align, burst_build


class TestBurst(unittest.TestCase):
    def setUp(self):
        prefix = "shogun-test-temp-"

        self.temp_dir = tempfile.TemporaryDirectory(prefix=prefix)

    def tearDown(self):
        self.temp_dir.cleanup()

    def test_burst_path(self):
        self.assertTrue(shutil.which("burst15") is not None)

    def test_burst_align(self):
        database = pkg_resources.resource_filename('shogun.tests', os.path.join('data', 'burst', 'genomes.small'))
        infile = pkg_resources.resource_filename('shogun.tests', os.path.join('data', 'combined_seqs.fna'))
        outfile = os.path.join(self.temp_dir.name, 'sims.b6')
        tax = pkg_resources.resource_filename('shogun.tests', os.path.join('data', 'genomes.small.tax'))
        self.assertTrue(burst_align(infile, outfile, database, tax=tax)[0] == 0)
        self.assertTrue(os.path.isfile(outfile) and os.path.getsize(outfile) > 0)

    def test_burst_build(self):
        fasta = pkg_resources.resource_filename('shogun.tests', os.path.join('data', 'genomes.small.fna'))
        outfile = os.path.join(self.temp_dir.name, 'genomes.small')
        print(burst_build(fasta, outfile, shell=False, clustradius=1050, shear=500))
        self.assertTrue(burst_build(fasta, outfile, shell=False, clustradius=0, shear=500)[0] == 0)
        # TODO: Proper Database file sniffing


if __name__ == '__main__':
    unittest.main()
