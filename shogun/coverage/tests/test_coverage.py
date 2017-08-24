"""
Copyright 2015-2017 Knights Lab, Regents of the University of Minnesota.

This software is released under the GNU Affero General Public License (AGPL) v3.0 License.
"""
import unittest
import pkg_resources
import os
import tempfile

from shogun.coverage import get_coverage_of_microbes
from shogun.redistribute import parse_bayes


class TestCoverage(unittest.TestCase):
    def setUp(self):
        prefix = 'shogun-temp-dir-'
        self.temp_dir = tempfile.TemporaryDirectory(prefix=prefix)

    def tearDown(self):
        self.temp_dir.cleanup()

    def test_coverage_report(self):
        bayes = pkg_resources.resource_filename('shogun.tests', os.path.join('data', 'sheared_bayes.32.txt'))
        df_bayes = parse_bayes(bayes)
        taxatable = pkg_resources.resource_filename('shogun.tests', os.path.join('data', 'results', 'burst_taxatable.txt'))
        infile = pkg_resources.resource_filename('shogun.tests', os.path.join('data', 'results', 'burst_results.b6'))

        self.assertTrue(get_coverage_of_microbes(infile, df_bayes, 6) is not None)
        self.assertTrue(get_coverage_of_microbes(infile, df_bayes, 7) is not None)
        self.assertTrue(get_coverage_of_microbes(infile, df_bayes, 8) is not None)

if __name__ == '__main__':
    unittest.main()
