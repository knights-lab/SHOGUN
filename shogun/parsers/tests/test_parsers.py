"""
Copyright 2015-2020 Knights Lab, Regents of the University of Minnesota.

This software is released under the GNU Affero General Public License (AGPL) v3.0 License.
"""

import unittest
import tempfile
import os
import pkg_resources

from shogun.parsers.sam import yield_alignments_from_sam_inf


class TestParsers(unittest.TestCase):
    def setUp(self):
        prefix = 'shogun-temp-dir-'
        self.temp_dir = tempfile.TemporaryDirectory(prefix=prefix)

    def tearDown(self):
        self.temp_dir.cleanup()

    def test_yield_alignments_from_samfile(self):
        inf_sam = pkg_resources.resource_filename('shogun.tests', os.path.join('data', 'results', 'bowtie2_results.sam'))
        gen = yield_alignments_from_sam_inf(inf_sam)
        i = len([1 for record in enumerate(gen)])
        assert i == 190
