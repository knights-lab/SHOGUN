"""
Copyright 2015-2020 Knights Lab, Regents of the University of Minnesota.

This software is released under the GNU Affero General Public License (AGPL) v3.0 License.
"""

import unittest
import shutil
import tempfile


class TestRefSeq(unittest.TestCase):
    def setUp(self):
        prefix = "shogun-test-temp-"
        self.temp_dir = tempfile.TemporaryDirectory(prefix=prefix)

    def tearDown(self):
        self.temp_dir.cleanup()

    def test_taxonkit_in_path(self):
        self.assertTrue(shutil.which("taxonkit") is not None)


if __name__ == '__main__':
    unittest.main()
