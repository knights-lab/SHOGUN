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


class TestRefSeq(unittest.TestCase):
    def setUp(self):
        prefix = "shogun-test-temp-"

        self.checksums = read_checksums(pkg_resources.resource_filename(
            'shogun.tests', os.path.join('data', 'bowtie2', 'checksums.txt')))

        self.temp_dir = tempfile.TemporaryDirectory(prefix=prefix)

    def tearDown(self):
        self.temp_dir.cleanup()

    def test_taxonkit_in_path(self):
        self.assertTrue(shutil.which("taxonkit") is not None)

    # def test_get_taxdmp(self):
    #     """
    #     This test take a long time because we must download the entire NCBI tax database
    #
    #     :return:
    #     """
    #     get_taxdump(self.temp_dir.name)




if __name__ == '__main__':
    unittest.main()
