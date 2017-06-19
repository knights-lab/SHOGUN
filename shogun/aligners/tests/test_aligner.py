"""
Copyright 2015-2017 Knights Lab, Regents of the University of Minnesota.

This software is released under the GNU Affero General Public License (AGPL) v3.0 License.
"""
import unittest
import shutil
import pkg_resources
import os
import tempfile

from shogun.aligners import EmbalmerAligner


class TestEmbalmer(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_embalmer_db(self):
        self.assertTrue(EmbalmerAligner.check_database(
            pkg_resources.resource_filename('shogun.wrappers.tests', os.path.join('data'))))

if __name__ == '__main__':
    unittest.main()
