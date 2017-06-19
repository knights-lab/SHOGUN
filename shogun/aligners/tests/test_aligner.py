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
        prefix = 'shogun-temp-dir-'
        self.temp_dir = tempfile.TemporaryDirectory(prefix=prefix)

    def tearDown(self):
        self.temp_dir.cleanup()

    def test_embalmer_db(self):
        self.assertTrue(EmbalmerAligner.check_database(
            pkg_resources.resource_filename('shogun.wrappers.tests', os.path.join('data'))))

    def test_embalmer_align(self):
        database = pkg_resources.resource_filename('shogun.wrappers.tests', os.path.join('data'))
        aligner = EmbalmerAligner(database, threads=1, shell=False)

        infile = pkg_resources.resource_filename('shogun.wrappers.tests', os.path.join('data', 'sims.fna'))
        outdir = os.path.join(self.temp_dir.name)

        self.assertTrue(aligner.align(infile, outdir)[0] == 0)



if __name__ == '__main__':
    unittest.main()
