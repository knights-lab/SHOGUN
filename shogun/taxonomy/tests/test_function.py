"""
Copyright 2015-2017 Knights Lab, Regents of the University of Minnesota.

This software is released under the GNU Affero General Public License (AGPL) v3.0 License.
"""

import unittest
import pkg_resources
import os
from yaml import load

from shogun.__main__ import _check_function_db

class TestFunctionCheck(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_read_taxonomy(self):
        database = pkg_resources.resource_filename('shogun.tests', os.path.join('data'))
        with open(os.path.join(database, 'metadata.yaml'), 'r') as stream:
            data_files = load(stream)
        results = _check_function_db(data_files, database)
        self.assertTrue(True)
        print(results)
