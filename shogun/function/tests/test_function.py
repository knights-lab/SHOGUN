"""
Copyright 2015-2017 Knights Lab, Regents of the University of Minnesota.

This software is released under the GNU Affero General Public License (AGPL) v3.0 License.
"""

import unittest
import pkg_resources
import os
from yaml import load

from shogun.__main__ import _parse_function_db, _parse_kegg_table

class TestFunctionCheck(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_check_database(self):
        database = pkg_resources.resource_filename('shogun.tests', os.path.join('data'))
        with open(os.path.join(database, 'metadata.yaml'), 'r') as stream:
            data_files = load(stream)
        results = _parse_function_db(data_files, database)
        self.assertTrue(results is not None)

    def test_function(self):
        taxatable = pkg_resources.resource_filename('shogun.tests', os.path.join('data', 'results', 'embalmer_taxatable.txt'))
        database =         database = pkg_resources.resource_filename('shogun.tests', os.path.join('data'))
        with open(os.path.join(database, 'metadata.yaml'), 'r') as stream:
            data_files = load(stream)
        function_db = _parse_function_db(data_files, database)

        kegg_table = _parse_kegg_table(function_db['strain'])
        print(kegg_table.head())


