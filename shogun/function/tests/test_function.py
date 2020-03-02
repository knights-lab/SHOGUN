"""
Copyright 2015-2020 Knights Lab, Regents of the University of Minnesota.

This software is released under the GNU Affero General Public License (AGPL) v3.0 License.
"""

import os
import tempfile
import unittest

import pkg_resources
import yaml

from shogun.__main__ import _function
from shogun.function import parse_function_db


class TestFunctionCheck(unittest.TestCase):
    def setUp(self):
        prefix = 'shogun-temp-dir-'
        self.temp_dir = tempfile.TemporaryDirectory(prefix=prefix)

    def tearDown(self):
        self.temp_dir.cleanup()

    def test_check_database(self):
        database = pkg_resources.resource_filename('shogun.tests', os.path.join('data'))
        with open(os.path.join(database, 'metadata.yaml'), 'r') as stream:
            data_files = yaml.load(stream, Loader=yaml.SafeLoader)
        results = parse_function_db(data_files, database)
        self.assertTrue(results is not None)

    def test_function(self):
        taxatable = pkg_resources.resource_filename('shogun.tests', os.path.join('data', 'results', 'burst_taxatable.txt'))
        database = pkg_resources.resource_filename('shogun.tests', os.path.join('data'))

        outdir = os.path.join(self.temp_dir.name)
        # Strain
        _function([taxatable, taxatable, taxatable], database, outdir, ['genus', 'species', 'strain'])
