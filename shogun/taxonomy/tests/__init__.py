"""
Copyright 2015-2017 Knights Lab, Regents of the University of Minnesota.

This software is released under the GNU Affero General Public License (AGPL) v3.0 License.
"""

import unittest
import pkg_resources
import os

from shogun.taxonomy import Taxonomy, parse_bayes

class TestTaxonomy(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_read_taxonomy(self):
        tax = pkg_resources.resource_filename('shogun.tests', os.path.join('data', 'genomes.small.tax'))
        d = Taxonomy.parse_taxonomy(tax)
        self.assertTrue(d['NC_002182.1'] == 'k__Bacteria;p__Chlamydiae;c__Chlamydiia;o__Chlamydiales;f__Chlamydiaceae;g__Chlamydia;s__Chlamydia_muridarum;t__Chlamydia_muridarum_str._Nigg')

    def test_read_bayes(self):
        bayes = pkg_resources.resource_filename('shogun.tests', os.path.join('data', 'sheared_bayes.txt'))
        parse_bayes(bayes)
