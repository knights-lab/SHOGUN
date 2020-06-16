"""
Copyright 2015-2020 Knights Lab, Regents of the University of Minnesota.

This software is released under the GNU Affero General Public License (AGPL) v3.0 License.
"""

import unittest
import pkg_resources
import os

from shogun.redistribute import parse_bayes, redistribute_taxatable
from shogun.utils.tree import Taxonomy


class TestRedistribute(unittest.TestCase):
    def test_read_taxonomy(self):
        tax = pkg_resources.resource_filename('shogun.tests', os.path.join('data', 'genomes.small.tax'))
        d = Taxonomy.parse_taxonomy(tax)
        self.assertTrue(d['NC_002182.1'] == 'k__Bacteria;p__Chlamydiae;c__Chlamydiia;o__Chlamydiales;f__Chlamydiaceae;g__Chlamydia;s__Chlamydia_muridarum;t__Chlamydia_muridarum_str._Nigg')

    def test_read_bayes(self):
        bayes = pkg_resources.resource_filename('shogun.tests', os.path.join('data', 'sheared_bayes.100.txt'))
        df = parse_bayes(bayes)

    def test_taxatable(self):
        bayes = pkg_resources.resource_filename('shogun.tests', os.path.join('data', 'sheared_bayes.32.txt'))
        df_bayes = parse_bayes(bayes)
        taxatable = pkg_resources.resource_filename('shogun.tests', os.path.join('data', 'results', 'burst_taxatable.txt'))
        taxatable_df_5 = redistribute_taxatable(taxatable, df_bayes, level=5)
        taxatable_df_6 = redistribute_taxatable(taxatable, df_bayes, level=6)
        taxatable_df_7 = redistribute_taxatable(taxatable, df_bayes, level=7)
        taxatable_df_8 = redistribute_taxatable(taxatable, df_bayes, level=8)
        taxatable_df_8.head()
