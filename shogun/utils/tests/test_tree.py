"""
Copyright 2015-2020 Knights Lab, Regents of the University of Minnesota.

This software is released under the GNU Affero General Public License (AGPL) v3.0 License.
"""
import unittest
import pkg_resources
import os

from shogun.utils.tree import build_tree_from_tax_file


class TestTree(unittest.TestCase):
    def test_build_tree_from_tax_file(self):
        tax = pkg_resources.resource_filename('shogun.tests', os.path.join('data', 'genomes.small.tax'))
        taxatree = build_tree_from_tax_file(tax)
        assert taxatree.num_nodes == len(taxatree.node_id_to_ancestors) == len(taxatree.node_id_to_taxa_name)
