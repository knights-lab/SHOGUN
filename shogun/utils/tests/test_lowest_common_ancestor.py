import unittest
import pkg_resources
import os

import pandas as pd

from shogun.utils.lowest_common_ancestor import gen_confidence_lowest_common_ancestor, build_lca_df
from shogun.utils.tree import build_tree_from_tax_file

class TestLowestCommonAncestor(unittest.TestCase):
    def test_num_samples_greater_than(self):
        tax = pkg_resources.resource_filename('shogun.tests', os.path.join('data', 'genomes.small.tax'))
        taxatree = build_tree_from_tax_file(tax)
        alignment_file = pkg_resources.resource_filename('shogun.tests', os.path.join('data', 'results', 'bowtie2_results.sam'))
        df_test_small_iter = build_lca_df(alignment_file, taxatree, samples_iter=2)
        df_test_large_iter = build_lca_df(alignment_file, taxatree, samples_iter=20)
        pd.testing.assert_frame_equal(df_test_large_iter, df_test_small_iter)


if __name__ == '__main__':
    unittest.main()
