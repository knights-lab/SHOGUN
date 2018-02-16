"""
Copyright 2015-2017 Knights Lab, Regents of the University of Minnesota.

This software is released under the GNU Affero General Public License (AGPL) v3.0 License.
"""
import unittest

import pandas as pd
import numpy as np

from shogun.utils import convert_to_relative_abundance

class TestRelativeAbundance(unittest.TestCase):
    def test_convert_relative_abundance(self):
        df = pd.DataFrame([[1, 1, 0], [1, 0, 0], [1, 0, 0]])
        df_expected = pd.DataFrame([[1./3, 1., np.nan], [1./3, 0., np.nan], [1./3, 0., np.nan]])
        df_ra = convert_to_relative_abundance(df)
        pd.testing.assert_frame_equal(df_ra, df_expected)

