"""
Copyright 2015-2020 Knights Lab, Regents of the University of Minnesota.

This software is released under the GNU Affero General Public License (AGPL) v3.0 License.
"""
import unittest

import pandas as pd
import numpy as np

from shogun.utils import convert_to_relative_abundance, least_common_ancestor


class TestRelativeAbundance(unittest.TestCase):
    def test_convert_relative_abundance(self):
        df = pd.DataFrame([[1, 1, 0], [1, 0, 0], [1, 0, 0]])
        df_expected = pd.DataFrame([[1./3, 1., np.nan], [1./3, 0., np.nan], [1./3, 0., np.nan]])
        df_ra = convert_to_relative_abundance(df)
        pd.testing.assert_frame_equal(df_ra, df_expected)


class TestLeastCommonAncestor(unittest.TestCase):
    def test_lca_family(self):
        # normal scenario
        taxa = [
            'k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Clostridiaceae;g__Clostridium;s__Clostridium acetobutylicum',
            'k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Clostridiaceae;g__Clostridium;s__Clostridium botulinum',
            'k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Clostridiaceae;g__Clostridium;s__Clostridium aldrichii',
            'k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Clostridiaceae;g__Caminicella;s__Caminicella sporogenes'
        ]
        obs = least_common_ancestor(taxa)
        exp = 'k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Clostridiaceae'
        self.assertEqual(obs, exp)

    def test_lca_kingdom(self):
        # scenario 1: lowest level is unclassified
        taxa = [
            'k__Bacteria;p__Deinococcus-Thermus;c__Deinococci;o__Thermales;f__Thermaceae;g__Thermus;s__Thermus thermophilus;t__',
            'k__Bacteria;p__Deinococcus-Thermus;c__Deinococci;o__Deinococcales;f__Deinococcaceae;g__Deinococcus;s__Deinococcus radiodurans;t__',
            'k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Clostridiaceae;g__Clostridium;s__Clostridium acetobutylicum;t__',
            'k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Clostridiaceae;g__Clostridium;s__Clostridium botulinum;t__'
        ]
        obs = least_common_ancestor(taxa)
        exp = 'k__Bacteria'
        self.assertEqual(obs, exp)

    def test_lca_phylum(self):
        # scenario 2: middle level is unclassified
        taxa = [
            'k__Viruses;p__ssDNA_viruses;c__Geminiviridae;o__;f__;g__Begomovirus;s__Sida_mosaic_Sinaloa_virus',
            'k__Viruses;p__ssDNA_viruses;c__Geminiviridae;o__;f__;g__Mastrevirus;s__Panicum_streak_virus',
            'k__Viruses;p__ssDNA_viruses;c__Geminiviridae;o__;f__;g__Begomovirus;s__Cotton_leaf_crumple_virus',
            'k__Viruses;p__ssDNA_viruses;c__Nanoviridae;o__;f__;g__Babuvirus;s__Banana_bunchy_top_virus'
        ]
        obs = least_common_ancestor(taxa)
        exp = 'k__Viruses;p__ssDNA_viruses'
        self.assertEqual(obs, exp)

    def test_lca_none(self):
        # scenario 3: top level is inconsistent
        taxa = [
            'k__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Vibrionales;f__Vibrionaceae;g__Vibrio;s__Vibrio_tasmaniensis;t__Vibrio_tasmaniensis_LGP32',
            'k__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Vibrionales;f__Vibrionaceae;g__Vibrio;s__Vibrio_maritimus;t__',
            'k__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Vibrionales;f__Vibrionaceae;g__Aliivibrio;s__Aliivibrio_wodanis;t__',
            'k__BacteriaPlasmid;p__Proteobacteria;c__Gammaproteobacteria;o__Vibrionales;f__Vibrionaceae;g__Aliivibrio;s__Aliivibrio_wodanis;t__'
        ]
        obs = least_common_ancestor(taxa)
        self.assertIsNone(obs)

    def test_blank_class(self):
        # normal scenario
        taxa = [
            'k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Clostridiaceae;g__Clostridium;s__Clostridium acetobutylicum',
            'k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Clostridiaceae;g__Clostridium;s__Clostridium botulinum',
            'k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Clostridiaceae;g__Clostridium;s__Clostridium aldrichii',
            'k__Bacteria;p__Firmicutes;c__;o__Clostridiales;f__Clostridiaceae;g__Caminicella;s__Caminicella sporogenes'
        ]
        obs = least_common_ancestor(taxa)
        exp = 'k__Bacteria;p__Firmicutes'
        self.assertEqual(obs, exp)

    def test_blanks_kingdom_class(self):
        # unclassified until family
        taxa = [
            'k__;p__Firmicutes;c__;o__Clostridiales;f__Clostridiaceae;g__Clostridium',
            'k__;p__Firmicutes;c__;o__Clostridiales;f__Clostridiaceae;g__Clostridium',
            'k__;p__Firmicutes;c__;o__Clostridiales;f__Clostridiaceae;g__Clostridium',
            'k__;p__Firmicutes;c__;o__Clostridiales;f__Clostridiaceae;g__Caminicella'
        ]
        obs = least_common_ancestor(taxa)
        exp = 'k__;p__Firmicutes;c__;o__Clostridiales;f__Clostridiaceae'
        self.assertEqual(obs, exp)

    def test_continued_blanks(self):
        # unclassified until phylum
        taxa = [
            'k__;p__Firmicutes;c__;o__;f__;g__Clostridium',
            'k__;p__Firmicutes;c__;o__;f__;g__Clostridium',
            'k__;p__Firmicutes;c__;o__;f__;g__Clostridium',
            'k__;p__Firmicutes;c__;o__;f__;g__Caminicella'
        ]
        obs = least_common_ancestor(taxa)
        exp = 'k__;p__Firmicutes'
        self.assertEqual(obs, exp)
