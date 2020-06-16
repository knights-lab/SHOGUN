import unittest
import pkg_resources
import os
import tempfile

import pandas as pd

from shogun.utils.lowest_common_ancestor import gen_confidence_lowest_common_ancestor, build_lca_df, gen_lowest_common_ancestor
from shogun.utils.tree import build_tree_from_tax_file


class TestLowestCommonAncestor(unittest.TestCase):
    def setUp(self):
        prefix = "shogun-test-temp-"
        self.temp_dir = tempfile.TemporaryDirectory(prefix=prefix)
        taxfile = pkg_resources.resource_filename('shogun.tests', os.path.join('data', 'tree.tax'))
        self.tree = build_tree_from_tax_file(taxfile)
        tax = pkg_resources.resource_filename('shogun.tests', os.path.join('data', 'genomes.small.tax'))
        self.genomes_small_taxatree = build_tree_from_tax_file(tax)

    def tearDown(self):
        self.temp_dir.cleanup()

    def test_num_samples_greater_than(self):
        alignment_file = pkg_resources.resource_filename('shogun.tests', os.path.join('data', 'results', 'bowtie2_results.sam'))
        df_test_small_iter = build_lca_df(alignment_file, self.genomes_small_taxatree, samples_iter=2)
        df_test_large_iter = build_lca_df(alignment_file, self.genomes_small_taxatree, samples_iter=20)
        pd.testing.assert_frame_equal(df_test_large_iter, df_test_small_iter)

    def yield_records(self, l):
        for ix, align in enumerate(l):
            yield [(ix, _) for _ in align]


    def test_lca_family(self):
        # normal scenario

        # taxa = [
        #     'k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Clostridiaceae;g__Clostridium;s__Clostridium acetobutylicum',
        #     'k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Clostridiaceae;g__Clostridium;s__Clostridium botulinum',
        #     'k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Clostridiaceae;g__Clostridium;s__Clostridium aldrichii',
        #     'k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Clostridiaceae;g__Caminicella;s__Caminicella sporogenes'
        # ]

        taxa = (('1', '2', '3', '4'),)
        obs = next(gen_lowest_common_ancestor(self.yield_records(taxa), self.tree))
        obs = self.tree.node_id_to_taxa_name[obs[1]]
        exp = 'k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Clostridiaceae'
        self.assertEqual(exp, obs)

    def test_lca_kingdom(self):
        # scenario 1: lowest level is unclassified

        # taxa = [
        #     'k__Bacteria;p__Deinococcus-Thermus;c__Deinococci;o__Thermales;f__Thermaceae;g__Thermus;s__Thermus thermophilus;t__',
        #     'k__Bacteria;p__Deinococcus-Thermus;c__Deinococci;o__Deinococcales;f__Deinococcaceae;g__Deinococcus;s__Deinococcus radiodurans;t__',
        #     'k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Clostridiaceae;g__Clostridium;s__Clostridium acetobutylicum;t__',
        #     'k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Clostridiaceae;g__Clostridium;s__Clostridium botulinum;t__'
        # ]

        taxa = (('5', '6', '7', '8'),)

        obs = next(gen_lowest_common_ancestor(self.yield_records(taxa), self.tree))
        obs = self.tree.node_id_to_taxa_name[obs[1]]
        exp = 'k__Bacteria'
        self.assertEqual(exp, obs)

    def test_lca_phylum(self):
        # scenario 2: middle level is unclassified
        # taxa = [
        #     'k__Viruses;p__ssDNA_viruses;c__Geminiviridae;o__;f__;g__Begomovirus;s__Sida_mosaic_Sinaloa_virus',
        #     'k__Viruses;p__ssDNA_viruses;c__Geminiviridae;o__;f__;g__Mastrevirus;s__Panicum_streak_virus',
        #     'k__Viruses;p__ssDNA_viruses;c__Geminiviridae;o__;f__;g__Begomovirus;s__Cotton_leaf_crumple_virus',
        #     'k__Viruses;p__ssDNA_viruses;c__Nanoviridae;o__;f__;g__Babuvirus;s__Banana_bunchy_top_virus'
        # ]

        taxa = (('9', '10', '11', '12'),)

        obs = next(gen_lowest_common_ancestor(self.yield_records(taxa), self.tree))
        obs = self.tree.node_id_to_taxa_name[obs[1]]
        exp = 'k__Viruses;p__ssDNA_viruses'
        self.assertEqual(exp, obs)

    def test_lca_none(self):
        # scenario 3: top level is inconsistent
        # taxa = [
        #     'k__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Vibrionales;f__Vibrionaceae;g__Vibrio;s__Vibrio_tasmaniensis;t__Vibrio_tasmaniensis_LGP32',
        #     'k__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Vibrionales;f__Vibrionaceae;g__Vibrio;s__Vibrio_maritimus;t__',
        #     'k__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Vibrionales;f__Vibrionaceae;g__Aliivibrio;s__Aliivibrio_wodanis;t__',
        #     'k__BacteriaPlasmid;p__Proteobacteria;c__Gammaproteobacteria;o__Vibrionales;f__Vibrionaceae;g__Aliivibrio;s__Aliivibrio_wodanis;t__'
        # ]

        taxa = (('13', '14', '15', '16'),)

        obs = next(gen_lowest_common_ancestor(self.yield_records(taxa), self.tree))
        obs = self.tree.node_id_to_taxa_name[obs[1]]
        self.assertEqual('root', obs)

    def test_blank_class(self):
        # normal scenario

        # taxa = [
        #     'k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Clostridiaceae;g__Clostridium;s__Clostridium acetobutylicum',
        #     'k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Clostridiaceae;g__Clostridium;s__Clostridium botulinum',
        #     'k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Clostridiaceae;g__Clostridium;s__Clostridium aldrichii',
        #     'k__Bacteria;p__Firmicutes;c__;o__Clostridiales;f__Clostridiaceae;g__Caminicella;s__Caminicella sporogenes'
        # ]

        taxa = (('17', '18', '19', '20'),)

        obs = next(gen_lowest_common_ancestor(self.yield_records(taxa), self.tree))
        obs = self.tree.node_id_to_taxa_name[obs[1]]
        exp = 'k__Bacteria;p__Firmicutes'
        self.assertEqual(exp, obs)

    def test_blanks_kingdom_class(self):
        # unclassified until family

        # taxa = [
        #     'k__;p__Firmicutes;c__;o__Clostridiales;f__Clostridiaceae;g__Clostridium',
        #     'k__;p__Firmicutes;c__;o__Clostridiales;f__Clostridiaceae;g__Clostridium',
        #     'k__;p__Firmicutes;c__;o__Clostridiales;f__Clostridiaceae;g__Clostridium',
        #     'k__;p__Firmicutes;c__;o__Clostridiales;f__Clostridiaceae;g__Caminicella'
        # ]

        taxa = (('21', '22', '23', '24'),)

        obs = next(gen_lowest_common_ancestor(self.yield_records(taxa), self.tree))
        obs = self.tree.node_id_to_taxa_name[obs[1]]
        exp = 'k__;p__Firmicutes;c__;o__Clostridiales;f__Clostridiaceae'
        self.assertEqual(obs, exp)

    def test_continued_blanks(self):
        # unclassified until phylum

        # taxa = [
        #     'k__;p__Firmicutes;c__;o__;f__;g__Clostridium',
        #     'k__;p__Firmicutes;c__;o__;f__;g__Clostridium',
        #     'k__;p__Firmicutes;c__;o__;f__;g__Clostridium',
        #     'k__;p__Firmicutes;c__;o__;f__;g__Caminicella'
        # ]

        taxa = (('25', '26', '27', '28'),)

        result = next(gen_lowest_common_ancestor(self.yield_records(taxa), self.tree))
        obs = self.tree.node_id_to_taxa_name[result[1]]
        exp = 'k__;p__Firmicutes'
        self.assertEqual(exp, obs)

    def test_all_match(self):
        # all match

        # taxa = [
        #     'k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Clostridiaceae;g__Clostridium;s__Clostridium acetobutylicum',
        #     'k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Clostridiaceae;g__Clostridium;s__Clostridium acetobutylicum',
        # ]

        taxa = (('29', '30'),)

        obs = next(gen_lowest_common_ancestor(self.yield_records(taxa), self.tree))
        obs = self.tree.node_id_to_taxa_name[obs[1]]
        exp = 'k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Clostridiaceae;g__Clostridium;s__Clostridium acetobutylicum'
        self.assertEqual(exp, obs)

    def test_all_match_different_length(self):
        # all match but different levels of the tree

        # taxa = [
        #     'k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Clostridiaceae;g__Clostridium;s__Clostridium acetobutylicum',
        #     'k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Clostridiaceae;g__Clostridium',
        #     'k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Clostridiaceae;g__Clostridium;s__Clostridium acetobutylicum',
        #     'k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Clostridiaceae;g__Clostridium;s__Clostridium acetobutylicum;t__'
        # ]

        taxa = (('31', '32', '33', '34'),)

        obs = next(gen_lowest_common_ancestor(self.yield_records(taxa), self.tree))
        obs = self.tree.node_id_to_taxa_name[obs[1]]
        exp = 'k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Clostridiaceae;g__Clostridium'
        self.assertEqual(exp, obs)

    def test_end_to_end_lca(self):
        alignment_file = pkg_resources.resource_filename('shogun.tests', os.path.join('data', 'truth.sam'))
        df_test = build_lca_df(alignment_file, self.genomes_small_taxatree, samples_iter=2)
        df_truth = pd.read_csv(pkg_resources.resource_filename('shogun.tests', os.path.join('data', 'truth.txt')), delimiter="\t", index_col=0)
        df_test.sort_index(inplace=True)
        df_truth.sort_index(inplace=True)
        df_test.index.name = "#OTU ID"
        pd.testing.assert_frame_equal(df_truth, df_test)


if __name__ == '__main__':
    unittest.main()
