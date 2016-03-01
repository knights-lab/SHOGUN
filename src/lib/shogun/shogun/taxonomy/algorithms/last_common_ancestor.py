#!/usr/bin/env python

from shogun.taxonomy.ncbi_tree import NCBITree


class LCA:
    def __init__(self, tree=NCBITree.load()):
        """

        :param tree: NCBITree
        :return:
        """
        self.tree_obj = tree

    def __call__(self, taxon_id_a, taxon_id_b):
        self.apply(taxon_id_a, taxon_id_b)

    def apply(self, taxon_id_a, taxon_id_b):
        # assumes that all nodes in the tree have a common root
        if taxon_id_a and taxon_id_b in self.tree_obj.tree.nodes():
            taxonomy_a = self.tree_obj.get_taxon_id_lineage_with_taxon_id(taxon_id_a)
            taxonomy_b = self.tree_obj.get_taxon_id_lineage_with_taxon_id(taxon_id_b)
            i = 0
            for taxon_ids in zip(reversed(taxonomy_a), reversed(taxonomy_b)):
                if taxon_ids[0] != taxon_ids[1]:
                        break
                i += 1
            return taxonomy_a[-i:]
        else:
            return []


def main():
    LCA().apply(1292, 106)

if __name__ == '__main__':
    main()
