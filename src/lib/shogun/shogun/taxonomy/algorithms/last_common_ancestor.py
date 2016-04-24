#!/usr/bin/env python

from shogun.taxonomy.ncbi_tree import NCBITree


class LCA:
    def __init__(self, tree, depth):
        """

        :param tree: NCBITree
        :return:
        """
        self.tree_obj = tree
        self.depth = depth

    def __call__(self, taxon_id_a, taxon_id_b):
        return self.apply(taxon_id_a, taxon_id_b)

    def apply(self, taxon_id_a, taxon_id_b):
        # assumes that all nodes in the tree have a common root
        # if taxon_id_a and taxon_id_b in self.tree_obj.tree.nodes():
        taxonomy_a = self.tree_obj.get_lineage_depth(taxon_id_a, self.depth)
        taxonomy_b = self.tree_obj.get_lineage_depth(taxon_id_b, self.depth)
        i = 0
        for taxon_ids in zip(reversed(taxonomy_a), reversed(taxonomy_b)):
            if taxon_ids[0] != taxon_ids[1]:
                    break
            i += 1
        if i == 0:
            return None
        else:
            return taxonomy_a[-i]

    def lca_mp(self, taxon_id, mp_taxonomy):
        taxonomy_1 = self.tree_obj.lca_mp(taxon_id)
        return lca_mp(taxonomy_1, mp_taxonomy)


def lca_mp(taxonomy_a, taxonomy_b):
    taxonomy_1 = taxonomy_a.split(';')
    taxonomy_2 = taxonomy_b.split(';')
    lca = []
    for i in zip(taxonomy_1, taxonomy_2):
        if i[0] == i[1]:
            lca.append(i[0])
        else:
            break
    if lca:
        return ';'.join(lca)


def main():
    tree = NCBITree()
    print(tree.mp_lineage(1292))
    print(tree.mp_lineage(106))
    print(LCA(tree, 4).apply(1292, 106))
    print(tree.mp_lineage(LCA(tree, 0).apply(1293, 106)))


if __name__ == '__main__':
    main()
