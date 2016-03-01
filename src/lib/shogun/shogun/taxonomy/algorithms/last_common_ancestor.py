#!/usr/bin/env python

from shogun.taxonomy.ncbi_tree import NCBITree


class LCA:
    def __init__(self, tree=NCBITree.load()):
        self.tree_obj = tree

    def __call__(self):
        pass

    def apply(self, taxon_id_a, taxon_id_b):
        current_node = taxon_id_a
        if current_node not in self.tree_obj.tree.nodes():
            return []
        while len(self.tree_obj.tree.successors(current_node)) > 0:
            current_node = self.tree_obj.tree.successors(current_node)


def main():
    LCA().apply(1292, 1290)

if __name__ == '__main__':
    main()
