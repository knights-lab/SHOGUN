from collections import defaultdict

import networkx as nx
import numpy as np
import csv


class Taxonomy:
    def __init__(self, filename: str):
        self.tax = self.parse_taxonomy(filename)

    @classmethod
    def parse_taxonomy(cls, filename: str) -> dict:
        with open(filename) as inf:
            csv_inf = csv.reader(inf, delimiter='\t')
            taxa_map = dict(csv_inf)
        return taxa_map

    def __call__(self, id: str):
        return self.tax[id]


class NXTaxonomy:
    def __init__(self,
                 # tree,
                 node_id_to_ancestors,
                 node_id_to_taxa_name,
                 ref_to_node_id,
                 ref_to_taxa_name
                 ):
        # if not nx.is_directed_acyclic_graph(tree):
        #     raise ValueError("Taxonomy file has to be acyclic and directed.")
        # self.tree = tree
        self.node_id_to_ancestors = node_id_to_ancestors
        self.node_id_to_taxa_name = node_id_to_taxa_name
        self.ref_to_node_id = ref_to_node_id
        self.ref_to_taxa_name = ref_to_taxa_name

    def lowest_common_ancestor(self, node_id_1, node_id_2):

        # get ancestors of both (intersection)
        #
        # lca = max(nx.ancestors(self.tree, node_id_1) & nx.ancestors(self.tree, node_id_2))
        # return lca
        lca = max(self.node_id_to_ancestors[node_id_1] & self.node_id_to_ancestors[node_id_2])
        return lca



TAX_LEVELS = ['k', 'p', 'c', 'o', 'f', 'g', 's', 't']


def tree(): return defaultdict(tree)


def add_tree(t, path):
  for node in path.split(';'):
    t = t[node]


# def build_tree_from_tax_file(filename: str) -> NXTaxonomy:
#     with open(filename) as inf:
#         csv_inf = csv.reader(inf, delimiter='\t')
#         ref_to_taxa_name = dict(csv_inf)
#     tree = nx.DiGraph()
#     tree.add_node(0, name="")
#     taxa_name_to_node_id = {"": 0}
#     current_node_id = 1
#     for ref, taxa_name in ref_to_taxa_name.items():
#         split = taxa_name.split(";")
#         for ix in range(len(split)):
#             taxa_name = ";".join(split[:ix+1])
#             if ix == 0:
#                 parent_name = "root"
#                 parent_node_id = -1
#             else:
#                 parent_name = ";".join(split[:ix])
#                 parent_node_id = taxa_name_to_node_id[parent_name]
#             if taxa_name in taxa_name_to_node_id:
#                 continue
#             else:
#                 taxa_name_to_node_id[taxa_name] = current_node_id
#                 tree.add_node(current_node_id, name=taxa_name)
#                 tree.add_edge(parent_node_id, current_node_id)
#                 current_node_id += 1
#
#     ref_to_node_id = {ref: taxa_name_to_node_id[taxa_name] for ref, taxa_name in ref_to_taxa_name.items()}
#     node_id_to_taxa_name = {node_id: taxa_name for taxa_name, node_id in taxa_name_to_node_id.items()}
#
#     return NXTaxonomy(
#         tree=tree,
#         node_id_to_taxa_name=node_id_to_taxa_name,
#         ref_to_node_id=ref_to_node_id,
#         ref_to_taxa_name=ref_to_taxa_name
#     )


def build_tree_from_tax_file(filename: str) -> NXTaxonomy:
    with open(filename) as inf:
        csv_inf = csv.reader(inf, delimiter='\t')
        ref_to_taxa_name = dict(csv_inf)
    tree = nx.DiGraph()
    tree.add_node(0, name="")
    taxa_name_to_node_id = {"": 0}
    current_node_id = 1
    node_id_to_ancestors = [set([0])]
    for ref, taxa_name in ref_to_taxa_name.items():
        split = taxa_name.split(";")
        for ix in range(len(split)):
            taxa_name = ";".join(split[:ix+1])
            if ix == 0:
                parent_name = "root"
                parent_node_id = 0
            else:
                parent_name = ";".join(split[:ix])
                parent_node_id = taxa_name_to_node_id[parent_name]
            if taxa_name in taxa_name_to_node_id:
                continue
            else:
                taxa_name_to_node_id[taxa_name] = current_node_id
                new_set = node_id_to_ancestors[parent_node_id].copy()
                new_set.add(parent_node_id)
                node_id_to_ancestors.append(new_set)
                # tree.add_node(current_node_id, name=taxa_name)
                # tree.add_edge(parent_node_id, current_node_id)
                current_node_id += 1

    ref_to_node_id = {ref: taxa_name_to_node_id[taxa_name] for ref, taxa_name in ref_to_taxa_name.items()}
    node_id_to_taxa_name = {node_id: taxa_name for taxa_name, node_id in taxa_name_to_node_id.items()}

    return NXTaxonomy(
        # tree=tree,
        node_id_to_ancestors=node_id_to_ancestors,
        node_id_to_taxa_name=node_id_to_taxa_name,
        ref_to_node_id=ref_to_node_id,
        ref_to_taxa_name=ref_to_taxa_name
    )
