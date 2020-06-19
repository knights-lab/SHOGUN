from collections import defaultdict

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


class LCATaxonomy:
    def __init__(self,
                 node_id_to_taxa_name: dict,
                 ref_to_node_id_ix_level: dict,
                 ref_to_taxa_name: dict,
                 node_id_to_ancestors
                 ):
        self.node_id_to_taxa_name = node_id_to_taxa_name
        self.ref_to_node_id_ix_level = ref_to_node_id_ix_level
        self.ref_to_taxa_name = ref_to_taxa_name
        self.num_nodes = len(self.node_id_to_taxa_name)
        self.node_id_to_ancestors = node_id_to_ancestors

TAX_LEVELS = ['k', 'p', 'c', 'o', 'f', 'g', 's', 't']


def tree(): return defaultdict(tree)


def add_tree(t, path):
  for node in path.split(';'):
    t = t[node]


def build_tree_from_tax_file(filename: str) -> LCATaxonomy:
    with open(filename) as inf:
        csv_inf = csv.reader(inf, delimiter='\t')
        ref_to_taxa_name = dict(csv_inf)
    taxa_name_to_node_id_ix_level = {"root": (0, 0, 0)}
    current_node_id = 1
    node_id_to_ancestors = [{0}]
    for ix, (ref, taxa_name) in enumerate(ref_to_taxa_name.items()):
        split = taxa_name.split(";")
        ancestors = [0]
        for level in range(len(split)):
            taxa_name = ";".join(split[:level+1])
            if taxa_name in taxa_name_to_node_id_ix_level:
                found_node_id, _, _ = taxa_name_to_node_id_ix_level[taxa_name]
                # Check if blank level
                if len(split[level]) > 3:
                    ancestors.append(found_node_id)
            else:
                taxa_name_to_node_id_ix_level[taxa_name] = (current_node_id, ix, level + 1)
                # Check if blank level
                if len(split[level]) > 3:
                    ancestors.append(current_node_id)
                current_node_id += 1
                node_id_to_ancestors.append(set(ancestors))

    ref_to_node_id_ix_level = {ref: taxa_name_to_node_id_ix_level[taxa_name] for ref, taxa_name in ref_to_taxa_name.items()}
    node_id_to_taxa_name = {node_id: taxa_name for taxa_name, (node_id, ix, level) in taxa_name_to_node_id_ix_level.items()}

    return LCATaxonomy(
        node_id_to_taxa_name=node_id_to_taxa_name,
        ref_to_node_id_ix_level=ref_to_node_id_ix_level,
        ref_to_taxa_name=ref_to_taxa_name,
        node_id_to_ancestors=node_id_to_ancestors
    )
