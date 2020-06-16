#!/usr/bin/env python
# simulate_test_sam_file.py <input_tax> <output_folder>
import csv
from collections import Counter, defaultdict
from functools import reduce
import sys
import os

from shogun.utils.lowest_common_ancestor import LCATaxonomy
from shogun.utils.tree import build_tree_from_tax_file

import pandas as pd
import numpy as np
import networkx as nx

np.random.seed(930525)


def build_nx_from_tax_file(filename: str) -> nx.Graph:
    with open(filename) as inf:
        csv_inf = csv.reader(inf, delimiter='\t')
        ref_to_taxa_name = dict(csv_inf)
    tree = nx.DiGraph()
    tree.add_node(0, name="root")
    taxa_name_to_node_id = {"root": 0}
    current_node_id = 1
    for ref, taxa_name in ref_to_taxa_name.items():
        split = taxa_name.split(";")
        for ix in range(len(split)):
            taxa_name = ";".join(split[:ix+1])
            if ix == 0:
                parent_node_id = 0
            else:
                parent_name = ";".join(split[:ix])
                parent_node_id = taxa_name_to_node_id[parent_name]
            if taxa_name in taxa_name_to_node_id:
                continue
            else:
                taxa_name_to_node_id[taxa_name] = current_node_id
                tree.add_node(current_node_id, name=taxa_name)
                tree.add_edge(parent_node_id, current_node_id)
                current_node_id += 1
    return tree


def gen_hits(tree: LCATaxonomy, nx_tree: nx.Graph, num_reads_to_simulate: int=1_000):
    num_hits = np.random.poisson(lam=1.0, size=num_reads_to_simulate) + 1
    tips = np.array(list(tree.ref_to_taxa_name.keys()))

    for num_hit in num_hits:
        choices = np.random.choice(tips, size=num_hit, replace=False)

        node_ids = [tree.ref_to_node_id_ix_level[c][0] for c in choices]

        if len(node_ids) == 1:
            lca = node_ids[0]
        else:
            lca = reduce(
                lambda x, y: nx.lowest_common_ancestor(nx_tree, x, y), node_ids)
        yield node_ids, lca, choices


def sam_file(outfile_folder: str, tree: LCATaxonomy, nx_tree: nx.Graph):
    lca_map = defaultdict(Counter)
    ix = 0

    with open(os.path.join(outfile_folder, "truth.sam"), "w") as outf:
        for s_id in range(5):
            gen = gen_hits(tree, nx_tree)
            for node_ids, lca, choices in gen:
                ix += 1
                for node_id, choice in zip(node_ids, choices):
                    outf.write("\t".join(
                        [f"S{s_id}_R{ix}", "10", choice, str(lca), str(tree.node_id_to_taxa_name[lca]), str(node_id),
                         str(tree.node_id_to_taxa_name[node_id])]) + "\n")
                if lca != 0:
                    lca_map[f"S{s_id}"].update([tree.node_id_to_taxa_name[lca]])
    df = pd.DataFrame(lca_map)
    df.drop("root", axis=0, errors="ignore", inplace=True)
    df.to_csv(os.path.join(outfile_folder, "truth.txt"), sep='\t', float_format="%d", na_rep=0, index_label="#OTU ID")


def main():
    tax_infile = sys.argv[1]
    outfile_folder = sys.argv[2]
    tree = build_tree_from_tax_file(tax_infile)
    nx_tree = build_nx_from_tax_file(tax_infile)
    sam_file(outfile_folder, tree, nx_tree)


if __name__ == '__main__':
    main()
