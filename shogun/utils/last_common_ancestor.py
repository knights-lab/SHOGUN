"""
Copyright 2015-2017 Knights Lab, Regents of the University of Minnesota.

This software is released under the GNU Affero General Public License (AGPL) v3.0 License.
"""

import typing

from shogun.redistribute import Taxonomy
from shogun.utils.tree import LCATaxonomy

import numpy as np


def build_lca_map(gen: typing.Iterator, tree: Taxonomy) -> dict:
    """
    Build a last common ancestor dictionary
    :param sam_inf: path to SAM infile
    :param extract_ncbi_tid: function to extract ncbi_tid
    :param tree: NCBITree
    :return: dict key (query name: string) value (ncbi_tid: int)
    """
    lca_map = {}
    for record in gen:
        taxs = (tree(rname) for qname, rname in record)
        # if qname in lca_map:
        #     current_tax = lca_map[qname]
        #     if current_tax:
        #         if current_tax != tax:
        #             lca_map[qname] = least_common_ancestor((tax, current_tax))
        # else:
        #     lca_map[qname] = tax
    return lca_map


def least_common_ancestor(taxa_set):
    lca = []
    taxa_set = [_.split(';') for _ in taxa_set]
    unclassified_flag = False
    lca_classified = []
    for level in zip(*taxa_set):
        if len(set(level)) > 1:
            if unclassified_flag:
                lca = lca_classified
            return ';'.join(lca) if lca else None
        elif len(level[0]) == 3 and level[0].endswith("__"):
            # hit an unclassified level
            if not unclassified_flag:
                lca_classified = lca.copy()
                unclassified_flag = True
        else:
            # reset classified flag
            unclassified_flag = False
        lca.append(level[0])


def build_lowest_common_ancestor_map(gen: typing.Iterator, tree: LCATaxonomy, confidence_threshold=1):
    lca_map = {}
    for ix, record in enumerate(gen):
        l_node_ids_ixs_levels = [tree.ref_to_node_id_ix_level[read[1]] for read in record]
        # fetch the ancestors
        ancestors = [tree.ix_to_ancestors[ix][:level] for node_id, ix, level in l_node_ids_ixs_levels]
        unique_elements, counts_elements = np.unique(np.vstack(ancestors), return_counts=True)
        counts_elements = counts_elements / len(record)
        mask = counts_elements >= confidence_threshold
        yield record[0][0], unique_elements[mask].max()
