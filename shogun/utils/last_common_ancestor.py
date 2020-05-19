"""
Copyright 2015-2017 Knights Lab, Regents of the University of Minnesota.

This software is released under the GNU Affero General Public License (AGPL) v3.0 License.
"""

import typing

from shogun.redistribute import Taxonomy


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
