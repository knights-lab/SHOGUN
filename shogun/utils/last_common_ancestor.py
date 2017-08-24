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
    for qname, rname in gen:
        tax = tree(rname)
        if qname in lca_map:
            current_tax = lca_map[qname]
            if current_tax:
                if current_tax != tax:
                    lca_map[qname] = least_common_ancestor((tax, current_tax))
        else:
            lca_map[qname] = tax
    return lca_map


def least_common_ancestor(taxa_set):
    taxa_set = [_.split(';') for _ in taxa_set]
    for i, level in enumerate(reversed(list(zip(*taxa_set)))):
        if len(set(level)) == 1:
            convert_i = lambda x: None if x == 0 else -x
            return ';'.join([taxon_level[0] for taxon_level in list(zip(*taxa_set))[:convert_i(i)]])
    return None
