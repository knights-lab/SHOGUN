"""
Copyright 2015-2017 Knights Lab, Regents of the University of Minnesota.

This software is released under the GNU Affero General Public License (AGPL) v3.0 License.
"""

from shogun.parsers import yield_alignments_from_sam_inf


def build_lca_map(sam_inf, extract_ncbi_tid, tree):
    """
    Build a last common ancestor dictionary
    :param sam_inf: path to SAM infile
    :param extract_ncbi_tid: function to extract ncbi_tid
    :param tree: NCBITree
    :return: dict key (query name: string) value (ncbi_tid: int)
    """
    lca_map = {}
    for qname, rname in yield_alignments_from_sam_inf(sam_inf):
        ncbi_tid = extract_ncbi_tid(rname)
        if qname in lca_map:
            current_ncbi_tid = lca_map[qname]
            if current_ncbi_tid:
                if current_ncbi_tid != ncbi_tid:
                    lca_map[qname] = tree.lowest_common_ancestor(ncbi_tid, current_ncbi_tid)
        else:
            lca_map[qname] = ncbi_tid
    return lca_map
