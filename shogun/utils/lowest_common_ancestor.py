"""
Copyright 2015-2020 Knights Lab, Regents of the University of Minnesota.

This software is released under the GNU Affero General Public License (AGPL) v3.0 License.
"""

import typing
from functools import reduce

from shogun.parsers import yield_alignments_from_sam_inf
from shogun.utils.tree import LCATaxonomy

import pandas as pd
import numpy as np


def build_lca_df(sam_file: str, tree: LCATaxonomy, confidence_threshold: float = 1.0, samples_iter: int = 50) -> pd.DataFrame:
    align_gen = yield_alignments_from_sam_inf(sam_file)
    if confidence_threshold == 1.0:
        lca_map_gen = gen_lowest_common_ancestor(align_gen, tree)
    else:
        lca_map_gen = gen_confidence_lowest_common_ancestor(align_gen, tree, confidence_threshold)

    sample_names_to_ix = dict()
    ix = 0
    mat_counts = np.zeros((tree.num_nodes, samples_iter), dtype=int)
    max_samples = samples_iter
    for rname, node_id in lca_map_gen:
        sample_name = rname.split('_')[0]
        if sample_name in sample_names_to_ix:
            c_ix = sample_names_to_ix[sample_name]
            mat_counts[node_id, c_ix] += 1
        else:
            if ix >= max_samples:
                b = np.zeros((tree.num_nodes, max_samples + samples_iter))
                b[:, :-samples_iter] = mat_counts
                mat_counts = b
                max_samples += samples_iter
            sample_names_to_ix[sample_name] = ix
            mat_counts[node_id, ix] += 1
            ix += 1

    sample_names = [k for k, v in sorted(sample_names_to_ix.items(), key=lambda item: item[1])]

    df = pd.DataFrame(mat_counts[:, :ix], dtype=int, columns=sample_names)
    # drop all node ids of all zeros
    df = df.loc[~(df == 0).all(axis=1)].copy()
    df.index = [tree.node_id_to_taxa_name[node_id] for node_id in df.index]
    return df


def gen_lowest_common_ancestor(gen: typing.Iterator, tree: LCATaxonomy):
    for ix, record in enumerate(gen):
        num_alignments = len(record)
        if num_alignments > 1:
            l_node_ids_ixs_levels = [tree.ref_to_node_id_ix_level[read[1]] for read in record]
            node_id = max(reduce(lambda x, y: x.intersection(y), (tree.node_id_to_ancestors[node_id] for node_id, ix, level in l_node_ids_ixs_levels)))
            yield record[0][0], node_id
        else:
            node_id, ix, level = tree.ref_to_node_id_ix_level[record[0][1]]
            yield record[0][0], node_id


def gen_confidence_lowest_common_ancestor(gen: typing.Iterator, tree: LCATaxonomy, confidence_threshold: float):
    for ix, record in enumerate(gen):
        num_alignments = len(record)
        if num_alignments > 1:
            l_node_ids_ixs_levels = [tree.ref_to_node_id_ix_level[read[1]] for read in record]
            # fetch the ancestors
            ancestors = [tree.ix_to_ancestors[ix][:level] for node_id, ix, level in l_node_ids_ixs_levels]
            unique_elements, counts_elements = np.unique(np.concatenate(ancestors), return_counts=True)
            counts_elements = counts_elements / num_alignments
            mask = counts_elements >= confidence_threshold
            yield record[0][0], unique_elements[mask].max()
        else:
            node_id, ix, level = tree.ref_to_node_id_ix_level[record[0][1]]
            yield record[0][0], node_id
