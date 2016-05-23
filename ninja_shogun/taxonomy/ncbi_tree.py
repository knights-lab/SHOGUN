#!/usr/bin/env python
"""
This library was adapted from:
        https://github.com/luo-chengwei/utilitomics

This library works with the taxonomy db download from NCBI FTP:

ftp://ftp.ncbi.nih.gov:/pub/taxonomy/taxdump.tar.gz

Init an taxonomy tree obj by:

t_tree = ncbi.NCBITree()
then you can get the taxonomy ranks by using:

path = t_tree.get_name_path_with_taxon_id(taxon_id)
"""

import os
import networkx as nx
import csv
import sys
from functools import lru_cache
from collections import defaultdict

from ninja_shogun.utilities.pickleable import Pickleable
from ninja_shogun.utilities.downloadable import download
from ninja_shogun.downloaders.download_ncbi_taxonomy import NCBITaxdmp


class NCBITree(Pickleable):
    def __init__(self, mp_ranks=None, _downloader=NCBITaxdmp()):
        # Private variables (should be set in settings)
        self._downloader = _downloader
        if mp_ranks is None:
            self.mp_ranks = dict(zip(('superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'),
                     ('k__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__')))
        else:
            self.mp_ranks = mp_ranks
        super().__init__()

    @download
    def _parse(self):
        # Initialize variables
        self.tree = nx.DiGraph()
        self.name2taxon_id = defaultdict(int)
        self.taxon_id2name = defaultdict(str)
        ncbi_taxdmp_dir = self._downloader.path
        csv.field_size_limit(sys.maxsize)

        with open(os.path.join(ncbi_taxdmp_dir, 'names.dmp'), 'r') as handle:
            csv_handle = csv.reader(handle, delimiter="\t")
            for cols in csv_handle:
                taxon_id = int(cols[0])
                name = cols[2]
                self.name2taxon_id[name] = taxon_id
                if cols[-2] == 'scientific name':
                    self.taxon_id2name[taxon_id] = name

        # construct node tree
        edges = []
        nodes = {}
        with open(os.path.join(ncbi_taxdmp_dir, 'nodes.dmp'), 'r') as handle:
            csv_handle = csv.reader(handle, delimiter="\t")
            for cols in csv_handle:
                parent_node = int(cols[2])
                child_node = int(cols[0])
                rank = cols[4]
                nodes[child_node] = rank
                if child_node != parent_node:
                    edges.append((child_node, parent_node))

        self.tree.add_edges_from(edges)
        nx.set_node_attributes(self.tree, 'rank', nodes)

    def get_taxon_id_lineage_with_taxon_id(self, taxon_id):
        if taxon_id in self.tree:
            for i in nx.dfs_preorder_nodes(self.tree, taxon_id):
                yield i

    def get_taxon_id_lineage_with_name(self, name):
        if name not in self.name2taxon_id:
            return []
        return self.get_taxon_id_lineage_with_taxon_id(self.name2taxon_id[name])

    def get_name_lineage_with_taxon_id(self, taxon_id):
        tid_lineage = self.get_taxon_id_lineage_with_taxon_id(taxon_id)
        name_lineage = []
        for x in tid_lineage:
            rank = self.tree.node[x]['rank']
            name = self.taxon_id2name[x]
            name_lineage.append((name, rank))
        return name_lineage

    def get_name_lineage_with_name(self, name):
        if name not in self.name2taxon_id:
            return []
        path = self.get_name_lineage_with_taxon_id(self.name2taxon_id[name])
        return path

    def get_rank_with_taxon_id(self, taxon_id, rank):
        taxon_id_lineage = self.get_taxon_id_lineage_with_taxon_id(taxon_id)
        name_lineage = self.get_name_lineage_with_taxon_id(taxon_id)
        for tid, (name, x) in zip(taxon_id_lineage, name_lineage):
            if x == rank:
                return tid, name
        return None, None

    def get_name_at_rank(self, name, rank):
        if name not in self.name2taxon_id:
            return None, None
        return self.get_rank_with_taxon_id(self.name2taxon_id[name], rank)

    def get_lineage(self, taxon_id, ranks={'superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'}):
        taxon_id_lineage = self.get_taxon_id_lineage_with_taxon_id(taxon_id)
        name_lineage = []
        for x in taxon_id_lineage:
            rank = self.tree.node[x]['rank']
            name = self.taxon_id2name[x]
            if rank in ranks:
                name_lineage.append((name, x, rank))
        return name_lineage

    def get_lineage_depth(self, taxon_id, depth, ranks=('superkingdom', 'phylum', 'class', 'order', 'family', 'genus',
                                                        'species')):
        taxon_id_lineage = self.get_taxon_id_lineage_with_taxon_id(taxon_id)
        ranks = set(ranks[depth:])
        lineage = []
        for x in taxon_id_lineage:
            rank = self.tree.node[x]['rank']
            if rank in ranks:
                lineage.append(x)
        return lineage

    # Note that this will create a global cache for all instances of NCBITree.
    # Will be fine unless if you want to compare trees.
    @lru_cache(maxsize=128)
    def mp_lineage(self, taxon_id):
        taxon_id_lineage = self.get_taxon_id_lineage_with_taxon_id(taxon_id)
        name_lineage = []
        for x in taxon_id_lineage:
            rank = self.tree.node[x]['rank']
            name = self.taxon_id2name[x]
            if rank in self.mp_ranks:
                name_lineage.append(self.mp_ranks[rank] + name.replace(' ', '_'))
        return ';'.join(reversed(name_lineage))


def main():
    ncbi_tree = NCBITree()
    print(ncbi_tree.get_taxon_id_lineage_with_taxon_id(640701851))

if __name__ == '__main__':
    main()
