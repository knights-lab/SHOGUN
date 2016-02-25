#!/usr/bin/env python

"""
This library was adapted from:
        https://github.com/luo-chengwei/utilitomics

This library works with the taxonomy db download from NCBI FTP:

ftp://ftp.ncbi.nih.gov:/pub/taxonomy/taxdump.tar.gz

basically you can init an taxonomy tree obj by:

t_tree = ncbi.NCBITree(path)   # tree being the tax dump db unzip dir
then you can basically get the taxonomy ranks using:

path = t_tree.get_name_path_with_taxon_id(taxon_id)

"""

import os
import networkx as nx
import pickle
import gzip
from shogun import SETTINGS


class NCBITree:
    def __init__(self, path=os.path.join(SETTINGS.data_dir, 'taxdmp')):
        self.tree = nx.DiGraph()

        # construct name -> taxon_id mapping
        self.name2taxon_id = {}
        self.taxon_id2name = {}
        with open(os.path.join(path, 'names.dmp'), 'rb') as handle:
            for line in handle:
                cols = line.split(b'\t')
                taxon_id = cols[0]
                name = cols[2]
                self.name2taxon_id[name] = taxon_id
                if cols[-2] == 'scientific name':
                    self.taxon_id2name[taxon_id] = name

        # construct node tree
        edges = []
        nodes = {}
        with open(os.path.join(path, 'nodes.dmp'), 'rb') as handle:
            for line in handle:
                cols = line.split(b'\t')
                parent_node = cols[2]
                children_node = cols[0]
                rank = cols[4]
                nodes[children_node] = rank
                if children_node != parent_node:
                    edges.append((children_node, parent_node))

        self.tree.add_edges_from(edges)
        nx.set_node_attributes(self.tree, 'rank', nodes)

    def save(self, _cache_dir=SETTINGS.cache_dir):
        self_dump = os.path.join(_cache_dir, "ncbi_taxon_tree.pklz")
        with gzip.open(self_dump, 'wb') as handle:
            pickle.dump(self, handle)

    @classmethod
    def load(cls, _cache_dir=SETTINGS.cache_dir):
        try:
            self_dump = os.path.join(_cache_dir, "ncbi_taxon_tree.pklz")
            with gzip.open(self_dump, 'rb') as handle:
                self = pickle.load(handle)
        except FileNotFoundError as error:
            raise

    def get_taxon_id_lineage_with_taxon_id(self, taxonID):
        path = [taxonID]
        current_node = taxonID
        if current_node not in self.tree.nodes():
            return []
        while len(self.tree.successors(current_node)) > 0:
            path.append(self.tree.successors(current_node)[0])
            current_node = self.tree.successors(current_node)[0]
        return path

    def get_taxon_id_lineage_with_name(self, name):
        if name not in self.name2taxon_id:
            return []
        path = self.get_taxon_id_lineage_with_taxon_id(self.name2taxon_id[name])
        return path

    def get_name_lineage_with_taxon_id(self, taxonID):
        taxon_id_lineage = self.get_taxon_id_lineage_with_taxon_id(taxonID)
        name_lineage = []
        for x in taxon_id_lineage:
            rank = self.tree.node[x]['rank']
            try:
                name = self.taxon_id2name[x]
            except KeyError:
                name = ''
            name_lineage.append((name, rank))
        return name_lineage

    def get_name_path_with_name(self, name):
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

    def get_rank_with_name(self, name, rank):
        if name not in self.name2taxon_id:
            return None, None
        return self.get_rank_with_taxon_id(self.name2taxon_id[name], rank)


def main():
    pass

if __name__ == '__main__':
    main()
