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
import csv

from shogun.utilities.pickle_class import PickleClass
from shogun import SETTINGS


class NCBITree(PickleClass):
    def __init__(self, _ncbi_taxdmp_dir=SETTINGS.ncbi_taxdmp_dir):
        super().__init__()
        self.tree = nx.DiGraph()
        # construct name -> taxon_id mapping
        self.name2taxon_id = {}
        self.taxon_id2name = {}

        # Private variables (should be set in settings)
        self._ncbi_taxdmp_dir = _ncbi_taxdmp_dir

        self._parse_ncbi_taxonomy()

    def _parse_ncbi_taxonomy(self):
        with open(os.path.join(self._ncbi_taxdmp_dir, 'names.dmp'), 'r') as handle:
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
        with open(os.path.join(self._ncbi_taxdmp_dir, 'nodes.dmp'), 'r') as handle:
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
        try:
            path = [taxon_id]
            current_node = taxon_id
            for node in self.tree.successors_iter(current_node):
                path.append(node)
            return path
        except nx.exception.NetworkXError:
            return []

    def get_taxon_id_lineage_with_name(self, name):
        if name not in self.name2taxon_id:
            return []
        path = self.get_taxon_id_lineage_with_taxon_id(self.name2taxon_id[name])
        return path

    def get_name_lineage_with_taxon_id(self, taxon_id):
        taxon_id_lineage = self.get_taxon_id_lineage_with_taxon_id(taxon_id)
        name_lineage = []
        for x in taxon_id_lineage:
            rank = self.tree.node[x]['rank']
            try:
                name = self.taxon_id2name[x]
            except KeyError:
                name = ''
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
            try:
                name = self.taxon_id2name[x]
            except KeyError:
                name = ''
            if rank in ranks:
                name_lineage.append((name, x, rank))
        return name_lineage

    def get_lineage_depth(self, taxon_id, depth, ranks=['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']):
        taxon_id_lineage = self.get_taxon_id_lineage_with_taxon_id(taxon_id)
        ranks = set(ranks[depth:])
        lineage = []
        for x in taxon_id_lineage:
            rank = self.tree.node[x]['rank']
            if rank in ranks:
                lineage.append(x)
        return lineage

    def mp_lineage(self, taxon_id, ranks={'superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'}
                   , nodes=('k__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__')):
        taxon_id_lineage = self.get_taxon_id_lineage_with_taxon_id(taxon_id)
        name_lineage = []
        for x in taxon_id_lineage:
            rank = self.tree.node[x]['rank']
            try:
                name = self.taxon_id2name[x]
            except KeyError:
                name = ''
            if rank in ranks:
                name_lineage.append(name)
        return mp_format(reversed(name_lineage), nodes=nodes)


def mp_format(taxa_array, nodes=('k__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__')):
    return ';'.join([i + j.replace(' ', '_') for i, j in zip(nodes, taxa_array)])


def main():
    ncbi_tree = NCBITree()
    ncbi_tree.save()
    # ncbi_tree = NCBITree.load()
    # print(ncbi_tree.get_taxon_id_lineage_with_taxon_id(640701851))

if __name__ == '__main__':
    main()
