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
import tarfile
import urllib.request
import csv

from shogun import SETTINGS
from shogun.utils.path import verify_make_path


class NCBITree:
    def __init__(self, cache=True, _ncbi_taxdmp_url=SETTINGS.ncbi_taxdmp_url,
                 _ncbi_taxdmp_path=os.path.join(SETTINGS.data_path, "ncbi_taxdmp"),
                 _pickle_path=os.path.join(SETTINGS.data_path, 'pickle')):
                self.tree = nx.DiGraph()
                # construct name -> taxon_id mapping
                self.name2taxon_id = {}
                self.taxon_id2name = {}
                self._pickle_path = _pickle_path
                verify_make_path(self._pickle_path)
                self._ncbi_taxdmp_url = _ncbi_taxdmp_url
                self._ncbi_taxdmp_path = _ncbi_taxdmp_path
                verify_make_path(self._ncbi_taxdmp_path)

                self._parse_taxonomy()

                if cache:
                    self.save()

    def _parse_taxonomy(self):
        with open(os.path.join(self._ncbi_taxdmp_path, 'names.dmp'), 'r') as handle:
            csv_handle = csv.reader(handle, delimiter="\t")
            for cols in csv_handle:
                taxon_id = cols[0]
                name = cols[2]
                self.name2taxon_id[name] = taxon_id
                if cols[-2] == 'scientific name':
                    self.taxon_id2name[taxon_id] = name

        # construct node tree
        edges = []
        nodes = {}
        with open(os.path.join(self._ncbi_taxdmp_path, 'nodes.dmp'), 'r') as handle:
            csv_handle = csv.reader(handle, delimiter="\t")
            for cols in csv_handle:
                parent_node = cols[2]
                children_node = cols[0]
                rank = cols[4]
                nodes[children_node] = rank
                if children_node != parent_node:
                    edges.append((children_node, parent_node))

        self.tree.add_edges_from(edges)
        nx.set_node_attributes(self.tree, 'rank', nodes)

    def save(self):
        self_dump = os.path.join(self._pickle_path, "ncbi_taxon_tree.pkl")
        with open(self_dump, 'wb') as handle:
            pickle.dump(self, handle)

    @classmethod
    def load(cls, _pickle_path=os.path.join(SETTINGS.data_path, 'pickle')):
        try:
            self_dump = os.path.join(_pickle_path, "ncbi_taxon_tree.pkl")
            with open(self_dump, 'rb') as handle:
                return pickle.load(handle)
        except FileNotFoundError as error:
            raise

    def get_taxon_id_lineage_with_taxon_id(self, taxon_id):
        path = [taxon_id]
        current_node = taxon_id
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

    def _download_ncbi_dump(self):
        req = urllib.request.Request(self._ncbi_taxdmp_url)
        with urllib.request.urlopen(req, 'rb') as ftp_stream:
            tfile = tarfile.open(fileobj=ftp_stream, mode='r|gz')
            tfile.extractall(self._ncbi_taxdmp_path)


def main():
    ncbi_tree = NCBITree.load()

if __name__ == '__main__':
    main()
