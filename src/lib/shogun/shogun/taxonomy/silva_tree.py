#!/usr/bin/env python
import csv
from glob import glob
import os

from shogun.utils.pickle_class import PickleClass
from shogun.taxonomy.ncbi_tree import NCBITree
from shogun import SETTINGS


class SilvaTree(PickleClass):
    def __init__(self,
                 ncbi_tree=NCBITree.load(),
                 _silva_taxdmp_dir=SETTINGS.silva_taxdmp_dir,
                 _silva_taxdmp_urls=SETTINGS.silva_taxdmp_urls):
        """

        :param ncbi_tree: NCBITree()
        :param _silva_taxdmp_dir:
        :param _silva_taxdmp_urls:
        """
        super().__init__()

        self.ncbi_tree = ncbi_tree

        self._silva_taxdmp_dir = _silva_taxdmp_dir
        self._silva_taxdmp_urls = _silva_taxdmp_urls

        self.silva_acc2taxon_id = self._parse_silva_taxonomy_file()

    def _parse_silva_taxonomy_file(self):
        silva_acc2taxon_id = {}
        files = glob(os.path.join(self._silva_taxdmp_dir, "*.txt"))
        for file in files:
            with open(file) as inf:
                csv_inf = csv.reader(inf, delimiter='\t')
                for row in csv_inf:
                    silva_acc = row[0].replace(" ", "")
                    ncbi_taxon_id = row[-1].replace(" ", "")
                    if ncbi_taxon_id in self.ncbi_tree.taxon_id2name:
                        silva_acc2taxon_id[silva_acc] = ncbi_taxon_id
        return silva_acc2taxon_id


def main():
    silva_tree = SilvaTree()
    silva_tree.save()

if __name__ == '__main__':
    main()
