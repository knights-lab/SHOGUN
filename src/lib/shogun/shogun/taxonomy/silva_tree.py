#!/usr/bin/env python
import csv
from glob import glob
import os

from shogun.utilities.pickle_class import PickleClass
from shogun.utilities.collections import reverse_dict
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

        self.silva_acc2taxon_id, self.taxon_id2silva_acc = self._parse_silva_taxonomy_file()

    def _parse_silva_taxonomy_file(self):
        silva_acc2taxon_id = {}
        files = glob(os.path.join(self._silva_taxdmp_dir, "*.txt"))
        for file in files:
            with open(file) as inf:
                csv_inf = csv.reader(inf, delimiter='\t')
                for row in csv_inf:
                    silva_acc = row[0].replace(" ", "")
                    ncbi_taxon_id = int(row[-1].replace(" ", ""))
                    if ncbi_taxon_id in self.ncbi_tree.taxon_id2name:
                        silva_acc2taxon_id[silva_acc] = ncbi_taxon_id
        return silva_acc2taxon_id, reverse_dict(silva_acc2taxon_id)

    def get_lineage(self, silva_acc, **kwargs):
        lineage = self.ncbi_tree.get_lineage(self.silva_acc2taxon_id[silva_acc], **kwargs)
        return lineage

    def __getstate__(self):
        d = dict(self.__dict__)
        del d['ncbi_tree']
        del d['taxon_id2silva_acc']
        return d

    def __setstate__(self, d):
        # TODO add try/except
        d['ncbi_tree'] = NCBITree.load()
        d['taxon_id2silva_acc'] = reverse_dict(d['silva_acc2taxon_id'])
        self.__dict__.update(d)


def main():
    # silva_tree = SilvaTree()
    # silva_tree.save()
    silva_tree = SilvaTree.load()
    print(silva_tree.get_lineage('AJGD01000050'))

if __name__ == '__main__':
    main()
