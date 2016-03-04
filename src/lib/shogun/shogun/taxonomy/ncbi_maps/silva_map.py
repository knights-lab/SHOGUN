#!/usr/bin/env python
import csv
from glob import glob
import os

from shogun.utilities.pickle_class import PickleClass
from shogun.utilities.collections import reverse_dict
from shogun import SETTINGS


class SilvaMap(PickleClass):
    def __init__(self, _silva_taxdmp_dir=SETTINGS.silva_taxdmp_dir):
        """

        :param ncbi_tree: NCBITree()
        :param _silva_taxdmp_dir:
        """
        super().__init__()

        self._silva_taxdmp_dir = _silva_taxdmp_dir

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
                    silva_acc2taxon_id[silva_acc] = ncbi_taxon_id
        return silva_acc2taxon_id, reverse_dict(silva_acc2taxon_id)

    def get(self, silva_acc):
        if silva_acc in self.silva_acc2taxon_id:
            return self.silva_acc2taxon_id[silva_acc]
        else:
            return None

    def __call__(self, silva_acc):
        self.get(silva_acc)

    def __getstate__(self):
        d = dict(self.__dict__)
        del d['taxon_id2silva_acc']
        return d

    def __setstate__(self, d):
        # TODO add try/except
        d['taxon_id2silva_acc'] = reverse_dict(d['silva_acc2taxon_id'])
        self.__dict__.update(d)


def main():
    silva_map = SilvaMap()
    silva_map.save()
    # silva_map = SilvaMap.load()
    # print(silva_map.get_lineage('AJGD01000050'))

if __name__ == '__main__':
    main()
