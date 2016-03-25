#!/usr/bin/env python
import csv
from glob import glob
import os

from shogun.utilities.pickleable import Pickleable
from shogun.utilities.collections import reverse_dict
from shogun.utilities.downloadable import download
from shogun.downloaders.download_silva2ncbi_taxonomy import SilvaMapping


class SilvaMap(Pickleable):
    def __init__(self, _downloader=SilvaMapping()):
        self._downloader = _downloader
        super().__init__()

    @download
    def _parse(self):
        self.silva_acc2taxon_id = {}

        silva_taxdmp_dir = self._downloader.path
        files = glob(os.path.join(silva_taxdmp_dir, "*.txt"))
        for file in files:
            with open(file) as inf:
                csv_inf = csv.reader(inf, delimiter='\t')
                for row in csv_inf:
                    silva_acc = row[0].replace(" ", "")
                    ncbi_taxon_id = int(row[-1].replace(" ", ""))
                    self.silva_acc2taxon_id[silva_acc] = ncbi_taxon_id
        self.taxon_id2silva_acc = reverse_dict(self.silva_acc2taxon_id)

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
    print(silva_map.get('AJGD01000050'))

if __name__ == '__main__':
    main()
