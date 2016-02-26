import io
import urllib
import os
import csv
from shogun.taxonomy.ncbi_tree import NCBITree

from shogun import SETTINGS


class SilvaTree:
    def __init__(self, silva_acc2taxon_id_path=os.path.join(SETTINGS.data_path, "silva_acc2taxon_id"),
                 _silva_taxdmp_urls=SETTINGS.silva_taxdmp_urls):
        # self.NCBI_tree = NCBITree.load()
        self._silva_acc2taxon_id_path = silva_acc2taxon_id_path
        self._silva_taxdmp_urls = _silva_taxdmp_urls

        self._silva_acc2taxon_id = {}
        self._download_silva_dump()
        self._parse_silva_taxonomy()

    def _parse_silva_taxonomy(self):
        pass

    def _download_silva_dump(self):
        for url in self._silva_taxdmp_urls:
            with urllib.request.urlopen(url) as stream:
                for line in stream:
                    cols = line.decode().strip().split('\t')
                    self._silva_acc2taxon_id[cols[0]] = cols[-1]

def main():
    ncbi_tree = SilvaTree()
    print("Hello, World")

if __name__ == '__main__':
    main()
