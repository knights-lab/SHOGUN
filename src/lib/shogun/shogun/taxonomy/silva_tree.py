import urllib
import os
import csv

from shogun.utils.pickle_class import PickleClass
from shogun import SETTINGS


class SilvaTree(PickleClass):
    def __init__(self, cache=False,
                 _silva_acc2taxon_id_path=os.path.join(SETTINGS.data_path, 'silva_acc2taxon_id.csv'),
                 _silva_taxdmp_urls=SETTINGS.silva_taxdmp_urls, _pickle_dir=SETTINGS.pickle_dir
                 ):

        self._silva_acc2taxon_id_path = _silva_acc2taxon_id_path
        self._silva_taxdmp_urls = _silva_taxdmp_urls
        self._pickle_dir = _pickle_dir

        self.silva_acc2taxon_id = {}
        self._download_parse_silva_url()
        self._parse_silva_taxonomy()

    def _parse_silva_taxonomy(self):
        pass

    def _download_parse_silva_url(self):
        with open(self._silva_acc2taxon_id_path, 'w') as outfile:
            csv_outfile = csv.writer(outfile, quotechar='"')
            csv_outfile.writerow(['silva_accession', 'ncbi_taxon_id'])
            for url in self._silva_taxdmp_urls:
                with urllib.request.urlopen(url) as stream:
                    for line in stream:
                        cols = line.decode().strip().split('\t')
                        csv_outfile.writerow([cols[0], cols[-1]])
                        self.silva_acc2taxon_id[cols[0]] = cols[-1]

def main():
    silva_tree = SilvaTree()
    print("Hello, World")

if __name__ == '__main__':
    main()
