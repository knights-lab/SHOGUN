import urllib
import os
import pickle
import csv
from shogun.taxonomy.ncbi_tree import NCBITree
from shogun.utils.path import verify_make_path

from shogun import SETTINGS


class SilvaTree:
    def __init__(self, cache=True,
                 _silva_acc2taxon_id_path=os.path.join(SETTINGS.data_path, 'silva_acc2taxon_id.csv'),
                 _silva_taxdmp_urls=SETTINGS.silva_taxdmp_urls, _pickle_path=os.path.join(SETTINGS.data_path, 'pickle')
                 ):
        # self.NCBI_tree = NCBITree.load()
        self._silva_acc2taxon_id_path = _silva_acc2taxon_id_path
        self._silva_taxdmp_urls = _silva_taxdmp_urls
        self._pickle_path = _pickle_path
        # verify the pickle path
        verify_make_path(self._pickle_path)

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

    def save(self):
        self_dump = os.path.join(self._pickle_path, "_silva_acc2taxon_id_path.pkl")
        with open(self_dump, 'wb') as handle:
            pickle.dump(self, handle)

    @classmethod
    def load(cls, _pickle_path=os.path.join(SETTINGS.data_path, 'pickle')):
        try:
            self_dump = os.path.join(_pickle_path, "_silva_acc2taxon_id_path.pkl")
            with open(self_dump, 'rb') as handle:
                return pickle.load(handle)
        except FileNotFoundError as error:
            raise

def main():
    silva_tree = SilvaTree()
    print("Hello, World")

if __name__ == '__main__':
    main()
