import urllib.request
import os
from shogun.taxonomy.ncbi_tree import NCBITree

from shogun import SETTINGS

class SilvaTree(NCBITree):
    def __init__(self, cache=True, _ncbi_taxdmp_url=SETTINGS.ncbi_taxdmp_url, _ncbi_taxdmp_path=os.path.join(SETTINGS.data_path, "ncbi_taxdmp"), _cache_path=SETTINGS.cache_path):
        super(NCBITree, self).__init__()
        self.silva_acc2taxon_id = {}
        self._parse_silva_taxonomy()

    def _parse_silva_taxonomy(self):
        pass

def main():
    ncbi_tree = NCBITree.load()
    print("Hello, World")

if __name__ == '__main__':
    main()
