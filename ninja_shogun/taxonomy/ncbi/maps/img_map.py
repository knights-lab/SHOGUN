#!/usr/bin/env python
import csv
import os

from ninja_shogun.utilities.pickleable import Pickleable
from ninja_shogun.utilities.utils import reverse_dict
from ninja_shogun.utilities.scroll import Scroll, scrolling
from ninja_shogun import SETTINGS


class IMGMap(Pickleable):
    def __init__(self, _scroll=Scroll(os.path.join(SETTINGS.img_taxdmp_dir, '00.taxon.tab.txt'))):
        """

        :param img_taxdmp_dir:
        """
        self._scroll = _scroll
        super().__init__()

    @scrolling
    def _parse(self):
        # init variables
        self.img2taxon_id = {}

        with open(self._scroll.path, 'r') as inf:
            csv_file = csv.reader(inf, delimiter='\t')
            # columns = ['taxon_oid', 'ncbi_taxon_id', 'domain', 'genus', 'species']
            row_1 = next(csv_file)

            # Checking for a header file
            if row_1[-1] == 'Finished':
                ncbi_taxon_id = int(row_1[1])
                self.img2taxon_id[int(row_1[0])] = ncbi_taxon_id

            for row in csv_file:
                try:
                    ncbi_taxon_id = int(row[1])
                    self.img2taxon_id[int(row[0])] = ncbi_taxon_id
                except ValueError as e:
                    continue

        self.taxon_id2img = reverse_dict(self.img2taxon_id)

    def get(self, img_id):
        if img_id in self.img2taxon_id:
            return self.img2taxon_id[img_id]
        else:
            return None

    def __call__(self, img_id):
        return self.get(img_id)

    def __getstate__(self):
        d = dict(self.__dict__)
        del d['taxon_id2img']
        return d

    def __setstate__(self, d):
        # TODO add try/except
        d['taxon_id2img'] = reverse_dict(d['img2taxon_id'])
        self.__dict__.update(d)


def main():
    img_map = IMGMap()
    print(img_map.get(637000263))

if __name__ == '__main__':
    main()
