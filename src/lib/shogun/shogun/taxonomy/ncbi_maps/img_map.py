#!/usr/bin/env python
import csv
from glob import glob
import os

from shogun.utilities.pickle_class import PickleClass
from shogun.utilities.collections import reverse_dict
from shogun import SETTINGS


class IMGMap(PickleClass):
    def __init__(self, _img_taxdmp_dir=SETTINGS.img_taxdmp_dir):
        """

        :param img_taxdmp_dir:
        """
        super().__init__()

        self._img_taxdmp_dir = _img_taxdmp_dir

        self.img2taxon_id, self.taxon_id2img = self._parse_silva_taxonomy_file()

    def _parse_silva_taxonomy_file(self):
        img2taxon_id = {}
        files = glob(os.path.join(self._img_taxdmp_dir, "*.txt"))
        for file in files:
            with open(file, 'r') as inf:
                csv_file = csv.reader(inf, delimiter='\t')
                # columns = ['taxon_oid', 'ncbi_taxon_id', 'domain', 'genus', 'species']
                row_1 = next(csv_file)

                # Checking for a header file
                if row_1[-1] == 'Finished':
                    ncbi_taxon_id = int(row_1[1])
                    img2taxon_id[int(row_1[0])] = ncbi_taxon_id

                for row in csv_file:
                    try:
                        ncbi_taxon_id = int(row[1])
                        img2taxon_id[int(row[0])] = ncbi_taxon_id
                    except ValueError as e:
                        continue

        return img2taxon_id, reverse_dict(img2taxon_id)

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
    img_map.save()
    # img_map = IMGMap.load()
    # print(img_map.get_lineage(637000263))

if __name__ == '__main__':
    main()
