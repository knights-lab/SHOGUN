"""
Copyright 2015-2017 Knights Lab, Regents of the University of Minnesota.

This software is released under the GNU Affero General Public License (AGPL) v3.0 License.
"""
import csv

class Taxonomy:
    def __init__(self, filename):
        self.tax = self.parse_taxonomy(filename)

    @classmethod
    def parse_taxonomy(cls, filename):
        with open(filename) as inf:
            csv_inf = csv.reader(inf, delimiter='\t')
            return dict(csv_inf)

    def __call__(self, id: str):
        return self.tax[id]

__all__ = ["Taxonomy"]
