"""
Copyright 2015-2017 Knights Lab, Regents of the University of Minnesota.

This software is released under the GNU Affero General Public License (AGPL) v3.0 License.
"""
import csv
import pandas as pd
import numpy as np

class Taxonomy:
    def __init__(self, filename: str):
        self.tax = self.parse_taxonomy(filename)

    @classmethod
    def parse_taxonomy(cls, filename: str) -> dict:
        with open(filename) as inf:
            csv_inf = csv.reader(inf, delimiter='\t')
            return dict(csv_inf)

    def __call__(self, id: str):
        return self.tax[id]

TAX_LEVELS = ['k', 'p', 'c', 'o', 'f', 'g', 's', 't']
class Tree:
    def __init__(self):
        pass

    def add_child(self, node):
        pass

def parse_bayes(filename: str) -> pd.DataFrame:
    columns = ["tax"] + TAX_LEVELS + ["genome_length"]
    df = pd.read_csv("shogun/tests/data/sheared_bayes.txt", sep="\t", header=None, names=columns, index_col = 0)
    df = df.sort_index()



    __all__ = ["Taxonomy", "parse_bayes"]
