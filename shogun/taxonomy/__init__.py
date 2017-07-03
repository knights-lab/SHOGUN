"""
Copyright 2015-2017 Knights Lab, Regents of the University of Minnesota.

This software is released under the GNU Affero General Public License (AGPL) v3.0 License.
"""
import csv
import pandas as pd
from collections import defaultdict
import numpy as np

def tree(): return defaultdict(tree)

def add_tree(t, path):
  for node in path.split(';'):
    t = t[node]

def longest_path_tree(t, path):
    s = []
    temp_spot = t
    for node in path.split(';'):
        if node in temp_spot:
            temp_spot = temp_spot[node]
            s.extend([node])
        else:
            break
    return ';'.join(s)

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

def parse_bayes(filename: str) -> pd.DataFrame:
    columns = ["tax"] + TAX_LEVELS + ["genome_length"]
    df = pd.read_csv(filename, sep="\t", header=None, names=columns, index_col = 0)
    return df.sort_index()

def pie_chart_taxatable(filename: str, counts_bayes: pd.DataFrame, level=8):
    df = pd.read_csv(filename, sep="\t", index_col=0)
    df = df[[type(_) == str for _ in df.index]]

    cb_index = tree()
    _ = [add_tree(cb_index, v) for v in counts_bayes.index]

    df['summary'] = [longest_path_tree(cb_index, v) for v in df.index]
    df = df.groupby('summary').sum()

    df['level'] = [_.count(';') + 1 if type(_) == str else 0 for _ in df.index]

    # summarize up
    below_level = df['level'] >= level
    leaf_counts_df = df[below_level]
    leaf_counts_df['summary'] = [';'.join(v.split(';')[:level]) for v in df[below_level].index]
    leaf_counts_df = leaf_counts_df.groupby('summary').sum()
    leaf_counts_df = leaf_counts_df.drop('level', axis=1)

    # summarize bayes to level
    counts_bayes_sum = _summarize_bayes_at_level(counts_bayes, leaf_counts_df.index, level=level)

    # summarize down
    for i, row in df[~below_level].sort_values('level', ascending=False).iterrows():
        # Get all children of item
        tmp_name = row.name
        leave_filter = _filter_leaves_for_tax(leaf_counts_df, tmp_name)
        num_leaves = np.sum(leave_filter)
        if num_leaves == 0 or num_leaves is None:
            print(tmp_name)
            # Filter back row names until in counts_bayes
            blank = ['k__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__', 't__']
            for i, _ in enumerate(row.name.split(';')):
                blank[i] = _
            tmp_counts_bayes_row = _summarize_bayes_at_level(counts_bayes, row.name, level=row.name.count(';') + 1)
            tmp_counts_bayes_row.name = ';'.join(blank[:level])
            row.name = tmp_counts_bayes_row.name
            leaf_counts_df = leaf_counts_df.append(row[:-1])
            if tmp_counts_bayes_row.name not in counts_bayes_sum.index:
                counts_bayes_sum = counts_bayes_sum.append(tmp_counts_bayes_row)
                counts_bayes_sum = counts_bayes_sum.fillna(0)
        elif num_leaves == 1:
            leaf_counts_df.loc[leave_filter] += row.values[:-1]
        elif num_leaves > 1:
            tmp_level = row.name.count(';')
            tmp_leaves = leaf_counts_df[leave_filter].sort_index()
            tmp_bayes = counts_bayes_sum.loc[tmp_leaves.index]
            # Series 1xn where n is the number of leave nodes below tax
            prob_tax_given_level = (tmp_bayes.iloc[:,tmp_level] + 1)/(tmp_bayes['genome_length'] + 1)
            prob_tax_given_level = prob_tax_given_level/np.sum(prob_tax_given_level)
            # Series 1xn where n is the number of unique reads for a given taxa
            uniqueness_per_genome = tmp_bayes.iloc[:,level-1]/tmp_bayes['genome_length']
            # Matrix divide each observed count by uniqueness
            counts_over_uniqueness = tmp_leaves.T / uniqueness_per_genome.values
            # Matrix divide each uniqueness count by sum of sample
            prob_tax = counts_over_uniqueness.T / counts_over_uniqueness.sum(axis=1)
            # Get the redistribution parameters
            # Should be taxa by samples same as the tmp_leaves
            # Each column should sum to 1
            redistribution_params = prob_tax.apply(lambda x: x*prob_tax_given_level.values, axis=0).apply(lambda x: x/x.sum(), axis=0)
            redistribution_numbers = (redistribution_params * row.values[:-1]).round()
            # Add the number back to the dataframe
            leaf_counts_df = leaf_counts_df.add(redistribution_numbers, fill_value=0)
    return leaf_counts_df


def _summarize_bayes_at_level(counts_bayes: pd.DataFrame, leave_names, level=7):
    if level < 8:
        counts_bayes['summary_taxa'] = [';'.join(_.split(';')[:level]) for _ in counts_bayes.index]
        counts_bayes = counts_bayes.groupby('summary_taxa').sum()
        counts = counts_bayes.iloc[:, level-1:8].sum(axis=1)
        counts_bayes.iloc[:, level-1] = counts
        counts_bayes = counts_bayes.drop(counts_bayes.columns[[level, 7]], axis=1)
        counts_bayes = counts_bayes.loc[leave_names]
    return counts_bayes


def _filter_leaves_for_tax(leaf_counts_df, taxa):
     return np.array([_.startswith(taxa + ';') for _ in leaf_counts_df.index])

def parse_taxatable2(filename: str, level: int = 7):
    with open(filename) as inf:
        csv_inf = csv.reader(inf, delimiter='\t')
        header = next(csv_inf)
        # sample_collection -> collection of samples
        num_samples = len(header)-1
        leaves = defaultdict(lambda: np.zeros(num_samples, dtype=int))
        inner_nodes = defaultdict(lambda: np.zeros(num_samples, dtype=int))
        for row in csv_inf:
            tax_row = row[0].split(';')
            if len(tax_row) >= level:
                tax = ';'.join(tax_row[:level])
                leaves[tax] += np.array(row[1:])
            else:
                tax = ';'.join(tax_row)
                inner_nodes[tax] = np.array(row[1:])
        df_leaves = pd.DataFrame(leaves)
        df_inner = pd.DataFrame(inner_nodes)

__all__ = ["Taxonomy", "parse_bayes", "pie_chart_taxatable"]
