"""
Copyright 2015-2017 Knights Lab, Regents of the University of Minnesota.

This software is released under the GNU Affero General Public License (AGPL) v3.0 License.
"""

import csv
import glob
import os
from collections import defaultdict, Counter

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix

from shogun import logger

TAXA = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'strain']

def function_run_and_save(input, func_db, output, level):
    prefix = ".".join(os.path.basename(input).split('.')[:-1])

    kegg_modules_df = func_db['modules']
    row_names = func_db['names']
    kegg_ids = func_db['kegg_ids']
    kegg_table_csr = func_db['csr']

    logger.debug("Level for summarization %d and starting summarizing KEGG Table at level with median." % level)
    if level < 8:
        kegg_table_csr, row_names = summarize_at_level(kegg_table_csr, row_names, kegg_ids, level)
    logger.debug("Head of row names %s" % str(list(row_names.keys())[:3]))

    if TAXA[level-1] not in prefix:
        prefix += "." + TAXA[level-1]

    logger.info("Reading in taxatable for functional prediction at %s." % os.path.abspath(input))
    taxatable_df = pd.read_csv(input, sep="\t", index_col=0)
    logger.debug("Taxatable for functional prediction shape %s" % str(taxatable_df.shape))
    taxatable_df = taxatable_df[[type(_) == str for _ in taxatable_df.index]]

    taxatable_df['summary'] = [';'.join(_.split(';')[:level]).replace(' ', '_') for _ in taxatable_df.index]
    taxatable_df = taxatable_df.groupby('summary').sum()

    # Normalizing for depth at median depth
    logger.debug("Normalizing to median depth")
    taxatable_df = taxatable_df.div(taxatable_df.sum(axis=0).div(taxatable_df.sum(axis=0).median()), axis=1).round().astype(int)

    logger.debug("Taxatable summarized shape %s" % str(taxatable_df.shape))

    logger.info("Starting functional prediction.")
    out_kegg_table_df, out_kegg_modules_df, out_kegg_modules_coverage = _do_function(taxatable_df, row_names, kegg_ids, kegg_table_csr, kegg_modules_df)
    out_kegg_table_df.to_csv(os.path.join(output, "%s.kegg.txt" % prefix), sep='\t', float_format="%d",na_rep=0, index_label="#KEGG ID")
    out_kegg_modules_df.to_csv(os.path.join(output, "%s.kegg.modules.txt" % prefix), sep='\t', float_format="%d",na_rep=0, index_label="#MODULE ID")
    out_kegg_modules_coverage.to_csv(os.path.join(output, "%s.kegg.modules.coverage.txt" % prefix), sep='\t', float_format="%f", na_rep=0, index_label="#MODULE ID")


def summarize_at_level(csr, names, kegg_ids, level):
    s = pd.DataFrame(sorted(names, key=names.get), columns=["names"])
    s['group'] = [";".join(_.split(';')[:level]) for _ in s['names']]
    indptr = [0]
    indices = []
    row_names = {}
    data = []
    for name, df in s.groupby("group"):
        _mat = csr[df.index, :].todense()
        # Threshold to over 80%
        _mat_thresh = np.divide((_mat > 0).sum(axis=0), _mat.shape[0]) > .8
        _mat_thresh = np.asarray(_mat_thresh).reshape(-1)
        if _mat_thresh.any():
            # Get the medians
            _medians = np.asarray(np.median(_mat[:, _mat_thresh], axis=0)).reshape(-1)
            indices.extend(np.where(_mat_thresh)[0])
            data.extend(_medians)
            row_names.setdefault(name, len(row_names))
            indptr.append(len(indices))
    return csr_matrix((data, indices, np.array(indptr)), dtype=np.int8), row_names


def _do_function(taxatable_df, row_names, column_names, kegg_table_csr, kegg_modules_df):
    _, num_kegg_ids = kegg_table_csr.shape
    logger.debug("Kegg table for functional prediction shape %s" % (str(kegg_table_csr.shape)))
    _, num_samples = taxatable_df.shape
    logger.debug("Taxatable for functional prediction shape %s" % (str(taxatable_df.shape)))

    kegg_table = np.zeros((num_samples, num_kegg_ids), dtype=np.int)

    for i, row in taxatable_df.iterrows():
        row_names_found = 0
        logger.debug(row.name)
        if row.name in row_names:
            row_names_found += 1
            idx = row_names[row.name]
            kegg_table += np.outer(row, kegg_table_csr.getrow(idx).todense())

    logger.debug("Row names found in taxatable %d" % row_names_found)

    out_kegg_table_df = pd.DataFrame(kegg_table, index=taxatable_df.columns, columns=sorted(column_names, key=column_names.get), dtype=np.int).T

    filtered_kegg_ids = out_kegg_table_df.loc[kegg_modules_df.columns].fillna(0)
    # kegg modules df
    out_kegg_modules_df = kegg_modules_df.dot(filtered_kegg_ids)

    out_kegg_modules_coverage = ((kegg_modules_df).dot(filtered_kegg_ids > 0).fillna(0)).div(kegg_modules_df.sum(axis=1), axis=0)

    # Filter out zeros
    out_kegg_modules_df = out_kegg_modules_df[(out_kegg_modules_df.T != 0).any()]

    # Filter out zeros
    out_kegg_table_df = out_kegg_table_df[(out_kegg_table_df.T != 0).any()]

    # Filter out zeros
    out_kegg_modules_coverage = out_kegg_modules_coverage[(out_kegg_modules_coverage.T != 0).any()]

    return out_kegg_table_df, out_kegg_modules_df, out_kegg_modules_coverage

def parse_function_db(metadata: dict, database: str) -> dict:
    if not 'function' in metadata:
        return {}
    else:
        file_set = set(glob.glob(os.path.join(database, metadata['function'] + '*')))
        suffices = ['module-annotations.txt', 'strain2ko.txt']
        files = ["%s-%s" % (os.path.join(database, metadata['function']), suffix) for suffix in suffices]
        for file in files:
            if file not in file_set:
                return {}

        modules_df = _parse_modules(files[0])
        #TODO: Implement the save csr, this works but requires write permissions to db folder
        #npz = os.path.join(database, metadata['function']) + "-strain2ko.npz"
        # if not npz in file_set:
        #    row_names, column_names, csr = _parse_kegg_table(files[1])
        #    save_csr_matrix(npz, csr, row_names, column_names)
        # else:
        #   row_names, column_names, csr = load_csr_matrix(npz)
        _strains = list(parse_kegg_table(files[1]))
        return dict(zip(('modules_file', 'file', 'names', 'kegg_ids', 'csr','modules'), files + _strains + [modules_df]))

def _parse_modules(infile):
    modules_keggs = defaultdict(Counter)
    with open(infile) as inf:
        csv_inf = csv.reader(inf, delimiter="\t")
        for row in csv_inf:
            modules_keggs[row[0]].update([row[-1][:7]])
    return pd.DataFrame(modules_keggs).fillna(0.0).astype(int)


def parse_kegg_table(infile):
    indptr = [0]
    indices = []
    kegg_ids = {}
    row_names = {}
    data = []
    with open(infile) as inf:
        csv_inf = csv.reader(inf, delimiter="\t")
        for line in csv_inf:
            row_names.setdefault(line[0], len(row_names))
            counts = Counter(line[1:])
            # Filter out the blank KEGG IDs
            for key, value in counts.items():
                if not key == '':
                    indices.append(kegg_ids.setdefault(key, len(kegg_ids)))
                    data.append(value)
            indptr.append(len(indices))
    return row_names, kegg_ids, csr_matrix((data, indices, indptr), dtype=np.int8)
