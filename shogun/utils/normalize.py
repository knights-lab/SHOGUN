"""
Copyright 2015-2020 Knights Lab, Regents of the University of Minnesota.

This software is released under the GNU Affero General Public License (AGPL) v3.0 License.
"""

from shogun import logger


def normalize_by_median_depth(df):
    logger.debug("Normalizing to median depth")
    return df.div(df.sum(axis=0).div(df.sum(axis=0).median()), axis=1).round().astype(int)
