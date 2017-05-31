"""
Copyright 2015-2017 Knights Lab, Regents of the University of Minnesota.

This software is released under the GNU Affero General Public License (AGPL) v3.0 License.
"""

import subprocess
import sys

from ..utils import run_command


def embalmer_search(input_fp, output_fp, embalmer_db, embalmer_tax, embalmer_acc, threads, pct_id, shell=True, taxa_ncbi=False):
    """

    :param input_fp:
    :param output_fp:
    :param embalmer_db:
    :param embalmer_tax:
    :param embalmer_acc:
    :param threads:
    :param pct_id:
    :param shell:
    :param taxa_ncbi:
    :return:
    """

    #TODO: Look up SOP

    cmd = [
        'e15cb',
        '--queries', input_fp,
        '--references', embalmer_db,
        '--taxonomy', embalmer_tax,
        '--accelerator', embalmer_acc,
        '--output', output_fp,
        '--threads', str(threads),
        '--mode', "CAPITALIST",
        '--id', str(pct_id),
        '--npenalize',
        '--taxasuppress',
        '--skipambig'
    ]

    if taxa_ncbi:
        cmd += ['--taxa_ncbi']

    return run_command(cmd, shell=shell)
