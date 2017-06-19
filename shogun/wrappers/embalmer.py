"""
Copyright 2015-2017 Knights Lab, Regents of the University of Minnesota.

This software is released under the GNU Affero General Public License (AGPL) v3.0 License.
"""

import os

from ..utils import run_command


def embalmer_align(input_fp, output_fp, embalmer_db_prefix, threads=1, pct_id=.98, tax=False,
                   accelerator=False, shell=False, taxa_ncbi=False):
    """

    :param input_fp:
    :param output_fp:
    :param embalmer_db_prefix:
    :param accelerator:
    :param tax:
    :param threads:
    :param pct_id:
    :param shell:
    :param taxa_ncbi:
    :return:
    """

    #TODO: Look up SOP

    cmd = [
        'emb15',
        '--queries', input_fp,
        '--references', embalmer_db_prefix + '.edb',
        '--output', output_fp,
        '--threads', str(threads),
        '--mode', 'CAPITALIST',
        '--id', str(pct_id),
        '--npenalize',
        '--taxasuppress',
        '--skipambig'
    ]

    if taxa_ncbi:
        cmd += ['--taxa_ncbi']

    if not accelerator:
        if os.path.exists(embalmer_db_prefix + '.acc'):
            cmd += ['--accelerator', embalmer_db_prefix + '.acc']
    else:
        cmd += ['--accelerator', accelerator]

    if not accelerator:
        if os.path.exists(embalmer_db_prefix + '.tax'):
            cmd += ['--taxaonomy', embalmer_db_prefix + '.tax']
    else:
        cmd += ['--taxonomy', tax]

    return run_command(cmd, shell=shell)


def embalmer_build(infile, outfile_prefix, accelerator=False, shell=False, cr=None, s=None):
    cmd = [
        'emb15',
        '-r', infile,
        '-o', outfile_prefix + ".edb",
        '-n',
        '-d',
        '-f',
    ]

    cmd += ['-s']

    if accelerator:
        cmd += ['--accelerator', accelerator]

    if s:
        cmd += [s]

    if cr:
        cmd += ['-cr', cr]

    return run_command(cmd, shell=shell)
