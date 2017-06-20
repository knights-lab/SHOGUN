"""
Copyright 2015-2017 Knights Lab, Regents of the University of Minnesota.

This software is released under the GNU Affero General Public License (AGPL) v3.0 License.
"""

import os

from shogun.utils import run_command


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
        '--references', embalmer_db_prefix + '.edx',
        '--output', output_fp,
        '--threads', str(threads),
        '--mode', 'CAPITALIST',
        '--id', str(pct_id),
        '--npenalize',
        '--taxasuppress',
        '--skipambig',
        '--forwardreverse'
    ]

    if taxa_ncbi:
        cmd += ['--taxa_ncbi']

    if accelerator:
        cmd += ['--accelerator', accelerator]

    if tax:
        cmd += ['--taxonomy', tax, '--taxasuppress']

    return run_command(cmd, shell=shell)


def embalmer_build(infile, outfile_prefix, accelerator=False, shell=False, cr=None, s=None):
    cmd = [
        'emb15',
        '--references', infile,
        '--output', outfile_prefix + ".edb",
        '--npenalize',
        '--makedb',
        '--fingerprint',
    ]

    if accelerator:
        cmd += ['--accelerator', accelerator]

    if s:
        cmd += ['--shear', str(s)]

    if cr:
        cmd += ['--clustradius', str(cr)]

    return run_command(cmd, shell=shell)


def embalmulate(infile, outdir, shell=False):
    cmd = [
        'embalmulate',
        infile,
        os.path.join(outdir, 'otutable.txt'),
        os.path.join(outdir, 'taxatable.txt'),
        'GGtrim'
    ]

    return run_command(cmd, shell=shell)
