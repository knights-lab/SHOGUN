"""
Copyright 2015-2017 Knights Lab, Regents of the University of Minnesota.

This software is released under the GNU Affero General Public License (AGPL) v3.0 License.
"""

import os

from shogun.utils import run_command


def burst_align(input_fp, output_fp, burst_db_prefix, threads=1, pct_id=.98, tax=False,
                   accelerator=False, taxacut=5, shell=False, taxa_ncbi=False):
    """

    :param input_fp:
    :param output_fp:
    :param burst_db_prefix:
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
        'burst15',
        '--queries', input_fp,
        '--references', burst_db_prefix + '.edx',
        '--output', output_fp,
        '--threads', str(threads),
        '--mode', 'CAPITALIST',
        '--id', str(pct_id),
        '--npenalize',
        '--skipambig',
        '--forwardreverse'
    ]

    if taxa_ncbi:
        cmd += ['--taxa_ncbi']

    if accelerator:
        cmd += ['--accelerator', accelerator]

    if tax:
        cmd += ['--taxonomy', tax, '--taxacut', taxacut]

    return run_command(cmd, shell=shell)


def burst_build(infile, outfile_prefix, accelerator=False, shell=False, clustradius=None, shear=None):
    cmd = [
        'burst15',
        '--references', infile,
        '--output', outfile_prefix + ".edb",
        '--npenalize',
        '--makedb',
        '--fingerprint',
    ]

    if accelerator:
        cmd += ['--accelerator', accelerator]

    if shear:
        cmd += ['--shear', str(shear)]

    if clustradius:
        cmd += ['--clustradius', str(clustradius)]

    return run_command(cmd, shell=shell)


def embalmulate(infile, outdir, shell=False):
    cmd = [
        'embalmulate',
        infile,
        os.path.join(outdir, 'burst_otutable.txt'),
        os.path.join(outdir, 'burst_taxatable.txt'),
        'GGtrim'
    ]

    return run_command(cmd, shell=shell)
