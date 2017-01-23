"""
Copyright 2015-2017 Knights Lab, Regents of the University of Minnesota.

This software is released under the GNU Affero General Public License (AGPL) v3.0 License.
"""

from ninja_utils.utils import run_command

from .. import SETTINGS


def embalmer_align(queries, references, outfile, mode='CAPITALIST2', num_threads=SETTINGS.N_jobs, shell=False):
    cmd = [
        'embalmer',
        '--queries', queries,
        '--references', references,
        '--output', outfile,
        '--threads', str(num_threads),
        '--mode', mode,
        '--shear',
        '--latency', 16,
        '--id', .985,
    ]
    return run_command(cmd, shell=shell)
