"""
Copyright 2015-2017 Knights Lab, Regents of the University of Minnesota.

This software is released under the GNU Affero General Public License (AGPL) v3.0 License.
"""

from ninja_utils.utils import run_command
import subprocess
import sys

from .. import SETTINGS

def embalmer_search(input_fp, output_fp, embalmer_db, embalmer_tax, embalmer_acc, threads, pct_id, shell=True):
    cmd = [
        'embalmer',
        '--queries', input_fp,
        '--references', embalmer_db,
        '--taxonomy', embalmer_tax,
        '--accelerator', embalmer_acc,
        '--output', output_fp,
        '--threads', str(threads),
        '--mode', "CAPITALIST",
        '--id', str(pct_id),
        '--npenalize',' '
        '--fingerprint', ' '
        '--taxasuppress', ' '
    ]

    print('Running: ' + ' '.join(cmd))

    process = subprocess.Popen(' '.join(cmd), shell=True, stderr=subprocess.PIPE)
    for line in iter(process.stderr.readline, b''):
        sys.stdout.write(line)
        sys.stdout.flush()
    stdout, stderr = process.communicate()
    return process.returncode
