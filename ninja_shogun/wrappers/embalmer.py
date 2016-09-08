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
        '--latency', 64,
        '--id', .985,
    ]
    return run_command(cmd, shell=shell)
