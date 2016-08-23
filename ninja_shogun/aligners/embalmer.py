from ninja_utils.utils import run_command

from .. import SETTINGS


def embalmer_align(queries, references, outfile, mode='CAPITALIST', num_threads=SETTINGS.N_jobs, shell=False):
    cmd = ['embalmer',
           '--queries', queries,
           '--references', references,
           '--output', outfile,
           '--threads', str(num_threads),
           '--mode', mode,
           '--shear',
           ]
    return run_command(cmd, shell=shell)
