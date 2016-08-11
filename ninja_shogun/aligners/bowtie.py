from ninja_utils.utils import run_command

from .. import SETTINGS


def bowtie2(infile, outfile, database, num_threads=SETTINGS.N_jobs, shell=False):
    cmd = ['bowtie2',
           '--no-unal',
           '-x', database,
           '-S', infile,
           '--np', '0',
           '--mp', '"1,1"',
           '--rdg', '"0,1"',
           '--rfg', '"0,1',
           '--score-min', '"L,0,-0.02"',
           '--norc',
           '-f', outfile,
           '--very-sensitive',
           '-k', num_threads,
           '-p', '24',
           '--no-hd']
    return run_command(cmd, shell=shell)


def bowtie2_build(infile, outfile, offrate=3, shell=False):
    cmd = ['bowtie2', '-f', '-o', offrate, infile, outfile]
    return run_command(cmd, shell=shell)
