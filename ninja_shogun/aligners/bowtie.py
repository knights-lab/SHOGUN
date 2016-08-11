from ninja_utils.utils import run_command

from .. import SETTINGS


def bowtie2(infile, outfile, database, alignments_to_report=32, num_threads=SETTINGS.N_jobs, shell=False):
    cmd = ['bowtie2',
           '--no-unal',
           '-x', database,
           '-S', outfile,
           '--np', '0',
           '--mp', '"1,1"',
           '--rdg', '"0,1"',
           '--rfg', '"0,1"',
           '--score-min', '"L,0,-0.02"',
           '--norc',
           '-f', infile,
           '--very-sensitive',
           '-k', str(alignments_to_report),
           '-p', str(num_threads),
           '--no-hd']
    return run_command(cmd, shell=shell)


def bowtie2_build(infile, outfile, offrate=3, shell=False):
    cmd = ['bowtie2-build', '-f', '-o', str(offrate), infile, outfile]
    return run_command(cmd, shell=shell)
