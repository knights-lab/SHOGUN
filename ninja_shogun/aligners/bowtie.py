from ninja_utils.utils import run_command

from .. import SETTINGS


def bowtie2(infile, outfile, database, alignments_to_report=32, num_threads=SETTINGS.N_jobs, shell=False):
    """
    Search a bowtie2 index with multiple alignment.
    :param infile: the query FASTA file
    :param outfile: the resulting SAM file
    :param database: path to the bowtie2 index
    :param alignments_to_report: the number of alignments to report (default=32)
    :param num_threads: the number of threads to use (default=SETTINGS)
    :param shell: whether to use the shell NOT RECOMMENDED (default=False)
    :return: the STDERR/STDOUT
    """
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
    """
    This function will build a bowtie2 index with a given infile and output to the outfile
    :param infile: the FASTA file to build the index with
    :param outfile: the prefix for the bowtie2 index
    :param offrate: offrate for the index (default=3)
    :param shell: whether to use the shell NOT RECOMMENDED (default=False)
    :return: teh STDERR/STDOUT
    """
    cmd = ['bowtie2-build', '-f', '-o', str(offrate), infile, outfile]
    return run_command(cmd, shell=shell)
