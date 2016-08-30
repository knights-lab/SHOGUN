from ninja_utils.utils import run_command

from .. import SETTINGS


def utree_build(input_fasta, input_fasta_labels, output_uncompressed_tree, threads=SETTINGS.N_jobs, shell=False):
    # usage: utree-build input_fasta.fa labels.map output.ubt [threads]
    cmd = [
        'utree-build',
       input_fasta,
       input_fasta_labels,
       output_uncompressed_tree,
       threads,
    ]
    cmd = [str(i) for i in cmd]
    return run_command(cmd, shell=shell)


def utree_build_gg(input_fasta, input_fasta_labels, output_uncompressed_tree, threads=SETTINGS.N_jobs, shell=False):
    # usage: utree-buildGG input_fasta.fa labels.map output.ubt [threads]
    cmd = [
        'utree-build_gg',
        input_fasta,
        input_fasta_labels,
        output_uncompressed_tree,
        threads,
    ]
    cmd = [str(i) for i in cmd]
    return run_command(cmd, shell=shell)


def utree_compress(input_uncompressed_tree, output_compressed_tree, shell=False):
    # usage: xtree-compress preTree.ubt compTree.ctr
    cmd = [
        'utree-compress',
        input_uncompressed_tree,
        output_compressed_tree
    ]
    cmd = [str(i) for i in cmd]
    return run_command(cmd, shell=shell)


def utree_search(input_compressed_tree, input_fasta_to_search, output, shell=False):
    # usage: xtree-search compTree.ctr fastaToSearch.fa output.txt
    cmd = [
        'utree-search',
        input_compressed_tree,
        input_fasta_to_search,
        output
    ]
    cmd = [str(i) for i in cmd]
    return run_command(cmd, shell=shell)


def utree_search_gg(input_compressed_tree, input_fasta_to_search, output, shell=False):
    # usage: xtree-searchGG compTreeGG.ctr fastaToSearch.fa output.txt
    cmd = [
        'utree-search_gg',
        input_compressed_tree,
        input_fasta_to_search,
        output
    ]
    cmd = [str(i) for i in cmd]
    return run_command(cmd, shell=shell)
