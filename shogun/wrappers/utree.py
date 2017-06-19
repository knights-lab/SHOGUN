"""
Copyright 2015-2017 Knights Lab, Regents of the University of Minnesota.

This software is released under the GNU Affero General Public License (AGPL) v3.0 License.
"""

from shogun.utils import run_command


def utree_build(input_fasta, input_fasta_labels, output_uncompressed_tree, threads=1, shell=False):
    # usage: utree-build input_fasta.fa labels.map output.ubt [threads]
    cmd = [
        'utree-build',
       input_fasta,
       input_fasta_labels,
       output_uncompressed_tree,
       threads,
    ]
    return run_command(cmd, shell=shell)


def utree_build_gg(input_fasta, input_fasta_labels, output_uncompressed_tree, threads=1, shell=False):
    # usage: utree-buildGG input_fasta.fa labels.map output.ubt [threads]
    cmd = [
        'utree-build_gg',
        input_fasta,
        input_fasta_labels,
        output_uncompressed_tree,
        threads,
    ]
    return run_command(cmd, shell=shell)


def utree_compress(input_uncompressed_tree, output_compressed_tree, shell=False):
    # usage: xtree-compress preTree.ubt compTree.ctr
    cmd = [
        'utree-compress',
        input_uncompressed_tree,
        output_compressed_tree
    ]
    return run_command(cmd, shell=shell)


def utree_search(input_compressed_tree, input_fasta_to_search, output, shell=False):
    # usage: xtree-search compTree.ctr fastaToSearch.fa output.txt
    cmd = [
        'utree-search',
        input_compressed_tree,
        input_fasta_to_search,
        output
    ]
    return run_command(cmd, shell=shell)


def utree_search_gg(input_compressed_tree, input_fasta_to_search, output, shell=False):
    # usage: xtree-searchGG compTreeGG.ctr fastaToSearch.fa output.txt
    cmd = [
        'utree-search_gg',
        input_compressed_tree,
        input_fasta_to_search,
        output
    ]
    return run_command(cmd, shell=shell)
