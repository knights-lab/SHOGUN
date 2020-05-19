"""
Copyright 2015-2020 Knights Lab, Regents of the University of Minnesota.

This software is released under the GNU Affero General Public License (AGPL) v3.0 License.
"""
import pysam


def yield_alignments_from_bam_inf(inf):
    fh = pysam.AlignmentFile(inf, "r", check_header=False, check_sq=False)
    for alignment in fh.fetch():
        yield alignment
