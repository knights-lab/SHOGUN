"""
Copyright 2015-2017 Knights Lab, Regents of the University of Minnesota.

This software is released under the GNU Affero General Public License (AGPL) v3.0 License.
"""


def yield_alignments_from_sam_inf(inf):
    with open(inf) as fh:
        for i in fh:
            line = i.split('\t')
            # this function yields qname, rname
            try:
                yield line[0], line[2]
            except IndexError:
                print('Incorrect SAM input %s' % inf)
