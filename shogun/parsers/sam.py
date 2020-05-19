"""
Copyright 2015-2020 Knights Lab, Regents of the University of Minnesota.

This software is released under the GNU Affero General Public License (AGPL) v3.0 License.
"""


def yield_alignments_from_sam_inf(inf):
    # this function yields qname, rname
    with open(inf) as fh:
        try:
            i = next(fh)
            line = i.split('\t')
            c_qname, rname = line[0], line[2]
            record = [[c_qname, rname]]
        except BaseException as e:
            print('Incorrect SAM input %s' % inf)
            raise e
        for i in fh:
            line = i.split('\t')
            try:
                qname, rname = line[0], line[2]
                if qname != c_qname:
                    yield record
                    record = []
                    c_qname = qname
                else:
                    record.append((qname, rname))
            except IndexError as e:
                print('Incorrect SAM input %s' % inf)
                raise e
        # final record
        yield record
