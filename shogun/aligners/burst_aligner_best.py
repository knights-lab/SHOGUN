"""
Copyright 2015-2017 Knights Lab, Regents of the University of Minnesota.

This software is released under the GNU Affero General Public License (AGPL) v3.0 License.
"""
import os
import csv

from .burst_aligner import  BurstAligner

from shogun.wrappers import burst_align_best
from shogun.utils import read_fasta
from shogun import logger


class BurstAlignerBest(BurstAligner):
    _name = 'filter'

    def align(self, infile, outdir, align=True):
        if not os.path.exists(outdir):
            os.makedirs(outdir)

        self.outfile = os.path.join(outdir, 'alignment.burst.best.b6')

        #TODO: pie chart and coverage
        if align:
            proc, out, err = burst_align_best(infile, self.outfile, self.database, accelerator=self.accelerator, shell=self.shell, threads=self.threads)
        else:
            proc, out, err = (None, None, None)
        if self.post_align:
            alignments = self._post_align(self.outfile)
            post_filtered = os.path.join(outdir, 'combined_seqs.filtered.fna')
            post_unfiltered = os.path.join(outdir, 'combined_seqs.fna')
            with open(infile) as input_fasta:
                fasta_gen = read_fasta(input_fasta)
                with open(post_filtered, 'w') as outf_post_filtered:
                    with open(post_unfiltered, 'w') as outf_post_unfiltered:
                        for title, seq in fasta_gen:
                            stripped_title = title.split()[0]
                            if stripped_title in alignments:
                                outf_post_filtered.write(">%s\n%s\n" % (title, seq))
                            else:
                                outf_post_unfiltered.write(">%s\n%s\n" % (title, seq))
        return proc, out, err

    def _post_align(self, outf):
        alignments = set()
        i = 0
        with open(outf) as alignment_file:
            alignment_gen = csv.reader(alignment_file, delimiter="\t")
            for row in alignment_gen:
                alignment_score = float(row[2])
                if alignment_score >= self.percent_id:
                    alignments.add(row[0])
                    i += 1
        logger.info("Human hits filter: %d" % i)
        return alignments

