"""
Copyright 2015-2017 Knights Lab, Regents of the University of Minnesota.

This software is released under the GNU Affero General Public License (AGPL) v3.0 License.
"""
import os
import csv

from .burst_aligner import  BurstAligner

from shogun.wrappers import burst_align_best
from shogun.utils import read_fasta

class BurstAlignerBest(BurstAligner):
    _name = 'burst_best'

    def align(self, infile, outdir):
        if not os.path.exists(outdir):
            os.makedirs(outdir)

        self.outfile = os.path.join(outdir, 'alignment.burst.best.b6')

        #TODO: pie chart and coverage
        proc, out, err = burst_align_best(infile, self.outfile, self.database, accelerator=self.accelerator, shell=self.shell, threads=self.threads)
        if self.post_align:
            post_align = self._post_align(infile, self.outfile)
            post_filtered = os.path.join(outdir, 'combined_seqs.filtered.fna')
            post_unfiltered = os.path.join(outdir, 'combined_seqs.fna')
            with open(post_filtered, 'w') as outf_post_filtered:
                with open(post_unfiltered, 'w') as outf_post_unfiltered:
                    for title, seq, alignment_score in post_align:
                        if alignment_score >= self.percent_id:
                            outf_post_filtered.write(">%s\n%s\n" % (title, seq))
                        else:
                            outf_post_unfiltered.write(">%s\n%s\n" % (title, seq))
        return proc, out, err

    def _post_align(self, inf, outf):
       with open(outf) as alignment_file:
           alignment_gen = csv.reader(alignment_file, delimiter="\t")
           with open(inf) as fasta_infile:
               fasta_gen = read_fasta(fasta_infile)
               for row, (title, seq) in zip(alignment_gen, fasta_gen):
                   alignment_score = float(row[2])
                   yield title, seq, alignment_score


