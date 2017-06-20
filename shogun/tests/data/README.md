Reads were simulated using the software package dwgsim from the genomes.small.fna.

We simulated only 100 reads at 50bps, with a 5% sequencing error, 2% substitution error rate, and a 2% insertion/deletion rate. The command used is below.

```
dwgsim -e 0.02 -B -f CTAG -c 2 -1 50 -2 0 -N 100 genomes.small.fna read
```
