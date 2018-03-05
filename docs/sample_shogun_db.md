# Make database for BURST
mkdir burst
burst15 --references genomes.small.fna --taxonomy genomes.small.tax --output burst/genomes.small.edx --npenalize --makedb --fingerprint --accelerator burst/genomes.small.acx

```
# Make database for Bowtie2
mkdir bowtie2
bowtie2-build -f genomes.small.fna bowtie2/genomes.small
```

```
# Make database for UTree
mkdir utree
utree-build_gg genomes.small.fna genomes.small.tax genomes.small.utr 8 0 RC

```

```
# Shear the database
```

```
# Prepare for Prokka
```
