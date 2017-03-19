#!/usr/bin/env Rscript
# usage
# Rscript *.r embalmer_taxon_table.tsv ref2taxa.txt > taxon_table_with_taxonomy_column.tsv
args <- commandArgs(trailing=TRUE)
ref2taxa <- read.table(args[2],sep='\t',head=T,row=1,check=F,comment='',quote='')
x <- read.table(args[1],sep='\t',head=T,check=F,row=1,comment='',quote='')
print(mean(rownames(x) %in% rownames(ref2taxa)))
x$taxonomy <- ref2taxa[match(rownames(x),rownames(ref2taxa)),1]
sink(args[1])
cat('#Taxon ID\t')
write.table(x,sep='\t',quote=F)
sink(NULL)
