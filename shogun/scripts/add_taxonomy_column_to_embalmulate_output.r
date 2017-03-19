#!/usr/bin/env Rscript
# usage
# Rscript *.r embalmer_taxon_table.tsv > taxon_table_with_taxonomy_column.tsv
args <- commandArgs(trailing=TRUE)
x <- read.table(args[1],sep='\t',head=T,check=F,row=1,comment='',quote='')
x$taxonomy <- rownames(x)
sink(args[1])
cat('#Taxon ID\t')
write.table(x,sep='\t',quote=F)
sink(NULL)
