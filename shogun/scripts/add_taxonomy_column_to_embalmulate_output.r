#!/usr/bin/env Rscript
# adds taxonomy column to tsv otu table
# if only two args provided, assumes rownames are the taxonomy names
#
# usage (including a mapping of refID to taxonomy)
# Rscript *.r embalmer_taxon_table.tsv ref2taxa.txt embalmer_taxon_table.tsv
#
# usage (not including a mapping of refID to taxonomy)
# Rscript *.r embalmer_taxon_table.tsv embalmer_taxon_table.tsv

args <- commandArgs(trailing=TRUE)

# skip first line if comment
line1 <- scan(args[1],w='c',sep='\t',nlines=1,quiet=TRUE)
line2 <- scan(args[1],w='c',sep='\t',nlines=1,skip=1,quiet=TRUE)
skip <- 0
if(length(line1) < length(line2)) skip <- 1

x <- read.table(args[1],sep='\t',head=T,check=F,row=1,comment='',quote='',skip=skip)
x <- x[rownames(x) != '',]

if(length(args) == 3){
    outfp <- args[3]
    ref2taxa <- read.table(args[2],sep='\t',head=T,row=1,check=F,comment='',quote='')
    x$taxonomy <- as.character(ref2taxa[match(rownames(x),rownames(ref2taxa)),1])
} else {
    outfp <- args[2]
    x$taxonomy <- rownames(x)
}

# add blanks for missing levels
for(i in 1:nrow(x)){
    taxonomy <- x$taxonomy[i]
    levels <- strsplit(taxonomy,';')[[1]]
    if(!grepl('p__',taxonomy)){
        taxonomy <- paste(taxonomy,'; p__; c__; o__; f__; g__; s__',sep='')
    } else if(!grepl('c__',taxonomy)){
        taxonomy <- paste(taxonomy,'; c__; o__; f__; g__; s__',sep='')
    } else if(!grepl('o__',taxonomy)){
        taxonomy <- paste(taxonomy,'; o__; f__; g__; s__',sep='')
    } else if(!grepl('f__',taxonomy)){
        taxonomy <- paste(taxonomy,'; f__; g__; s__',sep='')
    } else if(!grepl('g__',taxonomy)){
        taxonomy <- paste(taxonomy,'; g__; s__',sep='')
    } else if(!grepl('s__',taxonomy)){
        taxonomy <- paste(taxonomy,'; s__',sep='')
    }
    x$taxonomy[i] <- taxonomy
}


sink(outfp)
cat('#Taxon ID\t')
write.table(x,sep='\t',quote=F)
sink(NULL)
