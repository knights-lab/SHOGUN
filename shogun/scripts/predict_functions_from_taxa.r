#!/usr/bin/env Rscript
# takes mixed taxonomy abundance table
# and species- and strain-to-function mapping tables
# outputs a function table
#
# usage:
# Rscript *.r taxatable.txt species2ko.txt strain2ko.txt kotable.txt
#
# taxatable.txt is a QIIME-formatted species or mixed species/strain/other-level table
# such as would be put out by embalmer's embalmulate script
#
# species2ko.txt and strain2ko.txt are tab-delimited tables where the first column is
# taxon names (same format as those in taxatable.txt) and the rest of each line is
# a tab-delimited list of the KEGG KOs in that taxon. See shogun/databases/readme.txt
# for how these files would be generated. Current files of this type are:
# shogun/databases/megaGMG.species2ko.txt shogun/databases/megaGMG.strain2ko.txt
#
# kotable.txt is the output file name
#

args <- commandArgs(trailing=TRUE)
# skip first line if comment
line1 <- scan(args[1],w='c',sep='\t',nlines=1,quiet=TRUE)
line2 <- scan(args[1],w='c',sep='\t',nlines=1,skip=1,quiet=TRUE)
skip <- 0
if(length(line1) < length(line2)) skip <- 1

taxa <- read.table(args[1], sep='\t',head=T,row=1,check=F,comment='',quote='',skip=1)
taxa <- taxa[grep("s__",rownames(taxa)),]
rownames(taxa) <- gsub(' ','',rownames(taxa))
# drop last column if it is text-based
if(class(taxa[,ncol(taxa)]) != "numeric") taxa <- taxa[,-ncol(taxa)]
taxa <- as.matrix(taxa)

# load species 2 ko
cat('Loading species 2 KO map...\n')
lines <- readLines(args[2])
sp2ko <- strsplit(lines,'\t')
names(sp2ko) <- sapply(sp2ko,'[',1)
sp2ko <- sapply(sp2ko,'[',-1)

# load strain 2 ko
cat('Loading strain 2 KO map...\n')
lines <- readLines(args[3])
st2ko <- strsplit(lines,'\t')
names(st2ko) <- sapply(st2ko,'[',1)
st2ko <- sapply(st2ko,'[',-1)

tax2ko <- c(sp2ko,st2ko)
names(tax2ko) <- gsub(' ','',names(tax2ko))

cat('Identifying unique KOs...\n')
all.kos <- sort(unique(c(unlist(sp2ko),unlist(st2ko))))

kos <- matrix(0,nrow=length(all.kos), ncol=ncol(taxa))
rownames(kos) <- all.kos
colnames(kos) <- colnames(taxa)

taxa <- sweep(taxa,2,colSums(taxa),'/')

species.names <- sapply(strsplit(rownames(taxa),';'),function(xx) paste(xx[1:min(length(xx),7)],collapse=';'))

cat('Tallying KOs from taxa abundance...\n')
for(i in 1:nrow(taxa)){
    cat(i,' ')
    taxon <- rownames(taxa)[i]

    # if strain not in list, back off to species
    if(grepl('t__',taxon) && !(taxon %in% names(tax2ko))){
        cat('Strain',taxon,'not in ko mapping, backing off to species.\n')
        taxon <- species.names[i]
    }
    if(taxon %in% names(tax2ko)){
        taxon.kos <- tax2ko[[taxon]]
        for(ko in taxon.kos){
            # cat(ko,'...')
            # print(kos[ko,1:3])
            kos[ko,] <- kos[ko,] + taxa[i,]
        }
    }
}
cat('\n')

kos <- sweep(kos, 2, colSums(kos),'/')
kos[is.nan(kos)] <- 0
kos <- kos[rowSums(kos) > 0,]

sink(args[4])
cat('#OTU ID\t')
write.table(kos,sep='\t',quote=F)
sink(NULL)
