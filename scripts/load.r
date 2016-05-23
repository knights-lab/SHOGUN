DATADIR <- '~/drive/research/ibs/data/'
mapfp <- paste(DATADIR,'mappingfile_deidentified_nov2015_edit.txt',sep='')
taxfp <- paste(DATADIR,'taxonomy_table.txt',sep='')

m <- read.table(mapfp, sep='\t',head=T,row=1,check=F,comment='')
x <- t(read.table(taxfp, sep='\t',head=T,row=1,check=F,comment=''))
x <- x[rownames(m),]
x <- x[,colnames(x) != '_']

# species only
xsp <- x[,!grepl('_$',colnames(x))]
colnames(xsp) <- gsub('_',' ',colnames(xsp))

# genus only
genus <- sapply(strsplit(colnames(x),'_'),'[',1)
xgn <- t(apply(x,1, function(xx) sapply(split(xx,genus),sum)))

# drop low-abundance samples
xsp <- xsp[rowSums(xsp) > 1000,]
xgn <- xgn[rownames(xsp),]
m <- m[rownames(xsp),]

# normalize 
# xsp <- sweep(xsp, 1, rowSums(xgn), '/') # normalize species by genus-level assignments
xsp.raw <- xsp
xgn.raw <- xgn
xsp <- sweep(xsp, 1, rowSums(xsp), '/')
xgn <- sweep(xgn, 1, rowSums(xgn), '/')

# collapse by subject using mean
xspc <- apply(xsp,2,function(xx) sapply(split(xx,m$Patient_no.),mean))
rownames(xspc) <- sprintf('Subject_%03d',sapply(split(m$Patient_no.,m$Patient_no.),'[',1))
xgnc <- apply(xgn,2,function(xx) sapply(split(xx,m$Patient_no.),mean))
rownames(xgnc) <- sprintf('Subject_%03d',sapply(split(m$Patient_no.,m$Patient_no.),'[',1))
mc <- m[sapply(split(1:nrow(m),m$Patient_no.),'[',1),,drop=TRUE]
rownames(mc) <- sprintf('Subject_%03d',sapply(split(m$Patient_no.,m$Patient_no.),'[',1))

# collapse by subject using first sample, for example
xspc.first <- apply(xsp,2,function(xx) sapply(split(xx,m$Patient_no.),'[',1))
rownames(xspc.first) <- sprintf('Subject_%03d',sapply(split(m$Patient_no.,m$Patient_no.),'[',1))
xgnc.first <- apply(xgn,2,function(xx) sapply(split(xx,m$Patient_no.),'[',1))
rownames(xgnc.first) <- sprintf('Subject_%03d',sapply(split(m$Patient_no.,m$Patient_no.),'[',1))
mc <- m[sapply(split(1:nrow(m),m$Patient_no.),'[',1),,drop=TRUE]
rownames(mc) <- sprintf('Subject_%03d',sapply(split(m$Patient_no.,m$Patient_no.),'[',1))


xspc.raw <- apply(xsp.raw,2,function(xx) sapply(split(xx,m$Patient_no.),sum))
rownames(xspc.raw) <- sprintf('Subject_%03d',sapply(split(m$Patient_no.,m$Patient_no.),'[',1))
xgnc.raw <- apply(xgn.raw,2,function(xx) sapply(split(xx,m$Patient_no.),sum))
rownames(xgnc.raw) <- sprintf('Subject_%03d',sapply(split(m$Patient_no.,m$Patient_no.),'[',1))

# this was included only as an experiment to see if results could be retained using just one sample
# xspc <- xspc.first
# xgnc <- xgnc.first

# drop rare bugs (show up in less than 25% of subjects)
rare.sp.ix <- colMeans(xspc > 0) < .25
rare.gn.ix <- colMeans(xgnc > 0) < .25
xsp <- xsp[,!rare.sp.ix]
xgn <- xgn[,!rare.gn.ix]
xspc <- xspc[,!rare.sp.ix]
xgnc <- xgnc[,!rare.gn.ix]
xspc.first <- xspc.first[,!rare.sp.ix]
xgnc.first <- xgnc.first[,!rare.gn.ix]

xsp.raw <- xsp.raw[,!rare.sp.ix]
xgn.raw <- xgn.raw[,!rare.gn.ix]
xspc.raw <- xspc.raw[,!rare.sp.ix]
xgnc.raw <- xgnc.raw[,!rare.gn.ix]

# gather differentiation indicies
ixc.hc <- mc$Cohort == "Healthy"
ix.hc <- m$Cohort == "Healthy"
ixc.ibsc <- mc$Cohort == "IBS-C"
ix.ibsc <- m$Cohort == "IBS-C"
ixc.ibsd <- mc$Cohort == "IBS-D"
ix.ibsd <- m$Cohort == "IBS-D"
ixc.ibs <- mc$Cohort != "Healthy"
ix.ibs <- m$Cohort != "Healthy"

# difftest <- differentiation.test(xc, mc$Cohort)

sink('map.txt'); cat('#SampleID\t'); write.table(m, sep='\t',quote=F); sink(NULL)
sink('species.txt'); cat('Taxon\t'); write.table(t(xsp), sep='\t',quote=F); sink(NULL)
sink('genus.txt'); cat('Taxon\t'); write.table(t(xgn), sep='\t',quote=F); sink(NULL)
sink('map-collapsed.txt'); cat('#SampleID\t'); write.table(mc, sep='\t',quote=F); sink(NULL)
sink('species-collapsed.txt'); cat('Taxon\t'); write.table(t(xspc), sep='\t',quote=F); sink(NULL)
sink('genus-collapsed.txt'); cat('Taxon\t'); write.table(t(xgnc), sep='\t',quote=F); sink(NULL)
sink('map-collapsed-first-sample-per-subject.txt'); cat('#SampleID\t'); write.table(mc, sep='\t',quote=F); sink(NULL)
sink('species-collapsed-first-sample-per-subject.txt'); cat('Taxon\t'); write.table(t(xspc.first), sep='\t',quote=F); sink(NULL)
sink('genus-collapsed-first-sample-per-subject.txt'); cat('Taxon\t'); write.table(t(xgnc.first), sep='\t',quote=F); sink(NULL)

