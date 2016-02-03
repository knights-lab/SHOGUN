observed_staggered = read.table('./results/val-mock-community-2015-12-28/mock_community_staggered.csv', quote='"', header=T, sep=",", na.strings='')
expected_staggered = read.table('./results/val-mock-community-2015-12-28/mock_community_staggered_expected.csv', quote='"', header=T, sep=",")

observed_even = read.table('./results/val-mock-community-2015-12-28/mock_community_even.csv', quote='"', header=T, sep=",", na.strings='')
expected_even = read.table('./results/val-mock-community-2015-12-28/mock_community_even_expected.csv', quote='"', header=T, sep=",")

observed_simulated = read.table('./results/val-sim-12-30-2015/simulated_100000000.csv', quote='"', header=T, sep=",", na.strings='')
expected_simulated = read.table('./results/val-sim-12-30-2015/expected_counts.csv', quote='"', header=T, sep=",", na.strings='')

gene_lengths = read.table('./results/genome-lengths-12-31-2015/database_lengths.csv', quote='"', header=T, sep=",", na.strings='')

# columns to paste together
col_names <- c('domain', 'genus', 'species')

collapse = function(data, cols, name='name') {
    # create a new column `x` with the three columns collapsed together
    data$name <- apply( data[ , cols ] , 1 , paste , collapse = "-" )
    
    # remove the unnecessary rows
    data <- data[ , !( names( data ) %in% cols ) ]
}

make_correlation_plots = function(truth, observed, gene_lengths, name) {
    gene_lengths = collapse(gene_lengths, col_names)
    truth = collapse(truth, col_names)
    observed = collapse(observed, col_names)
    
    all = merge(truth, observed, by='name')
    all = merge(all, gene_lengths, by='name')
    
    # all$count.y = all$count.y/all$length
    
    all$relative.x = all$count.x/sum(all$count.x)
    all$relative.y = all$count.y/sum(all$count.y)
    
    pdf(paste('docs/', 'gene_length_', name, '.pdf', sep = ''))
    plot(log(all$length), (all$relative.x - all$relative.y), xlab='Number of BPs in Database', ylab='Error', main='Gene Length Bias')
    abline(0, 0, col='red')
    dev.off()
    
    pdf(paste('docs/', 'ratio_', name, '.pdf', sep = ''))
    plot(all$relative.x, all$relative.y, xlab='Expected', ylab='Observed', main='Ratio of Observed/Expected')
    abline(0, 1, col='red')
    dev.off()
    
    cor(all$relative.x, all$relative.y)
}

make_correlation_plots(expected_staggered, observed_staggered, gene_lengths, 'expected_staggered_species')
make_correlation_plots(expected_even, observed_even, gene_lengths, 'expected_even_species')
make_correlation_plots(observed_simulated, expected_simulated, gene_lengths, 'expected_simulated_species')

make_correlation_genus_plots = function(truth, observed, gene_lengths, name) {
    gene_lengths = tapply(gene_lengths$length, gene_lengths$genus, FUN=sum)
    truth = tapply(truth$count, truth$genus, FUN=sum)
    observed = tapply(observed$count, observed$genus, FUN=sum)
    
    gene_lengths = data.frame('length'=as.integer(gene_lengths), 'genus'=names(gene_lengths))
    truth = data.frame('count'=as.integer(truth), 'genus'=names(truth))
    observed = data.frame('count'=as.integer(observed), 'genus'=names(observed))
    
    all = merge(truth, observed, by='genus')
    all = merge(all, gene_lengths, by='genus')
    
    # all$count.y = all$count.y/all$length
    
    all$relative.x = all$count.x/sum(all$count.x)
    all$relative.y = all$count.y/sum(all$count.y)
    
    pdf(paste('docs/', 'gene_length_', name, '.pdf', sep = ''))
    plot(log(all$length), (all$relative.x - all$relative.y), xlab='Number of BPs in Database', ylab='Error', main='Genus Gene Length Bias')
    abline(0, 0, col='red')
    dev.off()
    pdf(paste('docs/', 'ratio_', name, '.pdf', sep = ''))
    plot(all$relative.x, all$relative.y, xlab='Expected', ylab='Observed', main='Genus ratio of Observed/Expected')
    abline(0, 1, col='red')
    dev.off()
    cor(all$relative.x, all$relative.y)
}

make_correlation_genus_plots(expected_staggered, observed_staggered, gene_lengths, 'expected_staggered_genus')
make_correlation_genus_plots(expected_even, observed_even, gene_lengths, 'expected_even_genus')
make_correlation_genus_plots(observed_simulated, expected_simulated, gene_lengths, 'expected_simulated_genus')
