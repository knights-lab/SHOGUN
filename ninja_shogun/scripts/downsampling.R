library(ggplot2)

k = 9

spearman = function(m) {
    counts = sapply(1:(length(m)-2), function(x) paste('count.', x, sep=''))
    vals = 1:length(counts)
    for (i in 1:length(counts)) {
        vals[i] = cor(m['count'], m[counts[i]])
    }
    return(vals)
}


precision = function(m) {
    counts = sapply(1:(length(m)-2), function(x) paste('count.', x, sep=''))
    vals = 1:length(counts)
    act = m['count'] > 0
    for (i in 1:length(counts)) {
        pred = m[counts[i]] > 0
        xTab = table(act, pred)
        print(xTab)
        vals[i] = xTab[2,2]/sum(xTab[,2]) # Hit Precision
    }
    return(vals)
}

recall = function(m) {
    counts = sapply(1:(length(m)-2), function(x) paste('count.', x, sep=''))
    vals = 1:length(counts)
    act = m['count'] > 0
    for (i in 1:length(counts)) {
        pred = m[counts[i]] > 0
        print(sum(m[counts[i]] > 0))
        xTab = table(act, pred)
        print(xTab)
        vals[i] <- xTab[2,2]/sum(xTab[2,]) # Hit Recall
    }
    return(vals)
}

data_dir = file.path('results', 'SKTSL-downsampling-2015-12-31-15')
full_stats = read.csv(file.path(data_dir, 'results_full.csv'))

# columns to paste together
col_names <- c('domain', 'genus', 'species')

collapse = function(data, cols, name='name') {
    # create a new column `x` with the three columns collapsed together
    data$name <- apply( data[ , cols ] , 1 , paste , collapse = "-" )
    
    # remove the unnecessary rows
    data <- data[ , !( names( data ) %in% cols ) ]
}

full_stats_species = collapse(full_stats, col_names)
full_stats_species[full_stats_species$count < k, 'count'] = 0

sizes = c(1000, 10000, 100000, 1000000, 10000000)

spearman_m = matrix(NA, nrow=length(sizes), ncol=3)
precision_m = matrix(NA, nrow=length(sizes), ncol=3)
recall_m = matrix(NA, nrow=length(sizes), ncol=3)

for (i in 1:length(sizes)) {
    size = sizes[i]
    temp_dir = sprintf('hits_%d', size)
    filenames = list.files(file.path(data_dir, temp_dir), pattern="*.csv", full.names=TRUE)
    experiments = lapply(filenames, read.csv)
    experiments = lapply(experiments, function(x) {collapse(x, col_names)})
    num = 1
    for (df in experiments) {
        colnames(df)[1] <- paste('count.', num, sep='')
        experiments[[num]] = df
        num = num + 1
    }
    experiments = Reduce(function(x,y) merge(x, y, by='name', all=TRUE), experiments)
    counts = sapply(1:(length(m)-2), function(x) paste('count.', x, sep=''))
    m = merge(experiments, full_stats_species, by='name', all=TRUE)
    m[is.na(m)] = 0
    m[which(rowSums(m[counts]) < log10(size)), counts] = 0
    temp = spearman(m)
    spearman_m[i,] = c(size, mean(temp), sd(temp))
    
    temp = precision(m)
    precision_m[i,] = c(size, mean(temp), sd(temp))
    
    temp = recall(m)
    recall_m[i,] = c(size, mean(temp), sd(temp))
}

col_names = cbind('depth', 'mean', 'sd')

precision_df = as.data.frame(precision_m)
colnames(precision_df) = col_names

recall_df = as.data.frame(recall_m)
colnames(recall_df) = col_names

spearman_df = as.data.frame(spearman_m)
colnames(spearman_df) = col_names

plot_error_bars = function(df, lab) {
    # Use geom_line()+geom_pointrange()
    plot1 = ggplot(df, aes(x=log10(depth), y=mean)) + 
        geom_line()+
        geom_pointrange(aes(ymin=mean-sd, ymax=mean+sd)) +
        labs(title=paste(lab, 'by Depth'), x='Log(Depth)', y=paste('Average', lab))

    ggsave(filename=file.path('docs', paste("log", "depth", lab,"pdf", sep='.')), plot=plot1)
        
    # Use geom_line()+geom_pointrange()
    plot2 = ggplot(df, aes(x=depth, y=mean)) + 
        geom_line()+
        geom_pointrange(aes(ymin=mean-sd, ymax=mean+sd)) +
        labs(title=paste(lab, 'by Depth'), x='Depth', y=paste('Average', lab))
    
    ggsave(filename=file.path('docs', paste("depth", lab,"pdf", sep='.')), plot=plot2)
}

plot_error_bars(precision_df, 'Precision')
plot_error_bars(recall_df, 'Recall')
plot_error_bars(spearman_df, 'Spearman Correlation')
