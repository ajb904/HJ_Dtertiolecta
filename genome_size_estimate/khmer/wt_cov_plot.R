library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

infile <- args[1]
k=20

print(paste0('File = ', infile))

summary_file <- gsub('.txt', '_hist.csv', infile, fixed = T)

covdata <- read.table(infile, header=F)
colnames(covdata) <- c('Read', 'Median', 'Mean', 'SD', 'ReadLength')

#Calculate estimated genome size from mode coverage and total sequence length
total.length <- sum(covdata$ReadLength)
print(paste0('Total sequence length = ', total.length, ' bp'))

total.length.kmers <- sum( covdata$ReadLength - k + 1)
print(paste0('Total sequence length in kmers = ', total.length.kmers))

#Based on previous histogram, the main peak is at <25 coverage. Use this to get number of 'single copy' reads
single.copy <- subset(covdata, Median < 25)
single.copy.length <- sum(single.copy$ReadLength)
print(paste0('Single copy sequence length = ', single.copy.length, ' bp'))

single.copy.length.kmers <- sum( single.copy$ReadLength - k + 1)
print(paste0('Single copy sequence length in kmers = ', single.copy.length.kmers))

breaks <- c(1:100, seq(200, 10000, 100), seq(20000,70000, 10000))
histdata <- hist(covdata$Median, plot = FALSE, breaks=breaks)
mode.cov <- histdata$breaks[order(histdata$counts, decreasing=T)][1]
print(paste0('Mode coverage = ', mode.cov))

est.genome.size <- (single.copy.length.kmers/mode.cov) / 1000000
print(paste0('Estimated genome size = ', round(est.genome.size, 2), ' Mb'))


#Plot Histogram (only for values < 100, otherwise plot is uninformative)
plotdata <- data.frame(breaks=histdata$breaks[1:length(breaks)-1], counts=histdata$counts, density=histdata$density)

pngfile <- gsub('txt', 'png', infile)

png(pngfile, width=1440, height=720)

p <- ggplot(plotdata, aes(x=breaks, y=counts)) + geom_bar(stat='identity') + scale_x_continuous(limits=c(0,100))
print(p)

dev.off()

#Write out histogram data in the same format as khmer's calc_median_dist.py script, so that we can check
#genome size estimate with khmer recipe 5
#Format is median, count, cumulative count, cumulative percentage

breaks=seq(1, 65535, 1)
fullhistdata <- hist(covdata$Median, plot = FALSE, breaks=breaks)
fullplotdata <- data.frame(med=fullhistdata$breaks[1:length(breaks)-1],
                           count=fullhistdata$counts)
fullplotdata['cumulative'] <- cumsum(fullplotdata$count)
fullplotdata['percent'] <- round(fullplotdata$cumulative / sum(fullplotdata$count) * 100, 3)
write.table(fullplotdata, summary_file, sep = ' ', row.names = F, col.names = F)