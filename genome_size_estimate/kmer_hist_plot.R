library(ggplot2)

genome.size <- function(jf.file, number.of.reads, k.size, read.length=250, plot.file = paste(jf.file, 'plot.png', sep='.')){
    
    df <- read.delim(jf.file, sep=' ', header=F, col.names = c('kmer.count', 'freq'))
    
    png(plot.file)
    p <- ggplot(df, aes(x=kmer.count, y=freq)) + geom_bar(stat='identity') + scale_x_continuous(limits=c(0,50))
    print(p)
    dev.off()
    
    peak.k <- df[order(df$freq, decreasing = T),][1,'kmer.count']
    
    coverage <- peak.k * read.length / (read.length - k.size + 1)
    
    genome.size <- number.of.reads * read.length / coverage
    
    return(genome.size)

    }


genomes <- data.frame(strain = rep(c('wt','L1','L2'),each=4),
                      num.reads = rep(c(10960620, 4975018, 5830192), each=4),
                      k.size = rep(c(17,21,25,31),3)
                      )

genomes['kmer.file'] <- paste(genomes$strain, '_kmer_counts_', genomes$k.size, '.hist', sep='')

genomes['genome.size'] <- mapply(genome.size, genomes$kmer.file, genomes$num.reads, genomes$k.size)

genomes