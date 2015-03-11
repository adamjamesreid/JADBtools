function(f)
    
    
    #! Rscript
    # Usage: Rscript -e "devtools::source_gist('c977a360a702e2b3ad10')"
    
    message('Running PE insert coverage - ', format(Sys.time(), "%a %b %d %X %Y"))
args <- commandArgs(TRUE)

quality_cutoff <- 0

require(Rsamtools)
require(rtracklayer)

what <- c("isize")
flag <- scanBamFlag(isUnmappedQuery = FALSE)
param <- ScanBamParam(what = what)


getisize <- function(f) {
    
    what <- c("isize")
    param <- ScanBamParam(what = what)
    
    message(basename(f))
    
    i <- scanBam(f, param=param)[[1]]$isize
    n <- na.omit(i)
    
    message('aligned: ', length(n)/length(i), '\n')
    
    p <- n[n>0]
    
    pdf(paste0('ReadDist_', basename(f), '.pdf'))
    hist(p, freq=FALSE, main=basename(f), xlab='')
    lines(density(p), col='red')
    tt <- summary(p)
    title(sub=paste(names(tt), tt, sep=': ', collapse = '; '), cex=0.8)
    dev.off()
    
    return(p)
} 

out <- sapply(dir(pattern='.bam$') , getisize)

q()    grng <- GRanges( seqnames =  rname, ranges=IRanges(start = pos, width = isize), strand = strand, seqinfo=SeqinfoForBSGenome('ce10') )
grng <- trim(grng)
message("All ranges: ", length(grng))
cover <- coverage(grng)
export.bw(cover, gsub('.bam$', '_insert_coverage.bw', basename(f)))

lg <- (a$mapq[ a$isize != 0 ] >= quality_cutoff)
grng <- grng[lg]
message("Mapq10 ranges: ", length(grng))
cover <- coverage(grng)
export.bw(cover, gsub('.bam$', '_mapq10_insert_coverage.bw', basename(f)))

grng <- unique(grng)
message("Unique and mapq10 ranges: ", length(grng))
cover <- coverage(grng)
export.bw(cover, gsub('.bam$', '_uniq_mapq10_insert_coverage.bw', basename(f)))