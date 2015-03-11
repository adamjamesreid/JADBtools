require('rtracklayer')
require('magrittr')

bed1 <- "/Volumes/data/_HOT_jointAWG_worm_fly_human/worm_Ce10_noDoubleCounting/CFP1_peakCalls/antiGFP^ab_rc03^F^Cfp1-GFP^L3_summits.bed.gz" 
peaks1 <- bed1 %>% file %>% import.bed

require(Rsamtools)
require(rtracklayer)
what <- c("rname", "pos", "isize", 'strand', 'mapq')
quality_cutoff <- 10
flag <- scanBamFlag(isUnmappedQuery = FALSE, isProperPair = TRUE, isMinusStrand = FALSE, isMateMinusStrand = TRUE)
param <- ScanBamParam(what = what, flag=flag)

require(parallel)

out <- sapply( dir(pattern='.bam$'), function(f) {
    message(basename(f), '\n')
    
    a <- scanBam(f, param=param)[[1]]
    isize <- a$isize[ a$isize != 0 ]
    pos <- a$pos[ a$isize != 0 ]
    rname <- a$rname[ a$isize != 0 ]
    strand <- a$strand[ a$isize != 0 ]
    
    pdf(paste0('ReadDist_', basename(f), '.pdf'))
    hist(isize, freq=FALSE, main=basename(f), xlab='')
    lines(density(isize), col='red')
    tt <- summary(isize)
    title(sub=paste(names(tt), tt, sep=': ', collapse = '; '), cex=0.8)
    dev.off()
    
    grng <- GRanges( seqnames =  rname, ranges=IRanges(start = pos, width = isize), strand = strand, seqinfo=SeqinfoForBSGenome('ce10'), mapq=a$mapq[ a$isize != 0 ] )
    grng <- trim(grng)
    
    out_temp <- grng
    
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
    
    return(out_temp)
})

ans <- sapply(out, subsetByOverlaps, peaks1)

M0 <- sapply(out, function(x) summary(width(x)))
M <- sapply(ans, function(x) summary(width(x)))

t(t( sapply(out, function(x) length(x))/10^6 ))

t(t(  sapply(ans, function(x) length(x))/10^6 ))

Map(function(gr, n) {
    message(n)
    pdf(paste0('CPF1peaks_rc03_ReadDist_', n, '.pdf'))
    hist(width(gr), freq=FALSE, main=n, xlab='')
    lines(density(width(gr)), col='red')
    tt <- summary(width(gr))
    title(sub=paste(names(tt), tt, sep=': ', collapse = '; '), cex=0.8)
    dev.off()
}, ans, names(ans))