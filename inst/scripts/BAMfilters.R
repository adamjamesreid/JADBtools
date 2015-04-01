require(Rsamtools)
require(rtracklayer)

#GenomicFiles


######### Get BAM of isize > 1kb ##################
what <- c("isize", 'mapq')
flag <- scanBamFlag(isUnmappedQuery = FALSE)
param <- ScanBamParam(what = what)


#f <- dir(pattern='.bam$')[3]

filter <- FilterRules(list(
    size = function(x) { 
        #message(x$isize, ' ', '\n')
        w <- abs(x$isize) >= 700
        w[is.na(w)] <- FALSE
        return(w)
    },
    qual = function(x) { 
        w <- x$mapq >= 10
        w[is.na(w)] <- FALSE
        return(w)
    }
))

flt <- function(f) {
    message('Started: ', f)
    filterBam(
        f,  gsub('.bam$', '_mapq10_isize_more700bp.bam', f), 
        filter=filter, indexDestination=TRUE,
        param=param
    )
    message('Done: ', f)
}

require(parallel)
outfiles1 <- mclapply(dir(pattern='*.bam$'), flt, mc.cores = 2)


outfiles1 <- dir(pattern='mapq10_isize_more700bp.bam$')
dir.create("longtags_mapq10")
file.rename(outfiles1, file.path('longtags_mapq10', outfiles1))
setwd("longtags_mapq10")

sapply(outfiles1, function(x) sortBam(x, gsub('.bam$', '_sorted', x)) )
sapply(dir(pattern='*.bam$'), indexBam )

require(GenomicAlignments)

# ex <- function(ff) {
#     message(ff)
#     gr <- readGAlignmentPairsFromBam(ff, use.names=TRUE)
#     message('Pairs: ', length(gr))
#     export.bw(coverage(gr), gsub('.bam$', '_tagsCoverage.bw', ff) )
#     export.bed(unlist(gr), gsub('.bam$', '_tags.bed', ff) )
#     export.bed( as(gr, 'GRanges'), gsub('.bam$', '_junctions.bed', ff), trackLine=new("TrackLine", name='junctions'))
# } 
# sapply(outfiles1, ex)


ex2 <- function(ff) {
    message(ff)
    gl <- readGAlignmentsList(ff, use.names=FALSE)
    gl <- gl[elementLengths(gl)>1]
    
    #inserts2 <- unique(unlist(range(grglist(gl, ignore.strand=TRUE)), use.names=FALSE))
    inserts <- unique(granges(gl, ignore.strand=TRUE))
    tags <- unique(granges(unlist(gl)))
    
    message('Pairs: ', length(inserts))
    export.bw(coverage(tags), gsub('.bam$', '_tagsCoverage.bw', ff) )
    export.bed(tags, gsub('.bam$', '_tags.bed', ff) )
    export.bed( inserts, gsub('.bam$', '_inserts.bed', ff))
    export.bed( inserts, gsub('.bam$', '_junctions.bed', ff), trackLine=new("TrackLine", name='junctions'))
    
    return(gl)
} 

out <- sapply(dir(patt='WQS.*sorted.bam$'), ex2)


for f in *sorted.bam; do 
echo "-> $f" >> flagstat.txt
samtools flagstat $f >> flagstat.txt; 
done;


######### Get BAM with SA tagse ##################


getSAgaps <- function(f) { 
    what <- c("isize", 'mapq')
    param <- ScanBamParam(what = what, tag = 'SA')
    #b <- scanBam(f, param=param)[[1]]
    
    message(f)
    
    #sap <- readGAlignmentPairsFromBam(f, param=param)
    ga <- readGAlignmentsFromBam(f, param=param)
    gr <- as(ga, 'GRanges')
    sa <- gr[!is.na(gr$SA)]
    seqinfo(sa) <- SeqinfoForBSGenome('ce10')
    
    
    split1 <- sapply(strsplit( sa$SA, ';'), '[[', 1)
    split <- do.call(rbind, strsplit( split1, ','))
    width <- as.numeric( gsub('.*([0-9]{2})M.*', '\\1', split[,4]) )
    grng2 <- GRanges( seqnames =  split[,1], 
                      ranges=IRanges(start = as.integer(split[,2]), width = width), 
                      strand = split[,3], seqinfo=SeqinfoForBSGenome('ce10'),
                      cigar=split[,4], mapq=as.integer(split[,5]))	
    
    flt <- (seqnames(grng2) == seqnames(sa)) & (strand(grng2) == strand(sa)) & (grng2$mapq >= 10) & (sa$mapq >= 10)
    #flt <- (seqnames(grng2) == seqnames(sa)) & (grng2$mapq >= 10) & (sa$mapq >= 10)
    
    left  <- resize(sa[flt], 1, fix='center'); strand(left) <- '*'
    right <- resize(grng2[flt], 1, fix='center'); strand(right) <- '*'
    
    end <- start(left) < start(right) 
    end(left[end]) <- start(right[end])
    start(left[!end]) <- start(right[!end])
    
    tags <- c(sa[flt,2][width(left) >= 1000],grng2[flt,2][width(left) >= 1000])
    out <-  left[width(left) >= 1000]
    
    pdf(gsub('.bam$', 'isize-to-SA.pdf', f)); 
    plot(x=width(out), y=abs(out$isize), xlab='split read distance (SA tag)', ylab='pair tag-to-tag distance (isize)'); 
    dev.off()
    
    export.bw( coverage(tags), gsub('.bam$', 'SAmore1kb_tagsCoverage.bw', f) )
    export.bed( tags, gsub('.bam$', 'SAmore1kb_tags.bed', f) )
    export.bed( out, gsub('.bam$', 'SAmore1kb_junctions_spider.bed', f), trackLine=new("TrackLine", name='junctions'))
    export.bed( out, gsub('.bam$', 'SAmore1kb_junctions_BED.bed', f))
    
    message(f, ': ', length(out), ' pairs')
    return(out)
}

require(parallel)
SA <- mclapply(dir(pattern='*.bam$'), getSAgaps, mc.cores = 8)
names(SA) <- dir(pattern='*.bam$')

ok <- SA[c(1, 4:7)]
for (i in 1:length(ok)) {
    f <- names(ok)[i]
    out <- ok[[i]]
    message(f)
    png(gsub('.bam$', 'isize-to-SA.png', f)); 
    plot(x=width(out), y=abs(out$isize), xlab='split read distance (SA tag)', 
         ylab='pair tag-to-tag distance (isize)', main=gsub('.bam$', '', f)); 
    dev.off()
}


return(sa)

# SA -  Other canonical alignments in a chimeric alignment, in the format of:
#       (rname,pos,strand,CIGAR,mapQ,NM ;)+. Each element in the semi-colon delimited
#       list represents a part of the chimeric alignment. Conventionally, at a supplementary
#       line, the first element points to the primary line.

filter <- FilterRules(list(
    sa = function(x) { 
        #message(x$isize, ' ', '\n')
        !is.na(x$tag$SA)
    }
))

flt <- function(f) {
    filterBam(
        f,  gsub('.bam$', '_SA_tags.bam', f), 
        filter=filter, indexDestination=TRUE,
        param=param
    )
}
outfiles2 <- sapply(dir(pattern='*2.bam$'), flt)
 





gl <- readGAlignmentsList(ff, use.names=TRUE)

RIPSeeker::galp2gal


W <- rbind( do.call(rbind, end(gl)), do.call(rbind, start(gl)) )

chr <- seqnames(unlist(gl)[rownames(W)])
apply(W, 1, max)


