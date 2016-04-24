QQ <- function(f, type=c('insert', 'read', 'tag'), stranded=TRUE, export=c('all', 'q10', 'uniq')) {

    type  <- match.arg(type)
    message('=> ', f)
    
    b <- BamFile(f[[1]])
    if(!file.exists(paste0(f, '.bai'))) {
        message('Indexing ', f)
        indexBam(b)
        message(f, ' indexed')
    }
    
    lst <- as.list(GRanges(seqinfo(b)))
    
    getGR <- function(chr, isFirstMateRead=TRUE) {
        if (type == 'insert') {
            flag <- scanBamFlag(isUnmappedQuery = FALSE, isProperPair = TRUE, isMinusStrand = FALSE, isMateMinusStrand = TRUE, isFirstMateRead=isFirstMateRead)
            what <- c("rname", "pos", "isize", 'strand', 'mapq')
            param <- ScanBamParam(what = what, flag=flag, which=chr)
            a <- scanBam(f, param=param)[[1]]
            isize <- a$isize[ a$isize > 0 ]
            pos <- a$pos[ a$isize > 0 ]
            rname <- a$rname[ a$isize > 0 ]
            strand <- a$strand[ a$isize > 0 ]
            mapq <- a$mapq[ a$isize > 0 ]
            grng <- GRanges( seqnames =  rname, ranges=IRanges(start = pos, width = isize), strand = strand )
            grng <- trim(grng)
            grng$mapq <- mapq
        }
        return(grng)
    }
    
   
    REV <- unlist(GRangesList(mclapply(lst, getGR, TRUE, mc.cores = parallel::detectCores())))
    FWD <- unlist(GRangesList(mclapply(lst, getGR, FALSE, mc.cores = parallel::detectCores())))
    
    nf <- (length(FWD) + length(REV))/1e6
    
    message(f, ' - mapped reads: ', nf*1e6)
    
    jobs <- list()
    jobs[[1]] <- mcparallel( export.bw(coverage(FWD), gsub('.bam$', '_FWD_insert_coverage.bw', basename(f))) )
    jobs[[2]] <- mcparallel( export.bw(coverage(REV), gsub('.bam$', '_REV_insert_coverage.bw', basename(f))) )
    jobs[[3]] <- mcparallel( export.bw(coverage(FWD)/nf, gsub('.bam$', '_FWD_read_coverage_readsPerMillion.bw', basename(f))) )
    jobs[[4]] <- mcparallel( export.bw(coverage(REV)/nf, gsub('.bam$', '_REV_read_coverage_readsPerMillion.bw', basename(f))) )
    mccollect(jobs)

}

# mclapply(dir(pattern='.bam$'), QQ, mc.cores = parallel::detectCores())