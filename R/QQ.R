QQ <- function(f, PE=TRUE, type=c('insert', 'read', 'tag'), stranded=TRUE, export=c('all', 'q10', 'uniq'), strandsRev=TRUE, stats=FALSE) {

    library(parallel)
    library(Rsamtools)
    library(rtracklayer)
    
    type  <- match.arg(type)
    #message('=> ', f)
    
    b <- BamFile(f[[1]])
    if(!file.exists(paste0(f, '.bai'))) {
        message('Indexing ', f)
        indexBam(b)
        message(f, ' indexed')
    }
    
    lst <- as.list(GRanges(seqinfo(b)))
    
    getGR <- function(chr, isFirstMateRead=TRUE) {
        if (PE & type == 'insert') {
            flag <- scanBamFlag(isUnmappedQuery = FALSE, isProperPair = TRUE, isMinusStrand = FALSE, isMateMinusStrand = TRUE, isFirstMateRead = isFirstMateRead)
            what <- c("rname", "pos", "isize", 'strand', 'mapq')
            param <- ScanBamParam(what = what, flag = flag, which = chr)
            a <- scanBam(f, param = param)[[1]]
            isize <- a$isize[ a$isize > 0 ]
            pos <- a$pos[ a$isize > 0 ]
            rname <- a$rname[ a$isize > 0 ]
            strand <- a$strand[ a$isize > 0 ]
            mapq <- a$mapq[ a$isize > 0 ]
            grng <- GRanges( seqnames =  rname, ranges = IRanges(start = pos, width = isize), strand = strand )
            grng <- trim(grng)
            grng$mapq <- mapq
        } else {
            what <- c("rname", "strand", "pos", "mapq", "qwidth")
            flag <- scanBamFlag(isUnmappedQuery = FALSE, isMinusStrand = isFirstMateRead)
            param <- ScanBamParam(flag = flag, simpleCigar = FALSE, what = what, which = chr)
            a <- scanBam(f, param=param)[[1]]
            grng <- GRanges( seqnames = a$rname, ranges = IRanges(a$pos, width = a$qwidth), strand = a$strand )
            grng <- trim(grng)
            grng$mapq <- a$mapq
        }
        return(grng)
    }
    
   
    REV <- unlist(GRangesList(mclapply(lst, getGR, strandsRev, mc.cores = parallel::detectCores())))
    FWD <- unlist(GRangesList(mclapply(lst, getGR, !strandsRev, mc.cores = parallel::detectCores())))
    
    nf <- (length(FWD) + length(REV))/1e6
    
    if(stats) {
        
        ALL <- c(FWD, REV)
        rRNA <- GRanges(
            c("chrI", "chrI", "chrI", "chrI", "chrM", "chrM", "chrV", "chrV", 
              "chrV", "chrV", "chrV", "chrV", "chrV", "chrV", "chrV", "chrV", 
              "chrV", "chrV", "chrV", "chrV", "chrV", "chrI"
             ), IRanges(c(15062072, 15069269, 15064290, 15064827, 898, 10403, 
                          17120944, 17115879, 17117839, 17118813, 17119989, 
                          17121920, 17122896, 17123872, 17124848, 17125824, 
                          17128380, 17129352, 17130338, 17131314, 17428854, 
                          6169086),
                        c(15063825, 15071022, 15064442, 15068335, 1593, 11354, 
                          17121062, 17115997, 17117957, 17118911, 17120087, 
                          17122038, 17123014, 17123990, 17124966, 17125942, 
                          17128498, 17129470, 17130456, 17131432, 17428971, 
                          6169197)
             )
        )
        
        jobs <- list()
        jobs[[1]] <- mcparallel( countBam(b)$records )
        jobs[[2]] <- mcparallel( length(ALL) )
        jobs[[3]] <- mcparallel( length(REV) )
        jobs[[4]] <- mcparallel( length(FWD) )
        jobs[[5]] <- mcparallel( sum(ALL$mapq >= 10) )
        jobs[[6]] <- mcparallel( length(unique(ALL)) )
        jobs[[7]] <- mcparallel( sum(unique(ALL)$mapq >= 10) )
        jobs[[8]] <- mcparallel( sum(duplicated(ALL)) )
        jobs[[9]] <- mcparallel( sum(ALL %over% rRNA) )
        jobs[[10]] <- mcparallel( sum(seqnames(ALL)=='chrM') )
        ans <- unlist(mccollect(jobs))
        names(ans) <- c('ALL', 'mapped', 'REV', 'FWD', 'q10', 'UNIQ', 'uniqQ10', 'dups', 'rRNA', 'chrM')

       
        st <- rbind(
            'reads [M]'=round(ans/1e6, 2),
            '% of ALL'=round(ans[]/ans[1]*100, 2)
        )
        sst <- paste0(c(round(ans[1]/1e6, 2), round(ans[c(2,5,7)]/ans[1]*100)), collapse='-')
        
        message(f, '\t[',sst,']:')
        print(st)
        
        # message(
        #     f, ' - all=', a,
        #     ' mapped=', nf*1e6,
        #     ' REV=', length(REV),
        #     ' FWD=', length(FWD),
        #     ' q10=', sum(ALL$mapq >= 10),
        #     ' UNIQ=', length(unique(ALL) >= 10),
        #     ' UNIQq10=', sum(unique(ALL)$mapq >= 10)
        # )
        return(st)
    }
    
    message(f, ' - mapped reads: ', nf*1e6)
    

    
    if(stranded == TRUE) {
        jobs <- list()
        jobs[[1]] <- mcparallel( export.bw(coverage(FWD), gsub('.bam$', paste0('_FWD_',type,'_coverage.bw'), basename(f))) )
        jobs[[2]] <- mcparallel( export.bw(coverage(REV), gsub('.bam$', paste0('_REV_',type,'_coverage.bw'), basename(f))) )
        jobs[[3]] <- mcparallel( export.bw(coverage(FWD)/nf, gsub('.bam$', paste0('_FWD_',type,'_coverage_readsPerMillion.bw'), basename(f))) )
        jobs[[4]] <- mcparallel( export.bw(coverage(REV)/nf, gsub('.bam$', paste0('_REV_',type,'_coverage_readsPerMillion.bw'), basename(f))) )
        mccollect(jobs)
    } else {
        ALL <- c(FWD, REV)
        jobs <- list()
        jobs[[1]] <- mcparallel( export.bw(coverage(ALL), gsub('.bam$', paste0('_', type,'_coverage.bw'), basename(f))) )
        jobs[[2]] <- mcparallel( export.bw(coverage(ALL)/nf, gsub('.bam$', paste0('_', type,'_coverage_readsPerMillion.bw'), basename(f))) )
        mccollect(jobs)
    }

}

# mclapply(dir(pattern='.bam$'), QQ, mc.cores = parallel::detectCores())

# library(parallel); mclapply(dir(pattern='.bam$'), QQ, PE=FALSE, type='read', stranded=TRUE, mc.cores = parallel::detectCores())
# library(parallel); mclapply(dir(pattern='.bam$'), QQ, PE=FALSE, type='read', stranded=FALSE, mc.cores = parallel::detectCores())

# rsync -av --progress *bw www-data@jadb:/mnt/ps/TRACKS
#   flag <- scanBamFlag(isUnmappedQuery = FALSE, isMinusStrand = isFirstMateRead)

# sum(unlist(mclapply(as.list(GRanges(seqinfo(b))), function(chr) countBam(b, param=ScanBamParam(which = chr))$records,  mc.cores = detectCores())))

# ST <- mclapply(dir(pattern='.bam$'), QQ, PE=FALSE, stats=TRUE, mc.cores = parallel::detectCores())
# names(ST) <- dir(pattern='.bam$')
# R <- do.call(rbind, lapply(st, function(stst) c(stat[1,], '%'=stat[2,-1])) )
