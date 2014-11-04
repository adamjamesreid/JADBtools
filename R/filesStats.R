#' Returms bam file stats
#' 
#' String is formated as: 
#' all=%d, aligned=%d[%.2f%%], mapq10=%d[%.2f%%], unique10=%d[%.2f%%]
#' 
#' @param f URL or path to BAM file  
#'   
#' @return stats string
#' 
#' @author Przemyslaw Stempor
#' 
#' @family stats
#' @export
#' 
bamStats <- function(f) {
    
    
    what <- c("rname", "strand", "pos", "mapq", "qwidth")
    flag <- scanBamFlag(isUnmappedQuery = FALSE)
    param <- ScanBamParam(flag = flag, simpleCigar = FALSE, what = what)
    
    a <- scanBam(f, param = param)[[1]]
    lg <- (a$mapq >= 10L)
    
    grng <- GRanges(
        seqnames = a$rname, ranges = IRanges(a$pos, width = a$qwidth), 
        strand = a$strand, seqinfo=SeqinfoForBSGenome('ce10')
    )
    
    all <- countBam(f)
    r <- all$records
    a <- length(grng)
    q <- sum(lg)
    u <- length(unique(grng[lg]))
    
    out <- sprintf(
        'all=%d, aligned=%d[%.2f%%], mapq10=%d[%.2f%%], unique10=%d[%.2f%%]',
        r, a, (a/r)*100, q, (q/r)*100, u, (u/r)*100
    )
    
    return(out)

}