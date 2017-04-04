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
bamStats <- function(f, aln=NULL) {
    
    data(rrnamodel)
    
    if(!is.list(a)) {
        what <- c("rname", "strand", "pos", "mapq", "qwidth")
        flag <- scanBamFlag(isUnmappedQuery = FALSE)
        param <- ScanBamParam(flag = flag, simpleCigar = FALSE, what = what)
        a <- scanBam(f, param = param)[[1]]
    } esle {
        a <- aln
    }
    
    lg <- (a$mapq >= 10L)
    
    grng <- GRanges(
        seqnames = a$rname, ranges = IRanges(a$pos, width = a$qwidth), 
        strand = a$strand
    )
    
    all <- countBam(f)
    r <- all$records
    a <- length(grng)
    q <- sum(lg)
    u <- length(unique(grng[lg]))
    rr <- sum(grng %over% rrnamodel)
    
    out <- sprintf(
        'all=%.2fM, aligned=%.2fM[%.0f%%], mapq10=%.2fM[%.0f%%], unique10=%.2fM[%.0f%%], rRNA=%.2fM[%.0f%%]',
        r/10^6, a/10^6, (a/r)*100, q/10^6, (q/r)*100, u/10^6, (u/r)*100, rr/10^6, (rr/r)*100
    )
    
    return(out)

}