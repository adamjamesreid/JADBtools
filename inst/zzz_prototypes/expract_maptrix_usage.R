#' Get data matrixes from BW collections
#' 
#' @param file the path to input FASTTQ
#'   
#' @return NULL; creates 2 FASTQ file in the dirname of file
#' 
#' @author Przemyslaw Stempor
#' 
#' @family shortRNA
#' @export
#' 
#' @example
#' fls <- dir('/Volumes/raid0/_Jcor/BWproject', pattern='*bw', full.names = TRUE)
#' 
get21Uand22G <- function(fsl) {
    
    M <- sapply(fls, function(x) {
        bwf <- BigWigFile(x)
        message(as.integer(seqlengths(bwf) / 1000))
#         z <- runValue(unlist( 
#             summary(bwf, as(seqinfo(bwf), 'GRanges')[1], size = as.integer(seqlengths(bwf)[1] / 1000), asRle = TRUE)
#         , use.names=FALSE))
        z <- unlist(summary(bwf, size = as.integer(seqlengths(bwf) / 1000))[1])$score
        return(z)
    })
    colnames(M) <- basename(colnames(M))
    #writeMat('Mtest1.mat', M=M)
    
}

ml <- lapply( fls, function(f, bin=1000) {
    cat(f, '\n')
    bw <- BigWigFile(f)
    return( as.numeric(unlist(summary(bw, size=seqlengths(bw)/bin), use.names=FALSE)) ) 
})