#' Extract 21U and 22G short RNA specied from FASTQ file
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
get21Uand22G <- function(file) {
    
    require(ShortRead)
    
    r1 <- readFastq(dirname(file), pattern = basename(file))
    seq <- r1@sread
    
    #21U
    #seq[subseq(seq,1,1) == 'T' & nchar(seq) >= 18 & nchar(seq) <= 30]
    r21U <- r1[subseq(seq,1,1) == 'T' & nchar(seq) >= 18 & nchar(seq) <= 30]
    #22G
    #seq[subseq(seq,1,1) == 'G' & nchar(seq) >= 18 & nchar(seq) <= 30]
    r22G <- r1[subseq(seq,1,1) == 'G' & nchar(seq) >= 18 & nchar(seq) <= 30]
    
    oldwd <- getwd(); setwd(dirname(file));
    writeFastq(r21U, gsub(".(fq|fastq.gz)", "_r21U.fastq.gz", basename(file)))
    writeFastq(r22G, gsub(".(fq|fastq.gz)", "_r22G.fastq.gz", basename(file)))
    setwd(oldwd)
    
    message('21U: ', length(r21U@sread), '\n22G: ', length(r22G@sread))

}

#' Extract 21U and 22G short RNA specied from FASTQ file
#' 
#' @param IDs IDs of trimmed files
#'   
#' @return NULL; creates 2 FASTQ file in the dirname of file
#' 
#' @author Przemyslaw Stempor
#' 
#' @family shortRNA
#' @export
#' 
get21Uand22GfromIDs <- function(IDs) {
    
    sapply(IDs, function(x) {
        fq <- getFilePath(x, processing = 'trimmed')
        download.file(fq, basename(fq))
        message(basename(fq))
        get21Uand22G(basename(fq))
    })
    
}
