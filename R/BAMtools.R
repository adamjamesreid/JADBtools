#' Calculates delation candidates and outputs them as BAM
#' 
#' @param ID the vector of ContactExpIDs
#' @param minsize minimum size of delation
#' @param download decides if file will be downloaded to local file system
#'   
#' @return GRanges with delation candidats
#' 
#' @author Przemyslaw Stempor
#' 
#' @family RNAseq
#' @export
#' 
#' @examples
#' #d1 <- bam2del('RC025')
#' #d2 <- bam2del('RC037')
#' #export.bed(intersect(d1, d2), 'RC025andRC037_del_more100)intersect.bed')

bam2del <- function(ID, minsize=100L, download=TRUE) {
    fp <- getFilePath(ID, 'bam', 'aligned')
    if(download) {
        tmp <- tempfile()
        download.file(fp, tmp)
        bam <- BamFile(tmp)
    } else {
        bam <- BamFile(fp)
    }
    
    message('Calculating coverage')
    c <- coverage(bam)
    
    message('Processing')
    gr <- as(c==0, 'GRanges')
    gr <- gr[gr$score]
    gr$score <- NULL
    
    out <- gr[width(gr) > minsize]
    out_name <- paste0(paste(rbeads:::ParseName(basename(fp))[c(
        'Factor', 'Antibody', 'ExtractID', 'Crosslinker', 'Strain',  'Stage', 'ContactExpID'
    )], collapse='_'), 'del_more', minsize, '.bed')
    export.bed(out, out_name)
    
    return(out)
}