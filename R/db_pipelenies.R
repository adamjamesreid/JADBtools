#' jadb_ChIPseq
#' 
#' @param IDs Vector of JADB ContactExpIDs
#'   
#' @return List 
#' 
#' @author Przemyslaw Stempor
#' 
#' @family pipelines
#' @export
#' 
#' @examples
#' #
jadb_ChIPseq <- function(ids) {
    
    setwd('/mnt/jadb/DBfile/DBfiles')
    library('rbeads')
    
    message('\t => \t Performing alignment')
    jadb_addAlignedBAM(ids)
    
    message('\t => \t Exporting tracks')
    jadb_addTracksFromBAM(ids)
    
    message('\t => \t Normalising NU')
    addNonUniqueQ10Beads(ids)
    
    message('\t => \t Adding zscored track')
    jadb_addScaledTrack(ids)
    
    message('\t => \t Ruinning FASTQC')
    jadb_processs_sceleton(
        ids, FUN = run_fastqc, 
        Processing = 'FastQC', Scale = 'NA', Resolution = 'NA', filetype_format = 'html',
        format = "bam", processing = "aligned", scale = "NA", resolution='NA'
    )
    
    message('\t => \t Ruinning MACS')
    callPeaksMACS(ids, local = FALSE)
}