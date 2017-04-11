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
jadb_ChIPseq <- function( ids, steps=c('aln', 'tracks', 'norm', 'fastqc', 'macs') ) {
    
    setwd('/mnt/jadb/DBfile/DBfiles')
    library('rbeads')
    
    if (getAnno(ids, anno = 'Factor', EXTABLE='labexperiment') == 'Input') {
        steps <- c('aln', 'tracks', 'fastqc')
    }
    
    if('aln' %in% steps) {
        message('\t => \t Performing alignment')
        jadb_addAlignedBAM(ids)
    }

    if('tracks' %in% steps) {
        message('\t => \t Exporting tracks')
        jadb_addTracksFromBAM(ids)
    }
    
    if('norm' %in% steps) {
        message('\t => \t Normalising NU')
        addNonUniqueQ10Beads(ids)
        
        message('\t => \t Adding zscored track')
        jadb_addScaledTrack(ids)
    }
    
    if('fastqc' %in% steps) {
        message('\t => \t Ruinning FASTQC')
        jadb_processs_sceleton(
            ids, FUN = run_fastqc, 
            Processing = 'FastQC', Scale = 'NA', Resolution = 'NA', filetype_format = 'html',
            format = "bam", processing = "aligned", scale = "NA", resolution='NA'
        )
    }
    
    if('macs' %in% steps) {
        message('\t => \t Ruinning MACS')
        callPeaksMACS(ids, local = FALSE)
    }

    #fastq Screen 
    
    #meme-chip
    

    

    

}


#' jadb_RNAseq
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
jadb_RNAseq <- function( ids, steps=c('trim', 'aln', 'tracks', 'rpm', 'fastqc') ) {
    
    setwd('/mnt/jadb/DBfile/DBfiles')
    
    if('trim' %in% steps) {
        message('\t => \t Trimming')
        jadb_trimm_fq(ids)
    }
    
    jadb_trimm_fq
    if('aln' %in% steps) {
        message('\t => \t Performing alignment')
        jadb_align_rnaseq(ids)
    }
    
    if('tracks' %in% steps) {
        message('\t => \t Exporting tracks')
        db_mran_seq_align_process(ids)
    }
    
    if('rpm' %in% steps) {
        db_mran_rpm_track(ids)
    }
    
    if('fastqc' %in% steps) {
        message('\t => \t Ruinning FASTQC')
        jadb_processs_sceleton(
            ids, FUN = run_fastqc, 
            Processing = 'FastQC', Scale = 'NA', Resolution = 'NA', filetype_format = 'html',
            format = "bam", processing = "aligned", scale = "NA", resolution='NA'
        )
    }
    
 
    
    #fastq Screen 
    
    # strand specific
    
    
    
    
    
    
}