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
jadb_ChIPseq <- function( 
    ids, 
    steps=c('aln', 'tracks', 'norm', 'log2_norm', 'map0_norm', 'log2_map0_norm', 'map10_uniq', 'fastqc', 'fastqscreen', 'macs', 'meme'), 
    genome='ce11', purge=TRUE, backup=FALSE 
) {
    
    owd <- getwd()
    on.exit(setwd(owd))
    
    setwd(MOUNT)
    library('rbeads')
    
    backupFiles <- function(id) {
        pth <- getFilePath(id, url = FALSE)
        pth <- gsub('/meme_chip', '', pth)
        dirname <- gsub('files', MOUNT, unique(dirname(pth)))
        message(paste(dirname, collapse='\n'))
        if(length(dirname) != 1) stop('mltp dir')
        
        bck <- gsub('DBfiles', 'DBfiles_legacy', dirname)
        dir.create(bck, recursive = TRUE)
        file.copy(file.path(dirname, '/'), bck, recursive = TRUE)
    }
    if(backup==TRUE) {
        message("backing up files")
        backupFiles(ids)
    }
    
    if(purge) jadb_purge_exp(ids)
    
    if (getAnno(ids, anno = 'Factor', EXTABLE='labexperiment') == 'Input') {
        steps <- c('aln', 'tracks', 'fastqc', 'fastqscreen')
    }
    
    
    
    ### ALN ###
    
    if('aln' %in% steps) {
        message('\t => \t Performing alignment')
        jadb_addAlignedBAM(ids, genome=genome)
    }

    if('tracks' %in% steps) {
        message('\t => \t Exporting tracks')
        jadb_addTracksFromBAM(ids, genome=genome)
    }
    
    ### QC ###
    
    if('fastqc' %in% steps) {
        message('\t => \t Ruinning FASTQC')
        jadb_processs_sceleton(
            ids, FUN = run_fastqc, 
            Processing = 'FastQC', Scale = 'NA', Resolution = 'NA', filetype_format = 'html',
            format = "bam", processing = "aligned", scale = "NA", resolution='NA', prefix_output = 'Q'
        )
    }
    
    if('fastqscreen' %in% steps) {
        message('\t => \t Ruinning fastqscreen')
        jadb_processs_sceleton(
            ids, FUN = run_fastq_screen, 
            Processing = 'fastqSCREEN', Scale = 'NA', Resolution = 'NA', filetype_format = 'html',
            format = "txt.gz", processing = "raw", scale = "NA", resolution='NA', prefix_output = 'Q'
        )
    }
    
    ### NORM ###
    
    if('norm' %in% steps) {
        message('\t => \t Normalising NU')
        addNonUniqueQ10Beads(ids, genome=genome)
        
        message('\t => \t Adding zscored track')
        jadb_addScaledTrack(ids, genome=genome)
    }
    
    if('log2_norm' %in% steps) {
        message('\t => \t Adding log2 track')
        jadb_addScaledTrack(ids, genome=genome, scale = 'log2')
    }
    
    ### NORM MAPQ0 ###
    
    if('map0_norm' %in% steps) {
        message('\t => \t Normalising NU MAPQ0')
        addNonUniqueNQBeads(ids, genome=genome)
        
        message('\t => \t Adding zscored track MAPQ0')
        jadb_addScaledTrack(ids, genome=genome, input = 'BEADSNQNU')
    }
    
    if('log2_map0_norm' %in% steps) {
        message('\t => \t Adding log2 track MAPQ0')
        jadb_addScaledTrack(ids, genome=genome, input = 'BEADSNQNU', scale = 'log2')
    }
    
    ### NORM UNIQ ###
    
    if('map10_uniq' %in% steps) {
        
        message('\t => \t Normalising UNIQ MAPQ10')
        addUniqueQ10Beads(ids, genome=genome)
        
        message('\t => \t Adding zscored track UNIQ MAPQ10')
        jadb_addScaledTrack(ids, genome=genome, input = 'BEADSQ10UNIQ')
        
        message('\t => \t Adding log2 track UNIQ MAPQ10')
        jadb_addScaledTrack(ids, genome=genome, input = 'BEADSQ10UNIQ', scale = 'log2')
        
    }
    
    
    ### TOOLS ###
    
    if('macs' %in% steps) {
        message('\t => \t Ruinning MACS')
        callPeaksMACS(ids, local = FALSE, genome=genome)
    }

    if('meme' %in% steps) {
        message('\t => \t Ruinning MEME')
        jadb_processs_sceleton(
            ids, FUN = run_meme_chip, 
            Processing = 'MEMEchip', Scale = 'NA', Resolution = '500bp', filetype_format = 'html',
            format = "bed", processing = "summits", scale = "MACS", resolution='q01',
            genome = genome, uniq='NA', rename_output=FALSE, prefix_output = 'M'
        )
    }
    
    message('Pipeline finished!')

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
jadb_RNAseq <- function( ids, steps=c('trim', 'aln', 'tracks', 'rpm', 'fastqc', 'star'), ... ) {
    
    require(magrittr)
    setwd(MOUNT)
    
    if('trim' %in% steps) {
        message('\t => \t Trimming')
        jadb_trimm_fq(ids)
    }
    
    #jadb_trimm_fq
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
    
    if('star' %in% steps) {
        message('\t => \t Ruinning STAR')
        jadb_processs_sceleton(
            ids, FUN = run_star, 
            Processing = 'STAR', Scale = 'NA', genome = 'ce11', filetype_format = 'bam',
            format = "txt.gz", processing = "raw", scale = "NA", resolution='NA', prefix_output = 'S'
        )
    }
    
    message('Pipeline finished!') 
    
    #fastq Screen 
    
    # strand specific
    
}