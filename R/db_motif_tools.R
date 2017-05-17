run_meme_chip <- function(file, interperor='bash') {
    library(BSgenome.Celegans.UCSC.ce11)
    
    summits <- import.bed(file)
    seqinfo(summits) <- seqinfo(Celegans)[seqlevels(summits)]
    
    seq <- getSeq(Celegans, trim(resize(summits, 500, fix = 'center')))
    names(seq) <- paste(summits)
    seq <- seq[width(seq)==500]
    writeXStringSet(seq, 'summit_sequence_500bp.fa')
    message('FASTA saved to: ', getwd(), '/summit_sequence_500bp.fa', ' [', length(seq), ' summits/sequences]')
    
    
    databases <- paste(
        "-db ~/meme/db/motif_databases/WORM/uniprobe_worm.meme",
        "-db ~/meme/db/motif_databases/EUKARYOTE/hallikas2006.meme",
        "-db ~/meme/db/motif_databases/EUKARYOTE/homeodomain.meme",
        "-db ~/meme/db/motif_databases/EUKARYOTE/jolma2013.meme",
        "-db ~/meme/db/motif_databases/EUKARYOTE/macisaac_theme.v1.meme",
        "-db ~/meme/db/motif_databases/EUKARYOTE/prodoric.meme", 
        "-db ~/meme/db/motif_databases/EUKARYOTE/regtransbase.meme",
        "-db ~/meme/db/motif_databases/EUKARYOTE/SwissRegulon_human_and_mouse.meme",
        "-db ~/meme/db/motif_databases/EUKARYOTE/wei2010_human_mws.meme",
        "-db ~/meme/db/motif_databases/EUKARYOTE/wei2010_mouse_mws.meme",
        "-db ~/meme/db/motif_databases/EUKARYOTE/wei2010_mouse_pbm.meme",
        "-db ~/meme/db/motif_databases/EUKARYOTE/zhao2011.meme",
        "-db ~/meme/db/motif_databases/MOUSE/uniprobe_mouse.meme"
    )
    
    cmd <- sprintf(
        "meme-chip -meme-p 8 -dreme-m 10 -spamo-skip -oc meme_chip %s %s",
        databases, 'summit_sequence_500bp.fa'
    )
    
    #output <- gsub('\\..+$', '', basename(file))
 
    cmd2 <- sprintf('echo "%s" | %s', cmd, interperor)
    message(cmd2)
    system(cmd2)
    
    if(file.exists('meme_chip/meme-chip.html')) message('MEME html exists!')
    #file.link('meme_chip/meme-chip.html', 'meme.html')
    
    return('meme_chip/meme-chip.html')
}


#' meme_chip_local
#' 
#' @param IDs Vector of JADB ContactExpIDs
#'   
#' @return List 
#' 
#' @author Przemyslaw Stempor
#' 
#' @family Peaks
#' @export
#' 
#' @examples
#' #meme_chip_local('AA001')
meme_chip_local <- function(ids) {
    library(BSgenome.Celegans.UCSC.ce11)
    library(rtracklayer)
    
    base_dir  <- getwd()
    on.exit(setwd(base_dir))
    
    fls <- sapply(ids, getFilePath, format = "bed", eq=TRUE, processing = "summits", scale = "MACS")
    if (length(fls) != 1) stop('No or more than 1 file.')
    
    prefix <- gsub("\\^[^\\^]+$", '', basename(fls))
    
    exp_dir <- gsub('files/', '', dirname(fls))
    
    con <- url(fls)
    summits <- import.bed(con)
    close(con)
    
    
    peaks_fls <- url(sapply(ids, getFilePath, format = "narrowPeak", processing = "PeakCalls", scale = "MACS"))
    peaks <- import(
        peaks_fls, format = "BED", 
        extraCols = c(singnalValue = "numeric", pValue = "numeric", qValue = "numeric", peak = "integer")
    )
    close(peaks_fls)
    
    seqinfo(summits) <- seqinfo(Celegans)[seqlevels(summits)]
    
    seq <- getSeq(Celegans, trim(resize(summits, 100, fix = 'center')))
    names(seq) <- paste(summits)
    seq <- seq[width(seq)==100]

    writeXStringSet(seq, 'summit_sequence.fa')
    #writeXStringSet(seq[1:25], 'summit_sequence_25_100bp.fa')
    
    output <- gsub('\\..+$', '', basename(fls))
    cmd <- sprintf(
        "meme-chip -meme-p 8 -oc meme_test_100_r -db ~/meme/db/motif_databases/WORM/uniprobe_worm.meme -db ~/meme/db/motif_databases/EUKARYOTE/hallikas2006.meme -db ~/meme/db/motif_databases/EUKARYOTE/homeodomain.meme -db ~/meme/db/motif_databases/EUKARYOTE/jolma2010.meme -db ~/meme/db/motif_databases/EUKARYOTE/jolma2013.meme -db ~/meme/db/motif_databases/EUKARYOTE/macisaac_theme.v1.meme -db ~/meme/db/motif_databases/EUKARYOTE/prodoric.meme -db ~/meme/db/motif_databases/EUKARYOTE/regtransbase.meme -db ~/meme/db/motif_databases/EUKARYOTE/SwissRegulon_human_and_mouse.meme -db ~/meme/db/motif_databases/EUKARYOTE/wei2010_human_mws.meme -db ~/meme/db/motif_databases/EUKARYOTE/wei2010_mouse_mws.meme -db ~/meme/db/motif_databases/EUKARYOTE/wei2010_mouse_pbm.meme -db ~/meme/db/motif_databases/EUKARYOTE/zhao2011.meme summit_sequence_25_100bp.fa"
    )
    cmd2 <- sprintf('echo "%s" | %s', cmd, interperor)
    message(cmd2)
    system(cmd2)
    
    
}

