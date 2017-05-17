#' jadb_trimm_fq
#' 
#' @param IDs Vector of JADB ContactExpIDs
#'   
#' @return null 
#' 
#' @author Przemyslaw Stempor
#' 
#' @family db_rnaseq
#' @export
#' 
jadb_trimm_fq <- function(ids) {
    
    library(JADBtools)
    library(RJSONIO)
    library(rtracklayer)
    
    base_dir  <- getwd()
    on.exit(setwd(base_dir))
    
    ids  %>% sapply(getFilePath, format = "txt.gz", eq=TRUE, processing = "raw", scale = "NA", url = FALSE)  -> fls
    if (length(ids) != 1) stop('No or more than 1 BAM files.')
    
    exp_dir <- gsub('files/', '', dirname(fls))
    message(exp_dir); setwd(exp_dir)
    
    fq.uid <- getFileUID(ids, format='txt.gz', processing='raw', scale='NA', eq=TRUE)
    ref.genome.version <- 'WS220/ce10'
    
    ##### Trimming #####
    
    trimmed <- addGenericFile( 
        ContactExpID=ids, Processing='trimmed36bp', prefix='U', filetype_format="fq.gz", 
        parent1_uid=fq.uid, genome =  NA,
        path=file.path(dirname(fls), formGenericPath(ids, 'labrnaseq', 'trimmed36bp', NA, NA))
    )
    
    Sys.setenv( 
        input=basename(fls),
        out=basename(trimmed$path),
        trimmer='fastx_trimmer'
    )
    cmd <- system('echo "zcat $input | $trimmer -z -l 36 -Q 33 -o $out"', intern=TRUE)
    message(cmd, '\n'); system(paste('echo "', cmd, '" | bash'), wait = TRUE)
    
    message('Trimming done!')
}

#' jadb_align_rnaseq
#' 
#' @param IDs Vector of JADB ContactExpIDs
#'   
#' @return null 
#' 
#' @author Przemyslaw Stempor
#' 
#' @family db_rnaseq
#' @export
#' 
jadb_align_rnaseq <- function(ids) {
    
    library(JADBtools)
    library(RJSONIO)
    library(rtracklayer)
    
    base_dir  <- getwd()
    on.exit(setwd(base_dir))
    
    ids  %>% sapply(getFilePath, format = "fq.gz", eq=TRUE, processing = "trimmed36bp", scale = "NA", url = FALSE)  -> fls
    if (length(ids) != 1) stop('No or more than 1 BAM files.')
    
    exp_dir <- gsub('files/', '', dirname(fls))
    message(exp_dir); setwd(exp_dir)
    
    parent.uid <- getFileUID(ids, format='fq.gz', processing = "trimmed36bp", scale='NA', eq=TRUE)
    ref.genome.version <- 'WS220/ce10'
    
    message('Aligning file: ', basename(fls))
    message(getwd())
    
    ALN <- JADBtools:::run_bwa(
        basename(fls), 
        file.path(base_dir, '_ref_genomes_/ce10.fa')
    )
    message('Aligned!')
    
    stats <- JADBtools::bamStats(ALN)
    
    aligned <- addGenericFile(
        ids,
        path = file.path(dirname(fls), formGenericPath(ids, 'labrnaseq', 'aligned', NA, NA)), 
        Processing = 'aligned', 
        Scale = 'NA', 
        Resolution = 'NA',
        filetype_format = 'bam', 
        prefix = 'R',
        comments = stats,
        parent1_uid=parent.uid,
        genome =  ref.genome.version
    )
    out <- file.rename(basename(ALN), basename(aligned$path))
    
    message("Indexing")
    indexBam(basename(aligned$path))
    
    message("Alignment done")
}

#' db_mran_seq_align_process
#' 
#' @param IDs Vector of JADB ContactExpIDs
#'   
#' @return null 
#' 
#' @author Przemyslaw Stempor
#' 
#' @family db_rnaseq
#' @export
#' 
db_mran_seq_align_process <- function(ids) {
    
    base_dir  <- getwd()
    on.exit(setwd(base_dir))
    
    ids  %>% sapply(getFilePath, format = "bam", eq=TRUE, processing = "aligned", scale = "NA", url = FALSE)  -> fls
    if (length(ids) != 1) stop('No or more than 1 BAM files.')
    
    exp_dir <- gsub('files/', '', dirname(fls))
    message(exp_dir); setwd(exp_dir)
    
    parent.uid <- getFileUID(ids, format='bam', processing = "aligned", scale='NA', eq=TRUE)
    ref.genome.version <- 'WS220/ce10'
    
    message('Adding tracks for: ', basename(fls))
    
    f <- basename(basename(fls))
    ID <- ids
    fq.dir <- file.path('files', exp_dir)
    bam_stats <- JADBtools::bamStats(f)
    
    require(Rsamtools)
    require(rtracklayer)
    
    what <- c("rname", "strand", "pos", "qwidth", "mapq")
    flag <- scanBamFlag(isUnmappedQuery = FALSE)
    param <- ScanBamParam(flag = flag, simpleCigar = FALSE, what = what)
    
    message('Processing file: ', basename(f))
    
    a <- scanBam(f, param = param)[[1]]
    grng <- GRanges(
        seqnames = a$rname, ranges = IRanges(a$pos, width = a$qwidth), 
        strand = a$strand, seqinfo=SeqinfoForBSGenome('ce10')
    )       
    message("All ranges: ", length(grng))
    
    addLinearBW <- function(track, proc) {
        outEntry <- addGenericFile( 
            ContactExpID=ID, Processing=proc, prefix='U', filetype_format="bw", 
            parent1_uid=parent.uid , genome =  ref.genome.version,
            Resolution='1bp', Scale='linear', comments = bam_stats,
            path=file.path(fq.dir, formGenericPath(ID, 'labrnaseq', proc, '1bp', 'linear'))
        )
        export.bw(track, basename(outEntry$path))
    }
    
    addLinearBW(coverage(grng), 'readCoverage')
    #addLinearBW(coverage(grng[strand(grng)=="+"]), 'FWDtagCoverage.bw')
    #addLinearBW(coverage(grng[strand(grng)=="-"]), 'REVtagCoverage.bw')
    
    quality_cutoff = 10L
    lg <- (a$mapq >= quality_cutoff)
    grng2 <- grng[lg]
    message("Mapq10 ranges: ", length(grng2))
    
    addLinearBW(coverage(grng2), 'mapq10readCoverage')
    #addLinearBW(coverage(grng2[strand(grng2)=="+"]), 'FWDmapq10tagCoverage.bw')
    #addLinearBW(coverage(grng2[strand(grng2)=="-"]), 'REVmapq10tagCoverage.bw')
    
    message('Creating summarized experiment')
    
    e <- summarizeBAMs(f)
    message(colnames(e))
    colData(e)$strain <- factor( unlist(getStrain(ID)) )
    colData(e)$stage <- factor( unlist(getStage(ID)) )
    
    #message('Saving summarized experiment')
    entry <- addGenericFile( 
        ContactExpID=ID, Processing='TagCounts', prefix='U', filetype_format="Rdata", 
        genome =  ref.genome.version, 
        path=file.path(fq.dir, formGenericPath(ID, 'labrnaseq', 'TagCounts', NA, NA))
    )
    save(e, file=basename(entry$path) ) 
    
    message('Exporting CSV')
    M <- assays(e)$counts
    RPKM <- as.data.frame(
        M / ((sum(width(gnmodel)) / 1e3) %*% t(colSums(M) / 1e6))
    )
    colnames(RPKM) <- paste0('RPKM_', colnames(RPKM))
    
    out <- cbind(
        elementMetadata(gnmodel), 
        data.frame(exon_width_kb=(sum(width(gnmodel)) / 1e3)),
        M,
        RPKM
    )
    
    entry <- addGenericFile( 
        ContactExpID=ID, Processing='RPKM', prefix='U', filetype_format="csv", 
        genome =  ref.genome.version, 
        path=file.path(fq.dir, formGenericPath(ID, 'labrnaseq', 'RPKM', NA, NA))
    )
    
    write.csv(out, basename(entry$path))
    
    message('Pipeline complete.')
}


#' db_mran_rpm_track
#' 
#' @param IDs Vector of JADB ContactExpIDs
#'   
#' @return null 
#' 
#' @author Przemyslaw Stempor
#' 
#' @family db_rnaseq
#' @export
#' 
db_mran_rpm_track <- function(ids) {
    require(magrittr)
    require(Rsamtools)
    require(BSgenome)
    
    base_dir  <- getwd()
    on.exit(setwd(base_dir))
    
    ids  %>% sapply(getFilePath, format = "bam", eq=TRUE, processing = "aligned", scale = "NA", url = FALSE)  -> fls
    if (length(ids) != 1) stop('No or more than 1 BAM files.')
    
    
    outnames_both <- sub('\\^.+', '', sub(
        '([^_]+)_([^_]+)_([^_]+)_([^_]+)_([^_]+)_([^_]+)_([^_]+)', 
        '\\1_\\2_\\3_\\4_\\5_readCoverage_RPM_1bp', 
        basename(fls)
    ))
    
    
    prefix <- basename(fls) %>% substr(start=0, stop=nchar(.)-13)
    
    exp_dir <- gsub('files/', '', dirname(fls))
    message(exp_dir)
    setwd(exp_dir)
    
    ref.genome.version <- 'WS220/ce10'
    
    
    what <- c("rname", "strand", "pos", "qwidth")
    flag <- scanBamFlag(isUnmappedQuery = FALSE)
    param <- ScanBamParam(flag = flag, simpleCigar = FALSE, what = what)
    
    message('File: ', basename(fls))
    
    a <- scanBam(basename(fls), param = param)[[1]]
    grng <- GRanges(
        seqnames = a$rname, ranges = IRanges(a$pos, width = a$qwidth), 
        strand = a$strand, seqinfo=SeqinfoForBSGenome('ce10')
    )
    
    nf <- (length(grng)/1e6)
    
    Entry <- addGenericFile(
        ids,
        path = file.path('files', exp_dir, outnames_both), 
        Processing = 'readCoverage', 
        Scale = 'RPM', 
        Resolution = '1bp',
        filetype_format = 'bw', 
        prefix = 'R',
        comments = JADBtools::bamStats(basename(fls)),
        genome =  ref.genome.version
    )
    message('Exporting ', basename(Entry$path))
    export.bw(coverage(grng)/nf, basename(Entry$path))
    message('Done!')
}