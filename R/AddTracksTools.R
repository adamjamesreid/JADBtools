#' addMapq0Track
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
#' #addMapq0Track('AA001')
addMapq0Track <- function(ids) {
    require(magrittr)
    require(Rsamtools)
    require(BSgenome)
    
    base_dir  <- getwd()
    on.exit(setwd(base_dir))
    
    ids  %>% sapply(getFilePath, format = "bam", eq=TRUE, processing = "aligned", scale = "NA", url = FALSE)  -> fls
    if (length(ids) != 1) stop('No or more than 1 BAM files.')
    
    outnames <- sapply(basename(fls), rbeads:::reName, proccesing = 'mapq0', scale = 'linear', resolution = '1bp', ext='.bw')
    prefix <- basename(fls) %>% substr(start=0, stop=nchar(.)-13)
    
    exp_dir <- gsub('files/', '', dirname(fls))
    message(exp_dir)
    setwd(exp_dir)
    
    
    what <- c("rname", "strand", "pos", "qwidth")
    flag <- scanBamFlag(isUnmappedQuery = FALSE)
    param <- ScanBamParam(flag = flag, simpleCigar = FALSE, what = what)
    
    message('File: ', basename(fls))
    
    a <- scanBam(basename(fls), param = param)[[1]]
    grng <- GRanges(
        seqnames = a$rname, ranges = IRanges(a$pos, width = a$qwidth), 
        strand = a$strand, seqinfo=SeqinfoForBSGenome('ce10')
    )
    
    grng <- trim(resize(grng, 200L))
    
    Entry <- addGenericFile(
        ids,
        path = file.path('files', exp_dir, gsub('aligned\\^NA\\^NA', 'mapq0^linear^1bp', prefix)), 
        Processing = 'mapq0', 
        Scale = 'linear', 
        Resolution = '1bp',
        filetype_format = 'bw', 
        prefix = 'P',
        comments = JADBtools::bamStats(basename(fls))
    )
    
    message("All ranges (200bp): ", length(grng))
    export.bw(coverage(grng), basename(Entry$path))
    
}

#' addStrandedRNAseq
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
#' #addMapq0Track('AA001')
addStrandedRNAseq <- function(ids) {
    require(magrittr)
    require(Rsamtools)
    require(BSgenome)
    
    base_dir  <- getwd()
    on.exit(setwd(base_dir))
    
    ids  %>% sapply(getFilePath, format = "bam", eq=TRUE, processing = "aligned", scale = "NA", url = FALSE)  -> fls
    if (length(ids) != 1) stop('No or more than 1 BAM files.')
    
    outnames_fwd <- sub(
        '([^_]+)_([^_]+)_([^_]+)_([^_]+)_([^_]+)_([^_]+)_([^_]+)_([^_]+)', 
        '\\1_\\2_\\3_\\4_\\5_FWDreadCoverage_RPM_1bp', 
        basename(fls)
    )
    
    outnames_rev <- sub(
        '([^_]+)_([^_]+)_([^_]+)_([^_]+)_([^_]+)_([^_]+)_([^_]+)_([^_]+)', 
        '\\1_\\2_\\3_\\4_\\5_REVreadCoverage_RPM_1bp', 
        basename(fls)
    )
    

    prefix <- basename(fls) %>% substr(start=0, stop=nchar(.)-13)
    
    exp_dir <- gsub('files/', '', dirname(fls))
    message(exp_dir)
    setwd(exp_dir)
    
    
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
    export.bw(coverage(grng[strand(grng)=="-"])/nf, outnames_fwd)
    export.bw(coverage(grng[strand(grng)=="+"])/nf, outnames_rev)
    
    Entry <- addGenericFile(
        ids,
        path = outnames_fwd, 
        Processing = 'FWDreadCoverage', 
        Scale = 'RPM', 
        Resolution = '1bp',
        filetype_format = 'bw', 
        prefix = 'P',
        comments = JADBtools::bamStats(basename(fls))
    )

    Entry <- addGenericFile(
        ids,
        path = outnames_fwd, 
        Processing = 'REVreadCoverage', 
        Scale = 'RPM', 
        Resolution = '1bp',
        filetype_format = 'bw', 
        prefix = 'P',
        comments = JADBtools::bamStats(basename(fls))
    )
    
    message("All ranges (200bp): ", length(grng))
    export.bw(coverage(grng), basename(Entry$path))
    
}