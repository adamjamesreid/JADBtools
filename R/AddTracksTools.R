#' addMapq0Beads
#' 
#' @param IDs Vector of JADB ContactExpIDs
#'   
#' @return List 
#' 
#' @author Przemyslaw Stempor
#' 
#' @family tracks
#' @export
#' 
#' @examples
#' #addMapq0Beads('AA001')
addMapq0Beads <- function(ids) {
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

    crosslink <- JADBtools::getAnno(ids, anno = 'Crosslinker', EXTABLE = 'labexperimentview')
    if(grepl('^e', crosslink, ignore.case = TRUE)) {
        input <- file.path(base_dir, 'Input/SummedInputs/HISeq_EGS_map0_NOTuniq_SummedInput_linear_25bp_SummedInput_linear_25bp.bw')
    }  else {
        input <- file.path(base_dir, 'Input/SummedInputs/HISeq_FRM_map0_NOTuniq_SummedInput_linear_25bp_SummedInput_linear_25bp.bw')
    }

    message('File: ', basename(fls), '\n vs. ', basename(input))
    
    message(getwd())
    NRM0 <- beads(
        basename(fls), 
        input, 
        file.path(base_dir, 'Input/SummedInputs/map1.bw'), 'ce10', 
        uniq = FALSE, insert = 200L, mapq_cutoff = 0, export = "BEADS", 
        rdata = FALSE, export_er = TRUE, quickMap = TRUE
    )
    
    
    final.path <- file.path('files', exp_dir, gsub('aligned\\^NA\\^NA', 'BEADSmapq0^linear^1bp', prefix))
    out <- file.rename(basename(path(NRM0)), basename(final.path))
    
    Entry <- addGenericFile(
        ids,
        path = final.path, 
        Processing = 'BEADSmapq0', 
        Scale = 'linear', 
        Resolution = '1bp',
        filetype_format = 'bw', 
        prefix = 'P',
        comments = JADBtools::bamStats(basename(fls))
    )
    
    message("Done!")
    
}


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

#' addMapq0Track and z-score it
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
addMapq0TrackZcs <- function(ids) {
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
        path = file.path('files', exp_dir, gsub('aligned\\^NA\\^NA', 'mapq0^zscore^1bp', prefix)), 
        Processing = 'mapq0', 
        Scale = 'zscore', 
        Resolution = '1bp',
        filetype_format = 'bw', 
        prefix = 'P',
        comments = JADBtools::bamStats(basename(fls))
    )
    
    cov <- coverage(grng)
    ucov <- unlist(cov)
    
    mi <- mean(ucov)
    mu <- sd(ucov)
    zsc <- (cov-mi)/mu
    
    
    message("All ranges (200bp): ", length(grng))
    export.bw(zsc, basename(Entry$path))
    
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
#' #parallel::mclapply(sprintf('rAM%.3i', 47:56), addStrandedRNAseq, mc.cores = 8)
addStrandedRNAseq <- function(ids) {
    require(magrittr)
    require(Rsamtools)
    require(BSgenome)
    
    base_dir  <- getwd()
    on.exit(setwd(base_dir))
    
    ids  %>% sapply(getFilePath, format = "bam", eq=TRUE, processing = "aligned", scale = "NA", url = FALSE)  -> fls
    if (length(ids) != 1) stop('No or more than 1 BAM files.')
    
    outnames_fwd <- sub('\\^.+', '', sub(
        '([^_]+)_([^_]+)_([^_]+)_([^_]+)_([^_]+)_([^_]+)_([^_]+)', 
        '\\1_\\2_\\3_\\4_\\5_FWDreadCoverage_RPM_1bp', 
        basename(fls)
    ))
    
    outnames_rev <- sub('\\^.+', '', sub(
        '([^_]+)_([^_]+)_([^_]+)_([^_]+)_([^_]+)_([^_]+)_([^_]+)', 
        '\\1_\\2_\\3_\\4_\\5_REVreadCoverage_RPM_1bp', 
        basename(fls)
    ))
    
    outnames_both <- sub('\\^.+', '', sub(
        '([^_]+)_([^_]+)_([^_]+)_([^_]+)_([^_]+)_([^_]+)_([^_]+)', 
        '\\1_\\2_\\3_\\4_\\5_readCoverage_RPM_1bp', 
        basename(fls)
    ))
    

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
    
    Entry <- addGenericFile(
        ids,
        path = file.path('files', exp_dir, outnames_fwd), 
        Processing = 'FWDreadCoverage', 
        Scale = 'RPM', 
        Resolution = '1bp',
        filetype_format = 'bw', 
        prefix = 'P',
        comments = JADBtools::bamStats(basename(fls))
    )
    message('Exporting ', basename(Entry$path))
    export.bw(coverage(grng[strand(grng)=="-"])/nf, basename(Entry$path))

    Entry <- addGenericFile(
        ids,
        path = file.path('files', exp_dir, outnames_rev), 
        Processing = 'REVreadCoverage', 
        Scale = 'RPM', 
        Resolution = '1bp',
        filetype_format = 'bw', 
        prefix = 'P',
        comments = JADBtools::bamStats(basename(fls))
    )
    message('Exporting ', basename(Entry$path))
    export.bw(coverage(grng[strand(grng)=="+"])/nf, basename(Entry$path))
    
    Entry <- addGenericFile(
        ids,
        path = file.path('files', exp_dir, outnames_both), 
        Processing = 'readCoverage', 
        Scale = 'RPM', 
        Resolution = '1bp',
        filetype_format = 'bw', 
        prefix = 'P',
        comments = JADBtools::bamStats(basename(fls))
    )
    message('Exporting ', basename(Entry$path))
    export.bw(coverage(grng)/nf, basename(Entry$path))
}