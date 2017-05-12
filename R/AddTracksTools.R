#' addMapq0BeadsRep
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
addMapq0BeadsRep <- function(ids) {
    require(magrittr)
    
    base_dir  <- getwd()
    on.exit(setwd(base_dir))
    
    ids  %>% sapply(getFilePath, format = "bw", eq=TRUE, processing = "NORM", scale = "linear", url = FALSE)  -> fls
    sub('\\^', '|', strsplit(fls, '_')[[1]][[7]]) -> ind_ids
    
    ind_ids  %>% sapply(getFilePath, format = "bw", processing = "BEADSmapq0", scale = "linear", url = FALSE)  -> ind_fls
    message(ind_fls[[1]], '\n', ind_fls[[2]])
    
    ##anno <- as.data.frame(t(sapply(basename(ind_fls[[1]]), rbeads:::ParseName)))
    out <- combineReps(strsplit(ind_ids, '\\|')[[1]], processing = 'BEADSmapq0', outdir = gsub('files/', '', dirname(fls)), res = 100)
    addGenericFile(
        ids, 
        path = file.path('files', out$out), 
        Processing = 'BEADSmapq0', 
        Resolution = '1bp', Scale = 'linear', filetype_format = 'bw', 
        prefix = 'R', repPath = TRUE
    )
    
    message('Done :)')
    
}



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
    
    
    final.path <- file.path('files', exp_dir, gsub('aligned\\^[A-Za-z0-9]+\\^[A-Za-z0-9]+', 'BEADSmapq0^linear^1bp', prefix))
    
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
    
    out <- file.rename(basename(path(NRM0)), basename(Entry$path))
    
    
    message("Done!")
    
}

#' addNonUniqueQ10Beads
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
addNonUniqueQ10Beads <- function(ids, genome='ce11') {
    require(magrittr)
    require(Rsamtools)
    require(BSgenome)
    
    base_dir  <- getwd()
    on.exit(setwd(base_dir))
    
    ids  %>% sapply(getFilePath, format = "bam", eq=TRUE, processing = "aligned", scale = "NA", url = FALSE)  -> fls
    if (length(ids) != 1) stop('No or more than 1 BAM files.')
    
    outnames <- sapply(basename(fls), rbeads:::reName, proccesing = 'mapq0', scale = 'linear', resolution = genome, ext='.bw')
    prefix <- basename(fls) %>% substr(start=0, stop=nchar(.)-13)
    
    exp_dir <- gsub('files/', '', dirname(fls))
    message(exp_dir)
    setwd(exp_dir)
    
    crosslink <- JADBtools::getAnno(ids, anno = 'Crosslinker', EXTABLE = 'labexperimentview')
    input_suffix <- '_HiSeq_nonUNIQ_MAPQ10_200bp_SummedInput_bin25bp.bw'
    if(grepl('^e', crosslink, ignore.case = TRUE)) {
        input <- file.path(base_dir, 'Input/SummedInputs', genome, paste0(genome, '_EGS', input_suffix))
    }  else {
        input <- file.path(base_dir, 'Input/SummedInputs', base_dir, paste0(genome, '_FRM', input_suffix))
    }
    
    
    mappability <- file.path(base_dir, paste0("_mappability_files_/", genome, "_gem-mappability_36bp.bw"))
    
    message('File: ', basename(fls), '\n vs. ', basename(input))
    
    message(getwd())
    NRM0 <- beads(
        basename(fls), 
        input, 
        mappability,
        genome, 
        uniq = FALSE, insert = 200L, mapq_cutoff = 10, export = "BEADS", 
        rdata = FALSE, export_er = TRUE, quickMap = TRUE
    )
    
    
    final.path <- file.path('files', exp_dir, gsub('aligned\\^[A-Za-z0-9]+\\^[A-Za-z0-9]+', paste0('BEADSQ10NU^linear^', genome), prefix))
    
    Entry <- addGenericFile(
        ids,
        path = final.path, 
        Processing = 'BEADSQ10NU', 
        Scale = 'linear', 
        Resolution = '1bp',
        filetype_format = 'bw', 
        prefix = 'P',
        genome = genome,
        comments = JADBtools::bamStats(basename(fls)),
        uniq = FALSE
    )
    
    out <- file.rename(basename(path(NRM0)), basename(Entry$path))
    
    
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
        path = file.path('files', exp_dir, gsub('aligned\\^[A-Za-z0-9]+\\^[A-Za-z0-9]+', 'mapq0^linear^1bp', prefix)), 
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
        path = file.path('files', exp_dir, gsub('aligned\\^[A-Za-z0-9]+\\^[A-Za-z0-9]+', 'mapq0^zscore^1bp', prefix)), 
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



#' addBEADSmapq0Track and z-score it
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
#' #parallel::mclapply(sprintf('REP%.3i', 36:42)[-6], addBEADSmapq0TrackZcs), mc.cores = 8)
#' # sapply(sprintf('REP%.3i', 3:24), addBEADSmapq0TrackZcs)
addBEADSmapq0TrackZcs <- function(ids) {
    require(magrittr)
    require(Rsamtools)
    require(BSgenome)
    
    base_dir  <- getwd()
    on.exit(setwd(base_dir))
    
    fls <- getFilePath(ids, processing = 'BEADSmapq0', format = 'bw', scale = 'linear', url = FALSE, eq=TRUE)
    
    outnames <- sapply(basename(fls), rbeads:::reName, proccesing = 'BEADSmapq0', scale = 'zscore', resolution = '1bp', ext='.bw')
    prefix <- basename(fls) %>% substr(start=0, stop=nchar(.)-13)
    
    exp_dir <- gsub('files/', '', dirname(fls))
    message(exp_dir)
    setwd(exp_dir)
    
    
    Entry <- addGenericFile(
        ids,
        path = file.path('files', exp_dir, gsub('linear', 'zscore', prefix)), 
        Processing = 'BEADSmapq0', 
        Scale = 'zscore', 
        Resolution = '1bp',
        filetype_format = 'bw', 
        prefix = 'P',
        comments = ''
    )
    
    message(basename(fls))
    cov <- import.bw(basename(fls), as='RleList')
    ucov <- unlist(cov)
    
    mi <- mean(ucov)
    mu <- sd(ucov)
    zsc <- (cov-mi)/mu
    
    
    message(paste(names(mean(zsc)), '\t', mean(zsc), '\n' ))
    export.bw(zsc, basename(Entry$path))
    
    message('Done.')
    
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