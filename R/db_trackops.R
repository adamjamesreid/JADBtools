#' jadb_addTracksFromBAM
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
#' #jadb_addTracksFromBAM('AA001')
jadb_addTracksFromBAM <- function(ids, genome='ce11') {
    require(magrittr)
    require(Rsamtools)
    require(BSgenome)
    
    base_dir  <- getwd()
    on.exit(setwd(base_dir))
    
    ids  %>% sapply(getFilePath, format = "bam", eq=TRUE, processing = "aligned", scale = "NA", url = FALSE)  -> fls
    if (length(ids) != 1) stop('No or more than 1 BAM files.')

    prefix <- gsub("\\^[^\\^]+$", '', basename(fls))
    
    exp_dir <- gsub('files/', '', dirname(fls))
    message(exp_dir)
    setwd(exp_dir)
    
    
    what <- c("rname", "strand", "pos", "qwidth", "mapq")
    flag <- scanBamFlag(isUnmappedQuery = FALSE)
    param <- ScanBamParam(flag = flag, simpleCigar = FALSE, what = what)
    
    message('File: ', basename(fls))
    
    a <- scanBam(basename(fls), param = param)[[1]]
    grng <- GRanges(
        seqnames = a$rname, ranges = IRanges(a$pos, width = a$qwidth), 
        strand = a$strand, seqinfo=SeqinfoForBSGenome(genome), mapq=a$mapq
    )
    grng <- trim(resize(grng, 200L))
    message('Aligned sequences: ', length(grng))
    
    stats <- JADBtools::bamStats(basename(fls), a)
    message(stats)
    
    Entry <- addGenericFile(
        ids,
        path = file.path('files', exp_dir, gsub('aligned\\^NA\\^NA', paste0('alnNQNU^linear^', genome), prefix)), 
        Processing = 'alnNQNU', 
        Scale = 'linear', 
        Resolution = '1bp',
        filetype_format = 'bw', 
        prefix = 'P',
        comments = stats,
        genome = genome,
        uniq = FALSE
    )
    
    message("All ranges (200bp): ", length(grng))
    export.bw(coverage(grng), basename(Entry$path))
    
    
    EntryQ10 <- addGenericFile(
        ids,
        path = file.path('files', exp_dir, gsub('aligned\\^NA\\^NA', paste0('alnQ10NU^linear^', genome), prefix)), 
        Processing = 'alnQ10NU', 
        Scale = 'linear', 
        Resolution = '1bp',
        filetype_format = 'bw', 
        prefix = 'P',
        comments = stats,
        genome = genome,
        uniq = FALSE
    )
    
    grngQ10 <- grng[grng$mapq >= 10]
    message("All ranges (200bp) mapQ10: ", length(grngQ10))
    export.bw(coverage(grngQ10), basename(EntryQ10$path))
    
}

#' jadb_addScaledTrack
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
jadb_addScaledTrack <- function(ids, scale='zscore', input='BEADSQ10NU', genome='ce11') {
    require(magrittr)
    require(Rsamtools)
    require(BSgenome)
    
    base_dir  <- getwd()
    on.exit(setwd(base_dir))
    
    fls <- getFilePath(ids, processing = input, format = 'bw', scale = 'linear', url = FALSE, eq=TRUE)
    prefix <- gsub("\\^[^\\^]+$", '', basename(fls))
    
    exp_dir <- gsub('files/', '', dirname(fls))
    message(exp_dir)
    setwd(exp_dir)
    
    
    Entry <- addGenericFile(
        ids,
        path = file.path('files', exp_dir, gsub('linear', 'zscore', prefix)), 
        Processing = input, 
        Scale = scale, 
        Resolution = '1bp',
        filetype_format = 'bw', 
        prefix = 'P',
        genome = genome,
        uniq = FALSE
    )
    
    message(basename(fls))
    cov <- import.bw(basename(fls), as='RleList')

    
    if (scale=='zscore') {
        ucov <- unlist(cov)
        mi <- mean(ucov)
        mu <- sd(ucov)
        zsc <- (cov-mi)/mu
        message(paste(names(mean(zsc)), '\t', mean(zsc), '\n' ))
        export.bw(zsc, basename(Entry$path))
    } else {
        stop(scale, ' not yet supported!')
    }
    
    
    message('Done.')
    
}



