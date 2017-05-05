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
jadb_addTracksFromBAM <- function(ids) {
    library(BSgenome.Celegans.UCSC.ce11)
    library(rtracklayer)
    
    base_dir  <- getwd()
    on.exit(setwd(base_dir))
    
    fls <- sapply(ids, getFilePath, format = "bed", eq=TRUE, processing = "summits", scale = "MACS")
    if (length(fls) != 1) stop('No or more than 1 file.')
    
    prefix <- gsub("\\^[^\\^]+$", '', basename(fls))
    
    exp_dir <- gsub('files/', '', dirname(fls))
    message(exp_dir)
    setwd(exp_dir)
    
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
    
    seq <- getSeq(Celegans, trim(resize(summits, 1000, fix = 'center')))
    seq <- seq[width(seq)==1000]
    writeXStringSet(seq, 'summit_sequence.fa')
    

    
    what <- c("rname", "strand", "pos", "qwidth", "mapq")
    flag <- scanBamFlag(isUnmappedQuery = FALSE)
    param <- ScanBamParam(flag = flag, simpleCigar = FALSE, what = what)
    
    message('File: ', basename(fls))
    
    a <- scanBam(basename(fls), param = param)[[1]]
    grng <- GRanges(
        seqnames = a$rname, ranges = IRanges(a$pos, width = a$qwidth), 
        strand = a$strand, seqinfo=SeqinfoForBSGenome('ce10'), mapq=a$mapq
    )
    grng <- trim(resize(grng, 200L))
    message('Aligned sequences: ', length(grng))
    
    stats <- JADBtools::bamStats(basename(fls), a)
    message(stats)
    
    Entry <- addGenericFile(
        ids,
        path = file.path('files', exp_dir, gsub('aligned\\^NA\\^NA', 'alnNU^linear^1bp', prefix)), 
        Processing = 'alnNU', 
        Scale = 'linear', 
        Resolution = '1bp',
        filetype_format = 'bw', 
        prefix = 'P',
        comments = stats
    )
    
    message("All ranges (200bp): ", length(grng))
    export.bw(coverage(grng), basename(Entry$path))
    
    
    EntryQ10 <- addGenericFile(
        ids,
        path = file.path('files', exp_dir, gsub('aligned\\^NA\\^NA', 'alnQ10NU^linear^1bp', prefix)), 
        Processing = 'alnQ10NU', 
        Scale = 'linear', 
        Resolution = '1bp',
        filetype_format = 'bw', 
        prefix = 'P',
        comments = stats
    )
    
    grngQ10 <- grng[grng$mapq >= 10]
    message("All ranges (200bp) mapQ10: ", length(grngQ10))
    export.bw(coverage(grngQ10), basename(EntryQ10$path))
    
}

