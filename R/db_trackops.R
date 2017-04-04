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


# require(BSgenome)
# require(JADBtools)
# fls <- out.name.bam
# outnames <- sapply(basename(fls), rbeads:::reName, proccesing = "mapq0", 
#                    scale = "linear", resolution = "1bp", ext = ".bw")
# prefix <- basename(fls) %>% substr(start = 0, stop = nchar(.) - 
#                                        13)
# what <- c("rname", "strand", "pos", "qwidth")
# flag <- scanBamFlag(isUnmappedQuery = FALSE)
# param <- ScanBamParam(flag = flag, simpleCigar = FALSE, what = what)
# message("File: ", basename(fls))
# a <- scanBam(basename(fls), param = param)[[1]]
# grng <- GRanges(seqnames = a$rname, ranges = IRanges(a$pos, 
#                                                      width = a$qwidth), strand = a$strand, seqinfo = SeqinfoForBSGenome("ce10"))
# grng <- trim(resize(grng, 200L))
# Entry <- addGenericFile(ID, path = file.path(
#     fq.dir, 
#      gsub("aligned\\^NA\\^NA", "mapq0^linear^1bp", prefix)), 
#                         Processing = "mapq0", Scale = "linear", Resolution = "1bp", 
#                         filetype_format = "bw", prefix = "P", comments = JADBtools::bamStats(basename(fls)))
# export.bw(coverage(grng), basename(Entry$path))