#! Rscript
# Usage: Rscript -e "devtools::source_gist('8636a090e43ab4b9f035803eb774ee3c')"

message('Running PE stranded insert coverage - ', format(Sys.time(), "%a %b %d %X %Y"))
args <- commandArgs(TRUE)

quality_cutoff <- 10

require(Rsamtools)
require(rtracklayer)
require(parallel)
what <- c("rname", "pos", "isize", 'strand', 'mapq')

flag <- scanBamFlag(isUnmappedQuery = FALSE, isProperPair = TRUE, isMinusStrand = FALSE, isMateMinusStrand = TRUE)
param <- ScanBamParam(what = what, flag=flag)

library(parallel)
options(srapply_fapply="parallel", mc.cores=4)

strPEcov <- function(f) {
    message('=> ', basename(f))
    
    a <- scanBam(f, param=param)[[1]]
    isize <- a$isize[ a$isize > 0 ]
    pos <- a$pos[ a$isize > 0 ]
    rname <- a$rname[ a$isize > 0 ]
    strand <- a$strand[ a$isize > 0 ]
    
    #All reads
    grng <- GRanges( seqnames =  rname, ranges=IRanges(start = pos, width = isize), strand = strand )
    grng <- trim(grng)
    
    message(basename(f), " - aligned ranges: ", length(grng))
    
    cover <- coverage(grng)
    nf <- length(grng)/1e6
    
    export.bw(coverage(grng[strand(grng)=="+"]), gsub('.bam$', '_FWD_insert_coverage.bw', basename(f)))
    export.bw(coverage(grng[strand(grng)=="-"]), gsub('.bam$', '_REV_insert_coverage.bw', basename(f)))
    export.bw(coverage(grng[strand(grng)=="+"])/nf, gsub('.bam$', '_FWD_read_coverage_readsPerMillion.bw', basename(f)))
    export.bw(coverage(grng[strand(grng)=="-"])/nf, gsub('.bam$', '_REV_read_coverage_readsPerMillion.bw', basename(f)))
    
    #mapq10 reads
    lg <- (a$mapq[ a$isize > 0 ] >= quality_cutoff)
    grng <- grng[lg]
    
    message(basename(f), " - mapq10 ranges: ", length(grng))
    
    cover <- coverage(grng)
    nf <- length(grng)/1e6
    
    export.bw(coverage(grng[strand(grng)=="+"]), gsub('.bam$', '_FWD_mapq10_insert_coverage.bw', basename(f)))
    export.bw(coverage(grng[strand(grng)=="-"]), gsub('.bam$', '_REV_mapq10_insert_coverage.bw', basename(f)))
    export.bw(coverage(grng[strand(grng)=="+"])/nf, gsub('.bam$', '_FWD_mapq10_read_coverage_readsPerMillion.bw', basename(f)))
    export.bw(coverage(grng[strand(grng)=="-"])/nf, gsub('.bam$', '_REV_mapq10_read_coverage_readsPerMillion.bw', basename(f)))
    
}


mclapply(dir(pattern='.bam$'), strPEcov, mc.cores = parallel::detectCores())

q()


b <- BamFile(f[[1]])
lst <- as.list(GRanges(seqinfo(b)))

gc()
alns <- mclapply(lst, function(chr) {
    readGAlignmentPairs(b, param = ScanBamParam(which=chr, what = 'mapq'))
}, mc.cores = parallel::detectCores())


grs <- mclapply(alns, granges, mc.cores = parallel::detectCores())
GR <- unlist(GRangesList(grs))
FWD <- GR[strand(GR) == '+']
REV <- GR[strand(GR) == '-']


nf <- (length(FWD) + length(REV))/1e6

expr <- list(
    expression( export.bw(coverage(FWD), gsub('.bam$', '_FWD_insert_coverage_p2.bw', basename(f))) ),
    expression( export.bw(coverage(REV), gsub('.bam$', '_REV_insert_coverage_p2.bw', basename(f))) ), 
    expression( export.bw(coverage(FWD)/nf, gsub('.bam$', '_FWD_read_coverage_readsPerMillion_p2.bw', basename(f))) ),
    expression( export.bw(coverage(REV)/nf, gsub('.bam$', '_REV_read_coverage_readsPerMillion_p2.bw', basename(f))) )
)
mclapply(expr, eval, mc.cores = parallel::detectCores())

###### -------
b <- BamFile(f[[1]])
indexBam(b)
lst <- as.list(GRanges(seqinfo(b)))

getGR <- function(chr, isFirstMateRead=TRUE) {
    flag <- scanBamFlag(isUnmappedQuery = FALSE, isProperPair = TRUE, isMinusStrand = FALSE, isMateMinusStrand = TRUE, isFirstMateRead=isFirstMateRead)
    param <- ScanBamParam(what = what, flag=flag, which=chr)
    a <- scanBam(f, param=param)[[1]]
    isize <- a$isize[ a$isize > 0 ]
    pos <- a$pos[ a$isize > 0 ]
    rname <- a$rname[ a$isize > 0 ]
    strand <- a$strand[ a$isize > 0 ]
    
    #All reads
    grng <- GRanges( seqnames =  rname, ranges=IRanges(start = pos, width = isize), strand = strand )
    grng <- trim(grng)
    grng$mapq <- a$mapq[ a$isize > 0 ]
    grng
}

FWD <- unlist(GRangesList(mclapply(lst, getGR, TRUE, mc.cores = parallel::detectCores())))
REV <- unlist(GRangesList(mclapply(lst, getGR, FALSE, mc.cores = parallel::detectCores())))


nf <- (length(FWD) + length(REV))/1e6

expr <- list(
    expression( export.bw(coverage(FWD), gsub('.bam$', '_FWD_insert_coverage.bw', basename(f))) ),
    expression( export.bw(coverage(REV), gsub('.bam$', '_REV_insert_coverage.bw', basename(f))) ), 
    expression( export.bw(coverage(FWD)/nf, gsub('.bam$', '_FWD_read_coverage_readsPerMillion.bw', basename(f))) ),
    expression( export.bw(coverage(REV)/nf, gsub('.bam$', '_REV_read_coverage_readsPerMillion.bw', basename(f))) )
)
mclapply(expr, eval, mc.cores = parallel::detectCores())


flag <- scanBamFlag(isUnmappedQuery = FALSE, isProperPair = TRUE, isMinusStrand = FALSE, isMateMinusStrand = TRUE)
param <- ScanBamParam(what = what, flag=flag)

library(parallel)
options(srapply_fapply="parallel", mc.cores=4)

strPEcov <- function(f) {
    message('=> ', basename(f))
    
    a <- scanBam(f, param=param)[[1]]
    
    