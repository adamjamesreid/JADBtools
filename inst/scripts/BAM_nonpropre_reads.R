require(Rsamtools)
require(rtracklayer)
require(GenomicAlignments)
require(parallel)
require(BiocParallel)

fff <- dir(pattern = '*bam$')

ff <- fff[1]


bf <- BamFile(ff, asMates=TRUE)
param <- ScanBamParam(what = c('mapq', 'isize', 'flag'))
gl <- readGAlignmentPairsFromBam(bf, use.names=FALSE, param=param)

save(gl, file=paste0(ff, '.allPairs.Rdata'))

flt <- lapply(seqlevels(gl), function(x) {
    f <- all(seqnames(gl) == x)
    message(x, ' -> ', f)
    return(f)
})

starts <- start(gl)
st <- do.call(rbind, stars)
flt_dup <- duplicated(st)

flt <- all(seqnames())