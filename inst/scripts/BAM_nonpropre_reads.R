require(Rsamtools)
require(rtracklayer)
require(GenomicAlignments)
require(parallel)
require(BiocParallel)

fff <- dir(pattern = '*bam$')

ff <- fff[1]


bf <- BamFile(ff, asMates=TRUE)
param <- ScanBamParam(what = c('mapq', 'isize', 'flag'))
gl <- readGAlignmentsListFromBam(bf, use.names=FALSE, param=param)

save(gl, file=paste0(ff, '.allPairsList.Rdata'))

pair <- elementLengths(gl)>1
p <- gl[pair]

#save(GApairs, file=paste0(ff, '.GApairs.Rdata'))

flt <- lapply(seqlevels(p), function(x) {
    message(x)
    f <- all(seqnames(p) == x)
    message(x, ' -> ', sum(f))
    return(f)
})
same.chr <- Reduce('|', flt)


seq.pos <- all(strand(p) == '+')
seq.neg <- all(strand(p) == '-')
same.seq <- (seq.pos | seq.neg)


starts <- start(p)
st <- do.call(rbind, starts)
dup <- duplicated(st)

