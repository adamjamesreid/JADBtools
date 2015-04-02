require(Rsamtools)
require(rtracklayer)
require(GenomicAlignments)
require(parallel)
require(BiocParallel)

info <- function(z) {
    sig<-sum(z)
    len<-length(z)
    message(sys.call()[2], ': ',
        'F=', len-sig, ' [', round((len-sig)/len*100, 2), 
        '%]  T=', sig, ' [', round(sig/len*100, 2),
        '%]  (A=', len, ')')
}

fff <- dir(pattern = '*bam$')

ff <- fff[2]


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

#Long
starts <- start(p)
st <- do.call(rbind, starts)
dup <- duplicated(st)

#mapq10

mapq <- elementMetadata(p@unlistData)$mapq
m1 <- mapq[c(TRUE, FALSE)]
m2 <- mapq[c(FALSE, TRUE)]
mapq10 <- (m1>=10) & (m2>=10)

informMe <- function() {
    ##logicl vectors of length N (aligned reads)
    info(pair) #indicates paired reads, i.e. not singletons 
    
    ##logicl vectors of length N-singletons (!pairs)
    info(same.chr) #indicates reads mapped to same chromosme (proper)
    info(same.seq) #indicates reads mapped to same strand (not proper)
    info(dup) #indicates duplicates
    info(mapq10) #both reads have mapping quality >=10
}

infoList <- list(pair=pair, same.chr=same.chr, same.seq=same.seq, dup=dup, mapq10=mapq10)

#How many pairs where each read maps to same chromosome and same strand?
info(same.chr & same.seq & (!dup))

# > informMe()
# pair: F=2638436 [4.66%]  T=54016519 [95.34%]  (A=56654955)
# same.chr: F=660457 [1.22%]  T=53356062 [98.78%]  (A=54016519)
# same.seq: F=52746822 [97.65%]  T=1269697 [2.35%]  (A=54016519)
# dup: F=23388366 [43.3%]  T=30628153 [56.7%]  (A=54016519)
# mapq10: F=5108610 [9.46%]  T=48907909 [90.54%]  (A=54016519)
# > info(same.chr & same.seq & (!dup))
# same.chr & same.seq & (!dup): F=53416566 [98.89%]  T=599953 [1.11%]  (A=54016519)


