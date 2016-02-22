vmatchPattern(fa[[4]], Celegans)

lst <- pbapply::pblapply(fa[1:100], function(x) vmatchPattern(x, Celegans) )
aln <- unlist(GRangesList(lst)[1:30])

chr <- names(fa)  %>% gsub('.+_(.+)-(.+)_(.+)_.+', '\\1', .)
pos <- names(fa)  %>% gsub('.+_(.+)-(.+)_(.+)_.+', '\\2', .)
str <- names(fa)  %>% gsub('.+_(.+)-(.+)_(.+)_.+', '\\3', .)
wid <- width(fa)

gr <- GRanges(
    paste0('chr', gsub('MtDNA', 'M', chr)),
    IRanges(as.numeric(pos)+1, width = wid),
    strand = c('-', '+')[as.numeric(factor(str))]
)

#download.file('http://hgdownload.cse.ucsc.edu/goldenPath/ce6/liftOver/ce6ToCe10.over.chain.gz', 'ce6ToCe10.over.chain.gz')
#chain <- import.chain('ce6ToCe10.over.chain')

chain <- import.chain(system.file('anno/WBcel235toCe10.chain', package='JADBtools'))
ce10 <- unlist(liftOver(gr, chain))
start(aln) - start(ce10[1:30]) #sanity check

names(ce10) <- names(fa)
score(ce10) <- names(fa)  %>% gsub('.+_(.+)-(.+)_(.+)_(.+)', '\\4', .)  %>% as.numeric

export.bed(ce10, 'All_piRNAs_ce10.bed')

names(ce10) <- NULL
score(ce10) <- NULL
export.bed(ce10, 'All_piRNAs_ce10_noAnno.bed')

fa_out <- getSeq(Celegans, promoters(ce10, 200, 0))
writeXStringSet(fa_out, 'All_piRNAs_200bp_upstream_e10.fa')

peakloc <- resize(shift(ce10, -90), 40)
fa_out <- getSeq(Celegans, peakloc)
writeXStringSet(fa_out, 'All_piRNAs_70bp_upstream_40bp_centered_e10.fa')



fa_promoters <- getSeq(Celegans, promoters(ce10, 200, 0)[1:200])



allI <- unlist( intWBcel235 )
allIlo <- reduce(liftOver(gr, chain), min.gapwidth=50)
names(allIlo) <- names(allI)


#### general impl

bed2fa <- function(bed, ups=1000, dns=100, tss=TRUE, fa=gsub('bed', 'fa', bed)) {
    message(fa)
    reg <- import.bed(bed)
    fa_promoters <- getSeq(Celegans, promoters(reg, upstream = ups, downstream = dns))
    writeXStringSet(fa_promoters, fa)
    return(fa_promoters)
}