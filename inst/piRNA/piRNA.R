
require(magrittr)
fa <- readDNAStringSet('~/All_piRNAs_new_reference.fasta')

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

names(ce10) <- sapply(strsplit(names(fa), ' '), '[[', 1) 
names(ce10) <- NULL
score(ce10) <- names(fa)  %>% gsub('.+_(.+)-(.+)_(.+)_(.+)', '\\4', .)  %>% as.numeric

export.bed(ce10, 'All_piRNAs_ce10.bed')

indi <- ce10[ce10$score < 0]
dep <- ce10[ce10$score >= 0]

export.bed(indi, 'indi_piRNAs_ce10.bed')
export.bed(dep, 'dep_piRNAs_ce10.bed')

names(ce10) <- NULL
score(ce10) <- NULL
export.bed(ce10, 'All_piRNAs_ce10_noAnno.bed')

export.bed(ce10[order(ce10$score, decreasing = TRUE)], 'All_piRNAs_ordered_by_scorece10.bed')

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

motif_bed2fa <- function(bed, ups=1000, dns=100, tss=TRUE, fa=gsub('bed', 'fa', bed)) {
    message(fa)
    reg <- import.bed(bed)
    
    if(tss==TRUE) {
        fa_promoters <- getSeq(Celegans, promoters(reg, upstream = ups, downstream = dns))
    } else {
        rr <- promoters( resize(reg, 1, fix='end'),  upstream = ups, downstream = dns ) 
        fa_promoters <- getSeq(Celegans, rr)
    }
    
    names(fa_promoters) <- reg$name
    writeXStringSet(fa_promoters, fa)
    return(fa_promoters)
}

motif_peaks2fa <- function(peaks, ups=1000, dns=100, tss=FALSE, fa=gsub('narrowPeak', 'fa', peaks)) {
    message(fa)
    library(BSgenome.Celegans.UCSC.ce10)
    library(rtracklayer)
    library(Biostrings)
    
    reg <- ChIPseeker::readPeakFile(peaks, header=FALSE) #http://genome.ucsc.edu/FAQ/FAQformat.html#format12
 
    fa_reg <- getSeq(Celegans, reg)
    names(fa_reg) <- reg$V4
    writeXStringSet(fa_reg, fa)
    invisible( return(fa_reg) )
}

data("repeatModel")
grangeslist2fa <- function(gr) {
    library(BSgenome.Celegans.UCSC.ce10)
    library(rtracklayer)
    library(Biostrings)
    require(GenomicRanges)
    
    sapply(names(gr), function(nam) {
        message(nam)
        reg <- gr[[nam]]
        fa_reg <- getSeq(Celegans, reg)
        names(fa_reg) <- paste0(seqnames(reg),'-',start(reg),':',end(reg))
        writeXStringSet(fa_reg, paste0(nam, '.fa'))
    })
}

data("modencodetfbs")
names(modencodetfbs)  %>% gsub('spp.idrOptimal.bf.ce.[^_]+_', '', .)  %>%  gsub('.bam.unique.tagAlign.narrowPeak.gz', '', .)  -> nn
names(modencodetfbs) <- nn

elementLengths(modencodetfbs[grep('HPL-2', nn)])

#macs2 callpeak -t TOFU5^AB290_TOFU5-GFP^NA^NA^NA_aligned^NA^NA_AK004^F0f32716.bam -c antiGFP^GFP-trap_em01^E^N2^youngAdult_aligned^NA^NA_AA243^F1b05732.bam --nomodel --extsize 150 -f BAM -g ce -n TOFU5_AB290_TOFU5-GFP_AK004 -q 0.01

##PEAKS
peak <- lapply(dir(pattern = 'narrowPeak'), function(x) {
    message(x)
    ChIPseeker::readPeakFile(x, header=FALSE)
})
names(peak) <- gsub('.narrowPeak', '', dir(pattern = 'narrowPeak'))
peaks <- GRangesList(peak)
grangeslist2fa(peaks)

##SUMMITS
summits <- lapply(dir(pattern = 'bed'), function(x) {
    message(x); resize(import.bed(x), 200, fix = 'center')
})
names(summits) <- gsub('.bed', '', dir(pattern = 'bed'))
grangeslist2fa(summits)


oneLines = memeOutput[grep("^MOTIF \\d+", memeOutput)]

