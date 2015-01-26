require()

setwd('/Volumes/raid0/Rtemp/shortRNA')
p1 <- readGAlignmentPairs('Yan01_S5_L001andL002.bam')
cov <- coverage(p1)
export.bw(cov, 'Yan01_S5_L001andL002.bw')

pp <- granges(p1)
export.bed(pp, 'Yan01_S5_L001andL002.bed')


require(Rsamtools)
require(rtracklayer)
what <- c("rname", "pos", "isize", 'strand', 'mapq', 'seq')

flag <- scanBamFlag(isUnmappedQuery = FALSE, isProperPair = TRUE, isMinusStrand = FALSE, isMateMinusStrand = TRUE)
param <- ScanBamParam(what = what, flag=flag)