l <- as(sort(vmatchPattern('CGAGCAATGGTCCAAAACTT', Celegans))[c(2:6, 1)], 'GAlignments')
r <- as(sort(vmatchPattern('GGACAGTGTACATTCGGCCT', Celegans)), 'GAlignments')
r
GAlignmentPairs(l, r)
GAlignmentPairs(r, l)
pcr1 <- GAlignmentPairs(r, l)
pcr1
granges(pcr1)


l <- as(sort(vmatchPattern('gatcatctgggaatctgcgt', Celegans)), 'GAlignments')
r <- as(sort(vmatchPattern('aagcatcggtgaaatggaac', Celegans)), 'GAlignments')
r
GAlignmentPairs(l, r)
GAlignmentPairs(r, l)
pcr1 <- GAlignmentPairs(r, l)
pcr1
granges(pcr1)
unlistedRepeatModel[unlistedRepeatModel  %over% reduce(granges(pcr1), ignore.strand=T)]


doPCR <- function(l, r) {
    ll <- as(sort(vmatchPattern(l, Celegans)), 'GAlignments')
    rr <- as(sort(vmatchPattern(r, Celegans)), 'GAlignments')[c(2:6, 1)]
    pcr1 <- GAlignmentPairs(rr, ll)
    gr <- granges(pcr1)
    print(unlistedRepeatModel[unlistedRepeatModel  %over% reduce(granges(pcr1), ignore.strand=T)])
    return(gr)
}
