gistToFile <- function(id, ...) {
    library(httr) 
    tmp <- tempfile(...)
    txt <- httr::content(GET(
        sprintf("https://api.github.com/gists/%s", id)
    ))$files[[1]]$content
    cat(txt, file = tmp)
    return(tmp)
}

annotatePeaks <- function(bed) {
    
    require(ChIPseeker)
    
    db <- system.file('anno/WBcel235_EnsemblGenes77_TxDb.sqlite', package = 'JADBtools')
    txdb <- loadDb(db)
    seqlevelsStyle(txdb) <- 'UCSC'
    chain <- import.chain(gistToFile('8e4a0022eaa4a9d7f97a'))
    
    
    bed <- '/Volumes/raid0/Rtemp/AMA1^Q2357_aa63^NONE^N2^L4_aligned^NA^NA_AA486^F3931279_EnrichedRegions.bed'
    
    peak_ce10 <- readPeakFile(bed)
    
    peak <- unlist(reduce(liftOver(peak_ce10, chain), min.gapwidth=50))
    
    peakAnno <- annotatePeak(peak, tssRegion=c(-500, 500), TxDb=txdb)
    plotAvgProf2(peak, TxDb=txdb, upstream=3000, downstream=3000, xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
    peakHeatmap(peak, TxDb=txdb, upstream=3000, downstream=3000)
}