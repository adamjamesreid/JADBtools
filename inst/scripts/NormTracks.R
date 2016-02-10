Use short chromatin RNA from 2015-05-20:
    
JA26r - rYD014
RC8 - rYD015
RC10 - rYD016

Long chromatin RNA also from 2015-05-20
JA26r - rYD017
RC8 - rYD018
RC10 - rYD019

c()
paths <- 

require(magrittr)
require(Rsamtools)

bw <- paste0('rYD0', 14:19) %>% sapply(JADBtools::getFilePath, format = 'bw')

bam <- paste0('rYD0', 14:19) %>% sapply(JADBtools::getFilePath, format = 'bam')

countBam(file, index=file, ..., param=ScanBamParam())

library(GenomicFeatures)
require(rtracklayer)
gistToFile <- function(id, ...) {
    require(httr)
    tmp <- tempfile(...)
    txt <- httr::content(GET(
        sprintf("https://api.github.com/gists/%s", id)
    ))$files[[1]]$content
    cat(txt, file = tmp)
    return(tmp)
}
rdb2 <- makeTranscriptDbFromBiomart(dataset="celegans_gene_ensembl", 
                                    filters=list( "biotype"="rRNA" ), biomart="ENSEMBL_MART_ENSEMBL", 
                                    host='www.ensembl.org')
seqlevelsStyle(rdb2) <- 'UCSC'

#Extract genes
WBcel235 <- genes(rdb2)

#convert co ce10

chain <- import.chain(gistToFile('307189b979b7bcf68cb6'))
ce10 <- unlist(reduce(liftOver(WBcel235, chain), min.gapwidth=50))
elementMetadata(ce10) <- elementMetadata(WBcel235)
names(ce10) <- names(WBcel235)
seqinfo(ce10) <- SeqinfoForUCSCGenome('ce10')[seqlevels(ce10)]

require(GenomicAlignments)


alns <- sapply( bam, function(x) readGAlignments(BamFile(x), param = ScanBamParam(what='mapq')) )

normBW <- function(EXP_ID) {
    n_rrna <- sum(alns[[EXP_ID]] %over% ce10)
    n_norm <- length(alns[[EXP_ID]]) - n_rrna
    message(EXP_ID, ' => Total: ', length(alns[[EXP_ID]]), '; rRNA: ', n_rrna, '; non-rRNA: ', n_norm)
    
    bw[,EXP_ID][[1]]  %>%  BigWigFile  %>% import.bw(as='Rle') -> track
    qt <- quantile(unlist(track[track>0]), .99)
    message(basename(bw[,EXP_ID][[1]]), ': cutoff=', qt)
    track[track >= qt] <- 0
    out <- track / (n_norm/10^6)
    out_name <- paste0( paste(strsplit(basename(bw[,EXP_ID][[1]]), '_')[[1]][1:6], collapse = '_'), '_readsPerMillionNorRBA_top01pctMassked_', paste(strsplit(basename(bw[,EXP_ID][[1]]), '_')[[1]][9:10], collapse = '_'))
    export.bw(out, out_name)
    summary(unlist(out))
}


sapply(dir(pattern='bw'), function(x) { 
    tr <- import.bw(x, as='Rle')
    sum(sum(tr)) 
    })


getTop <- function(EXP_ID) {
    
    mask1 <- alns[[EXP_ID]] %over% ce10
    n_rrna <- sum(mask1)
    n_norm <- length(alns[[EXP_ID]]) - n_rrna
    message(EXP_ID, ' => Total: ', length(alns[[EXP_ID]]), '; rRNA: ', n_rrna, '; non-rRNA: ', n_norm)
    
    track <- coverage(alns[[EXP_ID]])
    track_nonzero <- unlist(track[track>0])
    qt <- quantile(track_nonzero, .99)
    CV <- sd(track_nonzero)/mean(track_nonzero)
    qq <- quantile(track_nonzero)
    QCD <- (qq[4]-qq[2])/(qq[4]+qq[2])

    
    message(EXP_ID, ': cutoff=', qt)
    more <- as(track >= qt, 'GRanges')
    more <- more[more$score]
    
    mask2 <- alns[[EXP_ID]] %over% more
    
    message(EXP_ID, ' => Total: ', length(alns[[EXP_ID]]), '; top 1% reads: ', sum(mask2), ' CV: ', round(CV,2), ' QCD: ',QCD)
    
    score(more) <- NULL
    export.bed(more, paste0(EXP_ID, '_top1%regions.bed'))
    
    out_name <- paste0( paste(strsplit(basename(bam[EXP_ID][[1]]), '_')[[1]][1:6], collapse = '_'), '_readsPerMillion_top01pctANDrRNAreadsMassked_', paste(strsplit(basename(bam[EXP_ID][[1]]), '_')[[1]][9:10], collapse = '_'))
    out_name <- gsub('bam', 'bw', out_name)
    mask <- (!mask1) & (!mask2)
    
    message('Non-masked reads: ',  (sum(mask)/1e6), ' million')
    out <- coverage(alns[[EXP_ID]][mask]) / (sum(mask)/1e6)
    export.bw(out, out_name)
    summary(unlist(out))

}

getAll <- function(EXP_ID) {
    
    message(EXP_ID, ': ', length(alns[[EXP_ID]])/1e6)
    out_name <- paste0( paste(strsplit(basename(bam[EXP_ID][[1]]), '_')[[1]][1:6], collapse = '_'), '_PILEUP_', paste(strsplit(basename(bam[EXP_ID][[1]]), '_')[[1]][9:10], collapse = '_'))
    out_name <- gsub('bam', 'bw', out_name)

    out <- coverage(alns[[EXP_ID]]) / (length(alns[[EXP_ID]])/1e6)
    export.bw(out, out_name)
    
}





