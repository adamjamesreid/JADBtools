filterNonMappabeAndBlacklistedPeaks <- function(peaks, scoreTreshold=100, pValueTreshold=20) {
    
    library(rtracklayer)
    
    blacklist <- import.bed('https://gist.githubusercontent.com/Przemol/ef62ac7ed41d3a84ad6c478132417770/raw/56e98b99e6188c8fb3dfb806ff6f382fe91c27fb/CombinedBlacklists.bed')
    #nonmappable <- import.bed('https://gist.githubusercontent.com/Przemol/ef62ac7ed41d3a84ad6c478132417770/raw/56e98b99e6188c8fb3dfb806ff6f382fe91c27fb/non_mappable.bed')
    #filter <- reduce(c(blacklist, nonmappable))
    
    filter <- blacklist
    
    
    
    
    out <- peaks[!peaks %over% filter]
    #out <- out[out$V5 >= scoreTreshold]
    #out <- out[out$V5 >= pValueTreshold]
    
    
    message('Peaks in input: \t', length(peaks))
    message('Peaks in output: \t', length(out))
    message('Peaks in filtered out: \t', length(peaks) - length(out))
    return(out)
    
}

concave_flt <- function(file_concave, file_idr, bw, curvature_idx=500L) {
    message(file_concave, ' -vs- ', file_idr)
    extraCols_narrowPeak <- c(
        signalValue = "numeric", pValue = "numeric", qValue = "numeric", peak = "integer",
        X11 = "character", X12 = "character", X13 = "character", X14 = "character", X15 = "character", 
        X16 = "character", X17 = "character", X18 = "character", X19 = "character", X20 = "character"
    )
    
    concave <- rtracklayer::import.bed(file_concave)
    idr <- rtracklayer::import.bed(file_idr, extraCols = extraCols_narrowPeak)[,1:6]
    ac <- length(concave)
    
    concave <- concave[concave %over% idr]
    
    blacklist <- import.bed('https://gist.githubusercontent.com/Przemol/ef62ac7ed41d3a84ad6c478132417770/raw/56e98b99e6188c8fb3dfb806ff6f382fe91c27fb/CombinedBlacklists.bed')
    bl <- concave %over% blacklist
    message(sum(bl), ' blacklisted')
    concave <- concave[!bl]
    
    message('Rejecting ' , sum(!as.numeric(concave$name) > curvature_idx), ' peaks based on curvature idx')
    concave <- concave[as.numeric(concave$name) > curvature_idx]
    
    #rtracklayer::import.bw(bw, selection=rtracklayer::BigWigSelection(r1[1:10]), as='RleList')[r1[1:10]] %>% max()
    score <- unlist(summary(rtracklayer::BigWigFile(basename(bw)), concave, type="max"))$score
    concave$score <- score
    
    neighborhood <- setdiff(flank(concave, 2000, both = TRUE), flank(concave, 1000, both=TRUE))
    
    
    neighborhood <- neighborhood[!neighborhood %over% blacklist]
    n_score <-  unlist(summary(rtracklayer::BigWigFile(basename(bw)), neighborhood, type="max"))$score
    #browser()
    
    neighborhood_idx <- quantile(n_score, 0.001, na.rm=TRUE)
    
    message('Neighborhood score 0.001 quantile: ', neighborhood_idx, ', rejecting ', sum(!concave$score > neighborhood_idx), ' peaks')
    #cat(quantile(n_score), '\n')
    concave <- concave[concave$score > neighborhood_idx]
    
    message('concave / idr / concave_filt: ', ac, ' / ', length(idr), ' / ', length(concave))
    message('cri: ', sum(idr %over% concave)/length(idr) )
    rtracklayer::export.bed(concave, gsub(
        'concave.bed', 
        paste0('concave_flt_ci=', curvature_idx , 'sc=', round(neighborhood_idx), '.bed'), 
        file_concave
    ))
    
    return(invisible(concave))
    
}

#' Single-pass IDR from replicate
#'
#' @param id ID of second experiment in replicate
#' @param processing 
#' @param curvature_idx 
#'
#' @return GRanges
#' @export
#'
idr_concave_on_replicate <- function(id, processing='BEADSQ10NU', curvature_idx=500L) {
    bw <- getFilePath(id, processing = processing, scale = 'lin')
    idr <- getFilePath(id, processing = 'IDR2')
    if( !file.exists(basename(bw)) ) download.file(bw, basename(bw))
    if( !file.exists(basename(idr)) ) download.file(idr, basename(idr))
    
    if( !file.exists(gsub('.bw$', '_concave.bed', basename(bw))) )
        system(sprintf(
            '/Users/przemol/miniconda2/bin/python /usr/local/bin/concave_regions %s > %s',
            basename(bw),
            gsub('.bw$', '_concave.bed', basename(bw))
        ))
    
    peak <- concave_flt(
        gsub('.bw$', '_concave.bed', basename(bw)), 
        basename(idr), bw=basename(bw),
        curvature_idx = curvature_idx
    )
    
    return(peak)
} 


#' Two-pass IDR procedure with oracle peak call set
#'
#' @param r1 ID of first experiment in replicate
#' @param r2 ID of second experiment in replicate
#' @param processing Which tracks to use (default: BEADSQ10NU)
#' @param curvature_idx - curvature index treshold (default: 100)
#'
#' @return NULL
#' @export
#'
idr_concave_two_pass_with_oracle <- function(r1, r2, processing='BEADSQ10NU', curvature_idx=100) {
    
    require(JADBtools)
    require(GenomicRanges)
    require(rtracklayer)
    require(magrittr)
    require(dplyr)
    require(readr)
    
    # Concave calls on BW
    
    bw1 <- getFilePath(r1, processing = processing, scale = 'lin')
    bw2 <- getFilePath(r2, processing = processing, scale = 'lin')
    
    if( !file.exists(basename(bw1)) ) download.file(bw1, basename(bw1))
    if( !file.exists(basename(bw2)) ) download.file(bw2, basename(bw2))
    
    if( !file.exists(gsub('.bw$', '_concave.bed', basename(bw1))) )
        system(sprintf(
            '/Users/przemol/miniconda2/bin/python /usr/local/bin/concave_regions %s > %s',
            basename(bw1),
            gsub('.bw$', '_concave.bed', basename(bw1))
        ))
    
    if( !file.exists(gsub('.bw$', '_concave.bed', basename(bw2))) )
        system(sprintf(
            '/Users/przemol/miniconda2/bin/python /usr/local/bin/concave_regions %s > %s',
            basename(bw2),
            gsub('.bw$', '_concave.bed', basename(bw2))
        ))
    
    
    # Concave calls on combined BW
    
    out_rep <- combineReps(r1 = bw1, r2 = bw2, outdir = '.')
    system(sprintf(
        '/Users/przemol/miniconda2/bin/python /usr/local/bin/concave_regions %s > %s',
        basename(out_rep$out),
        gsub('.bw$', '_concave.bed', basename(out_rep$out))
    ))
    
    
    
    # IDR on MACS
    p1 <- getFilePath(r1, processing = 'peakCalls')
    p2 <- getFilePath(r2, processing = 'peakCalls')
    if( !file.exists(basename(p1)) ) download.file(p1, basename(p1))
    if( !file.exists(basename(p2)) ) download.file(p2, basename(p2))
    
    
    if( !file.exists(gsub('.bw$', '_idr.bed', basename(bw2))) )
        system(sprintf(
            '/Users/przemol/miniconda2/envs/py36/bin/idr --samples %s %s --output-file %s',
            basename(p1), basename(p2),
            gsub('.bw$', '_idr.bed', basename(bw2))
        ))
    
    
    
    
    peak1 <- concave_flt(
        gsub('.bw$', '_concave.bed', basename(bw1)), 
        gsub('.bw$', '_idr.bed', basename(bw2)),
        bw = basename(bw1),
        curvature_idx = 0
    )
    
    peak2 <- concave_flt(
        gsub('.bw$', '_concave.bed', basename(bw2)), 
        gsub('.bw$', '_idr.bed', basename(bw2)),
        bw = basename(bw2),
        curvature_idx = 0
    )
    
    peakR <- concave_flt(
        gsub('.bw$', '_concave.bed', basename(out_rep$out)), 
        gsub('.bw$', '_idr.bed', basename(bw2)),
        bw = basename(out_rep$out),
        curvature_idx = 0
    )
    
    peakR$name <- peak1$name <- peak2$name <- NULL
    
    
    
    as.data.frame(peak1) %>% tbl_df %>% transmute(seqnames,  start,    end, st='.', score, '.', a=0, b=0, c=0)
    
  
    as.data.frame(peak1) %>% tbl_df %>% 
        transmute(seqnames, start, end, st='.', score, '.', a=0, b=0, c=0) %>% 
        write_delim('p1.bed', col_names = FALSE, delim = '\t')

    as.data.frame(peak2) %>% tbl_df %>% 
        transmute(seqnames, start, end, st='.', score, '.', a=0, b=0, c=0) %>% 
        write_delim('p2.bed', col_names = FALSE, delim = '\t')
    
    
    as.data.frame(peakR) %>% tbl_df %>% 
        transmute(seqnames, start, end, st='.', score, '.', a=0, b=0, c=0) %>% 
        write_delim('pR.bed', col_names = FALSE, delim = '\t')
    
    ann <- rbind(rbeads:::ParseName(basename(bw1)), rbeads:::ParseName(basename(bw2)))
    out <- ann[c(-3, -9:-13)][1,] %>% paste(collapse = '_') %>% paste0('_', paste(ann$ContactExpID, collapse = '^'), '_IDRconcave.narrowPeak')
    system(sprintf(
        '/Users/przemol/miniconda2/envs/py36/bin/idr --input-file-type bed --rank score --plot --samples %s %s --peak-list %s --output-file %s --idr-threshold 0.05',
        'p1.bed', 'p2.bed', 'pR.bed', 
        out
    ))
} 

