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
    file.remove(c("p1.bed", "p2.bed", "pR.bed"))
} 


#' Overlap intersect of two peak calls
#'
#' @param r1 ID of first experiment in replicate
#' @param r2 ID of second experiment in replicate
#' @param processing Which tracks to use (default: BEADSQ10NU)
#' @param dc Create new directory
#'
#' @return NULL
#' @export
#'
overlap_intersect <- function(r1, r2, processing='BEADSQ10NU', dc=TRUE) {
    
    require(JADBtools)
    require(GenomicRanges)
    require(rtracklayer)
    require(magrittr)
    require(dplyr)
    require(readr)
    require(tracktables)
    
    
    # IDR on MACS
    p1 <- getFilePath(r1, processing = 'peakCalls')
    p2 <- getFilePath(r2, processing = 'peakCalls')
    
    ll <- strsplit(c(basename(p1), basename(p2)), '_|\\^|\\.')
    repp <- !ll[[1]] == ll[[2]]
    construct <- ll[[1]]
    construct[repp] <- gsub(' ', '|', paste(ll[[1]][repp], ll[[2]][repp]))
    #out <- paste0(paste0(construct[-length(construct)], collapse='_'), '.', construct[length(construct)])
    out <- paste0(paste0(construct[-length(construct)], collapse='_'), '.', 'bed')
    
    
    if(dc) {
        dirn <- gsub('_(F|E)_.+', '', gsub('\\|', '_', out))
        dir.create(dirn)
        setwd(dirn)
    }
    
    if( !file.exists(basename(p1)) ) download.file(p1, basename(p1))
    if( !file.exists(basename(p2)) ) download.file(p2, basename(p2))
    
    # Get BW for browsing
    bw1 <- getFilePath(r1, processing = processing, scale = 'lin')
    bw2 <- getFilePath(r2, processing = processing, scale = 'lin')
    if( !file.exists(basename(bw1)) ) download.file(bw1, basename(bw1))
    if( !file.exists(basename(bw2)) ) download.file(bw2, basename(bw2))
    
    extraCols_narrowPeak <- c(
        signalValue = "numeric", pValue = "numeric", qValue = "numeric", peak = "integer"
    )
    g1 <- import.bed(basename(p1), extraCols=extraCols_narrowPeak)
    g2 <- import.bed(basename(p2), extraCols=extraCols_narrowPeak)
    
    int <- IRanges::intersect(g1, g2)
    
    #int2 <- as(int, 'GRangesList')
    #g22 <- as(g1, 'GRangesList')
    #g12 <- as(g2, 'GRangesList')
    
    #cbind(subsetByOverlaps(g12, int2), subsetByOverlaps(g22, int2))
    
    blacklist <- import.bed('https://gist.githubusercontent.com/Przemol/ef62ac7ed41d3a84ad6c478132417770/raw/56e98b99e6188c8fb3dfb806ff6f382fe91c27fb/CombinedBlacklists.bed')
    nonmappable <- import.bed('https://gist.githubusercontent.com/Przemol/ef62ac7ed41d3a84ad6c478132417770/raw/56e98b99e6188c8fb3dfb806ff6f382fe91c27fb/non_mappable.bed')
    filter <- reduce(c(blacklist, nonmappable))
    
    
    chain <- import.chain('/Users/przemol/Downloads/ce10ToCe11.over.chain')
    filter_ce11 <- liftOver(filter, chain) %>% reduce(min.gapwidth=50) %>%  unlist


    int2 <- int[!int %over% filter_ce11]
    
    message(sprintf(
        'P1=%s; P2=%s; I=%s; SF=%s', 
        length(g1), length(g2), length(int), length(int2)
    ))
    export.bed(int2,out)
    
    makeTT(bw1,bw2,p1,p2,out)
    
    if(dc) {
        setwd('..')
    }
} 


overlap_intersect_nd <- function(
    r1, r2, processing='BEADSQ10NU', dc=FALSE, filter_map=TRUE, dbname=TRUE, 
    mode='intersection', loadIGV=FALSE
) {
    
    require(JADBtools)
    require(GenomicRanges)
    require(rtracklayer)
    require(magrittr)
    require(dplyr)
    require(readr)
    require(tracktables)
    
    
    # IDR on MACS
    p1 <- getFilePath(r1, processing = 'peakCalls')
    p2 <- getFilePath(r2, processing = 'peakCalls')
    
    
    if(dbname) {
        anno <- as.data.frame(t(sapply(c(basename(p1), basename(p2)), rbeads:::ParseName)))
        same <- anno[1,c('Factor', 'Strain', 'Stage', 'Processing', 'Scale', 'Resolution')]
        out <- paste0(paste0(unlist(same), collapse='_'), '_', paste(anno$ExtractID, collapse = '|'), '_', paste(anno$ContactExpID, collapse = '|'), '_', mode, '.bed')
        
    } else {
        ll <- strsplit(c(basename(p1), basename(p2)), '_|\\^|\\.')
        repp <- !ll[[1]] == ll[[2]]
        construct <- ll[[1]]
        construct[repp] <- gsub(' ', '|', paste(ll[[1]][repp], ll[[2]][repp]))
        out <- paste0(paste0(construct[-length(construct)], collapse='_'), '.', 'bed')
        
    }
    
    
    if(dc) {
        dirn <- gsub('_(F|E)_.+', '', gsub('\\|', '_', out))
        dir.create(dirn)
        setwd(dirn)
    }
    
    #if( !file.exists(basename(p1)) ) download.file(p1, basename(p1))
    #if( !file.exists(basename(p2)) ) download.file(p2, basename(p2))
    
    # Get BW for browsing
    bw1 <- getFilePath(r1, processing = processing, scale = 'lin')
    bw2 <- getFilePath(r2, processing = processing, scale = 'lin')
    #if( !file.exists(basename(bw1)) ) download.file(bw1, basename(bw1))
    #if( !file.exists(basename(bw2)) ) download.file(bw2, basename(bw2))
    
    extraCols_narrowPeak <- c(
        signalValue = "numeric", pValue = "numeric", qValue = "numeric", peak = "integer"
    )
    g1 <- import.bed(url(p1), extraCols=extraCols_narrowPeak)
    g2 <- import.bed(url(p2), extraCols=extraCols_narrowPeak)
    
    if(mode=='intersection') {
        int <- IRanges::intersect(g1, g2)
    } else if(mode=='union') {
        int <- IRanges::union(g1, g2)
    }
    
    #int2 <- as(int, 'GRangesList')
    #g22 <- as(g1, 'GRangesList')
    #g12 <- as(g2, 'GRangesList')
    
    #cbind(subsetByOverlaps(g12, int2), subsetByOverlaps(g22, int2))
    
    blacklist <- import.bed('https://gist.githubusercontent.com/Przemol/ef62ac7ed41d3a84ad6c478132417770/raw/56e98b99e6188c8fb3dfb806ff6f382fe91c27fb/CombinedBlacklists.bed')
    nonmappable <- import.bed('https://gist.githubusercontent.com/Przemol/ef62ac7ed41d3a84ad6c478132417770/raw/56e98b99e6188c8fb3dfb806ff6f382fe91c27fb/non_mappable.bed')
    
    if(filter_map) {
        filter <- reduce(c(blacklist, nonmappable))
    } else {
        filter <- blacklist
    }
    
    download.file(
        'https://gist.githubusercontent.com/Przemol/e9f1a3a5053619e69fafbd46759a17e4/raw/c69209270843e2724dc8c6ea9715ddd52930f41c/ce10ToCe11.over.chain', 
        (tempfile() -> chainf), quiet = TRUE
    )
    chain_ce10ToCe11 <- import.chain(chainf)
    filter_ce11 <- liftOver(filter, chain_ce10ToCe11) %>% reduce(min.gapwidth=50) %>%  unlist
    
    
    int2 <- int[!int %over% filter_ce11]
    
    message(sprintf(
        'P1=%s; P2=%s; I=%s; SF=%s', 
        length(g1), length(g2), length(int), length(int2)
    ))
    export.bed(int2,out)
    
    makeTT(bw1,bw2,p1,p2,out, local = FALSE, loadIGV=loadIGV)
    
    if(dc) {
        setwd('..')
    }
    return(out)
} 


makeTT <- function(bw1,bw2,p1,p2,out, local=TRUE, loadIGV=FALSE) {
 
    if(!local) {
        bigwigs <- c(
            bw1,
            bw2
        )
        intervals <- c(
            p1,
            p2,
            out
        )
    } else {
        bigwigs <- c(
            file.path(getwd(),basename(bw1)),
            file.path(getwd(),basename(bw2))
        )
        intervals <- c(
            file.path(getwd(),basename(p1)),
            file.path(getwd(),basename(p2)),
            file.path(getwd(),out)
        )
    }
    
    
    bigWigMat <- cbind(
        gsub("\\^[^\\^]+\\^[^\\^]+\\^[^\\^]+_[^_]+_[^_]+$","",basename(bigwigs)),
        bigwigs
    )
    intervalsMat <- cbind(
        gsub("\\^[^\\^]+\\^[^\\^]+\\^[^\\^]+_[^_]+_[^_]+$","",basename(intervals)),
        intervals
    )
    
    FileSheet <- merge(bigWigMat,intervalsMat,all=TRUE)
    FileSheet <- as.matrix(cbind(FileSheet,NA))
    colnames(FileSheet) <- c("SampleName","bigwig","interval","bam")
    
    SampleSheet <- cbind(
        as.vector(FileSheet[,"SampleName"]),
        c("R1","R2","RC")
    )
    colnames(SampleSheet) <- c("SampleName","Rep")
    
    MakeIGVSampleMetadata(SampleSheet, FileSheet, getwd())
    sessionxml <- MakeIGVSessionXML(
        FileSheet, getwd(), out, "ce11", locusName = "All", 
        colourBy = apply(col2rgb(c('darkred', 'darkgreen', 'black')), 2, function(x) paste0(x, collapse = ",")), 
        igvParams = igvParam(
            bigwig.color = 'black', bigwig.altColor = 'red',
            interval.color = 'black', interval.altColor = 'red',
            bigwig.autoScale = "false", bigwig.minimum = 0, bigwig.maximum = 10
        )
    )

    if(loadIGV) httr::GET(sprintf('http://localhost:60151/load?file=%s', sessionxml))
    
    #browseURL(sprintf('http://localhost:60151/load?file=%s', sessionxml))
    
    # HTMLreport <- maketracktable(
    #     fileSheet=FileSheet,
    #     SampleSheet=SampleSheet,
    #     filename="IGVExample_ce11.html",
    #     basedirectory=getwd(),
    #     genome="ce11"
    # )
    
    
}


makeTTall <- function(bw1,bw2,p1,p2,out) {
    
    bigwig <- dir(patt='\\^(F|E)\\^.+bw')
    
    interval <- dir(patt='\\^(F|E)\\^.+narrowPeak$')
    rr <- dir(patt='_(F|E)_.+bed$')
    
    FileSheet <- cbind(
        SampleName=gsub("_[^_]+_[^_]+_[^_]+$","",basename(bigwigs)),
        bigwig=bigwig,
        interval=interval,
        bam=NA
    ) %>% tbl_df
    
    SampleSheet <- cbind(
        SampleName=as.vector(FileSheet[,"SampleName"]),
        Rep=c("R1","R2")
    )
    
    MakeIGVSampleMetadata(SampleSheet, FileSheet, getwd())
    sessionxml <- MakeIGVSessionXML(
        FileSheet, getwd(), 'ALL', "ce11", locusName = "All", 
        colourBy = apply(col2rgb( rep(c('darkred', 'darkgreen'), length(bigwigs)/2)), 2, function(x) paste0(x, collapse = ",")), 
        igvParams = igvParam(
            bigwig.color = 'black', bigwig.altColor = 'red',
            interval.color = 'black', interval.altColor = 'red',
            bigwig.autoScale = "false", bigwig.minimum = 0, bigwig.maximum = 10
        )
    )
    
    httr::GET(sprintf('http://localhost:60151/load?file=%s', sessionxml))
    
    #browseURL(sprintf('http://localhost:60151/load?file=%s', sessionxml))
    
    # HTMLreport <- maketracktable(
    #     fileSheet=FileSheet,
    #     SampleSheet=SampleSheet,
    #     filename="IGVExample_ce11.html",
    #     basedirectory=getwd(),
    #     genome="ce11"
    # )
    
    
}


# overlap_intersect('RC009', 'RC010')
# overlap_intersect('AA761', 'AA762')
# overlap_intersect('AA764', 'AA765')

