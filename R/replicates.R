#' Adds replicates to database
#' 
#' @param IDs Vector of JADB ContactExpIDs
#' @param res Resolution of tracks
#' @param outdir Directory to output BW file
#'   
#' @return NULL 
#' 
#' @author Przemyslaw Stempor
#' 
#' @family BW
#' @export
#' 
#' @examples
#' #combineReps(IDs)
jadb_replicates_chip <- function(IDs, res=100L, prefix='REP', extra_stats=NULL) {
    
    require(magrittr)
    require(parallel)
    
    message('Sumarizing BW tracks')
    allbw <- lapply(IDs, getFilePath, format='bw', url=FALSE, mount=TRUE) %>% lapply(sort)
    if(length(allbw[[1]]) != length(allbw[[2]])) stop('err_diff_proc_on_bw')
    outputbw <- mcMap(combineReps, allbw[[1]], allbw[[2]], mc.cores=11)
    
    
    ## get metatdata
    con <- dbConnect(dbDriver("MySQL"), group = GROUP, default.file='~/.my.cnf')
    REP <- dbReadTable(con, "labchipseqrep")
    CXID <- sprintf('REP%03.0f', max(as.numeric(gsub(prefix, '', REP$ContactExpID)))+1)
    
    setwd(MOUNT);
    outdir <- file.path('REPLICATES', JADBtools:::getAnno(IDs[[1]], EXTABLE = 'labexperiment'), CXID)
    dir.create(outdir, recursive = TRUE)
    
    BAM <- sapply(IDs, getFilePath, format='bam', url=FALSE)
    anno <- as.data.frame(t(sapply(basename(BAM), rbeads:::ParseName)))
    
    ## Add experiment to db
    INSERT <- anno[1,-c(3,7,8,9,10,11,12)]
    INSERT[['ContactExpID']]<- CXID
    INSERT[['Created']] <- paste(Sys.Date())
    INSERT[['Updated']] <- paste(Sys.Date())
    INSERT[['Comments']] <- paste0('corA=', round(outputbw[[2]]$cor, 3), '; corN=', round(outputbw[[6]]$cor, 3))
    INSERT[['ExtractID']] <- paste0(unlist(outputbw[[1]]$anno$ExtractID), collapse='|')
    INSERT[['Experiments']] <- paste0(unlist(outputbw[[1]]$anno$ContactExpID), collapse='|')
    
    TABLE <- 'labchipseqrep'
    fileds.def <- dbGetQuery(con, sprintf("SHOW FIELDS FROM %s", TABLE))
    PK <- dbGetQuery(con, sprintf("SHOW INDEX FROM %s WHERE Key_name = 'PRIMARY'", gsub('view$', '', TABLE) ))[['Column_name']]
    
    sql <- paste("INSERT INTO ", TABLE,"(", paste(colnames(INSERT), collapse=", "),") VALUES('", paste(as.character(unlist(INSERT)), collapse="', '"), "')", collapse=", ", sep="")
    
    rs <- dbSendQuery(con, sql )
    info <- dbGetInfo(rs)
    
    dbDisconnect(con)
    
    #outMapq0 <- combineReps(IDs, processing = 'mapq0', outdir = outdir, res = res)
    #outNormLog2 <- combineReps(IDs, processing = 'NORM', outdir = outdir, scale = 'log2$', res = res)
    #outNormLog2zsc <- combineReps(IDs, processing = 'NORM', outdir = outdir, scale = 'log2zsc', res = res)
    #outNormZscore <- combineReps(IDs, processing = 'NORM', outdir = outdir, scale = 'zscore', res = res)
    
    message('Sumarizing BW peak calls')
    
    oldwd <- getwd(); setwd(outdir)
    peaksU <- file.path(outdir, combinePeaksToBed(IDs, mode = 'union'))
    peaksI <- file.path(outdir, combinePeaksToBed(IDs, mode = 'intersection'))
    peaksIDR <- file.path(outdir, combinePeaksIDR(IDs))
    #enreg <- file.path(outdir, enrichedRegionsCall(
    #    basename(outNorm$out),
    #    basename(gsub('PeakCalls_MACS_q01(.+)_union', 'EnrichedRegions_a75_b9\\1', peaksU))
    #))
    comment1 <- paste0('n_peaks=', length( readLines(peaksU) ) ) 
    comment2 <- paste0('n_peaks=', length( readLines(peaksI) ) ) 
    comment3 <- paste0('n_peaks=', length( readLines(peaksIDR) ) ) 
    
    addGenericFile(CXID, path = file.path('files', peaksU), comments = comment1, Processing = 'PeakUnion',     Resolution = 'q01', Scale = 'MACS', filetype_format = 'bed', prefix = 'R', repPath = TRUE)
    addGenericFile(CXID, path = file.path('files', peaksI), comments = comment2, Processing = 'PeakIntersect', Resolution = 'q01', Scale = 'MACS', filetype_format = 'bed', prefix = 'R', repPath = TRUE)
    addGenericFile(CXID, path = file.path('files', peaksIDR), comments = comment3, Processing = 'IDR2', Resolution = 'q01', Scale = 'MACS', filetype_format = 'narrowPeak', prefix = 'R', repPath = TRUE)
    
    #addGenericFile(CXID, path = file.path('files', enreg), Processing = 'EnrichedRegions', Resolution = 'a75', Scale = 'q9', filetype_format = 'bed', prefix = 'R', repPath = TRUE)
    setwd(oldwd)
    
    
    lapply(outputbw, function(x) {
        message(basename(x$out))
        fp <- file.path(MOUNT, outdir, basename(x$out))
        file.copy(x$out, fp)
        
        addGenericFile(
            CXID, path = gsub(MOUNT, 'files', fp), 
            Processing = unlist(unique(x$anno$Processing)),  
            Resolution = unlist(unique(x$anno$Resolution)), 
            Scale = unlist(unique(x$anno$Scale)), 
            filetype_format = 'bw', 
            prefix = 'R', 
            repPath = TRUE,
            parent1_uid = x$anno$UID[[1]],
            parent1_uid = x$anno$UID[[2]]
        )
    })
    
    oldwd <- getwd(); setwd(outdir)
    
}