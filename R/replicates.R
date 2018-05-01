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
    require(rtracklayer)
    library(tracktables)
    
    message('Sumarizing BW tracks')
    allbw <- lapply(IDs, getFilePath, format='bw', url=FALSE, mount=TRUE) %>% lapply(sort)
    if(length(allbw[[1]]) != length(allbw[[2]])) stop('err_diff_proc_on_bw')
    outputbw <- mcMap(combineReps, allbw[[1]], allbw[[2]], outdir=file.path(MOUNT, 'temp'), mc.cores=11)
    
    
    ## get metatdata
    con <- dbConnect(dbDriver(DRIVER), group = GROUP, default.file='~/.my.cnf')
    REP <- dbReadTable(con, "labchipseqrep")
    
    if(any(grepl(prefix, REP$ContactExpID))) {
        nid <- max(as.numeric(gsub(prefix, '', REP$ContactExpID[grepl(prefix, REP$ContactExpID)])))+1
    } else {
        nid <- 1
    }
    CXID <- sprintf('%s%03.0f', prefix, nid)
    
    setwd(MOUNT);
    outdir <- file.path('REPLICATES', JADBtools:::getAnno(IDs[[1]], EXTABLE = 'labexperiment'), CXID)
    dir.create(outdir, recursive = TRUE)
    
    BAM <- sapply(IDs, getFilePath, format='bam', url=FALSE)
    anno <- as.data.frame(t(sapply(basename(BAM), rbeads:::ParseName)))
    comment <- paste0('corA=', round(outputbw[[2]]$cor, 3), '; corN=', round(outputbw[[6]]$cor, 3))
    if(!is.null(extra_stats)) comment <- paste0(extra_stats, '; ', comment)
    
    ## Add experiment to db
    INSERT <- anno[1,-c(3,7,8,9,10,11,12)]
    INSERT[['ContactExpID']]<- CXID
    INSERT[['Created']] <- paste(Sys.Date())
    INSERT[['Updated']] <- paste(Sys.Date())
    INSERT[['Comments']] <- comment
    INSERT[['ExtractID']] <- paste0(unlist(outputbw[[1]]$anno$ExtractID), collapse='|')
    INSERT[['Experiments']] <- paste0(unlist(outputbw[[1]]$anno$ContactExpID), collapse='|')
    
    TABLE <- 'labchipseqrep'
    fileds.def <- dbGetQuery(con, sprintf("SHOW FIELDS FROM %s", TABLE))
    PK <- dbGetQuery(con, sprintf("SHOW INDEX FROM %s WHERE Key_name = 'PRIMARY'", gsub('view$', '', TABLE) ))[['Column_name']]
    
    sql <- paste("INSERT INTO ", TABLE,"(", paste(colnames(INSERT), collapse=", "),") VALUES('", paste(as.character(unlist(INSERT)), collapse="', '"), "')", collapse=", ", sep="")
    
    rs <- dbSendQuery(con, sql )
    info <- dbGetInfo(rs)
    
    dbDisconnect(con)
    
    message('Sumarizing peak calls')
    
    oldwd <- getwd(); setwd(outdir)
    #peaksU <- file.path(outdir, combinePeaksToBed(IDs, mode = 'union'))
    #peaksI <- file.path(outdir, combinePeaksToBed(IDs, mode = 'intersection'))
    
    peaksU <- file.path(outdir, overlap_intersect_nd(IDs[[1]], IDs[[2]], mode = 'union', filter_map = FALSE))
    peaksI <- file.path(outdir, overlap_intersect_nd(IDs[[1]], IDs[[2]], mode = 'intersection', filter_map = FALSE))
    
    message('Running IDR2')
    #peaksIDR <- file.path(outdir, combinePeaksIDR(IDs))
    #enreg <- file.path(outdir, enrichedRegionsCall(
    #    basename(outNorm$out),
    #    basename(gsub('PeakCalls_MACS_q01(.+)_union', 'EnrichedRegions_a75_b9\\1', peaksU))
    #))
    comment1 <- paste0('n_peaks=', length( readLines(file.path(MOUNT, peaksU)) ) ) 
    comment2 <- paste0('n_peaks=', length( readLines(file.path(MOUNT,peaksI)) ) ) 
    #comment3 <- paste0('n_peaks=', length( readLines(file.path(MOUNT,peaksIDR)) ) ) 
    
    addGenericFile(CXID, path = file.path('files', peaksU), comments = comment1, Processing = 'PeakUnion',     Resolution = 'q01', Scale = 'MACS', filetype_format = 'bed', prefix = 'R', repPath = TRUE)
    addGenericFile(CXID, path = file.path('files', peaksI), comments = comment2, Processing = 'PeakIntersect', Resolution = 'q01', Scale = 'MACS', filetype_format = 'bed', prefix = 'R', repPath = TRUE)
    #addGenericFile(CXID, path = file.path('files', peaksIDR), comments = comment3, Processing = 'IDR2', Resolution = 'q01', Scale = 'MACS', filetype_format = 'narrowPeak', prefix = 'R', repPath = TRUE)
    
    addGenericFile(CXID, path = file.path('files', paste0(peaksU, '.xml')), comments = comment1, Processing = 'PeakUnion',     Resolution = 'q01', Scale = 'MACS', filetype_format = 'xml', prefix = 'R', repPath = TRUE)
    addGenericFile(CXID, path = file.path('files', paste0(peaksI, '.xml')), comments = comment2, Processing = 'PeakIntersect', Resolution = 'q01', Scale = 'MACS', filetype_format = 'xml', prefix = 'R', repPath = TRUE)
    
    #addGenericFile(CXID, path = file.path('files', enreg), Processing = 'EnrichedRegions', Resolution = 'a75', Scale = 'q9', filetype_format = 'bed', prefix = 'R', repPath = TRUE)
    setwd(oldwd)
    
    
    lapply(outputbw, function(x) {
        message(basename(x$out))
        fp <- file.path(MOUNT, outdir, basename(x$out))
        file.rename(x$out, fp)
        
        addGenericFile(
            CXID, path = gsub(MOUNT, 'files', fp), 
            Processing = unlist(unique(x$anno$Processing)),  
            Resolution = unlist(unique(x$anno$Resolution)), 
            Scale = unlist(unique(x$anno$Scale)), 
            filetype_format = 'bw', 
            prefix = 'R', 
            repPath = TRUE,
            parent1_uid = x$anno$UID[[1]],
            parent2_uid = x$anno$UID[[2]]
        )
    })
    
    oldwd <- getwd(); setwd(outdir)
    
}