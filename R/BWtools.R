getBWarray <- function(IDs, res=1000L, processing = 'aligned') {
    
    con <- dbConnect(dbDriver("MySQL"), group = GROUP, default.file='~/.my.cnf')
    dbListTables(con)
    all <- dbReadTable(con, "labexperimentview")
    
    require(dplyr)
    all  %>% filter(Factor == 'HPL2', Strain == 'N2')  %>% select(ContactExpID)  %>% unlist -> IDs
    names(IDs) <- IDs
    
    IDs %>% sapply(getFilePath, format = 'bw', processing=processing) -> paths
    require(rtracklayer)
    bwfl <- BigWigFileList(unlist(paths))
    
    lst <- lapply(bwfl, function(bwf) extarct_vector(bwf, size = seqlengths(bwf)/res) )
    M <- do.call(cbind, lst)
    C <- cor(M)
    
    library(d3heatmap)
    d3heatmap(C, 
              colors = colorFunc, 
              symm=TRUE, 
              scale = 'none', 
              theme='dark', 
              revC = TRUE,
              zlim=c(0,1)
              )
    #heatmap.2(C,col = rev(RColorBrewer::brewer.pal(11, 'RdYlBu')), scale = 'none')
    plot.new()
    fields::image.plot(C, smallplot= c(0,1,.5,1), 
                       legend.only = TRUE, horizontal=TRUE, col = rev(RColorBrewer::brewer.pal(11, 'Spectral')))
        
}

extarct_vector <- function(track, size, which = as(seqinfo(track), "GenomicRanges")) {
    sum <- .Call(
        'BWGFile_summary', path.expand(path(track)),
        as.character(seqnames(which)), ranges(which), 
        S4Vectors::recycleIntegerArg(size, "size", length(which)), 
        "mean", as.numeric(NA_real_), PACKAGE='rtracklayer'
    )
    return(unlist(sum))
#     M <- do.call( rbind, sum )
#     if (!ignore_strand) 
#         M[as.character(strand(gr))=='-',] <- M[
#             as.character(strand(gr))=='-', ncol(M):1]
#     return(M)
}

#' Combine big wiggle files in the database into replicate track
#' 
#' @param IDs Vector of JADB ContactExpIDs
#' @param processing Processing level of tracks
#' @param res Resolution of tracks
#' @param outdir Directory to output BW file
#'   
#' @return List 
#' 
#' @author Przemyslaw Stempor
#' 
#' @family BW
#' @export
#' 
#' @examples
#' #combineReps(IDs)

combineReps <- function(r1=NULL, r2=NULL, IDs=NULL, processing='aligned', res=100L, scale='linear', outdir=tempdir()) {
    
    if( !is.null(IDs) ) {
        message('Getting files')
        paths <- sapply(IDs, getFilePath, format = 'bw', processing=processing, scale=scale, url = FALSE) 
        paths <- gsub('files', MOUNT, paths)
        if(length(IDs) != length(paths)) stop('Many/missing tracks per ID, cannot combine.')
    } else {
        message('Using replicate with r1=', r1, ' and r2=', r2)
        paths <- c(r1, r2)
    }
    
    
    anno <- as.data.frame(t(sapply(basename(paths), rbeads:::ParseName)))
    names(paths) <- anno$ContactExpID
    
    validRep <- sapply(c('Factor', 'Strain', 'Stage', 'Processing', 'Scale', 'Resolution'), function(x) length(unique(anno[x]))==1)
    if(!all(validRep)) stop('Not a valid replicate')
    
    message('Calculating correlation')
    require(rtracklayer)   
    bwfl <- BigWigFileList(unlist(paths))
    lst <- lapply(bwfl, function(bwf) {
        extarct_vector(bwf, size = seqlengths(bwf)/res)
    })
    M <- do.call(cbind, lst)
    C <- cor(M, use = 'na.or.complete')
    C <- mean(C[upper.tri(C)])
    
    for(n in names(bwfl)) {
        message('Summarizing ', n)
        if(exists('bw')) {
            bw <- bw + import.bw(bwfl[[n]], as='Rle')
        } else {
            bw <- import.bw(bwfl[[n]], as='Rle')
        }
    }
    bw <- bw / length(bwfl)
    
    same <- anno[1,c('Factor', 'Strain', 'Stage', 'Processing', 'Scale', 'Resolution')]
    outname <- paste0(paste0(unlist(same), collapse='_'), '_', paste(anno$ContactExpID, collapse = '^'), '_', paste(anno$UID, collapse = '^'), '.bw')
    outname <- file.path(outdir, outname)
    message('Exporting: ', outname)
    export.bw(bw, outname)
    
    return( list(out=outname, cor=C, anno=anno) )
    
}

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
addRepToJADB <- function(IDs, res=100L) {
    
    require(magrittr)
    require(parallel)
    
    message('Sumarizing BW tracks')
    allbw <- lapply(IDs, getFilePath, format='bw', url=FALSE, mount=TRUE) %>% lapply(sort)
    if(length(allbw[[1]]) != length(allbw[[2]])) stop('err_diff_proc_on_bw')
    outputbw <- mcMap(combineReps, allbw[[1]], allbw[[2]], mc.cores=11)
    
    
    ## get metatdata
    con <- dbConnect(dbDriver("MySQL"), group = GROUP, default.file='~/.my.cnf')
    REP <- dbReadTable(con, "labchipseqrep")
    CXID <- sprintf('REP%03.0f', max(as.numeric(gsub('REP', '', REP$ContactExpID)))+1)
    
    setwd(MOUNT);
    outdir <- file.path('REPLICATES', JADBtools:::getAnno(IDs[[1]], EXTABLE = 'labexperiment'), CXID)
    dir.create(outdir, recursive = TRUE)
    
    BAM <- sapply(IDs, getFilePath, format='bam', url=FALSE)
    anno <- as.data.frame(t(sapply(basename(BAM), rbeads:::ParseName)))
    
    
    #out <- combineReps(IDs, processing = 'alnNQNU', outdir = outdir, res = res)
    #outNorm <- combineReps(IDs, processing = 'BEADSNQNU', outdir = outdir, res = res)
    
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
    addGenericFile(CXID, path = file.path('files', peaksU), Processing = 'PeakUnion',     Resolution = 'q01', Scale = 'MACS', filetype_format = 'bed', prefix = 'R', repPath = TRUE)
    addGenericFile(CXID, path = file.path('files', peaksI), Processing = 'PeakIntersect', Resolution = 'q01', Scale = 'MACS', filetype_format = 'bed', prefix = 'R', repPath = TRUE)
    addGenericFile(CXID, path = file.path('files', peaksIDR), Processing = 'IDR2', Resolution = 'q01', Scale = 'MACS', filetype_format = 'narrowPeak', prefix = 'R', repPath = TRUE)
    
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
            repPath = TRUE
        )
    })
   
    oldwd <- getwd(); setwd(outdir)
    #dir.create('IDR', recursive = TRUE)
    

    
    cmd <- sprintf(
        'export PATH=/mnt/home1/ahringer/ps562/software/bin:$PATH; /mnt/home1/ahringer/jarun/miniconda2/bin/ipython ~/TEST/macs2_idr.ipy -- %s %s -c %s -p ./IDR/idr',
        gsub('^files', MOUNT, BAM[1]), 
        gsub('^files', MOUNT, BAM[2]), 
        if( grepl('^E', anno$Crosslinker[[1]]) ) {
            file.path(MOUNT, 'Input/SummedInputs/ce11/ce11_EGS_HiSeq_input.bam')
        } else {
            file.path(MOUNT, 'Input/SummedInputs/ce11/ce11_FRM_HiSeq_input.bam')
        }
    )
    message('--> ', cmd)
    #system(cmd)
    
    
    #same <- anno[1,c('Factor', 'Antibody', 'Strain', 'Stage', 'Processing', 'Scale', 'Resolution')]
    #same$Processing <- 'PeakCalls'; same$Scale <- 'MACS'; same$Resolution <- 'q01'
    #outname <- paste0(paste0(unlist(same), collapse='_'), '_', paste(anno$ContactExpID, collapse = '^'), '_', 'IDR', '.bed')
    
    
    #ile.copy('IDR/idr_final_peaks_0.05.narrowPeak', outname)
    #addGenericFile(CXID, path = file.path('files', outdir, outname), Processing = 'IDR', Resolution = 'q01', Scale = 'MACS', filetype_format = 'bed', prefix = 'R', repPath = TRUE)
    setwd(oldwd)
    
    #UPDATE `mydb`.`labfiles` SET `filetype_format`='bwz' WHERE `UID`='R3e31188';
    
}
    

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
addIDRtoJADB <- function(IDs, res=100L) {
    
    
    ## get metatdata
    con <- dbConnect(dbDriver("MySQL"), group = GROUP, default.file='~/.my.cnf')
    T <- dbReadTable(con, "labchipseqrep")
    CXID <- sprintf('REP%03.0f', max(as.numeric(gsub('REP', '', T$ContactExpID)))+1)
    
    outdir <- file.path('REPLICATES', JADBtools:::getAnno(IDs[[1]], EXTABLE = 'labexperiment'), CXID)
    dir.create(outdir, recursive = TRUE)
    
    BAM <- sapply(IDs, getFilePath, format='bam', url=FALSE)
    anno <- as.data.frame(t(sapply(basename(BAM), rbeads:::ParseName)))
    
    out <- combineReps(IDs, processing = 'aligned', outdir = outdir, res = res)
    outNorm <- combineReps(IDs, processing = 'NORM', outdir = outdir, res = res)
    
    ## Add experiment to db
    INSERT <- out$anno[1,-c(3,7,8,9,10,11,12)]
    INSERT[['ContactExpID']]<- CXID
    INSERT[['Created']] <- paste(Sys.Date())
    INSERT[['Updated']] <- paste(Sys.Date())
    INSERT[['Comments']] <- paste0('corA=', round(out$cor, 3), '; corN=', round(outNorm$cor, 3))
    INSERT[['ExtractID']] <- paste0(unlist(out$anno$ExtractID), collapse='|')
    INSERT[['Experiments']] <- paste0(unlist(out$anno$ContactExpID), collapse='|')
    
    TABLE <- 'labchipseqrep'
    fileds.def <- dbGetQuery(con, sprintf("SHOW FIELDS FROM %s", TABLE))
    PK <- dbGetQuery(con, sprintf("SHOW INDEX FROM %s WHERE Key_name = 'PRIMARY'", gsub('view$', '', TABLE) ))[['Column_name']]
    
    sql <- paste("INSERT INTO ", TABLE,"(", paste(colnames(INSERT), collapse=", "),") VALUES('", paste(as.character(unlist(INSERT)), collapse="', '"), "')", collapse=", ", sep="")
    
    rs <- dbSendQuery(con, sql )
    info <- dbGetInfo(rs)
    
    dbDisconnect(con)
    
    outMapq0 <- combineReps(IDs, processing = 'mapq0', outdir = outdir, res = res)
    outNormLog2 <- combineReps(IDs, processing = 'NORM', outdir = outdir, scale = 'log2$', res = res)
    outNormLog2zsc <- combineReps(IDs, processing = 'NORM', outdir = outdir, scale = 'log2zsc', res = res)
    outNormZscore <- combineReps(IDs, processing = 'NORM', outdir = outdir, scale = 'zscore', res = res)
    
    oldwd <- getwd(); setwd(outdir)
    peaksU <- file.path(outdir, combinePeaksToBed(IDs, mode = 'union'))
    peaksI <- file.path(outdir, combinePeaksToBed(IDs, mode = 'intersection'))
    enreg <- file.path(outdir, enrichedRegionsCall(
        basename(outNorm$out),
        basename(gsub('PeakCalls_MACS_q01(.+)_union', 'EnrichedRegions_a75_b9\\1', peaksU))
    ))
    
    setwd(oldwd)
    
    
   
    
    
    addGenericFile(CXID, path = file.path('files', out$out), Processing = 'aligned',  Resolution = '1bp', Scale = 'linear', filetype_format = 'bw', prefix = 'R', repPath = TRUE)
    addGenericFile(CXID, path = file.path('files', outNorm$out), Processing = 'NORM', Resolution = '1bp', Scale = 'linear', filetype_format = 'bw', prefix = 'R', repPath = TRUE)
    addGenericFile(CXID, path = file.path('files', outMapq0$out), Processing = 'outMapq0', Resolution = '1bp', Scale = 'linear', filetype_format = 'bw', prefix = 'R', repPath = TRUE)
    
    
    addGenericFile(CXID, path = file.path('files', outNormLog2$out), Processing = 'NORM', Resolution = '1bp', Scale = 'log2', filetype_format = 'bw', prefix = 'R', repPath = TRUE)
    addGenericFile(CXID, path = file.path('files', outNormLog2zsc$out), Processing = 'NORM', Resolution = '1bp', Scale = 'log2zsc', filetype_format = 'bw', prefix = 'R', repPath = TRUE)
    addGenericFile(CXID, path = file.path('files', outNormZscore$out), Processing = 'NORM', Resolution = '1bp', Scale = 'zscore', filetype_format = 'bw', prefix = 'R', repPath = TRUE)
    
    addGenericFile(CXID, path = file.path('files', peaksU), Processing = 'PeakUnion',     Resolution = 'q01', Scale = 'MACS', filetype_format = 'bed', prefix = 'R', repPath = TRUE)
    addGenericFile(CXID, path = file.path('files', peaksI), Processing = 'PeakIntersect', Resolution = 'q01', Scale = 'MACS', filetype_format = 'bed', prefix = 'R', repPath = TRUE)
    addGenericFile(CXID, path = file.path('files', enreg), Processing = 'EnrichedRegions', Resolution = 'a75', Scale = 'q9', filetype_format = 'bed', prefix = 'R', repPath = TRUE)
    
    oldwd <- getwd(); setwd(outdir)
    dir.create('IDR', recursive = TRUE)
    
    
    
    cmd <- sprintf(
        'export PATH=/home/ps562/software/bin:$PATH; /home/ps562/anaconda2/bin/ipython ~/TEST/macs2_idr.ipy -- %s %s -c %s -p ./IDR/idr',
        gsub('^files', MOUNT, BAM[1]), 
        gsub('^files', MOUNT, BAM[2]), 
        if( grepl('^E', anno$Crosslinker[[1]]) ) {
            file.path(MOUNT, 'Input/SummedInputs/EGS_HiSeq_input.bam')
        } else {
            file.path(MOUNT, 'Input/SummedInputs/FRM_HiSeq_input.bam')
        }
    )
    message('--> ', cmd)
    system(cmd)
    
    
    same <- anno[1,c('Factor', 'Antibody', 'Strain', 'Stage', 'Processing', 'Scale', 'Resolution')]
    same$Processing <- 'PeakCalls'; same$Scale <- 'MACS'; same$Resolution <- 'q01'
    outname <- paste0(paste0(unlist(same), collapse='_'), '_', paste(anno$ContactExpID, collapse = '^'), '_', 'IDR', '.bed')
    
    
    file.copy('IDR/idr_final_peaks_0.05.narrowPeak', outname)
    addGenericFile(CXID, path = file.path('files', outdir, outname), Processing = 'IDR', Resolution = 'q01', Scale = 'MACS', filetype_format = 'bed', prefix = 'R', repPath = TRUE)
    setwd(oldwd)
    
    #UPDATE `mydb`.`labfiles` SET `filetype_format`='bwz' WHERE `UID`='R3e31188';
    
}