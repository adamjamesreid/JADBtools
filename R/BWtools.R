getBWarray <- function(IDs, res=1000L, processing = 'aligned') {
    
    con <- dbConnect(dbDriver("MySQL"), group = "jadb", default.file='~/.my.cnf')
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

combineReps <- function(IDs, processing='aligned', res=100L, scale='linear', outdir=tempdir()) {
    
    message('Getting files')
    paths <- sapply(IDs, getFilePath, format = 'bw', processing=processing, scale=scale) 
    if(length(IDs) != length(paths)) stop('Many/missing tracks per ID, cannot combine.')
    anno <- as.data.frame(t(sapply(basename(paths), rbeads:::ParseName)))
    
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
    
    con <- dbConnect(dbDriver("MySQL"), group = "jadb", default.file='~/.my.cnf')
    T <- dbReadTable(con, "labchipseqrep")
    CXID <- sprintf('REP%03.0f', max(as.numeric(gsub('REP', '', T$ContactExpID)))+1)
    
    outdir <- file.path('REPLICATES', JADBtools:::getAnno(IDs[[1]], EXTABLE = 'labexperiment'), CXID)
    dir.create(outdir, recursive = TRUE)
    
    out <- combineReps(IDs, processing = 'aligned', outdir = outdir, res = res)
    outNorm <- combineReps(IDs, processing = 'NORM', outdir = outdir, res = res)
    outNormLog2 <- combineReps(IDs, processing = 'NORM', outdir = outdir, scale = 'log2$', res = res)
    outNormLog2zsc <- combineReps(IDs, processing = 'NORM', outdir = outdir, scale = 'log2zsc', res = res)
    outNormZscore <- combineReps(IDs, processing = 'NORM', outdir = outdir, scale = 'zscore', res = res)
    
    peaksU <- combinePeaksToBed(IDs, mode = 'union')
    peaksI <- combinePeaksToBed(IDs, mode = 'intersection')
    
    con <- dbConnect(dbDriver("MySQL"), group = "jadb", default.file='~/.my.cnf')
    T <- dbReadTable(con, "labchipseqrep")
    
    INSERT <- out$anno[1,-c(2,3,7,8,9,10,11,12)]
    INSERT[['ContactExpID']]<- CXID
    INSERT[['dateCreated']] <- paste(Sys.Date())
    INSERT[['dateUpdated']] <- paste(Sys.Date())
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
    
    addGenericFile(CXID, path = file.path('files', out$out), Processing = 'aligned',  Resolution = '1bp', Scale = 'linear', filetype_format = 'bw', prefix = 'R', repPath = TRUE)
    addGenericFile(CXID, path = file.path('files', outNorm$out), Processing = 'NORM', Resolution = '1bp', Scale = 'linear', filetype_format = 'bw', prefix = 'R', repPath = TRUE)
    
    addGenericFile(CXID, path = file.path('files', outNormLog2$out), Processing = 'NORM', Resolution = '1bp', Scale = 'log2', filetype_format = 'bw', prefix = 'R', repPath = TRUE)
    addGenericFile(CXID, path = file.path('files', outNormLog2zsc$out), Processing = 'NORM', Resolution = '1bp', Scale = 'log2zsc', filetype_format = 'bw', prefix = 'R', repPath = TRUE)
    addGenericFile(CXID, path = file.path('files', outNormZscore$out), Processing = 'NORM', Resolution = '1bp', Scale = 'zscore', filetype_format = 'bw', prefix = 'R', repPath = TRUE)
    
    addGenericFile(CXID, path = file.path('files', peaksU), Processing = 'PeakUnion',     Resolution = 'q01', Scale = 'MACS', filetype_format = 'bed', prefix = 'R', repPath = TRUE)
    addGenericFile(CXID, path = file.path('files', peaksI), Processing = 'PeakIntersect', Resolution = 'q01', Scale = 'MACS', filetype_format = 'bed', prefix = 'R', repPath = TRUE)
    
    
    
    #UPDATE `mydb`.`labfiles` SET `filetype_format`='bwz' WHERE `UID`='R3e31188';
    
}
    