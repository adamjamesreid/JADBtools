#' Run full differential expression pipeline
#' 
#' @param ID the vector of ContactExpIDs
#' @param mnt_point give the mount point for filesystem, use URLs if NULL, defaults to NULL
#'   
#' @return stats string
#' 
#' @author Przemyslaw Stempor
#' 
#' @family RNAseq
#' @export
#' 
ids2diffExpr <- function(ID, mnt_point=NULL) {
    fls <- sapply(ID, getFilePath, processing='aligned', format='bam', scale='NA', url=is.null(mnt_point), eq=F)
    if(!is.null(mnt_point)) fls <- file.path(mnt_point, fls)
    
    message('Sumarizing the experiments')
    e <- summarizeBAMs(fls)
    
    message('Getting annotations')
    colData(e)$strain <- factor( unlist(getStrain(colnames(e))) )
    
    message('Running DEseq2')
    res <- doDiffExpr(e, design = ~ strain)
    
    return( res$outname )
    
}

#' Get counts from bam files using embended gene meodel
#' 
#' File URLs can be obtained like this:
#' fls <- sapply(c('rAM001', 'rAM002', 'rAM003', 'rAM004'), getFilePath, processing='aligned', format='bam')
#' 
#' @param fls file names vector
#'   
#' @return eset
#' 
#' @author Przemyslaw Stempor
#' 
#' @family RNAseq
#' @export
#' 
summarizeBAMs <-function(fls) {
    bfl <- BamFileList(fls)
    
    data(gnmodel)
    expset_repeats <- summarizeOverlaps(gnmodel, bfl)
    return(expset_repeats)
    
}

#' Get starnd from IDs
#' 
#' @param ContactExpID vector of IDs
#'   
#' @return eset
#' 
#' @author Przemyslaw Stempor
#' 
#' @family RNAseq
#' @export
#' 
getStrain <- function( ContactExpID, EXTABLE='labrnaseq'){
    
    con <- dbConnect(dbDriver("MySQL"), group = "jadb")
    PK <- dbGetQuery(con, sprintf("SHOW INDEX FROM %s WHERE Key_name = 'PRIMARY'", gsub('view$', '', EXTABLE) ))[['Column_name']]
    
    strain <- sapply(ContactExpID, function(x) dbGetQuery(
        con, sprintf('SELECT %s FROM %s WHERE %s = "%s"', 'Strain', 'labrnaseq', PK, x)
    ))
    
    dbDisconnect(con)
    return(strain)
}

#' Run DEseq2 and output results
#' 
#' @param e eset
#' @param design experiment design
#'   
#' @return list
#' 
#' @author Przemyslaw Stempor
#' 
#' @family RNAseq
#' @export
#' 
doDiffExpr <- function(e, design = ~ strain) {
    
    dds <- DESeqDataSet(e, design = design )
    dds <- DESeq(dds)
    
    res <- results(dds)
    out <- as.data.frame(res)
    colnames(out) <- elementMetadata(res)$description
    
    M <- assays(e)$counts
    RPKM <- as.data.frame(
        M / ((sum(width(gnmodel)) / 1e3) %*% t(colSums(M) / 1e6))
    )
    colnames(RPKM) <- paste0('RPKM_', colnames(RPKM))
    
    out <- cbind(
        elementMetadata(gnmodel), 
        data.frame(exon_width_kb=(sum(width(gnmodel)) / 1e3)),
        RPKM, 
        out
    )         

    outname <- elementMetadata(res)$description[[2]]
    outname <- gsub(' ', '_', gsub('log2 fold change \\(MAP\\): ', '', outname))
    outname <- paste(paste(colnames(e), collapse='-'), outname, 'RPKM_and_DESeq2.csv', sep='_')
    write.csv(out, outname)
    
    return(list(dds=dds, out=out, outname=outname))
    
}
