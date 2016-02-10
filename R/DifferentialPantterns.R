#' Get stage from IDs
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
getAnno <- function( ContactExpID, anno='Factor', EXTABLE='labrnaseq'){
    
    con <- dbConnect(dbDriver("MySQL"), group = "jadb", default.file='~/.my.cnf')
    PK <- dbGetQuery(con, sprintf("SHOW INDEX FROM %s WHERE Key_name = 'PRIMARY'", gsub('view$', '', EXTABLE) ))[['Column_name']]
    
    strain <- lapply(ContactExpID, function(x) {
        res <- dbGetQuery(con, sprintf('SELECT %s FROM %s WHERE %s = "%s"', anno, EXTABLE, PK, x))
        names(res) <- x
        return(res)
    })
    
    dbDisconnect(con)
    if( any(elementLengths(strain) != 1) ) stop('Multiple annotation for singe ID')
    return(unlist(strain))
}

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
#' @examples
#' 
#' ID=c("AM004","FB006","FB013merge","FB008","FB014merge","FB009","FB015merge","AA103p","AA076","AA186","AA378","BC137","BC131","AA077","AA078","BC138","BC132","AA089","AA079")
#' 
ids2SE <- function(ID, genes=unlist(range(get(data(gnmodel)))), mnt_point=NULL) {
    fls <- sapply(ID, getFilePath, processing='aligned', format='bam', scale='NA', url=is.null(mnt_point), eq=T)
    if(!is.null(mnt_point)) fls <- file.path(mnt_point, gsub('^files/', '', fls))
    
    if(!all(file.exists(fls))) stop('BAM missing')
    if(!all(file.exists(paste0(fls, '.bai')))) stop('index missing')
    
    library(Rsamtools); library(GenomicAlignments)
    message('Sumarizing the experiments')
    bfl <- BamFileList(fls)
    e <- summarizeOverlaps(genes, bfl)
    
    message('Getting annotations')
    SE <- e
    colData(SE)$factor <- factor(unlist(getAnno(ID, anno='factor', EXTABLE='labexperimentview'), use.names=FALSE))
    colData(SE)$stage <- factor(unlist(getAnno(ID, anno='stage', EXTABLE='labexperimentview'), use.names=TRUE))
    
    message('Running DEseq2')
    res <- doDiffExpr(e, design = ~ strain)
    
    return( res )
    
}


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
#' @examples
#' 
#' ID=c("AM004","FB006","FB013merge","FB008","FB014merge","FB009","FB015merge","AA103p","AA076","AA186","AA378","BC137","BC131","AA077","AA078","BC138","BC132","AA089","AA079")
#' 
wigs2SE <- function(ID, genes=unlist(range(get(data(gnmodel)))), proc='BEADS', mnt_point=NULL) {
    fls <- sapply(ID, getFilePath, processing='aligned', format='bam', scale='NA', url=is.null(mnt_point), eq=T)
    if(!is.null(mnt_point)) fls <- file.path(mnt_point, gsub('^files/', '', fls))
    
    if(!all(file.exists(fls))) stop('BAM missing')
    if(!all(file.exists(paste0(fls, '.bai')))) stop('index missing')
    
    library(Rsamtools); library(GenomicAlignments)
    message('Sumarizing the experiments')
    bfl <- BamFileList(fls)
    e <- summarizeOverlaps(genes, bfl)
    
    message('Getting annotations')
    colData(e)$strain <- factor( unlist(getStrain(colnames(e))) )
    colData(e)$strain <- factor( unlist(getStrain(colnames(e))) )
    
    message('Running DEseq2')
    res <- doDiffExpr(e, design = ~ strain)
    
    return( res )
    
}



