
testChromosomeNames <-  function(tss, gnm, ret=FALSE) {
    if( !all(seqlevels(tss) %in% seqlevels(gnm)) ) { 
        try( seqlevelsStyle(tss) <- seqlevelsStyle(gnm) )
        if( !all(seqlevels(tss) %in% seqlevels(gnm)) & ret ) {
            seqlevels(tss) <- as.character(as.roman( gsub('^chr', '', gsub('.*(M|m).*', 'M', seqlevels(tss)), ignore.case = TRUE) ))
            try( seqlevelsStyle(tss) <- seqlevelsStyle(gnm) )
        }
        if( !all(seqlevels(tss) %in% seqlevels(gnm)) ) 
            stop('Chromosome names provided in the file does not match ones defined in reference genome. \nINPUT: [', 
                 paste(seqlevels(tss)[!seqlevels(tss) %in% seqlevels(gnm)], collapse=', '), "]\nGENOME: [", paste(head(seqlevels(gnm), 5), collapse=', '), ', ...]', call. = FALSE) 
    }
    if(ret) return(tss)
}

#' Fixes overlaps in bed files
#' 
#' @param x path to bw
#' @param gnm BSgenome
#' @param pth output path
#' 
#' @author Przemyslaw Stempor
#' 
#' @family fix
#' @export
#' 
#' @examples
#' #sapply(dir(pattern='wig') , fiXWig)
#' 
fiXWig <- function(x, gnm = BSgenome.Celegans.UCSC.ce10::Celegans, pth = sub('wig', 'bw', x)) {
    message('fixing: ', x)
    fcon=file(x); wig <- rtracklayer::import.wig( fcon ); close(fcon);
    if( grepl('list', class(wig), ignore.case = TRUE) ) wig <- unlist(wig, use.names=FALSE)
    wig <- testChromosomeNames(wig , gnm, ret=TRUE)
    seqlengths(wig) <- seqlengths(gnm)[seqlevels(wig)];
    cov <- coverage(wig, weight='score')
    correction <- coverage(wig)
    correction[correction==0] <- 1
    cov <- cov/correction
    export.bw(cov, pth);
}    

#' Fixes names
#' 
#' @param ID ContactExpID of fixed entry
#' 
#' @author Przemyslaw Stempor
#' 
#' @family fix
#' @export
#' 
#' @examples
#' #sapply(sprintf('rML%.3i', 3:24) , fiXPath)
#' 
fiXPath <- function(ID) {
    message('fixing: ', ID)
    p <- dbGetQuery(con, sprintf("SELECT path, UID FROM labfiles WHERE ContactExpID = '%s'", ID))
    p$path <- file.path('files', p$path)
    rownames(p) <- p$UID
    sapply( rownames(p), function(x) {
        dbGetQuery(con, sprintf("UPDATE labfiles SET path = '%s' WHERE UID = '%s'", p[x,]$path, x))
    })
}