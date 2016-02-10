#' Get ratio between worm and 10x fly reads
#' 
#' This ratio can be used as global normalizing parameter for sipke in.
#' 
#' @param ID the database id of the experiment in "XX001" format
#'   
#' @return named \code{\link[base]{numeric}} value
#' 
#' @author Przemyslaw Stempor
#' 
#' @family spikeIn
#' @export
#' 
getRatio <- function(ID) {
    mysql <- dbDriver("MySQL")
    con <- dbConnect(mysql, group = "jadb", default.file='~/.my.cnf')
    
    exp_file 	<- unlist( dbGetQuery(con, paste("SELECT path FROM labfiles WHERE ContactExpID = '", ID, "' AND Filetype_format='bam'; ", collapse="", sep="") ) )
    addr <- file.path("http://jadb.gurdon.private.cam.ac.uk/db4", dirname(exp_file), "mga/results.xml")
    
    data <- xmlParse(addr)
    xml_data <- xmlToList(data)
    lst <- xml_data$MultiGenomeAlignmentSummary$AlignmentSummaries
    
    ce10 <- as.numeric( lst[sapply(lst, function(x) x[[1]][['id']]) == 'ce10'][[1]]$AssignedCount )
    dm3 <- as.numeric( lst[sapply(lst, function(x) x[[1]][['id']]) == 'dm3'][[1]]$AssignedCount )
    
    dbDisconnect(con)
    
    return( ce10/(10*dm3) )
    
    
}

#' Normalize track by spike in ratio value ratios produced by 
#' 
#' The function produces normalized track using \code{\link{getRatio}} value.
#' 
#' @param the named numeric representing normakization coefficient named with 
#' experiment ID in "XX001" format
#'   
#' @return normalized track in BW format (connection)
#' 
#' @author Przemyslaw Stempor
#' 
#' @family spikeIn
#' @export
#' 
normTrack <- function(nc) {
    mysql <- dbDriver("MySQL")
    con <- dbConnect(mysql, group = "jadb", default.file='~/.my.cnf')
    
    
    ID <- names(nc)
    message(ID)
    exp_file     <- unlist( dbGetQuery(
        con, paste(
            "SELECT path FROM labfiles WHERE ContactExpID = '", 
            ID, "' AND Filetype_format='bw' AND  Processing='NORM'",
            "AND Scale='linear'", collapse="", sep=""
        ) 
    ) )
    addr <- file.path("http://jadb.gurdon.private.cam.ac.uk/db4",  exp_file )
    message("Reading")
    bw <- import.bw(addr,as="RleList")
    out <- bw * nc
    gr <- as(out, 'GRanges')
    gr <- gr[!is.na(gr$score)]
    gr <- gr[!is.infinite(gr$score)]
    message("Writeing")
    export.bw(gr, rbeads:::reName(basename(addr), proccesing = 'spikeinBEADS', scale = 'zsc'))
}
