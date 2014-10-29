#' Get file URL
#' 
#' Accepts regular expressions.
#' 
#' @param ID ContactExperimetID 
#' @param scale filter regex
#' @param processing filter regex 
#' @param format filter regex
#'   
#' @return URL string
#' 
#' @author Przemyslaw Stempor
#' 
#' @family dbtools
#' @export
#' 
getFilePath <- function(ID, format='.', processing='.', scale='.') {

    con <- dbConnect(dbDriver("MySQL"), group = "jadb")
  
    exp_file <- unlist( dbGetQuery(
        con, paste(
            "SELECT path FROM labfiles WHERE ContactExpID REGEXP '", ID, 
            "' AND Filetype_format REGEXP '",format,"' AND  Processing REGEXP '", processing,
            "'", "AND Scale REGEXP '", scale, "'", collapse="", sep=""
        ) 
    ))
    addr <- file.path("http://jadb.gurdon.private.cam.ac.uk/db4",  exp_file )
    dbDisconnect(con)
    return(addr)
}

#' Make UID
#' 
#' @param string file patch or other string for md5 
#' @param prefix UID prefix letter, defaults to "X"
#'   
#' @return IUD string
#' 
#' @author Przemyslaw Stempor
#' 
#' @family dbtools
#' @export
#' 
makeFileUID <- function(string='', prefix='X') {
    con <- dbConnect(dbDriver("MySQL"), group = "jadb")
    md5 <- substr(digest::digest(string, algo='md5'), 1, 2)
    num <- max(as.numeric(substr(unlist(dbGetQuery(con, "SELECT UID FROM labfiles")), 4, 8)))+1
    dbDisconnect(con)
    return(sprintf('%s%s%05d', prefix, md5, num))
}


#' Add file to database
#' 
#' @param ContactExpID ID string
#' @param path prefix UID prefix letter, defaults to "X"
#'   
#' @return entry list
#' 
#' @author Przemyslaw Stempor
#' 
#' @family dbtools
#' @export
#' 
addGenericFile <- function(
    ContactExpID,
    path,
    Processing="raw",
    Resolution=NA,
    Scale=NA,
    genome=NA,
    prefix='X',
    filetype_format="txt.gz",
    comments='',
    userID=NA,
    parent1_uid=NA,
    parent2_uid=NA,
    parent3_uid=NA,
    uniq=NA
    ) {
    
    con <- dbConnect(dbDriver("MySQL"), group = "jadb")
    
    UID=makeFileUID(gsub('\\..+', '',  basename(path)), prefix)
    INSERT <- list(
        ContactExpID=ContactExpID,
        Processing=Processing,
        Resolution=Resolution,
        Scale=Scale,     
        UID=UID,
        filetype_format=filetype_format,
        path=sprintf(
            '%s_%s.%s', 
            file.path(dirname(path), gsub('\\..+$', '', basename(path))),
            UID, filetype_format
        ),     
        dateCreated=paste(Sys.Date()),
        dateUpdated=paste(Sys.Date()),
        comments=comments,
        userID=userID,
        parent1_uid=parent1_uid,
        parent2_uid=parent2_uid,
        parent3_uid=parent3_uid,
        uniq=uniq,
        genome=genome
    )
    dbSendQuery(con, paste(
        "INSERT INTO mydb.labfiles(", 
        paste(names(INSERT), collapse=", "),
        ") VALUES('", paste(INSERT, collapse="', '"), "')",
        collapse=", ", sep=""
    ) )
    dbDisconnect(con)
    return(INSERT)    
}

#' Make generic file name
#' 
#' @param ContactExpID experiment ID string
#' @param EXTABLE experiment table in db schema
#' @param scale file name parameter
#' @param processing file name parameter
#' @param format file name parameter
#'   
#'   
#' @return name string
#' 
#' @author Przemyslaw Stempor
#' 
#' @family dbtools
#' @export
#' 
formGenericPath <- function( 
    ContactExpID, EXTABLE='labrnaseq', Processing=NA, Resolution=NA, Scale=NA
){
    
    con <- dbConnect(dbDriver("MySQL"), group = "jadb")

    PK <- dbGetQuery(con, sprintf("SHOW INDEX FROM %s WHERE Key_name = 'PRIMARY'", gsub('view$', '', EXTABLE) ))[['Column_name']]
    fileds.def <- dbGetQuery(con, sprintf("SHOW FIELDS FROM %s", EXTABLE))
    
    fld <- fileds.def[fileds.def$Null=='NO' & fileds.def$Key!='PRI',]$Field
    EXPERIMENT <-   dbGetQuery(
        con, sprintf('SELECT %s FROM %s WHERE %s = "%s"', paste(fld, collapse = ', '), EXTABLE, PK, ContactExpID)
    )
    fileName=paste(c(EXPERIMENT, 'raw', 'NA', 'NA', ContactExpID), collapse = '_')
    
    dbDisconnect(con)
    return(fileName)
}

getTrackSummary <- function(ID, bin=10, processing='NORM', scale='linear') {
    
    addr <- getFilePath(ID, processing=processing, scale=scale)
}


# message("Reading")
# bw <- import.bw(addr,as="RleList")
# out <- bw * nc
# gr <- as(out, 'GRanges')
# gr <- gr[!is.na(gr$score)]
# gr <- gr[!is.infinite(gr$score)]
# message("Writeing")
# export.bw(gr, rbeads:::reName(basename(addr), proccesing = 'spikeinBEADS', scale = 'zsc'))
# }