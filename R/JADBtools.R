#' Get file URL
#' 
#' Accepts regular expressions.
#' 
#' @param ID ContactExperimetID 
#' @param scale filter regex
#' @param processing filter regex 
#' @param format filter regex
#' @param url return url or just the path to the file
#' @param eq t
#'   
#' @return URL string
#' 
#' @author Przemyslaw Stempor
#' 
#' @family dbtools
#' @export
#' 
getFilePath <- function(ID, format='.', processing='.', scale='.', url=TRUE, eq=FALSE, mount=FALSE) {
    
    if(length(ID) > 1) {
        lapply(ID, getFilePath, format=format, processing=processing, scale=scale, url=url, eq=eq, mount=mount)
    } else {
        if (eq) { R <- '=' } else { R <- 'REGEXP' } 
        con <- dbConnect(dbDriver(DRIVER), group = GROUP, default.file='~/.my.cnf')
        
        
        exp_file <- unlist( dbGetQuery(
            con, paste(
                "SELECT path FROM labfiles WHERE ContactExpID = '", ID, 
                "' AND Filetype_format ",R," '",format,"' AND  Processing ",R," '", processing,
                "'", "AND Scale ",R," '", scale, "'", collapse="", sep=""
            ) 
        ), use.names = FALSE)
        
        names(exp_file) <- NULL
        dbDisconnect(con)
        
        if(mount) 
            return( gsub('files', MOUNT,  exp_file) )
        else if(url)
            return( file.path("http://jadb.gurdon.private.cam.ac.uk/db4",  exp_file ) )
        else
            return(exp_file)
    }
    
    
 
}

#' Get file UID
#' 
#' Accepts regular expressions.
#' 
#' @param ID ContactExperimetID 
#' @param scale filter regex
#' @param processing filter regex 
#' @param format filter regex
#' @param url return url or just the path to the file
#' @param eq t
#'   
#' @return URL string
#' 
#' @author Przemyslaw Stempor
#' 
#' @family dbtools
#' @export
#'
getFileUID <- function(ID, format='.', processing='.', scale='.', eq=FALSE) {
    
    if (eq) { R <- '=' } else { R <- 'REGEXP' } 
    con <- dbConnect(dbDriver(DRIVER), group = GROUP, default.file='~/.my.cnf')
    
    uid <- unlist( dbGetQuery(
        con, paste(
            "SELECT UID FROM labfiles WHERE ContactExpID ",R," '", ID, 
            "' AND Filetype_format ",R," '",format,"' AND  Processing ",R," '", processing,
            "'", "AND Scale ",R," '", scale, "'", collapse="", sep=""
        ) 
    ))
    dbDisconnect(con)
    return(uid) 
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
    con <- dbConnect(dbDriver(DRIVER), group = GROUP, default.file='~/.my.cnf')
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
    uniq=NA,
    repPath=FALSE
    ) {
    
    con <- dbConnect(dbDriver(DRIVER), group = GROUP, default.file='~/.my.cnf')
    
    UID=makeFileUID(gsub('\\..+', '',  basename(path)), prefix)
    INSERT <- list(
        ContactExpID=ContactExpID,
        Processing=Processing,
        Resolution=Resolution,
        Scale=Scale,     
        UID=UID,
        filetype_format=filetype_format,
        path=if(repPath) {
            path
        } else {
            sprintf(
                '%s^%s.%s', 
                file.path(dirname(path), gsub('\\..+$', '', basename(path))),
                UID, filetype_format
            )
        },     
        Created=paste(Sys.Date()),
        Updated=paste(Sys.Date()),
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
    
    con <- dbConnect(dbDriver(DRIVER), group = GROUP, default.file='~/.my.cnf')

    PK <- dbGetQuery(con, sprintf("SHOW INDEX FROM %s WHERE Key_name = 'PRIMARY'", gsub('view$', '', EXTABLE) ))[['Column_name']]
    fileds.def <- dbGetQuery(con, sprintf("SHOW FIELDS FROM %s", EXTABLE))
    
    fld <- fileds.def[fileds.def$Null=='NO' & fileds.def$Key!='PRI',]$Field
    EXPERIMENT <-   dbGetQuery(
        con, sprintf('SELECT %s FROM %s WHERE %s = "%s"', paste(fld, collapse = ', '), EXTABLE, PK, ContactExpID)
    )
    fileName=paste(c(EXPERIMENT, Processing, Resolution, Scale, ContactExpID), collapse = '_')
    
    dbDisconnect(con)
    return(fileName)
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
correctName <- function( 
    ContactExpID, EXTABLE='labrnaseq', Processing=NA, Resolution=NA, Scale=NA
    
){
    oldpath <- getFilePath(ContactExpID, url = FALSE)
    
    getFileUID(ContactExpID)
    newpath <- paste0(
        formGenericPath(ContactExpID, EXTABLE='labrnaseq', Processing='raw'),
        substr(oldpath, nchar(oldpath)-15, nchar(oldpath))
    )
    newpath <- file.path(dirname(oldpath), newpath)
    
    con <- dbConnect(dbDriver(DRIVER), group = GROUP, default.file='~/.my.cnf')
    
    EXPERIMENT <-   dbGetQuery(
        con, sprintf('UPDATE labfiles SET path = "%s" WHERE UID = "%s"', newpath, getFileUID(ContactExpID))
    )
    dbDisconnect(con)
    cat('mv', sub('^files/', '', oldpath), sub('^files/', '',newpath), '\n\n')
 
}


getTrackSummary <- function(ID, bin=1000, processing='NORM', scale='linear', chr=NULL, numeric=TRUE) {
    
    addr <- getFilePath(ID, processing=processing, scale=scale, format = 'bw')
    bwf <- BigWigFile(addr)
    if( is.null(chr) ) {
        out <-  unlist(summary(bwf, size = seqlengths(bwf) / bin))
    } else {
        out <-  unlist(summary(bwf, which = as(seqinfo(bwf), 'GRanges')[chr], size = seqlengths(bwf)[chr] / bin))
    }
    seqlengths(out) <- seqlengths(bwf)
    if(numeric) {
        out <- as.numeric(unlist(out)$score)
    }
    return(out)
}

getTracksSummary <- function(IDs, bin=1000, processing='NORM', scale='linear', chr=NULL, numeric=TRUE) {
    
    out <- lapply(IDs, function(x) {
        cat('.')
        getTrackSummary(x, bin, processing, scale, chr, numeric)
    })
    names(out) <- IDs
    return(out)
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