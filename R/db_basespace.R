#' jagui_get_table
#'
#' @param tab name of MySQL table, e.g. 'labexperiment'
#'
#' @return tibble
#' @export
#'
#' @examples
#' # jagui_get_table('labexperiment')
jagui_get_table <- function(tab) { 
    library(dplyr)
    library(tibble)
    library(RMySQL)
    
    con <- dbConnect(dbDriver("MySQL"), group = GROUP, default.file='~/.my.cnf')
    
    fileds.def <- dbGetQuery(con, sprintf("SHOW FIELDS FROM %s", tab))
    PK <- dbGetQuery(con, sprintf("SHOW INDEX FROM %s WHERE Key_name = 'PRIMARY'", gsub('view$', '', tab) ))[['Column_name']]
    PF <- dbListFields(con, gsub('view', '', tab))
    
    
    TAB <- as_tibble(dbReadTable(con, tab))
    dbDisconnect(con)
    
    if('Comments' %in% colnames(TAB)) {
        TAB$Comments <- enc2utf8(TAB$Comments)
    }
    
    attr(TAB, 'fileds.def') <- fileds.def
    attr(TAB, 'PK') <- PK
    attr(TAB, 'PF') <- PF
    return(TAB)
}

#' Title
#'
#' @param addr goodgle spreadsheets adderss
#'
#' @return data.frame
#' @export
#'
get_tab_from_google <- function(addr) {
    
    library(RCurl)
    
    if (grepl("csv", addr)) {
        myCsv <- RCurl::getURL(addr)
    } else {
        library(stringr)
        key <- stringr::str_extract(addr, "[[:alnum:]_-]{30,}")
        address <- paste0("https://docs.google.com/spreadsheets/d/", 
                          key, "/export?exportFormat=", "csv")
        myCsv <- RCurl::getURL(address)
    }
    
    con <- textConnection(myCsv)
    data <- read.csv(con, stringsAsFactors = FALSE, colClasses = c("character"))
    close(con)
    
    return(data)
}


#' validate_jadb_submission_entry
#'
#' @param entry 
#' @param EXTABLE 
#' @param ignore.exist 
#'
#' @return data structuren (list)
#' @export
#'
validate_jadb_submission_entry <- function(
    entry, 
    EXTABLE = "mydb.labexperiment", 
    ignore.exist=FALSE
) {
    
    library(BaseSpaceR)
    app_access_token <- "f58a0ccf1599418d8b4b09034a56bdd5"
    aAuth <- AppAuth(access_token = app_access_token, scope = "browse Sample")
    myProj <- listProjects(aAuth, Limit = 1000)
    PrAno <- data.frame(Name = Name(myProj), Id = Id(myProj))
    
    ALL <- jagui_get_table(EXTABLE)
    extract <- jagui_get_table("labextract")
    
    
    ## FIX in mysql def
    fileds.def <- attr(ALL, "fileds.def")
    fld <- fileds.def[fileds.def$Null == "NO" & fileds.def$Key != "PRI", ]$Field
    
    message("Processing: ", entry$ContactExpID, ': ', entry$SampleID, ', ', entry$Factor, ' - ')
    
    
    if (!any(PrAno$Name == entry$ProjectID)) {
        message("[ERROR] Unable to mach BaseSpace project ID: ", 
                entry$ProjectID, ". Projects on BaseSpace:\n", paste(PrAno$Name, 
                                                                     collapse = ", "))
        return(FALSE)
    } 
    
    # Deffine ids
    DBinsert <- entry[names(entry) %in% colnames(ALL)]
    prID <- as.character(entry[["ProjectID"]])
    smplID <- as.character(entry[["SampleID"]])
    
    
    # Validation steps
    
    # check if names ok
    if (any(grepl("_| |:|\\^|\\/", DBinsert[fld]))) {
        message("[ERROR] Not allowed character [_], [ ], [:], [^] or [/] in file name fileds: ", 
                paste(fld, collapse = ", "))
        return(FALSE)
    }
    
    mySmpl <- listSamples(
        aAuth, 
        projectId = subset(PrAno, Name == prID, Id, drop = TRUE), 
        Limit = 1000
    )
    SmAno <- data.frame(Name = Name(mySmpl), Id = Id(mySmpl))
    files <- listFiles(
        aAuth, 
        sampleId = subset(SmAno, Name == smplID, Id, drop = TRUE)[[1]]
    )
    
    if (!any(SmAno$Name == smplID)) {
        warning("Experiment ID [", smplID, "] does not match one(s) in BaseSpace, allowed values are:\n", 
                paste(SmAno$Name, collapse = ", "))
    }
    
    if (length(files) == 2) 
        message("SE experiemt - ")
    if (length(files) == 4) 
        message("PE experiemt - ")
    if (length(files) == 0) {
        message("[ERROR] No files in BaseSpace")
        return(FALSE)
    }
    
    id <- entry["ContactExpID"]
    
    if (!grepl("^r*[A-Z]{2}[0-9]{3}$", id)) {
        message("[ERROR] Missformated experiemt ID")
        return(FALSE)
    }
    
    if (any(id == ALL$ContactExpID) & !ignore.exist) {
        message("[ERROR] ContactExperimentID ", id, " alredy exists in database")
        return(FALSE)
    }
    
    if (EXTABLE == "mydb.labexperiment") {
        if (!any(entry[["ExtractID"]] == extract$ExtractID)) {
            message("[", id, "]: extract ", entry[["ExtractID"]], 
                    " not found in extracts table.")
            return(FALSE)
        }
    }
    
    message("[OK] \n")
    return(list(
        files=files,
        insert=entry,
        DBinsert=DBinsert,
        aAuth=aAuth,
        EXTABLE=EXTABLE
    ))
}


#' download_fastq_from_basespace
#'
#' @param data 
#'
#' @return file path
#' @export
#'
download_fastq_from_basespace <- function(data) {
    
    insert <- data$insert
    files <- data$files
    aAuth <- data$aAuth
    
    prID <- as.character( insert[['ProjectID']] )
    smplID <- as.character( insert[['SampleID']] )
    temp_dir <- file.path('temp', prID, smplID)
    dir.create(temp_dir, recursive = TRUE)
    
    #if( insert[['SeqType']] == 'PE') {
    #  stop('Not yet implemented')      
    #} else {
    
    ##Single end force mode, TODO: expand on 2nd end
    files <- files[grepl('R1', Name(files))]
    
    ##stop(subset(SmAno, Name == smplID, Id, drop = TRUE))
    
    message('Procesing: ', paste(insert, collapse = ' | '))
    #File download, takes time
    getFiles(aAuth, id = Id(files), destDir = temp_dir, verbose = TRUE)
    
    #File size check
    if( !all(files$Size == file.info(file.path(temp_dir, files$Path))$size) ) {
        stop('Downloaded file sizes does not match the inforamtion from BaseSpace.')
    }
    
    #File joining, takes time
    system( sprintf('cat %s > %s', paste(file.path(temp_dir, files$Path), collapse=' '), gsub('L002', 'L001andL002', file.path(temp_dir, files$Path)[2])), intern=TRUE)
    finalFilePath <- gsub('L002', 'L001andL002', file.path(temp_dir, files$Path)[2])
    if( !file.exists(finalFilePath) ) stop('Joined file does not exist!')
    
    return(finalFilePath)
}


#' insert_entry_to_jadb
#'
#' @param data 
#' @param finalFilePath 
#'
#' @return NULL
#' @export
#'
insert_entry_to_jadb <- function(data, finalFilePath) {
    
    insert <- data$insert
    DBinsert <- data$DBinsert
    EXTABLE <- data$EXTABLE
    
    #DATABASE: Add experiment entry
    
    DBinsert[['Created']] 	<- paste(Sys.Date())
    DBinsert[['Updated']] 	<- paste(Sys.Date())
    DBinsert[['OryginalFileName']]   <- basename(finalFilePath)
    fields <- names(DBinsert)
    
    sql <- paste(
        sprintf("INSERT INTO %s(", EXTABLE),
        paste(fields, collapse=", "), ") VALUES('",
        paste(DBinsert, collapse="', '"), "')", collapse=", ", sep=""
    )
    
    
    mysql <- dbDriver("MySQL")
    con <- dbConnect(dbDriver("MySQL"), group = GROUP, default.file='~/.my.cnf')
    
    rs <- dbSendQuery(con,  sql)
    
    
    #DATABASE: Add file entry
    
    id <- insert['ContactExpID']
    TABLE <- 'labfiles'
    if( grepl('labexperiment', EXTABLE) ) {
        
        EXPERIMENT <- 	dbGetQuery(con, paste("SELECT Factor, Antibody, ExtractID, Crosslinker, Strain, Stage FROM mydb.labexperimentview WHERE ContactExpID = '", id, "'", collapse="", sep=""))
        dirPath <- file.path(EXPERIMENT[['Factor']], EXPERIMENT[['Strain']], paste(id, EXPERIMENT[['ExtractID']], EXPERIMENT[['Antibody']], sep='_'))
        fileName <- paste(	paste(EXPERIMENT[['Factor']], EXPERIMENT[['Antibody']], sep='^'), 
                           paste(EXPERIMENT[['ExtractID']], EXPERIMENT[['Crosslinker']], EXPERIMENT[['Strain']], EXPERIMENT[['Stage']], sep='^'), 
                           paste('raw', 'NA', 'NA', sep='^'), id, sep='_')
        if(! is.null(dbReadTable(con, TABLE)$UID) ) { 
            fileUID <- sprintf('%.1s%.2s%05d', 'C', digest(fileName, algo='md5', serialize = FALSE), max(as.numeric(substr(dbReadTable(con, TABLE)$UID, 4, 8)))+1) 
        } else { 
            fileUID <- sprintf('%.1s%.2s%05d', 'C', digest(fileName, algo='md5', serialize = FALSE), 1) 
        }
        fileName <- sprintf('%s^%s.%s', fileName, fileUID, 'txt.gz')
        
        message('ChIP dirPath: ', dirPath)
        dir.success <- dir.create(dirPath, recursive = TRUE)
        if ( !file.exists(finalFilePath) ) {dbDisconnect(con); stop(paste('Temp file do not exists', insert['OryginalFileName']))}
        file.copy( finalFilePath, file.path(dirPath, fileName) )
        
        fileds.def <- dbGetQuery(con, sprintf("SHOW FIELDS FROM %s", TABLE))
        PK <- dbGetQuery(con, sprintf("SHOW INDEX FROM %s WHERE Key_name = 'PRIMARY'", gsub('view$', '', TABLE) ))[['Column_name']]
        
        INSERT <- list()
        INSERT[['Processing']] 		<- as.character( 'raw' 		)
        INSERT[['Scale']] 			<- as.character( 'NA'		)
        INSERT[['ContactExpID']] 	<- as.character( id	)
        INSERT[['UID']] 			<- as.character( fileUID )
        INSERT[['filetype_format']] <- as.character( 'txt.gz'	)
        INSERT[['path']] 			<- as.character( file.path('files', dirPath, fileName) )
        INSERT[['Created']] 	<- paste(Sys.Date())
        INSERT[['Updated']] 	<- paste(Sys.Date())
        INSERT[['uniq']]			<- NA
        sql <- paste("INSERT INTO ", TABLE,"(", paste(names(INSERT), collapse=", "),") VALUES('", paste(INSERT, collapse="', '"), "')", collapse=", ", sep="")
        
        rs <- dbSendQuery(con, sql )
        
    } else {
        
        PK <- dbGetQuery(con, sprintf("SHOW INDEX FROM %s WHERE Key_name = 'PRIMARY'", gsub('view$', '', EXTABLE) ))[['Column_name']]
        fileds.def <- dbGetQuery(con, sprintf("SHOW FIELDS FROM %s", EXTABLE))
        
        fld <- fileds.def[fileds.def$Null=='NO' & fileds.def$Key!='PRI',]$Field
        EXPERIMENT <-   dbGetQuery(
            con, sprintf('SELECT %s FROM %s WHERE %s = "%s"', paste(fld, collapse = ', '), EXTABLE, PK, id)
        )
        fileName=paste(c(EXPERIMENT, 'raw', 'NA', 'NA', id), collapse = '_')
        
        if(toupper(gsub('mydb.lab', '', EXTABLE)) == "RNASEQ") {
            prefix <- 'Q' 
        } else if(toupper(gsub('mydb.lab', '', EXTABLE)) == "DNASE") {
            prefix <- 'A' 
        } else {
            prefix <- 'X'
        }
        if(! is.null(dbReadTable(con, TABLE)$UID) ) { 
            fileUID <- sprintf('%.1s%.2s%05d', prefix, digest(fileName, algo='md5', serialize = FALSE), max(as.numeric(substr(dbReadTable(con, TABLE)$UID, 4, 8)))+1) 
        } else { 
            fileUID <- sprintf('%.1s%.2s%05d', prefix, digest(fileName, algo='md5', serialize = FALSE), 1) 
        }
        fileName <- sprintf('%s_%s.%s', fileName, fileUID, 'txt.gz')
        dirPath <- file.path(
            toupper(gsub('mydb.lab', '', EXTABLE)), 
            EXPERIMENT[['CellFraction']], EXPERIMENT[['LibraryType']], id
        )
        
        message('RNA-seq dirPath: ', dirPath)
        dir.success <- dir.create(dirPath, recursive = TRUE)
        if ( !file.exists(finalFilePath) ) {dbDisconnect(con); stop(paste('Temp file do not exists', insert['OryginalFileName']))}
        file.copy( finalFilePath, file.path(dirPath, fileName) )
        
        INSERT <- list()
        INSERT[['Processing']]    <- as.character( 'raw' 		)
        INSERT[['Scale']] 			  <- as.character( 'NA'		)
        INSERT[['ContactExpID']] 	<- as.character( id	)
        INSERT[['UID']] 			    <- as.character( fileUID )
        INSERT[['filetype_format']] <- as.character( 'txt.gz'	)
        INSERT[['path']] 			    <- as.character( file.path('files', dirPath, fileName) )
        INSERT[['Created']] 	<- paste(Sys.Date())
        INSERT[['Updated']] 	<- paste(Sys.Date())
        INSERT[['uniq']]			    <- NA
        sql <- paste("INSERT INTO ", TABLE,"(", paste(names(INSERT), collapse=", "),") VALUES('", paste(INSERT, collapse="', '"), "')", collapse=", ", sep="")
        
        rs <- dbSendQuery(con, sql )
        
    } 
    
    dbDisconnect(con)
    message(' - OK')
}


#' Title
#'
#' @param addr 
#' @param EXTABLE 
#' @param ignore.exist 
#' @param select_id 
#'
#' @return result
#' @export
#'
jadb_basespace <- function(addr, EXTABLE='mydb.labexperiment', ignore.exist = FALSE, select_id='') {
    
    
    # JADBtools:::jadb_renove_exp("AA690")
    
    owd <- getwd()
    on.exit(setwd(owd))
    setwd(file.path(MOUNT))
    
    message(getwd())
    
    dat <- get_tab_from_google(addr)
    
    if(nchar(select_id)) {
        dat <-  dat[dat$ContactExpID==select_id,]
    }
    
    lst <- apply(dat, 1, as.list)
    out <- lapply(lst, validate_jadb_submission_entry, EXTABLE = EXTABLE, ignore.exist = ignore.exist)
    
    if(any(lengths(out)==1)) stop('Validation finished with error!')
    
    res <- lapply(out, function(x) {
        temp_file <- download_fastq_from_basespace(x)
        insert_entry_to_jadb(x,temp_file)
    })
    return(res)
}
    








    

