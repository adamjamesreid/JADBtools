


mergeLanesFQ <- function() {
    require(magrittr)
    cmds <- vector()
    dir(pattern='txt.gz') %>% gsub('.+(Index[0-9]+).+', '\\1', .) %>% unique -> idx
    for(i in 1:length(idx)){
        out <- unique(gsub('(.+)(s_[0-9].)(Index[0-9]+)(.+)', '\\1BothLanes.\\3\\4', dir(pattern=paste0(idx[i],"\\."))))
        cmd <- sprintf('cat %s > %s', paste(dir(pattern=paste0(idx[i],"\\.")), collapse=' '), out)
        
        cmds <- c(cmds, cmd)
    }
    sapply(cmds, function(x) { message('Joining: ', x); system(x)  })
}

#' Removes single file from JADB
#'
#' @param ids 
#'
#' @return NULL
#' @export
#'
jadb_rm_file <- function(ids) {
    
    con <- dbConnect(dbDriver(DRIVER), group = GROUP, default.file='~/.my.cnf')
    
    pth <- dbGetQuery(con, sprintf(
        "SELECT path FROM `labfiles` WHERE `UID` = '%s';"
    , ids))
    pth <- gsub('^files', MOUNT, pth)
    file.remove(pth)
    
    dbExecute(con, sprintf(
        "DELETE FROM `labfiles` WHERE `UID` = '%s';"
    , ids))
    
    dbDisconnect(con)
    
    message(ids, ' removed [', pth, ']')
    
}


#' Purge an experiment from all processed files
#'
#' @param ids 
#'
#' @return NULL
#' @export
#'
jadb_purge_exp <- function(ids) {
    
    sapply(ids, getFilePath, processing = 'raw', scale = "NA", url = FALSE)  -> fls
    fls <- fls[!fls %in% sapply(ids, getFilePath, format='Rdata', url = FALSE)]
    
    if (length(fls) != 1) stop('No or more than 1 FQ files.')
    
    rm_file <- dir(gsub('^files', MOUNT, dirname(fls)), full.names = TRUE)
    rm_file <- rm_file[!basename(rm_file) %in% basename(fls)]
    unlink(rm_file, recursive=TRUE)
    
    con <- dbConnect(dbDriver(DRIVER), group = GROUP, default.file='~/.my.cnf')
    dbGetQuery(
        con, paste0(
            "DELETE FROM labfiles WHERE ContactExpID = '", ids, "' AND  Processing != 'raw' "
        ) 
    )
    dbGetQuery(
        con, paste0(
            "DELETE FROM labfiles WHERE ContactExpID = '", ids, "' AND  filetype_format = 'Rdata' "
        ) 
    )
    dbDisconnect(con)
    
    message(ids, ' cleaned up!')
    
}


#' remove experiment completly
#'
#' @param ids 
#'
#' @return NULL
#' @export
#'
jadb_renove_exp <- function(ids, EXTABLE='labexperiment') {
    
    sapply(ids, getFilePath, processing = 'raw', scale = "NA", url = FALSE)  -> fls
    fls <- fls[!fls %in% sapply(ids, getFilePath, format='Rdata', url = FALSE)]
    
    if (length(fls)) {
        rm_file <- dir(gsub('^files', MOUNT, dirname(fls)), full.names = TRUE)
        unlink(rm_file, recursive=TRUE)
    }
    
    con <- dbConnect(dbDriver(DRIVER), group = GROUP, default.file='~/.my.cnf')
    dbGetQuery(
        con, paste0(
            "DELETE FROM labfiles WHERE ContactExpID = '", ids, "'"
        ) 
    )
    dbGetQuery(
        con, paste0(
            "DELETE FROM ", EXTABLE," WHERE ContactExpID = '", ids, "'"
        ) 
    )
    dbDisconnect(con)
    
    message('===>> ', ids, ' obliterated!', ' <<===')
    
}

#' Add files to database using local CSV or google spreadsheet
#' 
#' Accepts regular expressions.
#' 
#' @param csv path to CSV or url to gsheet
#' @param root root of the database FS
#' @param EXTABLE table to add the experiment to ("mydb.labrnaseq" for RNAseq)
#' @param gsheet is the source gsheet or local
#'   
#' @return NULL
#' 
#' @author Przemyslaw Stempor
#' 
#' @family FSops
#' @export
#'
addFilesFromCsv <- function(csv, root=MOUNT, EXTABLE='mydb.labexperiment', gsheet=TRUE) {
    if (!grepl(Sys.info()[['user']], "www-data|jarun")) stop('Run as web server user!', call. = TRUE)	
    
    library(DBI)
    library(RMySQL)
    library(digest)
    library(gsheet)
    
    if (gsheet == TRUE) {
        data <- gsheet2tbl(csv)
    } else {
        data <- read.csv(csv)
    }
    
    
 
    
    mysql <- dbDriver(DRIVER)
    con <- dbConnect(dbDriver(DRIVER), group = GROUP, default.file='~/.my.cnf')
    ALL <- dbReadTable(con, EXTABLE)
    
    
    records <- as.matrix(data)
    DBrecords <- records[,colnames(records) %in% colnames(ALL)]
    
    
    
    
    
    
    for(i in 1:nrow(records)) {
        
        message('Processing: ', i)
        
        #set yp variables
        insert <- records[i,]
        DBinsert <- DBrecords[i,]
        #prID <- as.character( insert[['ProjectID']] )
       # smplID <- as.character( insert[['SampleID']] )
        
        if (any( grepl('_| |:|\\^|\\/', DBinsert[c('Factor', 'Antibody','ExtractID', 'Crosslinker', 'Strain', 'Stage')]) )) {
            stop('Not allowed character "_" or " " or ":" or "^" or "/" in name fileds.')
        }
        
        
        finalFilePath <- insert['OryginalFileName']
        if( !file.exists(finalFilePath) ) stop('Joined file does not exist!')
        
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
        
        rs <- dbSendQuery(con,  sql)
        
        id <- insert['ContactExpID']
        TABLE <- 'labfiles'
        
        if( EXTABLE ==  'mydb.labexperiment' ) {
            
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
            
            dir.success <- dir.create(file.path(root, dirPath), recursive = TRUE)
            if ( !file.exists(finalFilePath) ) {dbDisconnect(con); stop(paste('Temp file do not exists', insert['OryginalFileName']))}
            file.copy( finalFilePath, file.path(file.path(root, dirPath), fileName) )
            
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
            
            dir.success <- dir.create(file.path(root, dirPath), recursive = TRUE)
            if ( !file.exists(finalFilePath) ) {dbDisconnect(con); stop(paste('Temp file do not exists', insert['OryginalFileName']))}
            file.copy( finalFilePath, file.path(file.path(root, dirPath), fileName) )
            
            INSERT <- list()
            INSERT[['Processing']]    <- as.character( 'raw' 		)
            INSERT[['Scale']] 			  <- as.character( 'NA'		)
            INSERT[['ContactExpID']] 	<- as.character( id	)
            INSERT[['UID']] 			    <- as.character( fileUID )
            INSERT[['filetype_format']] <- as.character( 'txt.gz'	)
            INSERT[['path']] 			    <- as.character( file.path("files", dirPath, fileName) )
            INSERT[['Created']] 	<- paste(Sys.Date())
            INSERT[['Updated']] 	<- paste(Sys.Date())
            INSERT[['uniq']]			    <- NA
            sql <- paste("INSERT INTO ", TABLE,"(", paste(names(INSERT), collapse=", "),") VALUES('", paste(INSERT, collapse="', '"), "')", collapse=", ", sep="")
            
            rs <- dbSendQuery(con, sql )
            
        } 
        
        
    }
    
    message("Success")
    
    
} 


#' Add files to database using local CSV or google spreadsheet
#' 
#' Accepts regular expressions.
#' 
#' @param csv path to CSV or url to gsheet
#' @param root root of the database FS
#' @param EXTABLE table to add the experiment to ("mydb.labrnaseq" for RNAseq)
#' @param gsheet is the source gsheet or local
#'   
#' @return NULL
#' 
#' @author Przemyslaw Stempor
#' 
#' @family FSops
#' @export
#'
addFilesFromBaseSpace <- function(csv, root=MOUNT, EXTABLE='mydb.labexperiment', gsheet=TRUE, subset=NULL) {
    
    if (!grepl(Sys.info()[['user']], "www-data|jarun")) stop('Run as web server user!', call. = TRUE)	
	
    #version 2.0
    
    setwd( root )	

    library(RJSONIO)
    library(DBI)
    library(RMySQL)
    library(digest)
    
        library(BaseSpaceR)
        app_access_token <- "f58a0ccf1599418d8b4b09034a56bdd5"
        aAuth <- AppAuth(access_token = app_access_token, scope = "browse Sample")
        myProj <- listProjects(aAuth, Limit = 1000)
        PrAno <- data.frame(Name = Name(myProj), Id = Id(myProj))
        
        mysql <- dbDriver(DRIVER)
        con <- dbConnect(dbDriver(DRIVER), group = GROUP, default.file='~/.my.cnf')
        
        ALL <- dbReadTable(con, EXTABLE)
        
        if (gsheet == TRUE) {
            if(grepl('csv', csv)) {
                require(RCurl)
                myCsv <- getURL(csv)
                conCsv <- textConnection(myCsv)
                data <- read.csv(conCsv)
                close(conCsv)
            } else {
                require(RCurl)
                library(stringr)
                key <- stringr::str_extract(csv, "[[:alnum:]_-]{30,}")
                address <- paste0("https://docs.google.com/spreadsheets/d/",key,"/export?exportFormat=", 'csv')
                myCsv <- RCurl::getURL(address)
                conCsv <- textConnection(myCsv)
                data <- read.csv(conCsv)
                close(conCsv)
            }
            
        } else {
            data <- read.csv(csv)
        }
        
        if(length(subset)) data <- data[data$ContactExpID == subset]
        
        records <- as.matrix(data)
        DBrecords <- records[,colnames(records) %in% colnames(ALL)]
        
        for(i in 1:nrow(records)) {
            
            #set yp variables
            insert <- records[i,]
            DBinsert <- DBrecords[i,]
            prID <- as.character( insert[['ProjectID']] )
            smplID <- as.character( insert[['SampleID']] )
            temp_dir <- file.path('temp', prID, smplID)
            dir.create(temp_dir, recursive = TRUE)
            
            if (any( grepl('_| |:|\\^|\\/', DBinsert[c('Factor', 'Antibody','ExtractID', 'Crosslinker', 'Strain', 'Stage')]) )) {
                stop('Not allowed character "_" or " " or ":" or "^" or "/" in name fileds.')
            }
            
            
            mySmpl <- listSamples(aAuth, projectId=subset(PrAno, Name == prID, Id, drop = TRUE), Limit=1000)
            SmAno <- data.frame(Name = Name(mySmpl), Id = Id(mySmpl))
            files <- listFiles(aAuth, sampleId = subset(SmAno, Name == smplID, Id, drop = TRUE))
            
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
            
            rs <- dbSendQuery(con,  sql)
            
            # }
            
            #DATABASE: Add file entry
            
            id <- insert['ContactExpID']
            TABLE <- 'labfiles'
            if( EXTABLE ==  'mydb.labexperiment' ) {
                
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
            
            message(' - OK')
            
            
        }
        
        dbDisconnect(con)
        
        message <- paste(	'<br /> <i>Experiments</i>: ',  paste(unique(data$ProjectID), collapse=' '),
                          '<br />	<i>SampleID:</i><br /> ', paste(data$SampleID, collapse=' '),
                          '<br />	<i>ContactExpID:</i><br /> ', paste(data$ContactExpID, collapse=' '),
                          '<br />	<i>fileUID:</i><br /> ', fileUID, 
                          '<br /> <i>Command</i>: ', '' , sep='')
        cat( toJSON(list( success=TRUE, status='<span style="color:green;">Done</span>', details=message )) )
        
        
        # } else {
        #   message <- sprintf('<i>Experiment</i>: %s<br />	<i>File</i>: %s<br />', POST$ContactExpID, FILES[[1]]$name)
        #  cat(toJSON(list( success=FALSE, file=FILES[[1]]$name, status='Failed', details=message )))
        # }
        

    
    message("Success")
    
    return(data$ContactExpID)
    
    
} 



#' Add files to database using local CSV or google spreadsheet
#' 
#' Accepts regular expressions.
#' 
#' @param csv path to CSV or url to gsheet
#' @param root root of the database FS
#' @param EXTABLE table to add the experiment to ("mydb.labrnaseq" for RNAseq)
#' @param gsheet is the source gsheet or local
#'   
#' @return NULL
#' 
#' @author Przemyslaw Stempor
#' 
#' @family FSops
#' @export
#'
validateFilesFromBaseSpace <- function(csv, EXTABLE='mydb.labexperiment', gsheet=TRUE) {
    
    library(DBI)
    library(RMySQL)
    library(digest)
    library(gsheet)
    
    if (gsheet == TRUE) {
        if(grepl('csv', csv)) {
            require(RCurl)
            myCsv <- getURL(csv)
            con <- textConnection(myCsv)
            data <- read.csv(con)
            close(con)
        } else {
            require(RCurl)
            library(stringr)
            key <- stringr::str_extract(csv, "[[:alnum:]_-]{30,}")
            address <- paste0("https://docs.google.com/spreadsheets/d/",key,"/export?exportFormat=", 'csv')
            myCsv <- RCurl::getURL(address)
            con <- textConnection(myCsv)
            data <- read.csv(con)
            close(con)
        }
        
    } else {
        data <- read.csv(csv)
    }
    
    if( length(unique(data$ProjectID)) < 1 ) stop('FATAL: No BaseSpace project ID found!')
    if( length(unique(data$ProjectID)) > 1 ) stop('FATAL: More than one BaseSpace project ID in form, split to multiple forms.')
    
    prID <- as.character(unique(data$ProjectID))
    
    library(BaseSpaceR)
    app_access_token <- "f58a0ccf1599418d8b4b09034a56bdd5"
    aAuth <- AppAuth(access_token = app_access_token, scope = "browse Sample")
    myProj <- listProjects(aAuth, Limit = 1000)
    PrAno <- data.frame(Name = Name(myProj), Id = Id(myProj))
    
    
    if(!any(PrAno$Name == prID)) stop('FATAL: Unable to mach BaseSpace project ID: ', prID, '. Projects on BaseSpace:\n', paste(PrAno$Name, collapse=', '))
    
    mysql <- dbDriver(DRIVER)
    con <- dbConnect(dbDriver(DRIVER), group = GROUP, default.file='~/.my.cnf')
    ALL <- dbReadTable(con, EXTABLE)
    extract <- dbReadTable(con, 'labextract')
    
    PK <- dbGetQuery(con, sprintf("SHOW INDEX FROM %s WHERE Key_name = 'PRIMARY'", gsub('view$', '', EXTABLE) ))[['Column_name']]
    fileds.def <- dbGetQuery(con, sprintf("SHOW FIELDS FROM %s", EXTABLE))
    dbDisconnect(con)
    
    records <- as.matrix(data)
    DBrecords <- records[,colnames(records) %in% colnames(ALL)]
    

    
    fld <- fileds.def[fileds.def$Null=='NO' & fileds.def$Key!='PRI',]$Field
    
    for(i in 1:nrow(records)) {
        
        cat('Processing: ', DBrecords[i,c(1, 4, 5)])
        
        #set yp variables
        insert <- records[i,]
        DBinsert <- DBrecords[i,]
        prID <- as.character( insert[['ProjectID']] )
        smplID <- as.character( insert[['SampleID']] )
        
        if (any( grepl('_| |:|\\^|\\/', DBinsert[fld]) )) {
            stop('FATAL: Not allowed character "_" or " " or ":" or "^" or "/" in file name fileds:\n', paste(fld, collapse=', '))
        }
        
        
        mySmpl <- listSamples(aAuth, projectId=subset(PrAno, Name == prID, Id, drop = TRUE), Limit=1000)
        SmAno <- data.frame(Name = Name(mySmpl), Id = Id(mySmpl))
        files <- listFiles(aAuth, sampleId = subset(SmAno, Name == smplID, Id, drop = TRUE))
        
        if(!any(SmAno$Name == smplID)) warning('Experiment ID [', smplID, '] does not match one(s) in BaseSpace, allowed values are:\n', paste(SmAno$Name, collapse=', '))
        
        if(length(files)==2) cat(' - SE experiemt')
        if(length(files)==4) cat(' - PE experiemt')
        if(length(files)==0) stop('No files in BaseSpace')
        #File joining, takes time
        #system( sprintf('cat %s > %s', paste(file.path(temp_dir, files$Path), collapse=' '), gsub('L002', 'L001andL002', file.path(temp_dir, files$Path)[2])), intern=TRUE)
        #finalFilePath <- gsub('L002', 'L001andL002', file.path(temp_dir, files$Path)[2])
        
        
        
        id <- insert['ContactExpID']
        if(!grepl('^r*[A-Z]{2}[0-9]{3}$', id)) stop('Missformated experiemt ID')
        
        #EXPERIMENT <-   dbGetQuery(
         #   con, sprintf('SELECT %s FROM %s WHERE %s = "%s"', paste(fld, collapse = ', '), EXTABLE, PK, id)
        #)
        if(any(id == ALL$ContactExpID)) stop('ContactExperimentID ', id, ' alredy exists in database')
        
        
        if (EXTABLE=='mydb.labexperiment') {
            if( !any(insert[['ExtractID']] == extract$ExtractID) ) stop('[', id, ']: extract ',insert[['ExtractID']], ' not found in extracts table.')
        }
        
    cat(" - [OK] \n")
    
    
    } 
}

#' Add extract to database using local CSV or google spreadsheet
#' 
#' Accepts regular expressions.
#' 
#' @param csv path to CSV or url to gsheet
#' @param gsheet is the source gsheet or local
#'   
#' @return NULL
#' 
#' @author Przemyslaw Stempor
#' 
#' @family FSops
#' @export
#'
addExtractIDs <- function(csv, gsheet=TRUE) {
    library(DBI)
    library(RMySQL)
    library(digest)
    library(gsheet)
    require(dplyr)
    
    #if (gsheet == TRUE) {
    #    data <- gsheet2tbl(csv)
    #} else {
    #    data <- read.csv(csv)
    #}
    
    data <- tbl_df(get_tab_from_google(csv))
    
    data %>% select(ExtractID, Crosslinking, Strain, Stage) %>% unique -> ext
    
    TABLE <- 'labextract'

    mysql <- dbDriver(DRIVER)
    con <- dbConnect(dbDriver(DRIVER), group = GROUP, default.file='~/.my.cnf')
    ALL <- dbReadTable(con, TABLE)
     
    fileds.def <- dbGetQuery(con, sprintf("SHOW FIELDS FROM %s", TABLE))
    PK <- dbGetQuery(con, sprintf("SHOW INDEX FROM %s WHERE Key_name = 'PRIMARY'", gsub('view$', '', TABLE) ))[['Column_name']]
    
    for(i in 1:nrow(ext)) {
        
        if(ext$ExtractID[i] %in% ALL$ExtractID) next()
        
        INSERT <- list()
        INSERT[['ExtractID']] 	<- as.character( ext$ExtractID[i] )
        INSERT[['Crosslinker']] <- as.character( ext$Crosslinking[i] )
        INSERT[['Strain']]      <- as.character( ext$Strain[i]	)
        INSERT[['Stage']] 		<- as.character( ext$Stage[i] )
        
        sql <- paste("INSERT INTO ", TABLE,"(", paste(names(INSERT), collapse=", "),") VALUES('", paste(INSERT, collapse="', '"), "')", collapse=", ", sep="")
        
        rs <- dbSendQuery(con, sql )
    }
    
    dbDisconnect(con)
}