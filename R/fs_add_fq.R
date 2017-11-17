#' Add fastq to existing experiment
#'
#' @param id 
#' @param finalFilePath 
#' @param EXTABLE 
#'
#' @return NULL
#' @export
#'
fs_insert_fastq_to_jadb <- function(id, finalFilePath, EXTABLE='labrnaseq') {
    

    library(digest)
    
    mysql <- dbDriver("MySQL")
    con <- dbConnect(dbDriver("MySQL"), group = GROUP, default.file='~/.my.cnf')
    
    
    #DATABASE: Add file entry
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
        message('ChIP fileName: ', fileName)
        
        stop('ZZ')
        
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
        message('RNA-seq fileName: ', fileName)
        
        
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


#' Add files from internal FS
#'
#' @param addr 
#' @param select_id 
#' @param EXTABLE 
#' @param ignore.exist 
#'
#' @return res object
#' @export 
#'
fs_add_internal <- function(addr, select_id='', EXTABLE='labexperiment', ignore.exist=FALSE) {
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
    out <- lapply(lst, validate_jadb_submission_entry, EXTABLE = EXTABLE, ignore.exist = ignore.exist, basespace=FALSE)
    
    if(any(lengths(out)==1)) stop('Validation finished with error!')
    
    catFQandOutput <- function(data) {
        
        insert <- data$insert
        if( length(insert$OryginalFileName_R1) ) {
            R1 <- insert$OryginalFileName_R1
            R2 <- insert$OryginalFileName_R2
            td <- tempdir()
            finalFilePath <- file.path(td, gsub('L002', 'L001andL002', basename(R2)))
            
            if( !file.exists(R1)  ) stop('File R1 does not exist! ', R1)
            if( !file.exists(R2)  ) stop('File R2 does not exist! ', R2)
            
            #File joining, takes time
            system( sprintf('cat "%s" "%s" > %s',  R1, R2, finalFilePath), intern=TRUE)
            
        } else {
            finalFilePath <- insert$OryginalFileName
        }
        
        if( !file.exists(finalFilePath) ) stop('Joined file does not exist!')
        return(finalFilePath)
    }
    
    ##addExtractIDs(addr)
    res <- lapply(out, function(x) {
        temp_file <- catFQandOutput(x)
        insert_entry_to_jadb(x,temp_file)
    })
    
    
    
    
    return(res)
}