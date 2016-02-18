rename.files <- function( id, fields ) {
    
    #setwd('/Library/WebServer/Documents/db4')
    setwd( file.path(www.path, 'db4') )
    old_id <- id
    
    mysql <- dbDriver("MySQL")
    con <- dbConnect(mysql, group = "jadb", default.file='~/.my.cnf')
    
    #Get old files from DB
    files <- dbGetQuery(con, paste("SELECT * FROM labfiles WHERE ContactExpID = '", id,"'", sep=''))
    if(length( files[['path']] ) == 0) {dbDisconnect(con); return(NULL)}
    
    if( GET$table ==  'mydb.labexperiment' ) {
        
        #Create new paths
        if(! is.na( fields['ContactExpID'] ) ) { id <- fields[['ContactExpID']] }
        EXPERIMENT <- 	dbGetQuery(con, paste("SELECT Factor, Antibody, ExtractID, Crosslinker, Strain, Stage FROM mydb.labexperimentview WHERE ContactExpID = '", id, "'", collapse="", sep=""))
        dirPath <- file.path('files', EXPERIMENT[['Factor']], EXPERIMENT[['Strain']], paste(id, EXPERIMENT[['ExtractID']], EXPERIMENT[['Antibody']], sep='_'))
        fileName <- paste(	paste(EXPERIMENT[['Factor']], EXPERIMENT[['Antibody']], sep='^'), 
                           paste(EXPERIMENT[['ExtractID']], EXPERIMENT[['Crosslinker']], EXPERIMENT[['Strain']], EXPERIMENT[['Stage']], sep='^'), 
                           paste(files[['Processing']], files[['Scale']], files[['Resolution']], sep='^'), id, sep='_')
        fileName <- sprintf('%s^%s.%s', fileName, files[['UID']], files[['filetype_format']])
        new.path <- as.character( file.path(dirPath, fileName) )
        
    } else {
        
        EXTABLE <- GET$table
        PK <- dbGetQuery(con, sprintf("SHOW INDEX FROM %s WHERE Key_name = 'PRIMARY'", gsub('view$', '', EXTABLE) ))[['Column_name']]
        fileds.def <- dbGetQuery(con, sprintf("SHOW FIELDS FROM %s", EXTABLE))
        
        fld <- fileds.def[fileds.def$Null=='NO' & fileds.def$Key!='PRI',]$Field
        EXPERIMENT <-   dbGetQuery(
            con, sprintf('SELECT %s FROM %s WHERE %s = "%s"', paste(fld, collapse = ', '), EXTABLE, PK, id)
        )
        
        fileName=paste(c(EXPERIMENT, files[['Processing']], files[['Scale']], files[['Resolution']], id), collapse = '_')
        
        if(toupper(gsub('mydb.lab', '', EXTABLE)) == "RNASEQ") {
            prefix <- 'Q' 
        } else if(toupper(gsub('mydb.lab', '', EXTABLE)) == "DNASE") {
            prefix <- 'A' 
        } else {
            prefix <- 'X'
        }
        if(! is.null(dbReadTable(con, TABLE)$UID) ) { 
            fileUID <- sprintf('%.1s%.2s%05d', prefix, digest::digest(fileName, algo='md5', serialize = FALSE), max(as.numeric(substr(dbReadTable(con, TABLE)$UID, 4, 8)))+1) 
        } else { 
            fileUID <- sprintf('%.1s%.2s%05d', prefix, digest::digest(fileName, algo='md5', serialize = FALSE), 1) 
        }
        fileName <- sprintf('%s_%s.%s', fileName, fileUID, files[['filetype_format']])
        dirPath <- file.path(
            'files', toupper(gsub('mydb.lab', '', EXTABLE)), 
            EXPERIMENT[['CellFraction']], EXPERIMENT[['LibraryType']], id
        )
        
        new.path <- as.character( file.path(dirPath, fileName) )
        
    }
    
    #Return error for invalid paths
    if ( length( unique( dirname(files[['path']]) ) ) != 1 ) {dbDisconnect(con); stop('[FS] Single experiment files in multipe folders: contact admin!')}
    
    #create new directory
    if ( unique( dirname(files[['path']]) ) != dirPath ) dir.create(dirPath, recursive = TRUE)
    
    #Rename files
    n.reg <- 0; n.other <- 0;
    for (i in 1:nrow(files) ) {
        if ( files[['path']][i] != new.path[i] ) {
            if( file.rename(files[['path']][i], new.path[i]) ) {
                fnames = c('path', 'ContactExpID')
                finsert = c( new.path[i], id )
                sql <- paste("UPDATE ", 'labfiles', " SET ", paste(fnames, "='", finsert, "'",  collapse=", ", sep=''), " WHERE ", 'UID', " = '", files[['UID']][i], "'", sep='')
                dbSendQuery(con, sql )
                n.reg <- n.reg + 1
            } else stop( paste( "[FS] Error while copping", warnings()) )
        }
    }	
    dbDisconnect(con)
    
    #Move all other files to new experiment folder
    if ( unique( dirname(files[['path']]) ) != dirPath ) {
        
        #mt <- do.call (rbind, strsplit(basename(files[['path']]), '') )
        #common <- substr( basename( files[['path']] )[[1]], 1, which( diff( apply(mt, 2, function(x) all(x==x[1]) ) ) == -1 )[1] )
        to.move <- dir( unique( dirname(files[['path']]) ), all.files = TRUE, full.names = TRUE, recursive = TRUE )
        if (length(to.move) > 0) {
            to.move <- to.move[!(to.move %in% files[['path']])]
            new.exp.part <- paste( paste(EXPERIMENT[['Factor']], EXPERIMENT[['Antibody']], sep='^'), paste(EXPERIMENT[['ExtractID']], EXPERIMENT[['Crosslinker']], EXPERIMENT[['Strain']], EXPERIMENT[['Stage']], sep='^'), sep='_')
            old.exp.part <- paste(strsplit(basename(files[['path']] )[[1]], '_')[[1]][1:2], collapse = '_')
            #bool.move <- to.move != file.path( dirPath, gsub(old.exp.part, new.exp.part, basename(to.move), fixed = T) )
            #if(new.exp.part != old.exp.part) 
            other.new.names <- file.path( dirPath, gsub(old.exp.part, new.exp.part, basename(to.move), fixed = T))
            other.new.names <- file.path( dirPath, gsub(old_id, id, basename(other.new.names), fixed = T))
            n.other <- file.rename( to.move, other.new.names )
        }
        if ( length( dir( unique( dirname(files[['path']]) ), all.files = TRUE, full.names = TRUE, recursive = TRUE ) ) != 0 ) {
            stop('[FS] Some files failed to move: cannot remove the folder!')
        }
        file.remove( unique( dirname(files[['path']]) ) )
        
        ##Clean-up emty directories
        dir.counts = 0
        stop.me = 0
        while (sum(dir.counts==0)) {
            dir.counts <- sapply(list.dirs('files'), function(x) length(list.files(x)))
            unlink( names(dir.counts)[dir.counts==0], recursive = TRUE)
            stop.me <- stop.me +1
            if(stop.me > 100 ) {stop(c('Infinite loop!', a))}
        }
        
    }	
    return( sprintf('Renamed and/or moved:<br /> - %i registered files <br /> - %i other files.', n.reg, sum(n.other)) )
    
    
    
    
    
}


addFilesFromCsv <- function(csv, root='/mnt/jadb/DBfile/DBfiles', EXTABLE='mydb.labexperiment') {
    if (system('whoami', intern = TRUE) != 'www-data') stop('Run as web server user!', call. = TRUE)	
    
    dir(pattern='txt.gz') %>% gsub('.+(Index[0-9]+).+', '\\1', .) %>% unique -> idx
    for(i in 1:length(idx)){
        out <- unique(gsub('(.+)(s_[0-9].)(Index[0-9]+)(.+)', '\\1BothLanes.\\3\\4', dir(pattern=paste0(idx[i],"\\."))))
        cmd <- sprintf('cat %s > %s', paste(dir(pattern=paste0(idx[i],"\\.")), collapse=' '), out)
        message(cmd)
        system(cmd)
    }
    


    
    
    
    
    library(DBI)
    library(RMySQL)
    library(digest)
    
    
    
    mysql <- dbDriver("MySQL")
    con <- dbConnect(dbDriver("MySQL"), group = "jadb", default.file='~/.my.cnf')
    ALL <- dbReadTable(con, EXTABLE)
    

    
    
    

        
        data <- read.csv(csv)
        records <- as.matrix(data)
        DBrecords <- records[,colnames(records) %in% colnames(ALL)]
        
        for(i in 1:nrow(records)) {
            
            #set yp variables
            insert <- records[i,]
            DBinsert <- DBrecords[i,]
            prID <- as.character( insert[['ProjectID']] )
            smplID <- as.character( insert[['SampleID']] )
            temp_dir <- file.path('files/temp', prID, smplID)
            dir.create(temp_dir, recursive = TRUE)
            
            if (any( grepl('_| |:|\\^|\\/', DBinsert[c('Factor', 'Antibody','ExtractID', 'Crosslinker', 'Strain', 'Stage')]) )) {
                stop('Not allowed character "_" or " " or ":" or "^" or "/" in name fileds.')
            }
            

            files <- ##from csv
        

            
            #File joining, takes time
            system( sprintf('cat %s > %s', paste(file.path(temp_dir, files$Path), collapse=' '), gsub('L002', 'L001andL002', file.path(temp_dir, files$Path)[2])), intern=TRUE)
            finalFilePath <- gsub('L002', 'L001andL002', file.path(temp_dir, files$Path)[2])
            if( !file.exists(finalFilePath) ) stop('Joined file does not exist!')
            
            
            #DATABASE: Add experiment entry
            
            DBinsert[['dateCreated']] 	<- paste(Sys.Date())
            DBinsert[['dateUpdated']] 	<- paste(Sys.Date())
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
                dirPath <- file.path('files', EXPERIMENT[['Factor']], EXPERIMENT[['Strain']], paste(id, EXPERIMENT[['ExtractID']], EXPERIMENT[['Antibody']], sep='_'))
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
                INSERT[['path']] 			<- as.character( file.path(dirPath, fileName) )
                INSERT[['dateCreated']] 	<- paste(Sys.Date())
                INSERT[['dateUpdated']] 	<- paste(Sys.Date())
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
                    'files', toupper(gsub('mydb.lab', '', EXTABLE)), 
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
                INSERT[['path']] 			    <- as.character( file.path(dirPath, fileName) )
                INSERT[['dateCreated']] 	<- paste(Sys.Date())
                INSERT[['dateUpdated']] 	<- paste(Sys.Date())
                INSERT[['uniq']]			    <- NA
                sql <- paste("INSERT INTO ", TABLE,"(", paste(names(INSERT), collapse=", "),") VALUES('", paste(INSERT, collapse="', '"), "')", collapse=", ", sep="")
                
                rs <- dbSendQuery(con, sql )
                
            } 
            
            pb <- list(success=TRUE, progress=i/nrow(records), text=insert[['SampleID']], toptext=paste('File entries created: ', unlist(dbGetInfo(rs, what = "rowsAffected"))))
            cat( RJSONIO::toJSON(pb), file=paste('files/temp/', UID, '.json', sep='') )
            
            Sys.sleep(5)
            
        }
        
        dbDisconnect(con)
        
        pb <- list(success=TRUE, progress=1, text=insert[['SampleID']], toptext='Finalizing')
        cat( RJSONIO::toJSON(pb), file=paste('files/temp/', UID, '.json', sep='') )
        
        message <- paste(	'<br /> <i>Experiments</i>: ',  paste(unique(data$ProjectID), collapse=' '),
                          '<br />	<i>SampleID:</i><br /> ', paste(data$SampleID, collapse=' '),
                          '<br />	<i>ContactExpID:</i><br /> ', paste(data$ContactExpID, collapse=' '),
                          '<br />	<i>fileUID:</i><br /> ', fileUID, 
                          '<br /> <i>Command</i>: ', '' , sep='')
        cat( toJSON(list( success=TRUE, status='<span style="color:green;">Done</span>', details=message )) )
        
        
    } else {
        message <- sprintf('<i>Experiment</i>: %s<br />	<i>File</i>: %s<br />', POST$ContactExpID, FILES[[1]]$name)
        cat(toJSON(list( success=FALSE, file=FILES[[1]]$name, status='Failed', details=message )))
    }
    
}