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