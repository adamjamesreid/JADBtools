callPeaksMACS <- function(ids, local=TRUE, outdir='.') {
    con <- dbConnect(dbDriver("MySQL"), group = "jadb", default.file='~/.my.cnf')
    all_experiments <<- dbReadTable(con, "labexperimentview")
    dbDisconnect(con)
    
    library(dplyr)
    
    getFileInfo <- function (ID, info, format = ".", processing = ".", scale = ".", url = TRUE,
                             eq = FALSE)
    {
        if (eq) {
            R <- "="
        }
        else {
            R <- "REGEXP"
        }
        con <- dbConnect(dbDriver("MySQL"), group = "jadb")
        exp_file <- unlist(dbGetQuery(con, paste("SELECT ", info, " FROM labfiles WHERE ContactExpID ",
                                                 R, " '", ID, "' AND Filetype_format ", R, " '", format,
                                                 "' AND  Processing ", R, " '", processing, "'", "AND Scale ",
                                                 R, " '", scale, "'", collapse = "", sep = "")), use.names = FALSE)
        names(exp_file) <- NULL
        dbDisconnect(con)
        return(exp_file)
    }
    
    getAlignedReads <- function(x) {
        info = getFileInfo(x, 'Comments', format = 'bam')
        as.numeric(gsub('.*INFO: ([0-9]+).+', '\\1', info))
    }
    
    
    inp <- lapply(ids, function(x) {
        extract <- JADBtools::getAnno(x, anno = 'ExtractID', EXTABLE = 'labexperiment')
        
        if(extract == 'mE12') extract <- 'mE11'
        
        inputs <- filter(all_experiments, ExtractID==extract, Factor=='Input')$ContactExpID
        
        if(length(inputs) > 1) {
            names(which.max(sapply(inputs[[1]], getAlignedReads)))
        } else {
            inputs
        }
    })
    
    ids  %>% sapply(getFilePath, format = 'bam', url = local)  -> fls
    unlist(inp) %>% sapply(getFilePath, format = 'bam', url = local)  -> inputs
    
    outnames <- sapply(basename(fls), rbeads:::reName, proccesing = 'PeakCalls', scale = 'MACS', resolution = 'q01', ext='.narowPeak')
    outnames_summits <- sapply(basename(fls), rbeads:::reName, proccesing = 'summits', scale = 'MACS', resolution = 'q01', ext='.bed')
    prefix <- basename(fls) %>% substr(start=0, stop=nchar(.)-13)
    
    

    
    if(local) {
        download.file(fls, basename(fls))
        download.file(inputs, basename(inputs))
        
        command <- '~/anaconda/envs/py27/bin/macs2 callpeak -t %s -c %s -f BAM -g ce -n %s -q 0.01 2>&1 | tee %s'
        
        cmd <- sprintf(
            command,
            basename(fls),
            basename(inputs),
            prefix,
            paste0(prefix, '_log.txt')
        )
        message(cmd)
        system(cmd)
        return(c(paste0(prefix, '_peaks.narrowPeak'), paste0(prefix, '_summits.bed')))
        
    } else {
        
        base_dir  <- getwd()
        exp_dir <- gsub('files/', '', dirname(fls))
        message(exp_dir)
        setwd(exp_dir)
        
        command <- '/home/ps562/anaconda2/bin/macs2 callpeak -t %s -c %s -f BAM -g ce -n %s -q 0.01 2>&1 | tee %s'
        
        cmd <- sprintf(
            command,
            gsub('files', '/mnt/jadb/DBfile/DBfiles', fls),
            gsub('files', '/mnt/jadb/DBfile/DBfiles', inputs),
            basename(fls) %>% substr(start=0, stop=nchar(.)-13),
            paste0(basename(fls) %>% substr(start=0, stop=nchar(.)-13), '_log.txt')
        )
        
        message(cmd)
        system(cmd)
        
        if(!file.exists(paste0(prefix, '_peaks.narrowPeak'))) stop('No peaks created, see log for details')
        
        peakEntry <- addGenericFile(
            ids, 
            path = file.path('files', exp_dir, gsub('aligned^NA^NA', 'PeakCalls^MACS^q01', prefix)), 
            Processing = 'PeakCalls', 
            Scale = 'MACS', 
            Resolution = 'q01', 
            filetype_format = 'narrowPeak', 
            prefix = 'P',
            comments=paste('input: ', inp)
        )
        
        summitEntry <- addGenericFile(
            ids, 
            path = file.path('files', exp_dir, gsub('aligned^NA^NA', 'summits^MACS^q01', prefix)), 
            Processing = 'summits', 
            Scale = 'MACS', 
            Resolution = 'q01', 
            filetype_format = 'bed', 
            prefix = 'P',
            comments=paste('input: ', inp)
        )
        
        file.rename(paste0(prefix, '_peaks.narrowPeak'),  basename(peakEntry$path))
        file.rename(paste0(prefix, '_summits.bed'), basename(summitEntry$path))
        
        setwd(base_dir)
    }

}