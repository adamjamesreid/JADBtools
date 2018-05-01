macs2 <- function(ids, local=TRUE, extsize=150, qval=.99, summedinput=TRUE, genome=NULL) {
    
    if(is.numeric(extsize)) {
        message('Extsize set to ', extsize)
        command <- paste0(
            'macs2 callpeak -t %s -c %s --nomodel --llocal 1000 --extsize ',
            extsize, ' -f BAM -g ce -B -n %s -q ', qval, ' 2>&1 | tee %s'
        )
    } else {
        command <- paste0('macs2 callpeak -t %s -c %s -f BAM -g ce -B -n %s -q ', qval, ' 2>&1 | tee %s')
    }
    
    library(dplyr)
    fls <- getFilePath(ids, format = 'bam', mount=TRUE)
    if(is.null(genome)) genome <- rbeads:::ParseName(basename(fls))$Resolution
    
  
    
    getFileInfo <- function (ID, info, format = ".", processing = ".", scale = ".", url = TRUE,
                             eq = FALSE)
    {
        if (eq) {
            R <- "="
        }
        else {
            R <- "REGEXP"
        }
        con <- dbConnect(dbDriver(DRIVER), group = GROUP)
        exp_file <- unlist(dbGetQuery(con, paste("SELECT ", info, " FROM labfiles WHERE ContactExpID ",
                                                 '=', " '", ID, "' AND Filetype_format ", R, " '", format,
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
    
    
    
    if(summedinput) {
        inp <- 'summed_input'
        inputs <- sapply(ids, function(x) {
            crosslink <- JADBtools::getAnno(x, anno = 'Crosslinker', EXTABLE = 'labexperimentview')
            if(grepl('^e', crosslink, ignore.case = TRUE)) {
                paste0(MOUNT, '/Input/SummedInputs/', genome, '/', genome, '_EGS_HiSeq_input_20M.bam')
            }  else {
                paste0(MOUNT, '/Input/SummedInputs/', genome, '/', genome, '_FRM_HiSeq_input_20M.bam')
            }
        })
    } else {
        inp <- lapply(ids, function(x) {
            
            extract <- JADBtools::getAnno(x, anno = 'ExtractID', EXTABLE = 'labexperiment')
            
            con <- dbConnect(dbDriver(DRIVER), group = GROUP, default.file='~/.my.cnf')
            all_experiments <- tbl(con, "labexperimentview")
            inputs <- all_experiments %>% filter(ExtractID=='mE20', Factor=='Input') %>% collect %>% pull(ContactExpID)
            dbDisconnect(con)
            
            
            
            if (length(inputs) > 1) {
                names(which.max(sapply(inputs, getAlignedReads)))
            } else if (length(inputs) == 1) {
                inputs
            } else {
                stop('No matching inputs found for: ', extract)
            }
        })
        unlist(inp) %>% sapply(getFilePath, format = 'bam', url = local)  -> inputs
    }
    
    
    qq <- paste0('q', format(qval, scientific = TRUE))
    outnames <- sapply(basename(fls), rbeads:::reName, proccesing = 'PeakCalls', scale = 'MACS', resolution = qq, ext='.narowPeak')
    outnames_summits <- sapply(basename(fls), rbeads:::reName, proccesing = 'summits', scale = 'MACS', resolution = qq, ext='.bed')
    prefix <- basename(fls) %>% substr(start=0, stop=nchar(.)-13)
    
    prefix <- paste0(prefix, '_llocal1000')
    
    cmd <- sprintf(
        command,
        (fls),
        (inputs),
        prefix,
        paste0(prefix, '_log.txt')
    )
    message(cmd)
    system(cmd)
    return(c(paste0(prefix, '_peaks.narrowPeak'), paste0(prefix, '_summits.bed')))
    
    
    
    
    
  
        # download.file(fls, basename(fls))
        # if(summedinput) {
        #     dwnlme <- !file.exists(basename(unique(inputs)))
        #     download.file(
        #         file.path('http://ja-db.gurdon.private.cam.ac.uk/db4/files', unique(inputs)[dwnlme]),
        #         basename(unique(inputs)[dwnlme])
        #     )
        # } else {
        #     download.file(inputs, basename(inputs))
        # }
   
    
}

idr <- function(x) {
    peaks <- dir(pattern='narrowPeak')
    ll <- strsplit(dir(pattern='narrowPeak'), '_|\\^|\\.')
    repp <- !ll[[1]] == ll[[2]]
    construct <- ll[[1]]
    construct[repp] <- gsub(' ', '|', paste(ll[[1]][repp], ll[[2]][repp]))
    out <- paste0(paste0(construct[-length(construct)], collapse='_'), '.', construct[length(construct)])
    
    cmd <- sprintf('idr --plot -i 0.05 -o "%s" -s %s %s', out, peaks[1], peaks[2])
    system(cmd)
}
