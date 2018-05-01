#' Run peak calls and optionally adds them to DB
#' 
#' @param IDs Vector of JADB ContactExpIDs
#' @param local should the job be run on local or remote FS
#'   
#' @return List 
#' 
#' @author Przemyslaw Stempor
#' 
#' @family Peaks
#' @export
#' 
#' @examples
#' #callPeaksMACS(IDs)
callPeaksMACS <- function(ids, local=TRUE, extsize=150, summedinput=TRUE, genome='ce11') {
    
    con <- dbConnect(dbDriver(DRIVER), group = GROUP, default.file='~/.my.cnf')
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
    
    ids  %>% sapply(getFilePath, format = 'bam', url = local)  -> fls
    
    if(summedinput) {
        inp <- 'summed_input'
        inputs <- sapply(ids, function(x) {
            crosslink <- JADBtools::getAnno(x, anno = 'Crosslinker', EXTABLE = 'labexperimentview')
            if(grepl('^e', crosslink, ignore.case = TRUE)) {
                paste0('Input/SummedInputs/', genome, '/', genome, '_EGS_HiSeq_input.bam')
            }  else {
                paste0('Input/SummedInputs/', genome, '/', genome, '_FRM_HiSeq_input.bam')
            }
        })
    } else {
        inp <- lapply(ids, function(x) {
            extract <- JADBtools::getAnno(x, anno = 'ExtractID', EXTABLE = 'labexperiment')
            
            if(extract == 'mE12') extract <- 'mE11'
            if(extract == 'ak02') extract <- 'ak01'
            if(extract == 'em01' | extract == 'em02') extract <- 'aa04'
            
            
            inputs <- filter(all_experiments, ExtractID==extract, Factor=='Input')$ContactExpID
            
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
    

    
    outnames <- sapply(basename(fls), rbeads:::reName, proccesing = 'PeakCalls', scale = 'MACS', resolution = 'q01', ext='.narowPeak')
    outnames_summits <- sapply(basename(fls), rbeads:::reName, proccesing = 'summits', scale = 'MACS', resolution = 'q01', ext='.bed')
    prefix <- basename(fls) %>% substr(start=0, stop=nchar(.)-13)
    
    

    
    if(local) {
        download.file(fls, basename(fls))
        if(summedinput) {
            dwnlme <- !file.exists(basename(unique(inputs)))
            download.file(
                file.path('http://ja-db.gurdon.private.cam.ac.uk/db4/files', unique(inputs)[dwnlme]),
                basename(unique(inputs)[dwnlme])
            )
        } else {
            download.file(inputs, basename(inputs))
        }
        
        command <- 'macs2 callpeak -t %s -c %s -f BAM -g ce -n %s -q 0.01 2>&1 | tee %s'
        
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
        
        if(is.numeric(extsize)) {
            message('Extsize set to ', extsize)
            command <- paste0(
                'macs2 callpeak -t %s -c %s --nomodel --extsize ',
                extsize, ' -f BAM -g ce -n %s -q 0.01 2>&1 | tee %s'
            )
        } else {
            command <- 'macs2 callpeak -t %s -c %s -f BAM -g ce -n %s -q 0.01 2>&1 | tee %s'
        }
        
        cmd <- sprintf(
            command,
            gsub('files', MOUNT, fls),
            if(summedinput) 
                file.path(MOUNT, inputs) 
            else
                gsub('files', MOUNT, inputs),
            basename(fls) %>% substr(start=0, stop=nchar(.)-13),
            paste0(basename(fls) %>% substr(start=0, stop=nchar(.)-13), '_log.txt')
        )
        
        message(cmd)
        system(cmd)
        
        if(!file.exists(paste0(prefix, '_peaks.narrowPeak'))) stop('No peaks created, see log for details')
        if(!length(readLines(paste0(prefix, '_peaks.narrowPeak')))) stop('No peaks created, [0 peaks in file], see log for details')
        
        peaks_created <- import(
            paste0(prefix, '_peaks.narrowPeak'), format = "BED", 
            extraCols = c(singnalValue = "numeric", pValue = "numeric", qValue = "numeric", peak = "integer")
        )
        npeaks <- length(peaks_created)
        medpeak <- median(width(peaks_created))
        
        peakEntry <- addGenericFile(
            ids, 
            path = file.path('files', exp_dir, gsub('aligned\\^[A-Za-z0-9]+\\^[A-Za-z0-9]+', 'PeakCalls^MACS^q01', prefix)), 
            Processing = 'PeakCalls', 
            Scale = 'MACS', 
            Resolution = 'q01', 
            filetype_format = 'narrowPeak', 
            prefix = 'M',
            genome = genome,
            comments=paste('Npeaks: ', npeaks, '; MEDwidth: ', medpeak, '; Input: ', inp)
        )
        
        summitEntry <- addGenericFile(
            ids, 
            path = file.path('files', exp_dir, gsub('aligned\\^[A-Za-z0-9]+\\^[A-Za-z0-9]+', 'summits^MACS^q01', prefix)), 
            Processing = 'summits', 
            Scale = 'MACS', 
            Resolution = 'q01', 
            filetype_format = 'bed', 
            prefix = 'M',
            genome = genome,
            comments=paste('Npeaks: ', npeaks, '; MEDwidth: ', medpeak, '; Input: ', inp)
        )
        
        file.rename(paste0(prefix, '_peaks.narrowPeak'),  basename(peakEntry$path))
        file.rename(paste0(prefix, '_summits.bed'), basename(summitEntry$path))
        
        setwd(base_dir)
    }

}

#' Combines peak in replicates
#' 
#' @param IDs Vector of JADB ContactExpIDs
#' @param mode u for union or i for intersection
#' @return GRanges 
#' 
#' @author Przemyslaw Stempor
#' 
#' @family Peaks
#' @export
#' 
#' @examples
#' #combinePeaks(IDs)
combinePeaks <- function(ids, mode='u') {
    files <- sapply(ids, getFilePath, format = 'narrowPeak', processing = 'peakCalls')
    peaks <- lapply(files, ChIPseeker::readPeakFile, header=F)
    if(mode=='u') {
        Reduce(GenomicRanges::union, peaks) 
    } else if (mode=='i') {
        Reduce(GenomicRanges::intersect, peaks)
    }
}

#' Combines peak in replicates and exports to bed file
#' 
#' @param IDs Vector of JADB ContactExpIDs
#' @param mode union or intersection
#' @return bed file path 
#' 
#' @author Przemyslaw Stempor
#' 
#' @family Peaks
#' @export
#' 
#' @examples
#' #combinePeaksToBed(IDs)
combinePeaksToBed <- function(ids, mode='union') {
    files <- sapply(ids, getFilePath, format = 'narrowPeak', processing = 'peakCalls')
    peaks <- lapply(files, ChIPseeker::readPeakFile, header=F)
    if(mode=='union') {
        out <- Reduce(GenomicRanges::union, peaks) 
    } else if (mode=='intersection') {
        out <- Reduce(GenomicRanges::intersect, peaks)
    }
    anno <- as.data.frame(t(sapply(basename(files), rbeads:::ParseName)))
    same <- anno[1,c('Factor', 'Strain', 'Stage', 'Processing', 'Scale', 'Resolution')]
    #outname <- paste0(paste0(unlist(same), collapse='_'), '_', paste(anno$ContactExpID, collapse = '^'), '_', mode, '.bed')
    outname <- paste0(paste0(unlist(same), collapse='_'), '_', paste(anno$ExtractID, collapse = '|'), '_', paste(anno$ContactExpID, collapse = '|'), '_', mode, '.bed')
    
    export.bed(out, outname)
    return(outname)
}


#' Combines peak in replicates and exports to bed file
#' 
#' @param IDs Vector of JADB ContactExpIDs
#' @param mode union or intersection
#' @return bed file path 
#' 
#' @author Przemyslaw Stempor
#' 
#' @family Peaks
#' @export
#' 
#' @examples
#' #combinePeaksToBed(IDs)
combinePeaksIDR <- function(ids, mode='union') {
    files <- sapply(ids, getFilePath, format = 'narrowPeak', processing = 'peakCalls', mount=TRUE)
    anno <- as.data.frame(t(sapply(basename(files), rbeads:::ParseName)))
    same <- anno[1,c('Factor', 'Strain', 'Stage', 'Processing', 'Scale', 'Resolution')]
    outname <- paste0(paste0(unlist(same), collapse='_'), '_', paste(anno$ContactExpID, collapse = '^'), '_', 'IDR2', '.narrowPeak')
    cmd <- paste(
        '~/miniconda2/envs/py3/bin/idr --samples', files[1], files[2], 
        '--output-file', outname, '--output-file-type narrowPeak'
    )
    message(cmd)
    system(cmd)
    
    return(outname)
}

#' Calls enriched regions from BigWig file
#' 
#' @param file path/url to BW
#' @param bedoutput path to output BED
#' @return GRanges or bed file path 
#' 
#' @author Przemyslaw Stempor
#' 
#' @family Peaks
#' @export
#' 
#' @examples
#' #combinePeaksToBed(IDs)
enrichedRegionsCall <- function(file, bedoutput=NULL) {

        bw <- BigWigFile(file)
        REF <- seqinfo(bw)
        covtrack <- import.bw(bw, as='RleList')
        a <- quantile(unlist(covtrack), 0.75)
        
        message("INFO: a = ", a)

        enriched_regions <- IRanges::slice(covtrack, lower = a)
        peakSumsRep1 <- viewSums(enriched_regions)
        enriched_regions <- RangedData(as(enriched_regions[peakSumsRep1 >= 
                                                               quantile(unlist(peakSumsRep1), 0.9)], "IRangesList"))
        enriched_regions <- as(enriched_regions, "GRanges")
        seqinfo(enriched_regions) <- REF[seqlevels(enriched_regions)]
        
        if(is.null(bedoutput)) {
            return(enriched_regions)
        } else {
            export.bed(enriched_regions, bedoutput)
            return(bedoutput)
        }
        
}

