spike_get_coef <- function(experiment, input) {
    what <- c("rname")
    flag <- scanBamFlag(isUnmappedQuery = FALSE)
    param <- ScanBamParam(what = what, flag = flag)
    
    inp <- gsub('^files', MOUNT, getFilePath(input, url=FALSE, processing = 'aligned', format = 'bam'))
    message('File: ', basename(inp))
    
    a <- scanBam(inp, param = param)[[1]]
    cb3i <- grepl('^cb3_', a[[1]])
    ei <- as.numeric(sum(!cb3i))
    bi <- as.numeric(sum(cb3i))
    message('#CeInput=', ei, '; #CbInput=', bi, '; ratio=', ei/bi)
    
    
    ex <- gsub('^files', MOUNT, getFilePath(experiment, url=FALSE, processing = 'aligned', format = 'bam'))
    message('File: ', basename(ex))
    
    e <- scanBam(ex, param = param)[[1]]
    cb3e <- grepl('^cb3_', e[[1]])
    ee <- as.numeric(sum(!cb3e))
    be <- as.numeric(sum(cb3e))
    message('#CeExp=', ee, '; #CbExp=', be, '; ratio=', ee/be)
    
    coef <- (be/ee)*(ei/bi)
    return(coef)
    
}



