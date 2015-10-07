getBWarray <- function(IDs, res=1000L, processing = 'aligned') {
    
    con <- dbConnect(dbDriver("MySQL"), group = "jadb")
    dbListTables(con)
    all <- dbReadTable(con, "labexperimentview")
    
    require(dplyr)
    all  %>% filter(Factor == 'HPL2', Strain == 'N2')  %>% select(ContactExpID)  %>% unlist -> IDs
    names(IDs) <- IDs
    
    IDs %>% sapply(getFilePath, format = 'bw', processing=processing) -> paths
    require(rtracklayer)
    bwfl <- BigWigFileList(unlist(paths))
    
    lst <- lapply(bwfl, function(bwf) extarct_vector(bwf, size = seqlengths(bwf)/res) )
    M <- do.call(cbind, lst)
    C <- cor(M)
    
    library(d3heatmap)
    d3heatmap(C, 
              colors = colorFunc, 
              symm=TRUE, 
              scale = 'none', 
              theme='dark', 
              revC = TRUE,
              zlim=c(0,1)
              )
    #heatmap.2(C,col = rev(RColorBrewer::brewer.pal(11, 'RdYlBu')), scale = 'none')
    plot.new()
    fields::image.plot(C, smallplot= c(0,1,.5,1), 
                       legend.only = TRUE, horizontal=TRUE, col = rev(RColorBrewer::brewer.pal(11, 'Spectral')))
        
}

extarct_vector <- function(track, size, which = as(seqinfo(track), "GenomicRanges")) {
    sum <- .Call(
        'BWGFile_summary', path.expand(path(track)),
        as.character(seqnames(which)), ranges(which), 
        S4Vectors::recycleIntegerArg(size, "size", length(which)), 
        "mean", as.numeric(NA_real_), PACKAGE='rtracklayer'
    )
    return(unlist(sum))
#     M <- do.call( rbind, sum )
#     if (!ignore_strand) 
#         M[as.character(strand(gr))=='-',] <- M[
#             as.character(strand(gr))=='-', ncol(M):1]
#     return(M)
}