#' Get BigWig tracks from GEO microarray experiment
#' 
#' @param ID list of GSM IDs 
#' @param platform GPL id of platform
#' @param liftover if TRUE lifts over from ce4 tp ce10
#'   
#' @return List of tracks as RleList 
#' 
#' @author Przemyslaw Stempor
#' 
#' @family GEO
#' @export
#' 
#' @examples
#' input <- list(
#'     "H3K36me3_mes-4 RNAi_rep1"="GSM938394",
#'     "H3K36me3_mes-4 RNAi_rep2"="GSM938395",
#'     "H3K27me3_mes-4 RNAi_rep1"="GSM938396",
#'     "H3K27me3_mes-4 RNAi_rep2"="GSM938397",
#'     "H3K36me3_N2_rep1"="GSM938401",	
#'     "H3K36me3_N2_rep2"="GSM938402",
#'     "H3K36me3_N2_rep3"="GSM938403",
#'     "H3K27me3_N2_rep1"="GSM938404",	
#'     "H3K27me3_N2_rep2"="GSM938405"
#' )	
#' 
#' ANS <- GEOarray2bw(input)
#' H3K36me3_mes4 <- (ANS$"H3K36me3_mes-4 RNAi_rep1" + ANS$"H3K36me3_mes-4 RNAi_rep2")/2
#' export.bw(H3K36me3_mes4 , 'H3K36me3_mes4.bw')
#' 
#' H3K27me3_mes4 <- (ANS$"H3K27me3_mes-4 RNAi_rep1" + ANS$"H3K27me3_mes-4 RNAi_rep2")/2
#' export.bw(H3K27me3_mes4, 'H3K27me3_mes4.bw')
#' 
#' H3K36me3_N2 <- (ANS$"H3K36me3_N2_rep1" + ANS$"H3K36me3_N2_rep2" + ANS$"H3K36me3_N2_rep3")/3
#' export.bw(H3K36me3_N2, 'H3K36me3_N2.bw')
#' 
#' H3K27me3_N2    <- (ANS$"H3K27me3_N2_rep1"+ ANS$"H3K27me3_N2_rep2")/2
#' export.bw(H3K27me3_N2, 'H3K27me3_N2.bw')
#' 
GEOarray2bw <- function(ID, platform="GPL8647", liftover=TRUE) {
    require(GEOquery)
    require(magrittr)
    
    gpl <- getGEO(platform)
    
    df <- Table(gpl)
    df$RANGE_START <- as.numeric(df$RANGE_START)
    df$RANGE_END <- as.numeric(df$RANGE_END)
    df$CHROMOSOME <- droplevels(df$CHROMOSOME)
    
    GR <- GRanges(df$CHROMOSOME, IRanges(df$RANGE_START, width = 50), id=df$ID)
    names(GR) <- GR$id
    GenomeInfoDb::seqlevelsStyle(GR) <- 'UCSC'
    
    if(liftover) {
        download.file('http://hgdownload-test.cse.ucsc.edu/goldenPath/ce4/liftOver/ce4ToCe10.over.chain.gz', tmp<-tempfile(fileext = '.chain.gz'))
        gunzip(tmp, tmp2<-tempfile(fileext = '.chain'))
        chain <- import.chain(tmp2)
    }


    geo2bw <- function(id) {
        
        message( id )
        gsm <- getGEO(id)
        message( Meta(gsm)$title )
            
        out <- GR[Table(gsm)$ID_REF]
        out$score <- as.numeric(Table(gsm)$VALUE)
        if(liftover) out <- liftOver(out, chain)
        if(liftover) out <- unlist(out)
        covr <- coverage(out, weight = 'score')
        
        export.bw(covr, paste0(id, '_', gsub(' ', '_', Meta(gsm)$title), '.bw'))
        return(covr)
    }
    
    ID %>% sapply(geo2bw) -> ANS
    return(ANS)
}
# RANGE_END[is.na(RANGE_END)] <- RANGE_START[is.na(RANGE_END)]
# gr <- makeGRangesFromDataFrame(
#     df, seqnames.field = 'CHROMOSOME', start.field = 'RANGE_START', 
#     end.field = "RANGE_START"
# )


GEOseqrch <- function(x) {
    res = queryGEO('sarcoma AND "cel"[Supplementary Files] AND 20:20000[Number of Samples]')
    
    #BEST
    require(rentrez)
    quer <- '(("genome binding/occupancy profiling by high throughput sequencing"[DataSet Type] AND ("2014/01/01"[PDAT] : "3000"[PDAT])) AND ("2014/01/01"[PDAT] : "2014/12/31"[PDAT])) AND "high throughput sequencing"[Platform Technology Type]'
    web_env_search <- entrez_search(db = "gds", term = quer, retmax = 10000, usehistory="y")
    cookie <- web_env_search$WebEnv
    qk <- web_env_search$QueryKey
    
    data_summaries <- entrez_summary(db = "gds", id = web_env_search$ids[1:100])
    data_summaries2 <- entrez_summary(db = "gds", WebEnv = cookie, query_key = qk)
    
    OUT <- sapply(web_env_search$ids, function(x) {
        message(x)
        entrez_summary(db = "gds", id = x)
    })
    
    
    LIST <- lapply(web_env_search$ids, function(x) {
        message(x)
        entrez_summary(db = "gds", id = x)
    })
    
    test <- sapply(LIST, function(x) {
        x <- as.list.List(x)
        el <- elementLengths(x) 
      
        x[sapply(x, class) != "character"] <- el[sapply(x, class) != "character"]
        return(x)
    })
    
    do.call(rbind, LIST)  -> TAB
    
    library ( XLConnect )
    
    wb <- loadWorkbook("GEOgenome_binding_occupancy_profiling_2014.xlsx", create = TRUE)
    createSheet(wb, name = "GSE")
    #createSheet(wb, name = "experiments")
    TAB <- as.data.frame(TAB) 
    
    # Write built-in data set 'CO2' to the worksheet created above;
    # offset from the top left corner and with default header = TRUE
    writeWorksheet(wb, as.data.frame(TAB), sheet = "GSE")
    #writeWorksheet(wb, t(OUT), sheet = "experiments")
    
    
    # Save workbook (this actually writes the file to disk)
    saveWorkbook(wb)
    
    
    wb <- loadWorkbook("GEOgenome_binding_occupancy_profiling_2014_experimentsCounts.xlsx", create = TRUE)
    createSheet(wb, name = "experiments")
    writeWorksheet(wb, test, sheet = "experiments")
    saveWorkbook(wb)
    
    
    
    
    library(dplyr)
    db = src_sqlite('GEOmetadb.sqlite')
    gds = tbl(db,'gds')
    filter(gse,gse=='GSE2553')
    
    library(GEOmetadb)
    #if(!file.exists('GEOmetadb.sqlite')) getSQLiteFile()
    con <- dbConnect(SQLite(), '/Volumes/raid0/_GEO/GEOmetadb.sqlite')
    
    geo_tables <- dbListTables(con)
    geo_tables  %>% sapply(function(x) dbListFields(con,x))
    
    z <-  dbGetQuery(con,'select * from gse where type LIKE "%genome binding/occupancy profiling by high throughput sequencing%" and status LIKE "Public % 2014" limit 100000')
    
    dbGetQuery(con, 'select gsm from gse_gsm where gse = "GSE62833"')
    
    w = dbGetQuery(con, 'select * from gsm where gsm = "GSM1534077"')
    
    }
