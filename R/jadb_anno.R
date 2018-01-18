#' Add annotation to JADB
#'
#' @param file - path to file or GR
#' @param Type Type
#' @param Version Version
#' @param Description Description
#' @param Source Source
#' @param Genome Genome
#' @param user user ID
#'
#' @return NULL
#' @export
#'
anno_add <- function(file, Type='Genes', Version='ensembl90', Description='tets add', Source='', Genome='ce11', user='ps562') {
    
    if (class(file)=="GRanges") {
        tmp <- tempfile()
        rtracklayer::export.bed(file, tmp)
        file <- tmp
    }
    
    #con <- dbConnect(dbDriver("MySQL"), group = GROUP, default.file='~/.my.cnf')
    #labanno <- con %>% tbl('labanno')
    maxid <- jagui_get_table('labanno')$ContactExpID %>% sort %>% tail(1)
    newid <- sprintf('ANN%03i', as.numeric(gsub('ANN', '', maxid))+1)
    
    GR <- import.bed(file)
    n <- length(GR)
    gc <- round(sum(width(reduce(GR, ignore.strand=TRUE)))/sum(width(GRangesForBSGenome('ce11'))), 3)
    rd <- round( sum(sum(coverage(GR)>1))/sum(sum(coverage(GR)>0)), 3)
    
    anno_add_file(newid, file, ref.genome.version=Genome) 
    record <- c(
        ContactEXpID=newid,
        Type=Type,
        Genome=Genome,
        Version=Version,
        Description=Description,
        Source=Source,
        Stats=glue::glue('n={n}; gc={gc}; rd={rd}'),
        userID=user
    )
    res <- jagui_add_recoerd(record, TABLE = 'labanno')
    
}


anno_add_file <- function(ID, file, ref.genome.version='ce11') {
    dir <- file.path(MOUNT, 'ANNO', ID)
    file_base <- formGenericPath('ANN001', 'labanno', Processing = 'source', Resolution = ref.genome.version, Scale = 'ranges')
    db_path <- sub(MOUNT, 'files', file.path(dir, file_base))
    
    db_entry <- addGenericFile(
        ID,
        path = db_path, 
        Processing = 'source', 
        Scale = 'ranges', 
        Resolution = 'ce11',
        filetype_format = 'bed', 
        prefix = 'A',
        comments = '',
        genome =  ref.genome.version
    )
    
    final_path <- sub('files', MOUNT, db_entry$path)
    #remote <- configr::read.config(file='~/.my.cnf')[[GROUP]]$host
    
    if(!grepl('Rtmp', file)) {
        temp <- tempfile()
        file.copy(file, temp)
    } else {
        temp <- file
    }
    
    ssh.utils::mkdir.remote(
        dirname(final_path), 
        remote = "jarun@cb-head2.gurdon.private.cam.ac.uk"
    )
    cp.remote(
        remote.src = "", 
        path.src = temp,
        remote.dest = "jarun@cb-head2.gurdon.private.cam.ac.uk", 
        path.dest = final_path, 
        verbose = TRUE
    )
}


jagui_add_recoerd <- function(RECORD, input=NULL, TABLE, verbose=FALSE, host=NULL) {
    
    nms <- names(RECORD)
    
    if(is.null(input)) {
        values <- RECORD
    } else {
        values <- unlist(lapply(nms, function(x) {
            input[[x]]
        }))
    }
    
    sql_string <- paste0(
        "INSERT INTO ", TABLE, " (", 
        paste0(nms, collapse = ", "), 
        ") VALUES ('", 
        paste0(values, collapse = "', '"),
        "');"
    ) 
    
    if(verbose) message(sql_string)
    
    con <- dbConnect(dbDriver("MySQL"), group = "jadb", default.file='~/.my.cnf', host=host)
    res <- dbSendStatement(con, sql_string)
    out_res <- dbGetInfo(res)
    dbClearResult(res)
    dbDisconnect(con)
    
    return(out_res)
}


