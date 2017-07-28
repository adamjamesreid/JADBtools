#' jadb_dc - direct submission to JADB
#'
#' @param gurl - Google Docs URL
#' @param genome - reference genome
#' @param legacy_ce10 - add entries to ce10 database
#'
#' @return NULL
#' @export
#'
jadb_dc <- function(gurl, genome=c('ce11', 'cb3ce11'), legacy_ce10=TRUE, wipeout = FALSE) {
    
    genome <- match.arg(genome)
    message('Processing experiments: ', genome, ' reference version', if(legacy_ce10) ' with legacy ce10 database processing')
    
    
    addExtractIDs(gurl)
    validateFilesFromBaseSpace(gurl)
    
    sapply(ids, jacl_send_to_cluster, genome=genome, basespace_addr=gurl, remote = 'jarun@cb-head2.gurdon.private.cam.ac.uk', wipeout = wipeout)
    
    message('Processing done for ', genome)
    

    if(legacy_ce10) {
        detach("package:JADBtools", unload=TRUE)
        Sys.setenv(JADB_GROUP="ja-db")
        Sys.setenv(JADB_MOUNT="/mnt/jadb2/DBfile/DBfiles")
        require(JADBtools)
        
        addExtractIDs(gurl)
        validateFilesFromBaseSpace(gurl)
        
        
        sapply(ids, jacl_send_to_cluster, genome='ce10', basespace_addr=gurl, remote = 'jarun@cb-head2.gurdon.private.cam.ac.uk', wipeout = wipeout)
        
        
        message('Processing done for ce10')
        
        detach("package:JADBtools", unload=TRUE)
        Sys.setenv(JADB_GROUP="")
        Sys.setenv(JADB_MOUNT="")
        require(JADBtools)
    }
}

#' jacl_send_to_cluster
#' 
#' @param ID
#'   
#' @return NULL
#' 
#' @author Przemyslaw Stempor
#' 
#' @family cluster
#' @export
#' 
#' @examples 
#' # sapply(JADB_mass_parallel(50), jacl_send_to_cluster)
#' 
#' # ids <- strsplit('AA652 AA653 AA654 AA655 AA656 AA657 AA658 AA659 AA660 AA661 AA662 AA663', ' ')[[1]]
#' # sapply(ids, jacl_send_to_cluster)
#' 
#' # ids <- jacl_mass_parallel(50)
#' # sapply(ids, jacl_send_to_cluster)
#' 
#' # jacl_send_to_cluster('AA007', ops=", steps=c(\"log2_norm\", \"map0_norm\", \"log2_map0_norm\"), purge=FALSE")
#' 
#' 
#' # jadb_renove_exp("AA691")
#' # jacl_send_to_cluster("AA691", basespace_addr="https://docs.google.com/spreadsheets/d/1QpWQxl3WDL1hRHfJLOlql5qsGWPyFTosCBmJQWTQiQU/edit?usp=sharing")
jacl_send_to_cluster <- function(ID, genome='ce11', ops='', out_sufix='chip', pipeline='jadb_ChIPseq', basespace_addr='', remote='', wipeout=FALSE) {
    message(ID)
    
    cmd_lst <- c(
        "echo '#!/usr/bin/Rscript",
        'Sys.info();',
        if(genome=='ce10') 'Sys.setenv(JADB_GROUP=\\"ja-db\\")' else '',
        if(genome=='ce10') 'Sys.setenv(JADB_MOUNT=\\"/mnt/jadb2/DBfile/DBfiles\\")' else '',
        'library(JADBtools);',
        'logdir <- file.path(MOUNT, \\"_log\\");',
        
        if(wipeout) sprintf('JADBtools:::jadb_renove_exp(\\"%s\\");', ID) else '',
        if(nchar(basespace_addr)) sprintf('jadb_basespace(\\"%s\\", select_id=\\"%s\\");', basespace_addr, ID) else '',
        sprintf('%s(\\"%s\\", genome=\\"%s\\"%s);', pipeline, ID, genome, ops),
        
        'setwd(logdir);',
        sprintf('file.rename(\\"%s.%s\\", \\"done/%s.%s\\");', ID, out_sufix, ID, out_sufix),
        "'"
    )
    
    

    
    
    cmd <- sprintf(
        "%s | sbatch --job-name=%s --output=%s/%s.%s --ntasks-per-node=8", #--exclude=node9
        paste0(cmd_lst, collapse = '\n'), ID, file.path(MOUNT, '_log'), ID, out_sufix
    )
    z <- ssh.utils::run.remote(cmd, remote, verbose = TRUE)
    message(z$cmd.out)
    
    
    
    
}

#' jacl_mass_parallel
#' 
#' @param ContactExpID vector of IDs
#'   
#' @return eset
#' 
#' @author Przemyslaw Stempor
#' 
#' @family cluster
#' @export
#' 
jacl_mass_parallel <- function(n=50) {
    
    library(dplyr)
    library(tibble)
    library(RMySQL)
    
    con <- dbConnect(dbDriver("MySQL"), group = "jadb", default.file='~/.my.cnf')
    message(paste0('    ', names(dbGetInfo(con))[1:3], ': ', dbGetInfo(con)[1:3], collapse = '\n'))
    dbListTables(con)
    
    dbReadTable(con, "labfiles") %>% as_tibble() -> files
    
    dbReadTable(con, "labexperiment") %>% as_tibble() -> chip
    dbDisconnect(con)
    
    efiles <- left_join(chip, files, by = 'ContactExpID')
    efiles %>% filter(Processing=='raw', filetype_format!='Rdata') %>%  select(ContactExpID, path) -> raw_paths
    
    efiles %>% filter(filetype_format=='bam') -> all_exp_to_proc
    
    has_no_fq <- chip$ContactExpID[!chip$ContactExpID %in%  raw_paths$ContactExpID]
    
    tbl <- readr::read_table(pipe('squeue'))
    under_processing <- tbl$NAME
    
    to_proc <- all_exp_to_proc %>% filter(genome!='ce11' | genome!='cb3ce11' | is.na(genome), !ContactExpID %in% has_no_fq, !ContactExpID %in% under_processing) %>% .$ContactExpID
    head(to_proc, n)
}


#' jacl_send_to_cluster_ce10 for ce10
#' 
#' @param ID
#'   
#' @return NULL
#' 
#' @author Przemyslaw Stempor
#' 
#' @family cluster
#' @export
#' 
#' @examples 
#' # sapply(ids, jacl_send_to_cluster_ce10)
jacl_send_to_cluster_ce10 <- function(ID) {
    message(ID)
    
 
 
    
    cmd <- sprintf(
        "%s | sbatch --job-name=%s --output=%s.out --ntasks-per-node=8", #--exclude=node9
        paste0(cmd_lst, collapse = '\n'), ID, ID
    )
    system(cmd)
}


# set nice
# for i in $(squeue -u jarun -h -t PD -o %i); do; scontrol update jobid=$i nice=1000; done;
# for i in $(squeue -u jarun -h -t R -o %i); do; scancel $i; done;


# for i in $(squeue -u jarun -h -t R -o %i); do
#    scancel $i
# done