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
jacl_send_to_cluster <- function(ID, genome='ce11') {
    message(ID)
    
    cmd_lst <- c(
        "echo '#!/usr/bin/Rscript",
        "logdir <- getwd()",
        "Sys.info();",
        "library(JADBtools);",
        sprintf("jadb_ChIPseq(\"%s\", genome=\"%s\");", ID, genome),
        "setwd(logdir)",
        sprintf("file.rename(\"%s.out\", \"done/%s.out\");", ID, ID),
        "'"
    )
    
    cmd <- sprintf(
        "%s | sbatch --job-name=%s --output=%s.out --ntasks-per-node=8", #--exclude=node9
        paste0(cmd_lst, collapse = '\n'), ID, ID
    )
    system(cmd)
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
    
    to_proc <- all_exp_to_proc %>% filter(genome!='ce11' | is.na(genome), !ContactExpID %in% has_no_fq, !ContactExpID %in% under_processing) %>% .$ContactExpID
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
    
    cmd_lst <- c(
        "echo '#!/usr/bin/Rscript",
        "logdir <- getwd()",
        "Sys.info();",
        "Sys.setenv(JADB_GROUP=\"ja-db\")",
        "Sys.setenv(JADB_MOUNT=\"/mnt/jadb2/DBfile/DBfiles\")",
        "library(JADBtools);",
        sprintf("jadb_ChIPseq(\"%s\", genome = \"ce10\");", ID),
        "setwd(logdir)",
        sprintf("file.rename(\"%s.out\", \"done/%s.out\");", ID, ID),
        "'"
    )
 
    
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