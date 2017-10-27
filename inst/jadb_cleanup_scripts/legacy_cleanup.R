jacl_send_to_cluster('CG133' , genome='cb3ce11', basespace_addr='https://docs.google.com/spreadsheets/d/1fEoOtEq8BEi5AyWujD3GKhnarXmOp7E17PoeW7B6wjQ/edit?usp=sharing', remote = 'jarun@cb-head2.gurdon.private.cam.ac.uk')

jacl_send_to_cluster('CG133' , genome='ce10', basespace_addr='https://docs.google.com/spreadsheets/d/1fEoOtEq8BEi5AyWujD3GKhnarXmOp7E17PoeW7B6wjQ/edit?usp=sharing', remote = 'jarun@cb-head2.gurdon.private.cam.ac.uk')


csenge <- 'https://docs.google.com/spreadsheets/d/1BTNlUPl007Y_SpFjEGa0BkO0_BcRky7hAiiluuanev0/edit?usp=sharing'
jacl_send_to_cluster('TG026' , genome='cb3ce11', basespace_addr=tess, remote = 'jarun@cb-head2.gurdon.private.cam.ac.uk')


jacl_send_to_cluster('CG133' , genome='cb3ce11', basespace_addr=csenge, remote = 'jarun@cb-head2.gurdon.private.cam.ac.uk')

jacl_send_to_cluster('CG133' , genome='ce10', basespace_addr=csenge, remote = 'jarun@cb-head2.gurdon.private.cam.ac.uk')


sapply(log[!log %in% jar_squeue()$NAME], jacl_send_to_cluster, genome='ce10', remote = 'jarun@cb-head2.gurdon.private.cam.ac.uk')


jar_log()



#' # jacl_send_to_cluster("AA691", basespace_addr="https://docs.google.com/spreadsheets/d/1QpWQxl3WDL1hRHfJLOlql5qsGWPyFTosCBmJQWTQiQU/edit?usp=sharing")
scl <- function(ID, genome='ce10', ops=', purge=TRUE, backup=TRUE', out_sufix=NULL, pipeline='jadb_ChIPseq') {
    message(ID)
    
    if(is.null(out_sufix)) {
        out_sufix <- pipeline
    }
    
    cmd_lst <- c(
        "echo '#!/usr/bin/Rscript",
        'Sys.info();',
        if(genome=='ce10') 'Sys.setenv(JADB_GROUP=\"ja-db\")' else '',
        if(genome=='ce10') 'Sys.setenv(JADB_MOUNT=\"/mnt/jadb2/DBfile/DBfiles\")' else '',
        'library(JADBtools);',
        'logdir <- file.path(MOUNT, \"_log\");',
        
        sprintf('%s(\"%s\", genome=\"%s\"%s);', pipeline, ID, genome, ops),
        
        'setwd(logdir);',
        sprintf('file.rename(\"%s.%s\", \"done/%s.%s\");', ID, out_sufix, ID, out_sufix),
        "'"
    )
    
    cmd <- sprintf(
        "%s | sbatch --job-name=%s --output=%s/%s.%s --ntasks-per-node=8 --nice=1000", #--exclude=node9
        paste0(cmd_lst, collapse = '\n'), ID, 
        file.path(MOUNT, '_log'), 
        ID, out_sufix
    )
    message(cmd)
    system(cmd)
    #z <- ssh.utils::run.remote(cmd, remote='localhost', verbose = TRUE)
    #message(z$cmd.out)
    
    
    
    
}

sc2l <- function(ID, genome='ce11', ops=', steps = \"map10_uniq\", purge = FALSE, backup = FALSE', out_sufix='uniq_beads', pipeline='jadb_ChIPseq') {
    message(ID)
    
    if(is.null(out_sufix)) {
        out_sufix <- pipeline
    }
    
    
    
    cmd_lst <- c(
        "echo '#!/usr/bin/Rscript",
        'Sys.info();',
        'library(JADBtools);',
        'logdir <- file.path(MOUNT, \"_log\");',
        
        sprintf('%s(\"%s\", genome=\"%s\"%s);', pipeline, ID, genome, ops),
        
        'setwd(logdir);',
        sprintf('file.rename(\"%s/%s.%s\", \"%s/done/%s.%s\");', file.path(MOUNT, '_log'), ID, out_sufix, file.path(MOUNT, '_log'), ID, out_sufix),
        "'"
    )
    
    cmd <- sprintf(
        "%s | sbatch --job-name=%s --output=%s/%s.%s --ntasks-per-node=8 --nice=10", #--exclude=node9
        paste0(cmd_lst, collapse = '\n'), ID, 
        file.path(MOUNT, '_log'), 
        ID, out_sufix
    )
    message(cmd)
    system(cmd)

    
    
    
    
}




log <- jar_nlog()
log[!log %in% jar_squeue()$NAME]


set_ce10 <- function(x) {
    detach("package:JADBtools", unload=TRUE)
    Sys.setenv(JADB_GROUP="ja-db")
    Sys.setenv(JADB_MOUNT="/mnt/jadb2/DBfile/DBfiles")
    require(JADBtools)
}

set_ce11 <- function(x) {
    detach("package:JADBtools", unload=TRUE)
    Sys.setenv(JADB_GROUP="jadb")
    Sys.setenv(JADB_MOUNT="/mnt/jadb/DBfile/DBfiles")
    require(JADBtools)
}

ja_switch <- function() {
    if(GROUP=="jadb") {
        GROUP <<- "ja-db"
        MOUNT <<- "/mnt/jadb2/DBfile/DBfiles"
    } else {
        GROUP <<- "jadb"
        MOUNT <<- "/mnt/jadb/DBfile/DBfiles"
    }
    message('Switched to: [', GROUP, '] mouned on [', MOUNT, ']')
    testConnection()
}




backupFiles <- function(id) {
    pth <- getFilePath(id, url = FALSE)
    dirname <- gsub('files', MOUNT, unique(dirname(pth)))
    
    bck <- gsub('DBfiles', 'DBfiles_legacy', dirname)
    dir.create(bck, recursive = TRUE)
    file.copy(file.path(dirname, '/'), bck, recursive = TRUE)
}


log <- jar_nlog(db='jadb2')
log[!log %in% jar_squeue()$NAME]

sapply(log[!log %in% jar_squeue()$NAME], jacl_send_to_cluster, genome='ce10', remote = 'jarun@cb-head2.gurdon.private.cam.ac.uk')

detach("package:JADBtools", unload=TRUE)
Sys.setenv(JADB_GROUP="")
Sys.setenv(JADB_MOUNT="")
require(JADBtools)



jadb_dc('https://docs.google.com/spreadsheets/d/1QpWQxl3WDL1hRHfJLOlql5qsGWPyFTosCBmJQWTQiQU/edit?usp=sharing', genome='legacy_only', cherry_pick = c('AA698', 'AA699', 'AA700'))


us <- c('AA688', make_jadb_ids(25:27, 'TG'))
jagui_get_table(tab = 'labexperiment') %>% filter(ContactExpID %in% us) -> tab
dir.create('submitt')
setwd('submitt')


jadb_download_files(us)
#rename files
export.csv(tab, 'submit.csv')

addFilesFromCsv('submit.csv', gsheet = FALSE)


