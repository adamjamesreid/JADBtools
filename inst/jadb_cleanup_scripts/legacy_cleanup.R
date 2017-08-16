jacl_send_to_cluster('CG133' , genome='cb3ce11', basespace_addr='https://docs.google.com/spreadsheets/d/1fEoOtEq8BEi5AyWujD3GKhnarXmOp7E17PoeW7B6wjQ/edit?usp=sharing', remote = 'jarun@cb-head2.gurdon.private.cam.ac.uk')

jacl_send_to_cluster('CG133' , genome='ce10', basespace_addr='https://docs.google.com/spreadsheets/d/1fEoOtEq8BEi5AyWujD3GKhnarXmOp7E17PoeW7B6wjQ/edit?usp=sharing', remote = 'jarun@cb-head2.gurdon.private.cam.ac.uk')


csenge <- 'https://docs.google.com/spreadsheets/d/1BTNlUPl007Y_SpFjEGa0BkO0_BcRky7hAiiluuanev0/edit?usp=sharing'
jacl_send_to_cluster('TG026' , genome='cb3ce11', basespace_addr=tess, remote = 'jarun@cb-head2.gurdon.private.cam.ac.uk')


jacl_send_to_cluster('CG133' , genome='cb3ce11', basespace_addr=csenge, remote = 'jarun@cb-head2.gurdon.private.cam.ac.uk')

jacl_send_to_cluster('CG133' , genome='ce10', basespace_addr=csenge, remote = 'jarun@cb-head2.gurdon.private.cam.ac.uk')


sapply(log[!log %in% jar_squeue()$NAME], jacl_send_to_cluster, genome='ce10', remote = 'jarun@cb-head2.gurdon.private.cam.ac.uk')


jar_log()




log <- jar_nlog()
log[!log %in% jar_squeue()$NAME]


detach("package:JADBtools", unload=TRUE)
Sys.setenv(JADB_GROUP="ja-db")
Sys.setenv(JADB_MOUNT="/mnt/jadb2/DBfile/DBfiles")
require(JADBtools)


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


