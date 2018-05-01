con <- dbConnect(dbDriver(DRIVER), group = GROUP, default.file='~/.my.cnf')

tab <- dbReadTable(con, 'labfiles') %>%  tbl_df()
pth <- filter(tab, Scale=='log2')$path

pth <- gsub('files/', '', pth)
x=pth[[1]]

require(rtracklayer)

sapply(pth, function(x) {
    message("-> processing: ", which(pth==x), '/2944 => ', x)
    log2sc_gr <- import.bw(x)

    finite_finitos <- is.finite(log2sc_gr$score)
    if( sum(!finite_finitos) > 0 ) {
        log2sc_gr <- log2sc_gr[is.finite(log2sc_gr$score)] 
        export.bw(log2sc_gr, x)
    } else {
        message(x, ' [OK] ')
    }
    message('Mean after log2: ', mean(log2sc_gr$score))
})

ex <- jagui_get_table('labexperimentview')
fls <- jagui_get_table('labfiles')
bt <- left_join(ex, fls, by='ContactExpID')

procme <- bt %>% filter(Processing=='BEADSNQNU', Scale=='linear') %>% pull(ContactExpID)
    

recalc_nqnu <- function(x) {
    uid <- bt %>% filter(Processing=='BEADSNQNU', ContactExpID==x) %>% pull(UID)
    genome <- bt %>% filter(Processing=='BEADSNQNU', ContactExpID==x) %>% pull(genome) %>% unique
    
    sapply(uid, jadb_rm_file)
    jadb_ChIPseq(x,steps=c('map0_norm', "log2_map0_norm"), genome=genome, purge = FALSE)
}
require(pbmcapply)

procme <- bt %>% filter(Processing=='BEADSNQNU', Scale=='linear') %>% filter(Updated.y!='2018-01-14') %>% filter(Updated.y!='2018-01-15') %>% pull(ContactExpID)
res <- pbmclapply(procme, recalc_nqnu, mc.cores = 32)


ex <- jagui_get_table('labexperimentview')
fls <- jagui_get_table('labfiles')
bt <- left_join(ex, fls, by='ContactExpID')

require(pbmcapply)
bfl <- bt %>% filter(Processing=='BEADSNQNU', Scale=='linear') %>% pull(path) %>% sub('^files', MOUNT, .) %>% BigWigFileList
ans2 <- pbmclapply(bfl, JADBtools:::extarct_vector, 1L, GRanges('chrII:7516896-7517981'), mc.cores = detectCores())

fl <- bt %>% filter(Processing=='BEADSNQNU', Scale=='linear') %>% pull(path) %>% sub('^files', MOUNT, .)
ans3 <- pbmclapply(fl, file.exists, mc.cores = detectCores())

fl <- bt %>% pull(path) %>% sub('^files', "http://jadb.gurdon.private.cam.ac.uk/db4", .)
ans4 <- pbmclapply(fl, file.exists, mc.cores = detectCores())

#------#

sapply(paste0('AA6', 52:63), jacl_send_to_cluster, genome='cb3ce11')
sapply(paste0('TG00', 3:8), jacl_send_to_cluster, genome='cb3ce11')


jacl_send_to_cluster('CG073', genome = 'cb3ce11')
jacl_send_to_cluster('CG074', genome = 'cb3ce11')
jacl_send_to_cluster('CG076', genome = 'cb3ce11')
jacl_send_to_cluster('CG077', genome = 'cb3ce11')
jacl_send_to_cluster('CG078', genome = 'cb3ce11')




fls <- jagui_get_table('labfiles')
pth <- filter(fls, Scale=='log2')$path
pth <-  gsub('^files', '/mnt/jadb/DBfile/DBfiles', pth)

fixlog2 <- function(fixme) {
    log2sc <- import.bw(fixme)
    log2sc_gr <- as(log2sc, 'GRanges')
    log2sc_gr <- log2sc_gr[is.finite(log2sc_gr$score)]
    export.bw(log2sc_gr, fixme)
}
require(pbmcapply)
require(rtracklayer)
pbmclapply(pth, fixlog2, mc.cores = 24)

 -> ids
sapply(as.character(b2), jacl_send_to_cluster)

sapply(c('TG010', 'TG012', 'TG015'), jacl_send_to_cluster)

sapply(c("TG016","TG017","TG018", "TG020"), jacl_send_to_cluster)



