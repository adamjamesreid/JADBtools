jacl_send_to_cluster('CG133' , genome='cb3ce11', basespace_addr='https://docs.google.com/spreadsheets/d/1fEoOtEq8BEi5AyWujD3GKhnarXmOp7E17PoeW7B6wjQ/edit?usp=sharing', remote = 'jarun@cb-head2.gurdon.private.cam.ac.uk')

jacl_send_to_cluster('CG133' , genome='ce10', basespace_addr='https://docs.google.com/spreadsheets/d/1fEoOtEq8BEi5AyWujD3GKhnarXmOp7E17PoeW7B6wjQ/edit?usp=sharing', remote = 'jarun@cb-head2.gurdon.private.cam.ac.uk')


csenge <- 'https://docs.google.com/spreadsheets/d/1BTNlUPl007Y_SpFjEGa0BkO0_BcRky7hAiiluuanev0/edit?usp=sharing'
jacl_send_to_cluster('TG026' , genome='cb3ce11', basespace_addr=tess, remote = 'jarun@cb-head2.gurdon.private.cam.ac.uk')


jacl_send_to_cluster('CG133' , genome='cb3ce11', basespace_addr=csenge, remote = 'jarun@cb-head2.gurdon.private.cam.ac.uk')

jacl_send_to_cluster('CG133' , genome='ce10', basespace_addr=csenge, remote = 'jarun@cb-head2.gurdon.private.cam.ac.uk')


sapply(log[!log %in% jar_squeue()$NAME], jacl_send_to_cluster, genome='ce10', remote = 'jarun@cb-head2.gurdon.private.cam.ac.uk')


jar_log()



####### RAN seq
require(GenomicAlignments)
require(Rsamtools)
id <- c("rYD077","rYD076","rYD101", "rYD078","rYD079","rYD102")
bam <- getFilePath(id, format = 'bam', processing = 'aligned', mount = TRUE)
bam <- BamFileList(unlist(bam))

data("gnmodel")
require(DESeq2)

SE <- summarizeOverlaps(gnmodel, bam, ignore.strand = TRUE)
colData(SE)$strain <- as.factor(unlist(getStrain(id)))
colnames(SE) <- id
dds <- DESeqDataSet(SE, ~strain)
dds <- DESeq(dds)

rl <- rlog(dds)
rownames(colData(rlog)) <- id
z <- plotPCA(rl, intgroup='strain', ntop = 200000)
require(ggrepel)

sub <- 
'The count data are transforemd to the log2 scale in a way which minimizes differences between samples
for rows with small counts, and which normalizes with respect to library size. The rlog transformation 
produces a similar variance stabilizing effect as varianceStabilizingTransformation,though rlog is more 
robust in the case when the size factors vary widely. The transformation is useful when checking for outliers.'
z +  ggtitle('PCA: rlog transformed counts on BWA aligend data', subtitle = sub) + geom_label_repel(aes(label=name)) + theme_classic(base_size = 12)


mySVD <- function(SE, type='center', sf=TRUE, hmap=FALSE , main='') {
    
    
    Strain <- colData(SE)$strain
    dds <- DESeqDataSet(SE, ~strain)
    dds <- estimateSizeFactors(dds)
    
    m <- counts(dds, normalized=sf)
    
    if(type=='center') {
        X <-  t( apply(m, 1, function(V){V-mean(V)} ) )
    } else if(type=='tags') {
        X <- m
    } else if(type=='rlog') {
        X <- assay(rlog(dds))
        #X <-  t( apply(X, 1, function(V){V-mean(V)} ) )
    }
    
    
    svd  =  svd(X)
    V  = svd$v
    S  = matrix(0, ncol(m), ncol(m))
    for (i in 1:ncol(m)) { S [i,i] =  svd$d[i] }
    
    pvc  <- round((svd$d^2/sum(svd$d^2))*100, 2)
    
    if(hmap) {
        image( 
            1:ncol(m), 1:ncol(m), S%*%t(V),
            main="Singular Value Decomposition", 
            #ylab = sprintf("Samples 1 - %i", ncol(m)), 
            ylab = "",
            xlab = sprintf("%i components", ncol(m)),
            yaxt="n"
        )
        axis(2, 1:ncol(m),labels=colnames(m), las=2)
    }
    #PCA plot
    #par( mfrow = c(1,2) )
    
    
    set.seed(1)
    '%p%' <- paste0
    require(ggrepel)
    require(ggplot2)
    ggplot(as.data.frame(V), aes(x=V1, y=V2, color=Strain, fill = Strain)) + 
        geom_point() +
        geom_label_repel(
            aes(label=colnames(X)),
            fontface = 'bold', color = 'white',
            box.padding = 0.35, point.padding = .5,
            segment.color = 'grey50'
        ) + theme_classic(base_size = 16) +
        xlab('PC1 [' %p% pvc[1] %p% '% variance]') + ylab('PC2 [' %p% pvc[2] %p% '% variance]') +
        ggtitle(main)
    
}
mySVD(SE)




results(dds, contrast = c('strain', 'cfp1', 'N2')) %>% as.data.frame() %>% 
    rownames_to_column('wb') %>% tbl_df() %>% arrange(padj) %>% filter(padj < 0.05) %>% .$wb -> r1

results(dds2, contrast = c('strain', 'cfp1', 'N2')) %>% as.data.frame() %>% 
    rownames_to_column('wb') %>% tbl_df() %>% arrange(padj) %>% filter(padj < 0.05) %>% .$wb -> r2

venn(list(OUR=r1, FP=r2))

id <- make_jadb_ids(1:12, 'rFP')
bam <- getFilePath(id, format = 'bam', processing = 'aligned', mount = TRUE)
bam <- BamFileList(unlist(bam))

SE2 <- summarizeOverlaps(gnmodel, bam, ignore.strand = TRUE)
colData(SE2)$strain <- as.factor(unlist(getStrain(id)))
colnames(SE2) <- id
#dds2 <- DESeqDataSet(SE2[,c(1:3, 7:9)], ~strain)
dds2 <- DESeqDataSet(SE2, ~strain)
dds2 <- DESeq(dds2)

rl2 <- rlog(dds2)
z <- plotPCA(rl2, intgroup='strain', ntop = 200000)
z +  ggtitle('PCA: rlog transformed counts on BWA aligend data', subtitle = sub) + geom_label_repel(aes(label=name)) + theme_classic(base_size = 12)

mySVD(SE2, sf=FALSE, type='rlog')



SEall <- cbind(SE, SE2)
colData(SEall)$source <- factor(gsub('[0-9]', '', rownames(colData(SEall))))

############

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


con <- dbConnect(dbDriver("MySQL"), group = "jadb", default.file='~/.my.cnf')
message(paste0('    ', names(dbGetInfo(con))[1:3], ': ', dbGetInfo(con)[1:3], collapse = '\n'))
dbListTables(con)

dbReadTable(con, "labfiles") %>% as_tibble() -> files

dbReadTable(con, "labexperiment") %>% as_tibble() -> chip
dbDisconnect(con)

efiles <- left_join(chip, files, by = 'ContactExpID')

efiles$ContactExpID %>% table -> zz
zz %>% table

zz2 <- zz[zz==15]

getFilePath('AA654', mount = TRUE) %>% file.exists() %>% '!'(.) %>%  getFilePath('AA654', mount = TRUE)[.]

getDup <- function(x) {
    getFilePath(x, mount = TRUE) %>% sub('\\^.[A-Za-z0-9.]+bw$', '', .) %>% duplicated() %>% getFilePath(x, mount = TRUE)[.]
}


getDup('AA659') %>% gsub('.+\\^|.bw', '', .) %>% sapply(jadb_rm_file)


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
        "%s | sbatch --job-name=%s --output=%s/%s.%s --ntasks-per-node=8", # --nice=10", #--exclude=node9
        paste0(cmd_lst, collapse = '\n'), ID, 
        file.path(MOUNT, '_log'), 
        ID, out_sufix
    )
    message(cmd)
    system(cmd)

    
    
    
    
}

sapply(names(zz2), sc2l, genome='cb3ce11')


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


