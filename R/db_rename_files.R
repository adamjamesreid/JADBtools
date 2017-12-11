jadb_rename_chip <- function(ID) {
    
    con <- dbConnect(dbDriver("MySQL"), group = GROUP, default.file = "~/.my.cnf")
    
    fileds.def <- dbGetQuery(con, sprintf("SHOW FIELDS FROM %s", 'labexperimentview'))
    
    tbl(con, 'labexperimentview') %>% filter(ContactExpID == ID) %>% collect() %>% 
        select(Factor, Antibody, ExtractID, Crosslinker, Strain, Stage, ContactExpID) -> ent
    
    
    tbl(con, 'labfiles') %>% filter(ContactExpID == ID) %>% collect() %>% 
        select(Processing, genome, Scale, ContactExpID, UID, filetype_format, path) %>% left_join(ent, .) -> fls
    
    dbDisconnect(con)
    
    fls %<>% mutate(path=gsub("files", MOUNT, path))
    
    fls %<>% mutate(pat2=sprintf(
        '%s^%s_%s^%s^%s^%s_%s^%s^%s_%s^%s.%s', 
        Factor, Antibody, ExtractID, Crosslinker, Strain, Stage, Processing,  Scale, genome, ContactExpID, UID, filetype_format
    )) 
    
    fls %<>% mutate(pat_refix=sprintf(
        '%s^%s_%s^%s^%s^%s', 
        Factor, Antibody, ExtractID, Crosslinker, Strain, Stage
    )) 
    
    
    fls %<>% mutate(dirname=sprintf(
        '%s/%s/%s/%s_%s_%s', 
        MOUNT, Factor, Strain, ContactExpID, ExtractID, Antibody
    ))
    

    
    #diff <- basename(fls$path) == fls$pat2
    
    dir <- unique( gsub('/meme_chip', '', dirname(fls$path)) )
    dir.create(unique(fls$dirname), recursive = TRUE)
    
    file.rename(
        fls$path, 
        file.path(unique(fls$dirname), fls$pat2)
    )
    
    opref <- gsub('([A-Za-z0-9\\^]_[A-Za-z0-9\\^]+).+', '\\1', basename(fls$path)[[1]])
    
    file.rename(
        dir(dir,full.names=TRUE),
        file.path(
            unique(fls$dirname), sub(opref, unique(fls$pat_refix), dir(dir), fixed = TRUE)
        )
    )
    file.remove(dir)
    
    
    SQL <- sprintf('UPDATE labfiles SET path = "%s" WHERE UID = "%s"', file.path(sub(MOUNT, 'files', fls$dirname), fls$pat2), fls$UID)
    con <- dbConnect(dbDriver("MySQL"), group = GROUP, default.file = "~/.my.cnf")
    sapply(SQL, function(x) dbGetQuery(con, x))
    dbDisconnect(con)
    
    path <- unique(dirname(regFLS))
    
    rbeads:::ParseName(regFLS[[2]])
    rbeads:::
    
    rbeads:::ParseName()
    regFLS <- JADBtools::getFilePath(ID, mount = TRUE)
    
}
    
    
    