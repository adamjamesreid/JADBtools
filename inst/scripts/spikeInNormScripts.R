
getRatio('RC019') #input N2
getRatio('RC025') #input cfp1

getRatio('RC024') #H3 N2
getRatio('RC030') #H3 cfp1

##Normalized on H3
norm_n2 <- sapply( paste0('RC0', 19:24), getRatio) / getRatio('RC024')
norm_cfp1 <- sapply( paste0('RC0', 25:30), getRatio) / getRatio('RC030')

all <- c(norm_n2[2:5], norm_cfp1[2:5])
for (f in 1:length(all)) {
    message(names(all[f]))
    normTrack(all[f])
}
normTrack(norm_n2[5])
normTrack(norm_cfp1[5])

##Normalized on input
inp_n2 <- sapply( paste0('RC0', 19:24), getRatio) / getRatio('RC019')
inp_cfp1 <- sapply( paste0('RC0', 25:30), getRatio) / getRatio('RC025')



exp_file     <- unlist( dbGetQuery(con, paste("SELECT path FROM labfiles WHERE ContactExpID = '", ID, "' AND Filetype_format='bam'; ", collapse="", sep="") ) )
M <- dbReadTable(con, 'labfiles')
N <- dbReadTable(con, 'labexperiment')

fls <- file.path("http://jadb.gurdon.private.cam.ac.uk/db4", grep('txt.gz$', M[M$ContactExpID %in% IDs,]$path, val=TRUE))
sapply(fls, function(f) download.file(f, basename(f)))
