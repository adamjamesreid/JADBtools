
urm <- unlistedRepeatModel
mcols(urm) <- NULL
urm <- as(urm, 'GRangesList')
mcols(urm) <- mcols(unlistedRepeatModel)
seqlevels(urm) <- seqlevels(urm)[-8]
seqinfo(urm) <- seqinfo(gnmodel)
names(urm) <- gsub('\\:(\\+|\\-)$', '', paste0(unlist(urm)))
genes_and_urepeats_model <- c(gnmodel, urm)