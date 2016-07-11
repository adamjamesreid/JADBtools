urm <- unlist(repeatModel)
positions <- paste0(seqnames(urm),'-',start(urm),':',end(urm))

data  %>% sapply(class) -> cl
experiment_names <- fls[cl=='numeric']

m <- do.call(cbind, data[cl=='numeric'])
colnames(m) <- basename(experiment_names)
rownames(m) <- positions
attr(m, 'id') <- urm$id
attr(m, 'cls') <- urm$classification_tags
attr(m, 'url') <- experiment_names
attr(m, 'strand') <- strand(urm)

m_with_reps <- m
m <- m[, -grep('REPLICATES', experiment_names)]
attr(m, 'url') <- experiment_names[-grep('REPLICATES', experiment_names)]
attr(m, 'id') <- urm$id
attr(m, 'cls') <- urm$classification_tags
attr(m, 'strand') <- strand(urm)


means <- colMeans(m, na.rm = TRUE)
rank <- order(means, decreasing = TRUE)
en <- experiment_names[rank][means[rank] > 2]
table(sapply(strsplit(basename(en), '\\^'), '[[', 1))

require(matrixStats)
medians <- colMedians(m, na.rm = TRUE)
rank <- order(medians, decreasing = TRUE)
en <- experiment_names[rank][medians[rank] > 1.2]

t(t(table(sapply(strsplit(basename(en), '\\^'), '[[', 1))))

by_cls <- sapply( unique(attr(m, "cls")), function(x) {
    colMeans(m[attr(m, "cls")==x,], na.rm = TRUE)
})

by_cls_names <- sapply( unique(attr(m, "cls")), function(x) {
    rank <- order(by_cls[,x], decreasing = TRUE) 
    enriched <- colnames(m)[by_cls[,x] > 3]
    enriched_vals <- by_cls[,x][by_cls[,x] > 3]
    if(length(enriched)) {
        out <- table(sapply(strsplit(enriched, '\\^'), '[[', 1))
        mvals <- sapply(names(out), function(z) {
            median(enriched_vals[grep(z, names(enriched_vals))])
        })
        out <- data.frame(out, mvals, row.names = NULL)
        out <- out[order(mvals, decreasing = TRUE),]
        out
    }
})


