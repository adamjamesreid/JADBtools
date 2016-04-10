
big <- tail(sort(lengths(repeatModel)), 10)
rm <- repeatModel[names(big)]

dat <- data.frame( nam=unlist(rm)$id, wd=width(unlist(rm)))

p <- ggplot(dat, aes(factor(nam), wd))
p + geom_violin(aes(fill = factor(nam)), scale = "width") + scale_y_log10()



rep <- unlist(repeatModel)
peak <- unlist(reduce(liftOver(rep, chain), min.gapwidth=50))
peakAnno <- annotatePeak(peak, tssRegion=c(-500, 500), TxDb=txdb)

#families
unique(unlist(repeatModel)$classification_tags)  %>%  length()

#number
length(unlist(repeatModel))


#number
length(unlist(repeatModel))
