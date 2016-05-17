
big <- tail(sort(lengths(repeatModel)), 10)
rm <- repeatModel[names(big)]
rm <- repeatModel[heli]

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

#msa

load(repeatModel)
heli <- grep('Helitron', names(repeatModel), value = TRUE)
rep <- repeatModel[[heli[[1]]]]
paste(rep)  %>%  gsub('\\:(\\-|\\+)', '', .)  -> names(rep)
seq <- getSeq(Celegans, rep)
require(msa)
myFirstAlignment <- msa(seq)

msaPrettyPrint(myFirstAlignment, output="pdf", showNames="none", showLogo="none", askForOverwrite=FALSE, verbose=FALSE)