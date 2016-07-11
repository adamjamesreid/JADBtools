
big <- tail(sort(lengths(repeatModel)), 10)
rm <- repeatModel[names(big)]
rm <- repeatModel[heli]

dat <- data.frame( Class=factor(unlist(rm)$id), wd=width(unlist(rm)))

p <- ggplot(dat, aes(Class, wd)) + xlab('Class') + ylab('repeat widts')
p + geom_violin(aes(fill = Class), scale = "width") + scale_y_log10()



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

data(repeatModel)
require(magrittr)
require(BSgenome.Celegans.UCSC.ce10)
heli <- grep('Helitron', names(repeatModel), value = TRUE)
rep <- repeatModel[[heli[[1]]]]
paste(rep)  %>%  gsub('\\:(\\-|\\+)', '', .)  -> names(rep)
seq <- getSeq(Celegans, rep)
require(msa)
myFirstAlignment <- msa(seq)
conMat <- consensusMatrix(myFirstAlignment)


seq2 <- getSeq(Celegans, repeatModel)
cc <- sapply(seq2, function(seq) mean(vcountPattern('TACBGTA', seq, fixed = FALSE)))
cc2 <- sapply(seq2, function(seq) sum(vcountPattern('TACBGTA', seq, fixed = FALSE))/(sum(width(seq))/1e3)  )




hpl2val <- sapply(levels(M$class), function(n) mean(M[M$class==n,2], na.rm=T) )

msaPrettyPrint(myFirstAlignment, output="pdf", showNames="none", showLogo="none", askForOverwrite=FALSE, verbose=FALSE)

