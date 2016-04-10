require(magrittr)

files <- dir('/Volumes/data/_HOT_jointAWG_worm_fly_human/worm_Ce10_noDoubleCounting/source_data', full.names = TRUE)
stage <- strsplit(dir('/Volumes/data/_HOT_jointAWG_worm_fly_human/worm_Ce10_noDoubleCounting/source_data'),'_')  %>% sapply('[[', 3)
factor <- strsplit(dir('/Volumes/data/_HOT_jointAWG_worm_fly_human/worm_Ce10_noDoubleCounting/source_data'),'_')  %>% sapply('[[', 2)

lst <- lapply(files, ChIPseeker::readPeakFile, header=FALSE)
names(lst) <- basename(files)
modencodetfbs <- GRangesList(lst)

let <- lapply(dir('inst/let418/', pattern = 'bed', full.names = TRUE)[-1], import.bed)
let <- GRangesList(let)
names(let) <- dir('inst/let418/', pattern = 'bed')
int <- Reduce(intersect, let)
uni <- unique(unlist(let))

ovr <- findOverlaps(promoters(uni, upstream = 500, downstream = 100), modencodetfbs)

#% of genes having annotated TF
length(unique(queryHits(ovr))) / length(uni)

#% of genes having annotated TF
length(unique(queryHits(ovr))) / length(uni)

sort(factor[subjectHits(ovr)]  %>% table, decreasing = TRUE)

sort(factor[subjectHits(ovr)]  %>% table, decreasing = TRUE) / length(uni)



o2 <- findOverlaps(promoters(int, upstream = 500, downstream = 100), modencodetfbs)
sort(factor[subjectHits(o2)]  %>% table, decreasing = TRUE) / length(int)


hi <- import.bed('http://ws190.gurdon.private.cam.ac.uk:3838/seqplots/files/TSS_top_20pct.bed')
ovrR <- findOverlaps(promoters(sample(hi, 10), upstream = 500, downstream = 100), modencodetfbs)
sort(factor[subjectHits(ovrR)]  %>% table, decreasing = TRUE) / length(uni)


unique(factor)  %>%  sapply(function(x) length(unlist(modencodetfbs[factor==x])) ) -> n

