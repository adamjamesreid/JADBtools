require(magrittr)
require(Rsamtools)
require(GenomicAlignments)

load("~/code/DEview/data/AM_RNAseq_repeats.Rdata")

ids <- colnames(repeatCounts)  %>% strsplit('_')  %>% do.call(rbind, .)  %>% .[,9]

ids %>% lapply(JADBtools::getFilePath, format = 'bam') -> bam
names(bam) <- ids
if (all(bam  %>% elementLengths == 1)) bam <- unlist(bam)
bam  %>% sapply(function(x) download.file(x, basename(x)))




SE <- summarizeBAMs(fls = basename(bam), model = "fullmodel")


bfl <- BamFileList(fls)
expset_repeats <- summarizeOverlaps( get(data(list = model)), bfl[1] )
return(expset_repeats)