require(magrittr)
require(Rsamtools)
require(GenomicAlignments)

#### Build unlisted repeats model ####

#fullmodel[mcols(fullmodel)$type == 'repeat_dfam']  %>% mcols  %>% str
data(repeatModel)
unlistedRepeatModel <- unlist(repeatModel)
names(unlistedRepeatModel)  <- unlistedRepeatModel$name
meta <- DataFrame(geneName=mcols(unlistedRepeatModel)$description, desc=mcols(unlistedRepeatModel)$all_text, seqID=mcols(unlistedRepeatModel)$id, type='repeat_dfam')    
mcols(unlistedRepeatModel) <- meta


#### Sumarize ####

load("~/code/DEview/data/AM_RNAseq_repeats.Rdata")




ids <- colnames(repeatCounts)  %>% strsplit('_')  %>% do.call(rbind, .)  %>% .[,9]

ids %>% lapply(JADBtools::getFilePath, format = 'bam') -> bam
names(bam) <- ids
if (all(bam  %>% elementLengths == 1)) bam <- unlist(bam)
bam  %>% sapply(function(x) download.file(x, basename(x)))



##FullModle
SE <- summarizeBAMs(fls = basename(bam), model = "fullmodel")
colData(SE) <-  colData(repeatCounts)[colnames(SE),]


### Combinatory 

mcols(fullmodel)$type  %>% table()
prep <- mcols(fullmodel)$type %in% c('protein_coding', 'repeat_dfam')
gr <- fullmodel[prep]
seqlevels(gr) <- seqlevelsInUse(gr)

bfl <- BamFileList(dir(pattern='bam$'))
SE  <- summarizeOverlaps(gr, bfl, mode='IntersectionNotEmpty', ignore.strand=FALSE, singleEnd=FALSE)
colData(SE)$strain <- factor(c(rep('N2', 2), rep('let418ts', 2)))
colData(SE)$stage <- factor(rep('LE', 4))

system('rsync -av --progress /mnt/home3/ahringer/ps562/tests/kallisto*Rdata ps562@ws190.gurdon.private.cam.ac.uk:/var/shiny-server/www/DEview/data')


save(SE, file='AM_nuclearRNA_genes_and_repeats_PEcounts_IntersectionNotEmpty.Rdata')
save(SE, file='AM_individual_repeat_counts_mapq0.Rdata')



m1 <- unlistedRepeatModel; seqlevels(m1) <- seqlevelsInUse(m2)
elementMetadata(m1) <- NULL
grlist <- as(m1, "GRangesList")
mcols(grlist) <- mcols(unlistedRepeatModel)
mcols(grlist)$geneName <- gsub(':\\+|:\\-', '', paste0(unlistedRepeatModel))

m2 <- fullmodel[mcols(fullmodel)$type == 'protein_coding']; seqlevels(m2) <- seqlevelsInUse(m2)


mod2 <- c(grlist, m2)
mcols(grlist)$geneName


SEunlisted  <- summarizeOverlaps(mod2, bfl, mode='IntersectionNotEmpty', ignore.strand=FALSE, singleEnd=FALSE)
colData(SEunlisted)$strain <- factor(c(rep('N2', 2), rep('let418ts', 2)))
colData(SEunlisted)$stage <- factor(rep('LE', 4))
save(SEunlisted, file='AM_nuclearRNA_genes_and_repeats_PEcounts_IntersectionNotEmpty_unlisted_repeats.Rdata')

##unlistedRepeatsModel

### mapq0
SE <- summarizeBAMs(fls = basename(bam), model = "unlistedRepeatModel")
attr(SE, 'annosource') <- 'http://www.dfam.org/entry/'
colData(SE) <-  colData(repeatCounts)[colnames(SE),]

rr <- rowRanges(SE)
rowRanges(SE)$geneName <-  paste0(seqnames(rr), ':', start(rr), '-', end(rr))

save(SE, file='/Users/przemol/code/DEview/data/AM_individual_repeat_counts_mapq0.Rdata')

### mapq10
mapq_filter <- function(features, reads, algorithm,
                        ignore.strand, inter.feature)
{ 
    require(GenomicAlignments) # needed for parallel evaluation
    Union(features, reads[mcols(reads)$mapq >= 10], algorithm,
          ignore.strand, inter.feature) 
}
param <- ScanBamParam(what="mapq")
bfl <- BamFileList( basename(bam))

SEuniq <- summarizeOverlaps(unlistedRepeatModel, bfl, mode=mapq_filter, param=param)

attr(SEuniq, 'annosource') <- 'http://www.dfam.org/entry/'
colData(SEuniq) <-  colData(repeatCounts)[colnames(SEuniq),]


cd <- DataFrame(strain=unlist(JADBtools::getStrain(ids)), stage=unlist(JADBtools::getStage(ids)), row.names = ids)

rr <- rowRanges(SEuniq)
rowRanges(SEuniq)$geneName <-  paste0(seqnames(rr), ':', start(rr), '-', end(rr))

save(SEuniq, file='/Users/przemol/code/DEview/data/AM_individual_repeat_counts_mapq10.Rdata')






