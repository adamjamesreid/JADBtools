require(rtracklayer)
library(GenomicFeatures)
devtools::load_all(".")
data("gnmodel")

chain <- import.chain(system.file('anno/WBcel235toCe10.chain', package='JADBtools'))
rdb <- loadDb('inst/anno/WBcel235_EnsemblGenes77_TxDb.sqlite')
seqlevelsStyle(rdb) <- 'UCSC'

exon <- unlist(gnmodel)
strand(exon) <- '*'

genes <- genes(rdb)
genes <- unlist(reduce(liftOver(genes, chain), min.gapwidth=45))
strand(genes) <- '*'


uRepModelNoExonNS <- unlistedRepeatModel[!(unlistedRepeatModel  %over% exon)]
uRepModelNoGenesNS <- uRepModelNoExon[!(uRepModelNoExon  %over% genes)]

uRepModelNoExonNS  %>% export.bed('uRepModelNoExonNS.bed')
uRepModelNoGenesNS %>% export.bed('uRepModelNoGenesNS.bed')
uRepModelAllNS <- unlistedRepeatModel

repgenes <- genes[overlapsAny(genes, unlistedRepeatModel, type='within')]

uRepModelNoExon <- c(uRepModelNoExon, unlistedRepeatModel[unlistedRepeatModel  %over%  repgenes])
uRepModelNoGenes <- c(uRepModelNoGenes, unlistedRepeatModel[unlistedRepeatModel  %over%  repgenes])
require(gProfileR)

uRepModelNoExonNS <- unique(uRepModelNoExon)
uRepModelNoGenesNS <- uRepModelNoGenes
uRepModelAllNS <- unlistedRepeatModel


uRepModelNoExonNS  %>% export.bed('uRepModelNoExonNS.bed')
uRepModelNoGenesNS %>% export.bed('uRepModelNoGenesNS.bed')

save(uRepModelNoExonNS, file='data/uRepModelNoExonNS.rda')
save(uRepModelNoGenesNS, file='data/uRepModelNoGenesNS.rda')
save(uRepModelAllNS, file='data/uRepModelAllNS.rda')

                                
                                