require(GenomicFeatures)
require(biomaRt)

###### TranscriptDb ######
rdb2 <- makeTranscriptDbFromBiomart(dataset="celegans_gene_ensembl", 
                                    filters=list( "biotype"="protein_coding" ), biomart="ENSEMBL_MART_ENSEMBL", 
                                    host='www.ensembl.org')
seqlevelsStyle(rdb2) <- 'UCSC'


download.file('https://gist.githubusercontent.com/Przemol/307189b979b7bcf68cb6/raw/6e91f187a56920c6461c9e32a35ad7e2c26c08ee/WBcel235toCe10.chain', destfile = 'inst/anno/WBcel235toCe10.chain', method="curl")
chain <- import.chain(system.file('anno/WBcel235toCe10.chain', package='JADBtools'))
intWBcel235 <- exonsBy(rdb2, by="gene")

allI <- unlist( intWBcel235 )
allIlo <- reduce(liftOver(allI, chain), min.gapwidth=50)
names(allIlo) <- names(allI)

gnmodel <- relist(unlist(allIlo), intWBcel235)

map <- read.csv(system.file('anno/WBgene2geneName.csv.gz', package='JADBtools'))
map <- data.frame(geneName=map$Associated.Gene.Name, desc=map$Description, row.names=map$Ensembl.Gene.ID) 
elementMetadata(gnmodel) <- map[names(gnmodel),]


# map2 <- read.csv(system.file('anno/WBgene2geneNameAndSeqName.csv.txt', package='JADBtools'))
# map2$Ensembl.Transcript.ID <- sub('(^.+\\.[0-9]+)(.+)', '\\1', map2$Ensembl.Transcript.ID); map2 <- unique(map2)
# map2 <- data.frame(geneName=map2$Associated.Gene.Name, desc=map2$Description, row.names=map2$Ensembl.Gene.ID, seqID=) 
# elementMetadata(gnmodel) <- map2[names(gnmodel),]


map2 <- read.csv(system.file('anno/WBgene2geneSeqName.csv.gz', package='JADBtools'))
map2 <- data.frame(seqID=map2$WormBase.Gene.Sequence.name.Accession, row.names=map2$Ensembl.Gene.ID) 
elementMetadata(gnmodel)$seqID <- map2[names(gnmodel),]

export.bed(asBED(gnmodel), 'inst/anno/gnmodel.bed')
seqinfo(gnmodel) <- SeqinfoForBSGenome('ce10')[seqlevels(gnmodel)]

spikein <- import.gff('/Users/przemol/code/JADBtools/inst/anno/ERCC92.gtf')
g <- elementMetadata(spikein)$group
n <- as.character(seqnames(spikein))
elementMetadata(spikein) <- NULL
spikein <-as(spikein, 'GRangesList')

elementMetadata(spikein)$geneName <- n
elementMetadata(spikein)$desc <- g
elementMetadata(spikein)$seqID <- n


spikeinmodel <- c(gnmodel, spikein)

### rRAN
require(GenomicFeatures)
require(biomaRt)
require(rtracklayer)
rdb2 <- makeTxDbFromBiomart(dataset="celegans_gene_ensembl", 
                                    filters=list( "biotype"="rRNA" ), biomart="ENSEMBL_MART_ENSEMBL", 
                                    host='www.ensembl.org')
seqlevelsStyle(rdb2) <- 'UCSC'


chain <- import.chain(system.file('anno/WBcel235toCe10.chain', package='JADBtools'))
intWBcel235 <- exonsBy(rdb2, by="gene")

allI <- unlist( intWBcel235 )
allIlo <- unlist( liftOver(allI, chain) )
names(allIlo) <- names(allI)

rrnamodel <- relist(unlist(allIlo), intWBcel235)



require(AnnotationDbi)
require(rtracklayer)
require(XLConnect)

db <- loadDb('inst/anno/WBcel235_EnsemblGenes77_TxDb.sqlite')
seqlevelsStyle(db) <- 'UCSC'
chain <- import.chain(system.file('anno/WBcel235toCe10.chain', package='JADBtools'))
g0 <- genes(db)

genes <- unlist(reduce(liftOver(g0, chain), min.gapwidth=50))
names(genes) <- names(g0)

map <- read.csv(system.file('anno/WBgene2geneName.csv.gz', package='JADBtools'))
map <- data.frame(geneName=map$Associated.Gene.Name, desc=map$Description, row.names=map$Ensembl.Gene.ID) 
map2 <- read.csv(system.file('anno/WBgene2geneSeqName.csv.gz', package='JADBtools'))
map2 <- data.frame(seqID=map2$WormBase.Gene.Sequence.name.Accession, row.names=map2$Ensembl.Gene.ID) 

elementMetadata(genes) <- map[names(genes),]
elementMetadata(genes)$seqID <- map2[names(genes),]


###ChromatinStatesGeneAssignment_unique

setwd("/Volumes/raid0/_ActiveInactive/")
a <- import.bed("S3_Table_activeL3.bed"); seqlevelsStyle(a) <- 'UCSC'
i <- import.bed("S4_Table_inactiveL3.bed"); seqlevelsStyle(i) <- 'UCSC'
m <- import.bed("S5_Table_marked-borders.bed"); seqlevelsStyle(m) <- 'UCSC'
r <- import.bed("S6_Table_regulatory-borders.bed"); seqlevelsStyle(r) <- 'UCSC'

out <- list()
out$active <- genes[genes %over% a]; gg <- genes[!genes %over% a]
out$inactive <- gg[gg %over% i]; gg <- gg[!gg %over% i]
out$marked <- gg[gg %over% m]; gg <- gg[!gg %over% m]
out$regulatory <- gg[gg %over% r]; gg <- gg[!gg %over% r]


wb <- loadWorkbook("ChromatinStatesGeneAssignment_unique.xlsx", create = TRUE)
lapply(names(out), function(x) {
    createSheet(wb, name = x)
    export.bed(out[[x]], paste0(x, '.bed'))
    writeWorksheet(wb, as.data.frame(out[[x]]), sheet = x) 
})
saveWorkbook(wb)


out <- list()
out$active_NONunique <- genes[genes %over% a]
out$inactive_NONunique <- genes[genes %over% i]
out$marked_NONunique <- genes[genes %over% m]
out$regulatory_NONunique <- genes[genes %over% r]


wb <- loadWorkbook("ChromatinStatesGeneAssignment_NONunique.xlsx", create = TRUE)
lapply(names(out), function(x) {
    createSheet(wb, name = x)
    export.bed(out[[x]], paste0(x, '.bed'))
    writeWorksheet(wb, as.data.frame(out[[x]]), sheet = x) 
})
saveWorkbook(wb)

spaned <- genes[genes %over% reduce(md, ignore.strand=TRUE) & genes %over% a & genes %over% i]
spanedm <- genes[genes %over% m & genes %over% a & genes %over% i]
spanedr <- genes[genes %over% r & genes %over% a & genes %over% i]
export.bed(c(spanedr, spanedm), 'spaned.bed')

gg <- genes[!(genes %over% m & genes %over% a & genes %over% i) | (genes %over% r & genes %over% a & genes %over% i)]

out <- list()
out$marked_edge_inactive <- gg[gg %over% m & gg %over% i]
out$marked_edge_active <- gg[gg %over% m & gg %over% a]
out$regulatory_edge_inactive <- gg[gg %over% r & gg %over% i]
out$regulatory_edge_active <- gg[gg %over% r & gg %over% a]

out <- list()
out$within_inactive <- gg[gg %within% i]
out$within_active <- gg[gg %within% a]
out$within_marked <- gg[gg %within% m]
out$within_regulatory <- gg[gg %within% r]

md <- import.bed("ans-g8aw2-worm-WS220-r37-seed-oC-gaps-with-bstate-verNz.bed"); seqlevelsStyle(md) <- 'UCSC'
rd <- import.bed("ans-g8aw2-worm-WS220-r37-seed-oC-gaps-without-bstate-verNz.bed"); seqlevelsStyle(rd) <- 'UCSC'

out <- list()
out$marked_edge_inactive <- genes[genes %over% reduce(resize(md, 1, fix="start"), ignore.strand=TRUE)]
out$marked_edge_active <- genes[genes %over% reduce(resize(md, 1, fix="end"), ignore.strand=TRUE)]
out$regulatory_edge_inactive <- genes[genes %over% reduce(resize(rd, 1, fix="start"), ignore.strand=TRUE)]
out$regulatory_edge_active <- genes[genes %over% reduce(resize(rd, 1, fix="end"), ignore.strand=TRUE)]

wb <- loadWorkbook("ChromatinStatesGeneAssignmentWithin.xlsx", create = TRUE)
lapply(names(out), function(x) {
    createSheet(wb, name = x)
    writeWorksheet(wb, as.data.frame(out[[x]]), sheet = x) 
})
saveWorkbook(wb)

lapply(names(out), function(x) {
    export.bed(out[[x]], paste0(x, '.bed'))
})

gr <- GRangesList(i=reduce(i),a=reduce(a),m=reduce(m),r=reduce(r))
d <- unlist(width(gr))
dat <- data.frame(cond=names(d), d=d)
ggplot(dat, aes(x=d, fill=cond)) + geom_density(alpha=.3) + xlim(0, 10000) + ggtitle("Width of [a]ctive, [i]nactive, [m]arked and [r]egulatory states\n")
ggplot(dat, aes(x=d, fill=cond)) + geom_density(alpha=.3, position="stack") + xlim(0, 10000) + ggtitle("Width of [a]ctive, [i]nactive, [m]arked and [r]egulatory states\n")
ggplot(dat, aes(x=d, fill=cond)) + geom_density(alpha=.3, position="fill") + xlim(0, 10000) + ggtitle("Width of [a]ctive, [i]nactive, [m]arked and [r]egulatory states\n")


out <- list()
genes[genes %over% i]  %>% .[nearest(resize(md, 1, fix="start"), ., ignore.strand=TRUE)] %>% unique -> out$marked_edge_inactive
genes[genes %over% a]  %>% .[nearest(resize(md, 1, fix="end"),   ., ignore.strand=TRUE)] %>% unique -> out$marked_edge_active
genes[genes %over% i]  %>% .[nearest(resize(rd, 1, fix="start"), ., ignore.strand=TRUE)] %>% unique -> out$regulatory_edge_inactive
genes[genes %over% a]  %>% .[nearest(resize(rd, 1, fix="end"),   ., ignore.strand=TRUE)] %>% unique -> out$regulatory_edge_active


wb <- loadWorkbook("ChromatinStatesGeneAssignment_gaps_active_inactive.xlsx", create = TRUE)
lapply(names(out), function(x) {
    createSheet(wb, name = x)
    export.bed(out[[x]], paste0(x, '.bed'))
    writeWorksheet(wb, as.data.frame(out[[x]]), sheet = x) 
})
saveWorkbook(wb)

gi <- genes[genes %over% i]
ga <- genes[genes %over% a] 
ans <- lapply(md, function(x){
    unique(c(
        gi[nearest(resize(x, 1, fix="start"), gi)],
        ga[nearest(resize(x, 1, fix="end"  ), ga)]
    ))
    
})
grl <- GRangesList(ans)

GetGenePairs <- function(gaps, desc){
    genes[genes %over% i]  %>% distanceToNearest(resize(gaps, 1, fix="start"), .,  ignore.strand=TRUE) -> d1
    genes[genes %over% a]  %>% distanceToNearest(resize(gaps, 1, fix="end"),   .,  ignore.strand=TRUE) -> d2
    
    findDups <- function(d) {
        dup <- d@subjectHits[duplicated(d@subjectHits)]
        unlist(sapply(dup, function(x) {
            o <- order(elementMetadata(d[d@subjectHits == x])$distance)
            d[d@subjectHits == x]@queryHits[o][-1]
        }))
    }
    dups <- c(findDups(d1), findDups(d2))
    d1 <- d1[-sort(dups)]
    d2 <- d2[-sort(dups)]
    
    require(GenomicAlignments)
    t <- GAlignmentPairs(as(gi[d1@subjectHits], 'GAlignments'), as(ga[d2@subjectHits], 'GAlignments'), rep(TRUE, length(d1@subjectHits)))
    t <- t[granges(first(t)) != granges(last(t))]
    t <- t[!overlapsAny(first(t), last(t), ignore.strand=T)]
    seqinfo(t) <- SeqinfoForBSGenome('ce10')
    export(t, paste0(desc, ".bam"))
    
    g <- grglist(t); strand(g) <- '*'; export.bed(asBED(g), paste0(desc, '.bed'))
    return(t)
}
mp <- GetGenePairs(md, "GenePairs_marked")
rp <- GetGenePairs(rd, "GenePairs_regulatory")


table(strand(first(mp)) == strand(last(mp)))
table(strand(first(rp)) == strand(last(rp)))

###### FULL TranscriptDb ######

require(GenomicFeatures)
require(biomaRt)
rdb2 <- makeTxDbFromBiomart(dataset="celegans_gene_ensembl", biomart="ENSEMBL_MART_ENSEMBL", host='www.ensembl.org')
seqlevelsStyle(rdb2) <- 'UCSC'

mart <- useMart(dataset="celegans_gene_ensembl",biomart="ENSEMBL_MART_ENSEMBL", host='www.ensembl.org')
anno <- getBM(c("ensembl_gene_id", "wormbase_gene_seq_name", 'external_gene_name', "gene_biotype", "description"), mart=mart)

download.file('https://gist.githubusercontent.com/Przemol/307189b979b7bcf68cb6/raw/6e91f187a56920c6461c9e32a35ad7e2c26c08ee/WBcel235toCe10.chain', destfile = 'inst/anno/WBcel235toCe10.chain', method="curl")
chain <- import.chain(system.file('anno/WBcel235toCe10.chain', package='JADBtools'))

intWBcel235 <- exonsBy(rdb2, by="gene")
allI <- unlist( intWBcel235 )
allIlo <- reduce(liftOver(allI, chain), min.gapwidth=50)
names(allIlo) <- names(allI)

model <- relist(unlist(allIlo), intWBcel235)

map <- DataFrame(geneName=anno$external_gene_name, desc=anno$description, seqID=anno$wormbase_gene_seq_name, type=anno$gene_biotype, row.names=anno$ensembl_gene_id) 
elementMetadata(model) <- map[names(model),]
seqinfo(model) <- SeqinfoForBSGenome('ce10')[seqlevels(model)]


repeatModel %>% sapply(function(x) x[1])  %>% GRangesList   %>% unlist  %>% mcols -> ranno
rownames(ranno) <- ranno$id
names(rm) <- ranno[names(rm),]$name

rmap <- DataFrame(geneName=ranno$description, desc=ranno$all_text, seqID=ranno$id, type='repeat_dfam', row.names=ranno$name) 

elementMetadata(rm) <- rmap[names(rm),]
seqinfo(rm) <- SeqinfoForBSGenome('ce10')[seqlevels(rm)]

fullmodel <- c(model, rm, ercc)

