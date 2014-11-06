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


export.bed(asBED(gnmodel), 'inst/anno/gnmodel.bed')
seqinfo(gnmodel) <- SeqinfoForBSGenome('ce10')[seqlevels(gnmodel)]



