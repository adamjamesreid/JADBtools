"+" = function(x,y) {
    if(is.character(x) || is.character(y)) {
        return(paste(x , y, sep=""))
    } else {
        .Primitive("+")(x,y)
    }
}

####### RAN seq
require(GenomicAlignments)
require(Rsamtools)
require(ggrepel)
require(DESeq2)
require(tidyverse)
devtools::load_all("/mnt/home1/ahringer/jarun//JADBtools")

data("gnmodel")
sub <- 
    'The count data are transforemd to the log2 scale in a way which minimizes differences between samples
for rows with small counts, and which normalizes with respect to library size. The rlog transformation 
produces a similar variance stabilizing effect as varianceStabilizingTransformation,though rlog is more 
robust in the case when the size factors vary widely. The transformation is useful when checking for outliers.'


id <- c("rYD077","rYD076","rYD101", "rYD078","rYD079","rYD102")
bam <- getFilePath(id, format = 'bam', processing = 'aligned', mount = TRUE)
bam <- BamFileList(unlist(bam))



SEce10 <- summarizeOverlaps(gnmodel, bam, ignore.strand = TRUE)
colData(SEce10)$strain <- as.factor(unlist(getStrain(id)))
colnames(SEce10) <- id

ddsce10 <- DESeqDataSet(SEce10, ~strain)
ddsce10 <- DESeq(ddsce10)

RES_ce10 <- results(ddsce10, contrast = c('strain', 'cfp1', 'N2')) %>% as.data.frame() %>% 
    rownames_to_column('wb') %>% tbl_df %>% arrange(padj)
write_csv(left_join(RES,ann), 'DE_cfp1_vs_N2_rYD077_rYD076_rYD101_rYD078_rYD079_rYD102.csv')

rl <- rlog(dds)
rownames(colData(rlog)) <- id
z <- plotPCA(rl, intgroup='strain', ntop = 200000)



z +  ggtitle('PCA: rlog transformed counts on BWA aligend data', subtitle = sub) + geom_label_repel(aes(label=name)) + theme_classic(base_size = 12)
z <- plotPCA(rl2, intgroup='strain', ntop = 200000)
z +  ggtitle('PCA: rlog transformed counts on BWA aligend data', subtitle = sub) + geom_label_repel(aes(label=name)) + theme_classic(base_size = 12)


mySVD <- function(SE, type='center', sf=TRUE, hmap=FALSE , main='') {
    
    
    Strain <- colData(SE)$strain
    dds <- DESeqDataSet(SE, ~strain)
    dds <- estimateSizeFactors(dds)
    
    m <- counts(dds, normalized=sf)
    
    if(type=='center') {
        X <-  t( apply(m, 1, function(V){V-mean(V)} ) )
    } else if(type=='tags') {
        X <- m
    } else if(type=='rlog') {
        X <- assay(rlog(dds))
        #X <-  t( apply(X, 1, function(V){V-mean(V)} ) )
    }
    
    
    svd  =  svd(X)
    V  = svd$v
    S  = matrix(0, ncol(m), ncol(m))
    for (i in 1:ncol(m)) { S [i,i] =  svd$d[i] }
    
    pvc  <- round((svd$d^2/sum(svd$d^2))*100, 2)
    
    if(hmap) {
        image( 
            1:ncol(m), 1:ncol(m), S%*%t(V),
            main="Singular Value Decomposition", 
            #ylab = sprintf("Samples 1 - %i", ncol(m)), 
            ylab = "",
            xlab = sprintf("%i components", ncol(m)),
            yaxt="n"
        )
        axis(2, 1:ncol(m),labels=colnames(m), las=2)
    }
    #PCA plot
    #par( mfrow = c(1,2) )
    
    
    set.seed(1)
    '%p%' <- paste0
    require(ggrepel)
    require(ggplot2)
    ggplot(as.data.frame(V), aes(x=V1, y=V2, color=Strain, fill = Strain)) + 
        geom_point() +
        geom_label_repel(
            aes(label=colnames(X)),
            fontface = 'bold', color = 'white',
            box.padding = 0.35, point.padding = .5,
            segment.color = 'grey50'
        ) + theme_classic(base_size = 16) +
        xlab('PC1 [' %p% pvc[1] %p% '% variance]') + ylab('PC2 [' %p% pvc[2] %p% '% variance]') +
        ggtitle(main)
    
}
mySVD(SE)


### SE2
id <- make_jadb_ids(1:12, 'rFP')
bam <- getFilePath(id, format = 'bam', processing = 'aligned', mount = TRUE)
bam <- BamFileList(unlist(bam))

SE2 <- summarizeOverlaps(gnmodel, bam, ignore.strand = TRUE)
colData(SE2)$strain <- as.factor(unlist(getStrain(id)))
colnames(SE2) <- id
#dds2 <- DESeqDataSet(SE2[,c(1:3, 7:9)], ~strain)
dds2 <- DESeqDataSet(SE2, ~strain)
dds2 <- DESeq(dds2)

rl2 <- rlog(dds2)
z <- plotPCA(rl2, intgroup='strain', ntop = 200000)
z +  ggtitle('PCA: rlog transformed counts on BWA aligend data', subtitle = sub) + geom_label_repel(aes(label=name)) + theme_classic(base_size = 12)

mySVD(SE2, sf=FALSE, type='rlog')


#### res
results(dds, contrast = c('strain', 'cfp1', 'N2')) %>% as.data.frame() %>% 
    rownames_to_column('wb') %>% tbl_df() %>% arrange(padj) %>% filter(padj < 0.05) %>% .$wb -> r1

results(dds2, contrast = c('strain', 'cfp1', 'N2')) %>% as.data.frame() %>% 
    rownames_to_column('wb') %>% tbl_df() %>% arrange(padj) %>% filter(padj < 0.05) %>% .$wb -> r2

venn(list(JA=r1, FP=r2))

#### all
SEall <- cbind(SE, SE2)
colData(SEall)$source <- factor(gsub('[0-9]', '', rownames(colData(SEall))))










ddsA <- DESeqDataSet(SEall, ~strain)
ddsA <- DESeq(ddsA)
fp <- fpkm(ddsA)
ann <- cbind(as.data.frame(mcols(gnmodel)), wb=names(gnmodel)) %>% tbl_df()

rownames(fp) <- mcols(gnmodel[rownames(fp)])$seqID

require(readr)
time.data <- read_delim('~/analises/embryo_unified_dcpm.140308a', delim = ' ')
td <- time.data  %>% select(-Transcripts) %>% as.matrix()
rownames(td) <- time.data$Transcripts
common <- base::intersect(rownames(fp), time.data$Transcripts)

corel_p <- cor(td[common,], fp[common,], method = 'pearson')
corel_k <- cor(td[common,], fp[common,], method = 'kendall')
corel_s <- cor(td[common,], fp[common,], method = 'spearman')

ggplot( melt(corel_p), aes(x=Var1, color=Var2, y=value, group=Var2)) + geom_line()

ggplot( melt(corel_k), aes(x=Var1, color=Var2, y=value, group=Var2)) + geom_line()
ggplot( melt(corel_s), aes(x=Var1, color=Var2, y=value, group=Var2)) + geom_line()



p <- cbind(id=colnames(fp), stage_pearson=colnames(td)[apply(corel_p, 2, which.max)]) %>% tbl_df()
k <- cbind(id=colnames(fp), stage_kendall=colnames(td)[apply(corel_k, 2, which.max)]) %>% tbl_df()
s <- cbind(id=colnames(fp), stage_spearman=colnames(td)[apply(corel_s, 2, which.max)]) %>% tbl_df()
kable(left_join(left_join(p, k), s))

ans <- matrix(nrow = ncol(td), ncol = ncol(fp), dimnames = list(colnames(td), colnames(fp)))
for(i in 1:ncol(td)) {
    for(j in 1:ncol(fp)) {
        ans[i,j] <- sum((td[common,i] - fp[common,j])^2)
    }
}

stage <- colnames(td)[apply(ans, 2, which.min)]
stage <- colnames(td)[apply(ans, 2, which.max)]
cbind(colnames(fp), stage) %>% tbl_df()

ggplot( melt(log2(1/ans[,7:8])), aes(x=Var1, color=Var2, y=value, group=Var2)) + geom_line()
ggplot( melt(cbind(td[common,], fp[common,])), aes(x=Var2, color=Var2, y=value)) + geom_violin()

boxplot(cbind(td[common,], fp[common,]), las=2, ylim=c(1,100))

rl3 <- rlog(ddsA)
require(ggplot2)
z <- plotPCA(rl3, intgroup='strain', ntop = 200000)
z +  ggtitle('PCA: rlog transformed counts on BWA aligend data', subtitle = sub) + geom_label_repel(aes(label=name)) + theme_classic(base_size = 12)


    

SEwt <- SEall[,colData(SEall)$strain =='N2']
dds3 <- DESeqDataSet(SEwt, ~source)
dds3 <- DESeq(dds3)
plotMA(dds3)


RES3 <- results(dds3, contrast = c('source', 'rYD', 'rFP')) %>% as.data.frame() %>% 
    rownames_to_column('wb') %>% tbl_df %>% arrange(padj)
write_csv(left_join(RES3,ann), 'DE_YDN2_vs_FPN2.csv')


results(dds3, contrast = c('source', 'rYD', 'rFP')) %>% as.data.frame() %>% 
    rownames_to_column('wb') %>% tbl_df() %>% arrange(padj) %>% filter(padj < 0.05) %>% .$wb -> r1

#############
# [X] Staging based on Waterston RNAseq data
# [ ] Boxplots of expressed genes that are cfp-1 pos vs negative (log2FC of cfp-1/N2)
# [X] Francesca N2 vs our N2

cfp1 <- import.bed('~/analises/Chen2013_Supplemental_Information_S4_cfp1_p100.bed')
require(tidyverse)
require(GenomicRanges)

outron <- read_tsv('https://raw.githubusercontent.com/jurgjn/relmapping/master/annot/S2_regulatory_annotation/S2_outron-extended_genes.bed', col_names = FALSE)
outron %<>% makeGRangesFromDataFrame(seqnames.field = 'X1', start.field = 'X2', end.field = 'X3', keep.extra.columns = TRUE)
    
reg2 <- outron[outron %over% cfp1 ]$X4

jjano <- read_tsv('https://raw.githubusercontent.com/jurgjn/relmapping/master/annot/S2_regulatory_annotation/S2_regulatory_annotation.tsv')
jjf <- jjano %>% filter(promoter_gene_id_fwd!='.' & annot_fwd=='coding_promoter')  %>% select(chrom, start, end, promoter_gene_id_fwd, promoter_gene_id_rev) %>% makeGRangesFromDataFrame(keep.extra.columns = TRUE)
jjr <- jjano %>% filter(promoter_gene_id_rev!='.' & annot_rev=='coding_promoter')  %>% select(chrom, start, end, promoter_gene_id_fwd, promoter_gene_id_rev) %>% makeGRangesFromDataFrame(keep.extra.columns = TRUE)

rr <- unique(c(jjf, jjr)) 
tmp_all <- unique(c(rr$promoter_gene_id_fwd, rr$promoter_gene_id_rev))

cfp1_rr <- rr[rr%over% cfp1 ]
reg <- unique(c(cfp1_rr$promoter_gene_id_fwd, cfp1_rr$promoter_gene_id_rev))

par(mfcol=c(1,2))
venn(list(JJanno=tmp_all, reg=reg, ce11_genes=rownames(dds)))
venn(list(JJanno=tmp_all, reg=reg, ce10_genes=rownames(SEce10)))






dds <- DESeqDataSet(SE, ~strain)
dds <- DESeq(dds)
fpkm(dds)[, 1:3] %>% rowMeans() %>% cut_number(n=5) -> div
levels(div) <- paste('bin', 1:5, levels(div), sep='_')

bins <- sapply(levels(div), function(x) {
    rownames(fpkm(dds))[div==x]
})
bins <- cbind(wb=rownames(fpkm(dds)), bin=as.character(div)) %>% tbl_df()


RES <- results(dds, contrast = c('strain', 'cfp1', 'N2')) %>% as.data.frame() %>% 
    rownames_to_column('wb') %>% tbl_df %>% arrange(padj)

list(
    no_overlap = RES %>% filter(!wb %in% reg) %>% filter( padj < 0.05) %>% .$log2FoldChange,
    cfp1_reg_overlap = RES %>% filter(wb %in% reg) %>% filter( padj < 0.05) %>% .$log2FoldChange
) %>% boxplot(varwidth=TRUE, notch=TRUE, main='Significantly missregulated genes (padj < 0.05) \nLFC of cfp-1 regulated genes')

RES %>% mutate(cfp=wb %in% reg) %>% filter( padj < 0.05) %>% ggplot(aes(x=cfp, y=log2FoldChange)) + 
    geom_jitter() +geom_violin(alpha=0.8) + geom_boxplot(alpha=0.5,notch = TRUE) +
    ggtitle('Significantly missregulated genes (padj < 0.05) \nLFC of cfp-1 regulated genes')+ xlab('Marked by CFP1')

RES %>% mutate(cfp=wb %in% reg) %>% filter( baseMean > 1000 & baseMean < 10000) %>% ggplot(aes(x=cfp, y=log2FoldChange, fill=cfp)) + 
    geom_violin(alpha=0.6) + geom_boxplot(alpha=0.3,notch = TRUE) +
    ggtitle('Base mean 1000-10000 - LFC of cfp-1 regulated genes') + 
    xlab(sprintf('Marked by CFP1'))

RES %>% left_join(bins) %>% mutate(cfp=wb %in% reg) %>% ggplot(aes(x=cfp, y=log2FoldChange, fill=cfp)) + 
    geom_violin(alpha=0.6) + geom_boxplot(alpha=0.3,notch = TRUE, outlier.shape = NA) +
    ggtitle('LFC of cfp-1 regulated genes - 5 bins based on WT expression', paste(sum(RES$wb %in% reg), 'genes marked with CFP1')) + 
    xlab(sprintf('Marked by CFP1')) + facet_wrap(~bin, nrow = 1) +
    coord_cartesian(ylim = c(-1,1))

RES %>% mutate(bin=cut_number(baseMean, n=5)) %>% mutate(cfp=wb %in% reg) %>% ggplot(aes(x=cfp, y=log2FoldChange, fill=cfp)) + 
    geom_violin(alpha=0.6) + geom_boxplot(alpha=0.3,notch = TRUE, outlier.shape = NA) +
    ggtitle('LFC of cfp-1 regulated genes - 5 bins based on baseMean', paste(sum(RES$wb %in% reg), 'genes marked with CFP1')) + 
    xlab(sprintf('Marked by CFP1')) + facet_wrap(~bin, nrow = 1) +
    coord_cartesian(ylim = c(-1,1))



list(
    cfp1_overlap =RES %>% filter(wb %in% reg2) %>% filter( padj < 0.05) %>% .$log2FoldChange,
    no_overlap = RES %>% filter(!wb %in% reg2) %>% filter( padj < 0.05) %>% .$log2FoldChange
) %>% boxplot(varwidth=TRUE, notch=TRUE)

write_csv(left_join(RES,ann), 'DE_cfp1_vs_N2_rYD077_rYD076_rYD101_rYD078_rYD079_rYD102.csv')


############

pp <- getFilePath(c("AA758","AA759","AA760"), processing = 'raw', mount = TRUE)
require(Biostrings)
require(ShortRead)

r1 <- readFastq(pp[[1]])@sread
r2 <- readFastq(pp[[2]])@sread
r3 <- readFastq(pp[[3]])@sread

index <- read_delim('~/analises/SNAPchip.txt', delim = ' ')
pd <- PDict(index$sequence)
cc1 <- vcountPDict(pd,  r1)
cc2 <- vcountPDict(pd,  r2)
cc3 <- vcountPDict(pd,  r3)

tabb <- cbind(index, H3K27me2=rowSums(cc1), H3K9me3=rowSums(cc2), Input=rowSums(cc3))


bf <- BamFileList(unlist(getFilePath(c("AA758","AA759","AA760"), processing = 'align', mount = TRUE)))

unalnseq_r1 <- scanBam(bf[[1]], param = ScanBamParam(what = "seq", flag = scanBamFlag(isUnmappedQuery = TRUE)))
cc1 <- vcountPDict(pd,  unalnseq_r1[[1]]$seq)

unalnseq_r2 <- scanBam(bf[[1]], param = ScanBamParam(what = "seq", flag = scanBamFlag(isUnmappedQuery = TRUE)))
cc2 <- vcountPDict(pd,  unalnseq_r2[[1]]$seq)


unalnseq_r3 <- scanBam(bf[[3]], param = ScanBamParam(what = "seq", flag = scanBamFlag(isUnmappedQuery = TRUE)))
cc3 <- vcountPDict(pd,  unalnseq_r3[[1]]$seq)

tabb2 <- cbind(index, H3K27me2=rowSums(cc1), H3K9me3=rowSums(cc2), Input=rowSums(cc3))

#p1 <- scanBam(pp[[1]], param = ScanBamParam(what="seq"))

#cc <- vcountPattern(DNAString('TTCGCGCGTAACGACGTACCGT'), p1[[1]]$seq)

#cc <- vcountPattern(reverseComplement(DNAString('TTCGCGCGTAACGACGTACCGT')), p1[[1]]$seq)



cc <- vcountPDict(PDict(index$sequence),  p1[[1]]$seq)

cc2 <- vcountPDict(reverseComplement(DNAStringSet(index$sequence)),  p1[[1]]$seq)

Rsamt
reads <- readLines(pp[[1]])
############