####### RAN seq
require(JADBtools)
require(GenomicAlignments)
require(Rsamtools)

data("gnmodel")
require(DESeq2)


id <- c("rYD077","rYD076","rYD101", "rYD078","rYD079","rYD102")
bam <- getFilePath(id, format = 'bam', processing = 'aligned', mount = TRUE)
bam <- BamFileList(unlist(bam))



SE <- summarizeOverlaps(gnmodel, bam, ignore.strand = TRUE)
colData(SE)$strain <- as.factor(unlist(getStrain(id)))
colnames(SE) <- id
dds <- DESeqDataSet(SE, ~strain)
dds <- DESeq(dds)

rl <- rlog(dds)
rownames(colData(rlog)) <- id
z <- plotPCA(rl, intgroup='strain', ntop = 200000)
require(ggrepel)

sub <- 
    'The count data are transforemd to the log2 scale in a way which minimizes differences between samples
for rows with small counts, and which normalizes with respect to library size. The rlog transformation 
produces a similar variance stabilizing effect as varianceStabilizingTransformation,though rlog is more 
robust in the case when the size factors vary widely. The transformation is useful when checking for outliers.'
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




results(dds, contrast = c('strain', 'cfp1', 'N2')) %>% as.data.frame() %>% 
    rownames_to_column('wb') %>% tbl_df() %>% arrange(padj) %>% filter(padj < 0.05) %>% .$wb -> r1

results(dds2, contrast = c('strain', 'cfp1', 'N2')) %>% as.data.frame() %>% 
    rownames_to_column('wb') %>% tbl_df() %>% arrange(padj) %>% filter(padj < 0.05) %>% .$wb -> r2

venn(list(OUR=r1, FP=r2))

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



SEall <- cbind(SE, SE2)
colData(SEall)$source <- factor(gsub('[0-9]', '', rownames(colData(SEall))))
