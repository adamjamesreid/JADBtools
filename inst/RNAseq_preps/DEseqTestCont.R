

require(GenomicFeatures)
require(rtracklayer)
txdb_fls <- file.path(MOUNT, "_ref_genomes_/ce11_star/Caenorhabditis_elegans.WBcel235.90.gtf")
txdb_gtf <- import.gff(txdb_fls, format='gtf')
txdb <- makeTxDbFromGRanges(txdb_gtf[txdb_gtf$gene_biotype == 'protein_coding'])

elementMetadata(txdb_gtf) %>% as.data.frame %>% tbl_df %>% 
    filter(gene_biotype == 'protein_coding', type=='exon') %>% 
    mutate(seq_name=(gsub('^([A-Za-z0-9]+\\.[0-9+]).+', '\\1', exon_id))) %>% 
    select(wb=gene_id, gene_name, seq_name) %>% unique -> ANNO

ce11model <- exonsBy(txdb, 'gene')

## Exon vs. intron
introns <- intronsByTranscript(txdb, use.names=TRUE)
SEintrons <- summarizeOverlaps(introns, bam, ignore.strand = TRUE)
require(DESeq2)
ddsi <- DESeqDataSet(SEintrons, ~1)
fi <- fpkm(ddsi, robust = TRUE)

exons <- exonsBy(txdb, 'tx', use.names=TRUE)
SEexons <- summarizeOverlaps(exons, bam, ignore.strand = TRUE)
ddse <- DESeqDataSet(SEexons, ~1)
fe <- fpkm(ddse, robust = TRUE)


log2(fi[,3:4]/fi[,1:2]) %>% head

SEtest <- cbind(SEexons[,1:2], SEexons[,1:2])
colData(SEtest)$cnd <- factor(c('exon', 'exon', 'intron', 'intron'))
ddst <- DESeqDataSet(SEtest, ~cnd)
ddst <- DESeq(ddst)
res <- results(ddst, contrast = c('cnd', 'exon', 'intron'))



require(ggrepel)
data.frame(exon=rowMeans(fe[,1:2]), intron=rowMeans(fi[,1:2])) %>% 
    rownames_to_column('wb') %>% tbl_df %>% 
    ggplot(aes(x=exon, y=intron, label=wb)) + 
    geom_point() + scale_x_log10() + scale_y_log10() +
    ggtitle('WT')

data.frame(exon=rowMeans(fe[,3:4]), intron=rowMeans(fi[,3:4])) %>% 
    rownames_to_column('wb') %>% tbl_df %>% 
    ggplot(aes(x=exon, y=intron, label=wb)) + 
    geom_point() + scale_x_log10() + scale_y_log10()+
    ggtitle('hpl-2')

+
    ggrepel::geom_label_repel(
        data = . %>% arrange(padj) %>% head(10), color='black'
    ) + theme_bw(base_size = 16) +
    ggtitle('Means of WT (rYD077,rYD076,rYD101) vs. cfp1 (rYD078,rYD079,rYD102) \nLog10 of counts, size and color denote significance')



id <- 

de_getSTAR_dds <- function(id) {
    require(BiocParallel)
    bpparam()
    bam <- getFilePath(id, format = 'bam', processing = 'STAR', mount = TRUE)
    bam <- BamFileList(unlist(bam))
    
    message('Counting ', length(path(bam)), ' bam files')
    SE <- summarizeOverlaps(ce11model, bam, ignore.strand = TRUE)
    colData(SE)$strain <- as.factor(unlist(getStrain(id)))
    colnames(SE) <- id
    
    dds <- DESeqDataSet(SE, ~strain)
    dds <- DESeq(dds)
    return(dds)
}

dds_ja_ce11 <- de_getSTAR_dds(c("rYD077","rYD076","rYD101", "rYD078","rYD079","rYD102"))
dds_ja_set2_ce11 <- de_getSTAR_dds(c("rYD077","rYD076","rYD101", "rYD078","rYD079","rYD102", 'rYD103', 'rYD104'))
dds_fp_ce11 <- de_getSTAR_dds(make_jadb_ids(1:12, 'rFP'))

de_get_res <- function(x, contrast) {
    RES <- results(x, contrast = c('strain', contrast, 'N2')) %>% as.data.frame() %>% 
        rownames_to_column('wb') %>% tbl_df %>% left_join(ANNO, .) %>% arrange(padj) 
    return(RES)
}

export_dds <- function(dds, name='test1') {
    
    fpkm(dds) %>% as.data.frame %>% rownames_to_column('wb') %>% tbl_df %>% 
        write_xlsx(path = paste0('fpkm_', name, '_roboust.xlsx'))
    assay(dds) %>% as.data.frame %>% rownames_to_column('wb') %>% tbl_df %>% 
        write_xlsx(path = paste0('counts_', name, '.xlsx'))
    colData(dds) %>% as.data.frame %>% rownames_to_column('ID') %>% tbl_df %>% 
        write_xlsx(path = paste0('coldata_', name, '.xlsx'))
    save(dds, file=paste0('dds_', name, '.Rdata'))
    
}
export_dds(dds_ja_set2_ce11, 'ja_set2_cfp1_ce11')


#dds_fp_ce11 %>% colData %>% .$strain %>%  levels



fpkm(dds_fp_ce11) %>% as.data.frame %>% rownames_to_column('wb') %>% tbl_df %>% write_xlsx(path = 'fpkm_sin3_set2_cfp1_fp_ce11_roboust.xlsx')
assay(dds_fp_ce11) %>% as.data.frame %>% rownames_to_column('wb') %>% tbl_df %>% write_xlsx(path = 'counts_sin3_set2_cfp1_pf_ce11_roboust.xlsx')
colData(dds_fp_ce11) %>% as.data.frame %>% rownames_to_column('ID') %>% tbl_df %>% write_xlsx(path = 'coldata_sin3_set2_cfp1_fp_ce11_roboust.xlsx')




all_res <- list(
    ja_cfp1=dds_ja_ce11 %>%  de_get_res('cfp1'),
    fp_cfp1=dds_fp_ce11 %>%  de_get_res('cfp1'),
    fp_set2=dds_fp_ce11 %>%  de_get_res('set2'),
    fp_sin3=dds_fp_ce11 %>%  de_get_res('sin3')
)
sig_res <- lapply(all_res, dplyr::filter, padj<0.05)


pp <- promoters(unlist(range(ce11model)), upstream = 500, downstream = 0)
seqlevelsStyle(pp) <- 'UCSC'

wb_prom <- data_frame(wb=names(pp), promoter=paste(pp))
wb_prom$promoter %>% GRanges()

tab <- read_tsv('https://raw.githubusercontent.com/jurgjn/relmapping/master/annot/S2_regulatory_annotation/S2b_promoter_annotation_6Dec17.tsv')
jj <- tab %>% makeGRangesFromDataFrame(keep.extra.columns = TRUE)

sig <- unlist(summary(BigWigFile(getFilePath('RC009', processing = 'SQ10N', scale = 'lin', mount = TRUE)), pp))
sig <- unlist(summary(BigWigFile(getFilePath('RC010', processing = 'SQ10N', scale = 'lin', mount = TRUE)), pp))

sig <- data_frame(promoter=paste(sig), val=sig$score) 

sig <-left_join(wb_prom, sig)
ja_cfp1_ex <- left_join(ja_cfp1, sig)
ja_cfp1_ex %<>% mutate(signif=padj<0.05)

ja_cfp1_ex

give.n <- 
gg0 <- ja_cfp1_ex %>% filter(baseMean > 5, !is.na(signif)) %>% select(val, signif) %>% 
    ggplot(aes(x=signif, y=val, fill=signif)) + scale_y_log10() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(legend.position="bottom") +
    guides(fill=guide_legend(nrow=1,byrow=TRUE, title="")) 

gg1 <- gg0 + geom_boxplot() + geom_violin(alpha=0.3) #+ geom_jitter(alpha=0.2) #
gg1 <- gg1 +  stat_summary(fun.data = function(x){return(c(y = min(x)*1.1, label = length(x)))}, geom = "text")

ggslackr(gg1, width = 16, height = 10)

gg2 <- gg1 + facet_grid(~ cl)
ggslackr(gg2, width = 16, height = 10)


gg3 <- gg0 +  geom_boxplot(aes(x=gene, y=e, fill=cl)) 
ggslackr(gg3, width = 16, height = 10)



library(writexl)
library(googledrive)
xlsx <- writexl::write_xlsx(all_res, path = "cfp1_sin3_set2_de_res_star_WBcel235e90.xlsx")
#options(httr_oob_default=TRUE)
gfile <- drive_upload(
    xlsx, 
    'DE/cfp1_sin3_set2_de_res_star_WBcel235e90_auto_clust', 
    type='spreadsheet'
)
drive_browse(gfile)
    


str(sig, 1)

# Pull treatment and control genes out of population
gnT <- unlist(range(ce11model[ sig$ja_cfp$wb ]))
sig$ja_cfp$baseMean %>% min
gnC <- unlist(range(ce11model[ all_res$ja_cfp %>% dplyr::filter(baseMean>8.489622) %>% dplyr::filter(!wb %in% sig$ja_cfp$wb) %>% pull(wb) ]))
require(JADBtools)
seqlevelsStyle(gnT) <- 'UCSC'
d1 <- unlist(summary(
    BigWigFile(getFilePath('RC009', processing = 'SQ10N', scale = 'lin', mount = TRUE))
    , gnT)
)$score

seqlevelsStyle(gnC) <- 'UCSC'
d2 <- unlist(summary(
    BigWigFile(getFilePath('RC009', processing = 'SQ10N', scale = 'lin', mount = TRUE))
    , gnC)
)$score

boxplot(list(d1, d2), ylim=c(0,5))
title('')




RES <- results(dds_fp_ce11, contrast = c('strain', 'cfp1', 'N2')) %>% as.data.frame() %>% 
    rownames_to_column('wb') %>% tbl_df %>% left_join(ANNO, .) %>% arrange(padj)




require("org.Ce.eg.db")
require(clusterProfiler)
rownames(dds) %>% head %>% bitr("WORMBASE", "SYMBOL", "org.Ce.eg.db") %>% tbl_df

ce11model



plotMA(dds, ylim=c(-4, 4))

assay(SE) %>% as.data.frame() %>% rownames_to_column('wb') %>% tbl_df %>% 
    mutate(N2=(rYD077+rYD076+rYD101)/3, cfp1=(rYD078+rYD079+rYD102)/3) %>% 
    left_join(RES) %>% mutate(col=-10*log10(padj), significant=padj < 0.05) %>% 
    left_join(ANNO) %>% 
    ggplot(aes(x=N2, y=cfp1, color=significant, size=col, label=gene_name)) + 
    geom_point(aes(alpha=col)) + scale_x_log10() + scale_y_log10() +
    ggrepel::geom_label_repel(
        data = . %>% arrange(padj) %>% head(10), color='black'
    ) + theme_bw(base_size = 16) +
    ggtitle('Means of WT (rYD077,rYD076,rYD101) vs. cfp1 (rYD078,rYD079,rYD102) \nLog10 of counts, size and color denote significance')


require(matrixStats)
assay(SE) %>% as.data.frame() %>% rownames_to_column('wb') %>% tbl_df %>% 
    mutate( N2=rowMedians(cbind(rYD077,rYD076,rYD101)), cfp1=rowMedians(cbind(rYD078,rYD079,rYD102)) ) %>% 
    left_join(RES) %>% mutate(col=-10*log10(padj), significant=padj < 0.05) %>%
    left_join(ANNO) %>% 
    ggplot(aes(x=N2, y=cfp1, color=significant, size=col, label=gene_name)) + 
    geom_point(aes(alpha=col)) + scale_x_log10() + scale_y_log10() +
    ggrepel::geom_label_repel(
        data = . %>% arrange(padj) %>% head(10), color='black'
    ) + theme_bw(base_size = 14) +
    ggtitle('Medians of WT (rYD077,rYD076,rYD101) vs. cfp1 (rYD078,rYD079,rYD102) \nLog10 of counts, size and color denote significance')



RES %>% left_join(ANNO) %>% mutate(Significant=padj < 0.05) %>% 
    ggplot(aes(x = log2FoldChange, y = -log10(pvalue))) +
        geom_point(aes(color = Significant)) +
        scale_color_manual(values = c("grey", "red")) +
        theme_bw(base_size = 16) +
        geom_text_repel(
            data = . %>% arrange(padj) %>% head(25),
            aes(label = gene_name),
            size = 3,
            box.padding = 0.25,
            point.padding = 0.3
        )





paste(id, collapse='|')



ce11_genes <- import.gff()
gemes makeTxDbFromGFF(txdb)

gui_('labfiles')

### Junctions and coverage 
a1 <- readGAlignments(bam[[1]])



#### anno 2
gn <- genes(txdb)
gnp <- resize(gn, width(gn+500), fix="end", use.names=TRUE,ignore.strand=FALSE)
##gnp <- promoters(gn)
seqlevelsStyle(gnp) <- 'UCSC'
reg2 <- names(gnp[gnp %over% cfp1])
venn(list(JJpromoters=reg, WBptomGB=reg2))


promoters(, upstream = 500, downstream = width(genes(txdb)))


venn(list(JJpromoters=reg, WBptomGB=reg2))

## diff splice
ids <- c('rAM057',  'rAM058', 'rAM061',  'rAM062')
bamlst <- getFilePath(ids, format = 'bam', processing = 'STAR', mount = TRUE)
bamlst <- BamFileList(unlist(bamlst))
exonicParts <- disjointExons( txdb, aggregateGenes=FALSE )

require(GenomicAlignments)
SEdex <- summarizeOverlaps( 
    exonicParts, bamlst, mode="Union",
    singleEnd=FALSE, ignore.strand=TRUE, inter.feature=FALSE,
    fragments=TRUE 
)
#SEdex_SE <- summarizeOverlaps( exonicParts, bamlst, mode="Union",singleEnd=TRUE, ignore.strand=TRUE, inter.feature=FALSE)
colData(SEdex)$strain <- as.factor(unlist(getStrain(ids)))
colnames(SEdex) <- ids

require(DEXSeq)
dxd <- DEXSeqDataSetFromSE( SEdex, design= ~ sample + exon + strain:exon )





