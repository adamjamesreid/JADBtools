require(GenomicFeatures)
require(rtracklayer)
require(JADBtools)
require(magrittr)
require(tidyverse)

txdb_fls <- file.path(MOUNT, "_ref_genomes_/ce11_star/Caenorhabditis_elegans.WBcel235.90.gtf")
txdb_gtf <- import.gff(txdb_fls, format='gtf')
txdb <- makeTxDbFromGRanges(txdb_gtf[txdb_gtf$gene_biotype == 'protein_coding'])

elementMetadata(txdb_gtf) %>% as.data.frame %>% tbl_df %>% 
    filter(gene_biotype == 'protein_coding', type=='exon') %>% 
    mutate(seq_name=(gsub('^([A-Za-z0-9]+\\.[0-9+]).+', '\\1', exon_id))) %>% 
    select(wb=gene_id, gene_name, seq_name) %>% unique -> ANNO


ce11model <- exonsBy(txdb, 'gene')
ce11genes <- genes(txdb)



intron <- IRanges::setdiff(ce11genes, unlist(ce11model))
ovr <- findOverlaps(ce11genes, intron)
as.list(ovr) %>% str



require(pbapply)
#Takes ~2min on single core
intron_lst <- pblapply(as.list(ovr), function(x) intron[x])
intron_grl <- GRangesList(intron_lst)
names(intron_grl) <- names(ce11genes)

export.bed(unlist(intron_grl), 'intron_grl_v2.bed')

ce11introns <- intron_grl

## Exon vs. intron

de_exons_introns <- function(id) {
    require(BiocParallel)
    require(GenomicAlignments)
    bpparam()
    bam <- getFilePath(id, format = 'bam', processing = 'STAR', mount = TRUE)
    bam <- BamFileList(unlist(bam))
    
    message('Counting exons ', length(path(bam)), ' bam files')
    SEexons <- summarizeOverlaps(ce11model, bam, ignore.strand = TRUE)
    colData(SEexons)$strain <- as.factor(unlist(getStrain(id)))
    colnames(SEexons) <- id
    
    message('Counting introns ', length(path(bam)), ' bam files')
    SEintron <- summarizeOverlaps(ce11introns, bam, ignore.strand = TRUE)
    colData(SEintron)$strain <- as.factor(unlist(getStrain(id)))
    colnames(SEintron) <- id
    
    #dds <- DESeqDataSet(SE, ~strain)
    #dds <- DESeq(dds)
    return(list(SEexons=SEexons, SEintron=SEintron))
}

#rAM062

##HPL2 in starvedL1
#rAM051|rAM052|rAM047|rAM048
IDS <- c('rAM051', 'rAM052', 'rAM047', 'rAM048')

#set25met2 in YA
IDS <- c("rAM086","rAM085", "rAM058","rAM057")


getDS <- function(IDS) {
    se_lst <- de_exons_introns(IDS)
    require(DESeq2)
    ddse <- DESeqDataSet(se_lst$SEexons, ~1)
    ddsi <- DESeqDataSet(se_lst$SEintron, ~1)
    
    fi <- fpkm(ddsi, robust = TRUE)
    colnames(fi) <- paste0('intron_', colnames(fi))
    fi %<>% as.data.frame %>% rownames_to_column('wb') %>% tbl_df
    fi %<>% dplyr::mutate( intron_hpl2 = (intron_rAM051+intron_rAM052)/2 ) %>% dplyr::mutate(intron_wt=(intron_rAM047+intron_rAM048)/2)
    
    
    
}

#attach(se_lst)





fe <- fpkm(ddse, robust = TRUE)
colnames(fe) <- paste0('exon_', colnames(fe))
fe %<>% as.data.frame %>% rownames_to_column('wb') %>% tbl_df
fe %<>% dplyr::mutate( exon_hpl2 = (exon_rAM051+exon_rAM052)/2 ) %>% dplyr::mutate(exon_wt=(exon_rAM047+exon_rAM048)/2)

fi %<>% select(-intron_rAM051, -intron_rAM052, -intron_rAM047, -intron_rAM048)
fe %<>% select(-exon_rAM051, -exon_rAM052, -exon_rAM047, -exon_rAM048)

plot_data <- fe %>% left_join(fi)
plot_data <- plot_data[rowSums(is.na(plot_data[,-1]))==0,]
plot_data <- plot_data[rowSums(plot_data[,-1]==0)==0,]



plot_data %<>% mutate(exon_ratio=exon_hpl2/exon_wt, intron_ratio=intron_hpl2/intron_wt)
plot_data %<>% mutate(FC=abs(log2(exon_ratio/intron_ratio)))
plot_data %<>% mutate(Significant=FC > 2.5)
plot_data <- left_join(plot_data, ANNO)


gg0 <- plot_data %>% 
    ggplot(aes(x=exon_ratio, y=intron_ratio, label=gene_name)) + 
    geom_abline(slope = 1, intercept = 0) +
    #geom_abline(slope = 1, intercept = c(-log10(2.5)*2, sqrt(2.5)*2), col='grey') +
    geom_point(aes(alpha=FC, size=FC, col=Significant)) + 
    scale_x_log10() + scale_y_log10() +
    ggrepel::geom_label_repel(
        data = . %>% arrange(desc(FC)) %>% head(10), color='black'
    ) + theme_bw(base_size = 14) +
    ggtitle('hpl-2/WT ratios of exonin to intronic reads')

ggsave('hpl2WT_ratios_of_exonin_to_intronic_reads.pdf')

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



########

exons <- exonsBy(txdb, 'tx', use.names=TRUE)
SEexons <- summarizeOverlaps(exons, bam, ignore.strand = TRUE)
ddse <- DESeqDataSet(SEexons, ~1)
fe <- fpkm(ddse, robust = TRUE)


log2(fi[,3:4]/fi[,1:2]) %>% head
