#source('https://gist.githubusercontent.com/andrewparkermorgan/286f245ef5f5701a8b6e/raw/8b9eedcb294f1aebcce0d550736e86a89dee11ea/rnaseq.R')

#download.file('https://gist.githubusercontent.com/andrewparkermorgan/286f245ef5f5701a8b6e/raw/8b9eedcb294f1aebcce0d550736e86a89dee11ea/rnaseq.R', 'inst/spliceing/sashimiplot_fun.R')




de_exons_introns <- function(id) {
    require(BiocParallel)
    require(GenomicAlignments)
    bpparam()
    bam <- getFilePath(id, format = 'bam', processing = 'STAR', mount = TRUE)
    if(all(file.exists(unlist(bam)))) {
        bam <- BamFileList(unlist(bam))
    } else {
        jadb_download_files(IDS, processing = 'STAR', auto = TRUE)
        bam <- BamFileList(basename( unlist(bam) ))
    }
    
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

getDS <- function(se_lst) {
    #se_lst <- de_exons_introns(IDS)
    require(DESeq2)
    library(tibble)
    library(dplyr)
    library(magrittr)
    ddse <- DESeqDataSet(se_lst$SEexons, ~1)
    ddsi <- DESeqDataSet(se_lst$SEintron, ~1)
    
    
    
    fi <- fpkm(ddsi, robust = TRUE)
    mi <- cbind(
        int_mut=rowMeans(fi[,colData(ddsi)$strain!='N2']),
        int_wt=rowMeans(fi[,colData(ddsi)$strain=='N2'])
    ) %>% as.data.frame %>% rownames_to_column('wb') %>% tbl_df
    
    fe <- fpkm(ddse, robust = TRUE)
    me <- cbind(
        exn_mut=rowMeans(fe[,colData(ddsi)$strain!='N2']),
        exn_wt=rowMeans(fe[,colData(ddsi)$strain=='N2'])
    ) %>% as.data.frame %>% rownames_to_column('wb') %>% tbl_df
    
    # Filter missing val
    plot_data <- me %>% left_join(mi)
    plot_data <- plot_data[rowSums(is.na(plot_data[,-1]))==0,]
    plot_data <- plot_data[rowSums(plot_data[,-1]==0)==0,]
    
    # Calculate stats
    plot_data %<>% mutate(exon_ratio = exn_mut / exn_wt)
    plot_data %<>% mutate(intron_ratio= int_mut / int_wt)
    
    plot_data %<>% mutate(FC = abs(log2(exon_ratio/intron_ratio)) )
    plot_data %<>% mutate(Significant=FC > 2.5)
    plot_data <- left_join(plot_data, ANNO)
    
    #colnames(fe) <- paste0('exon_', colnames(fe))
    #fe %<>% as.data.frame %>% rownames_to_column('wb') %>% tbl_df
    #fe %<>% dplyr::mutate( exon_hpl2 = (exon_rAM051+exon_rAM052)/2 ) %>% dplyr::mutate(exon_wt=(exon_rAM047+exon_rAM048)/2)
    
    
    #colnames(fi) <- paste0('intron_', colnames(fi))
    #fi %<>% as.data.frame %>% rownames_to_column('wb') %>% tbl_df
    #fi %<>% dplyr::mutate( intron_hpl2 = (intron_rAM051+intron_rAM052)/2 ) %>% dplyr::mutate(intron_wt=(intron_rAM047+intron_rAM048)/2)
    
    
    
}


plt <- function(genename, IDS) {
    fls <- lapply(getFilePath(IDS, format = 'bam', processing = 'STAR', mount = TRUE), basename)
    names(fls) <- sapply(fls, strsplit, '_') %>% lapply('[', 4:3) %>% sapply(paste0, collapse='\n')
    
    gene <- txdb_gtf[txdb_gtf$gene_name==genename & !is.na(txdb_gtf$gene_name) ]
    genef <- reduce(c(flank(range(gene), 500,start = FALSE), range(gene),flank(range(gene), 500,start = TRUE)))
    LST <- sapply(fls, readGAlignments, param = ScanBamParam(which = genef))
    
    splot <- sashimiplot(
        aln = LST, tx = gene, min.splices = 2, smooth=5, splice.scale = c(-0.1,-50),
        colours = rev(c('red', 'red', 'grey', 'grey')), title = genename
    )
    
    ggsave( paste0(paste(IDS, collapse = '_'), '__', genename, '.pdf'),splot)
    return(splot)
}





require(rtracklayer)
source('./inst/spliceing/sashimiplot_fun.R')


txdb_fls <- file.path('/Users/przemol/Downloads/Caenorhabditis_elegans.WBcel235.90.gtf.gz')
txdb_gtf <- import.gff(txdb_fls, format='gtf')

txdb <- GenomicFeatures::makeTxDbFromGRanges(txdb_gtf[txdb_gtf$gene_biotype == 'protein_coding'])
ce11model <- GenomicFeatures::exonsBy(txdb, 'gene')
ce11genes <- GenomicFeatures::genes(txdb)

intron <- IRanges::setdiff(ce11genes, unlist(ce11model))
ovr <- findOverlaps(ce11genes, intron)
as.list(ovr) %>% str



require(pbapply)
#Takes ~2min on single core
intron_lst <- pblapply(as.list(ovr), function(x) intron[x])
intron_grl <- GRangesList(intron_lst)
names(intron_grl) <- names(ce11genes)
ce11introns <- intron_grl
# 
# 
elementMetadata(txdb_gtf) %>% as.data.frame %>% tbl_df %>%
    filter(gene_biotype == 'protein_coding', type=='exon') %>%
    mutate(seq_name=(gsub('^([A-Za-z0-9_]+\\.[0-9+]).+', '\\1', exon_id))) %>%
    select(wb=gene_id, gene_name, seq_name) %>% unique -> ANNO
# 
# 
# ce11model <- exonsBy(txdb, 'gene')
# ce11genes <- genes(txdb)

#set25met2 in YA
IDS <- c("rAM086","rAM085", "rAM058","rAM057")
ds <- de_exons_introns(IDS)
pd <- getDS(ds)
top25 <- pd %>% arrange(desc(FC)) %>% select(gene_name, FC, exon_ratio, intron_ratio  ) %>% head(25) #%>% mutate(DS='-          NO          -') %>% knitr::kable()

gg0 <- pd %>% 
    ggplot(aes(x=exon_ratio, y=intron_ratio, label=gene_name)) + 
    geom_abline(slope = 1, intercept = 0) +
    geom_abline(slope = 1, intercept = c(-log10(2.5)*2, log10(2.5)*2), col='grey') +
    geom_point(aes(alpha=FC, size=FC, col=Significant)) + 
    scale_x_log10() + scale_y_log10() +
    ggrepel::geom_label_repel(
        data = . %>% arrange(desc(FC)) %>% head(10), color='black'
    ) + theme_bw(base_size = 14) +
    ggtitle('set25met2 /WT ratios of exonin to intronic reads ->NEW<- ')

ggsave(paste0(paste(IDS, collapse = '_'), '_splice.pdf'), gg0)
suppressMessages(pblapply(top25$gene_name, plt, IDS))

require(biobroom)

de <- DESeq(DESeqDataSet(ds$SEexons, ~strain))
library("IHW")
ihwRes <- ihw(pvalue ~ baseMean,  data = deRes, alpha = 0.1)

sig <- tidy(results(de)) %>% arrange(p.adjusted) %>% left_join(mutate(ANNO, gene=wb), .) %>% filter(p.adjusted < 0.01)
pd %>% arrange(desc(FC)) %>% filter(wb %in% sig$gene)

pd %>% left_join(sig, .) %>%  arrange(desc(FC)) %>% filter(!is.na(FC)) %>% arrange(p.adjusted) %>% 
    select(gene_name, p.adjusted, FC) %>% filter(FC>4)

'plk-3' %>% plt(IDS) %>% plot




gene <- txdb_gtf[txdb_gtf$gene_name=='set-25' & !is.na(txdb_gtf$gene_name) ]

genef <- reduce(c(flank(range(gene), 500,start = FALSE), range(gene),flank(range(gene), 500,start = TRUE)))

LST <- sapply(c(
    N2_r1='/Users/przemol/BIGfiles/HCF/leafcutter/test/totalRNA_random_rAM057_N2_YA_STAR^NA^ce11^Sb862645.bam',
    N2_r2='/Users/przemol/BIGfiles/HCF/leafcutter/test/totalRNA_random_rAM058_N2_YA_STAR^NA^ce11^S4762644.bam',
    mt_r1='/Users/przemol/BIGfiles/HCF/leafcutter/test/totalRNA_random_rAM086_set25met2_YA_STAR^NA^ce11^S2762653.bam',
    mt_r2='/Users/przemol/BIGfiles/HCF/leafcutter/test/totalRNA_random_rAM085_set25met2_YA_STAR^NA^ce11^S4262652.bam'
), readGAlignments, param = ScanBamParam(which = genef))


   
    


ggsave('test.pdf', sashimiplot(
aln = LST, tx = gene, min.splices = 2, smooth=25, splice.scale = c(-0.1,-50),
colours = c('red', 'red', 'grey', 'grey')
))


require(plyr)
require(ggrepel)
#require(argyle)



plot(sashimiplot(
    aln = LST, tx = gene, min.splices = 2, smooth=5, splice.scale = c(-0.1,-50), 
    colours = c('darkred', 'darkred', 'darkgreen', 'darkgreen')
))


require(plyr)
#devtools::install_github("andrewparkermorgan/argyle")
require(argyle)

PLOT <- sashimiplot(aln = LST, tx = gene)

GA <- list(
    N2_r1=ga_N2_r1,
    N2_r2=ga_N2_r2,
    mt_r1=ga_mt_r1, 
    mt_r2=ga_mt_r2
)




meta = NULL; smooth = 0; min.coverage = 0; min.splices = 0; max.coverage = Inf; log.coverage = FALSE;
splice.scale = c(1.5e3, 5); colours = NULL; colour.by = NULL

aln = list(ga); tx = gene




