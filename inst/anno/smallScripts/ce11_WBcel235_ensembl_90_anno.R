require(GenomicFeatures)
require(rtracklayer)
txdb_fls <- 'ftp://ftp.ensembl.org/pub/release-90/gtf/caenorhabditis_elegans/Caenorhabditis_elegans.WBcel235.90.gtf.gz'
txdb_gtf <- import.gff(txdb_fls, format='gtf')
txdb <- makeTxDbFromGRanges(txdb_gtf[txdb_gtf$gene_biotype == 'protein_coding'])

elementMetadata(txdb_gtf) %>% as.data.frame %>% tbl_df %>% 
    filter(gene_biotype == 'protein_coding', type=='exon') %>% 
    mutate(seq_name=(gsub('^([A-Za-z0-9]+\\.[0-9+]).+', '\\1', exon_id))) %>% 
    select(wb=gene_id, gene_name, seq_name) %>% unique -> ANNO

ce11model <- exonsBy(txdb, 'gene')
ce11gene <- genes(txdb)


pp <- promoters(ce11gene, upstream = 500, downstream = 0)
seqlevelsStyle(pp) <- 'UCSC'

wb_prom <- data_frame(wb=names(pp), promoter=paste(pp))
#wb_prom$promoter %>% GRanges()

tab <- read_tsv('https://raw.githubusercontent.com/jurgjn/relmapping/master/annot/S2_regulatory_annotation/S2b_promoter_annotation_6Dec17.tsv')
jj <- tab %>% makeGRangesFromDataFrame(keep.extra.columns = TRUE)
jj_prom <- data_frame(wb=jj$gene_id, jjpromoter=paste(jj))

all_prom <- left_join(wb_prom, jj_prom)
all_prom <- left_join(ANNO, all_prom)

all_prom %>% filter(!is.na(jjpromoter)) -> z
    distance(resize(GRanges(z$promoter), 1, fix = 'end'), GRanges(z$jjpromoter))
z$d <- distance(resize(GRanges(z$promoter), 1, fix = 'end'), GRanges(z$jjpromoter))

all_prom %>% pull(jjpromoter) %>% na.omit %>% GRanges()

get_signal <- function(prom, ID) {
    file <- BigWigFile(getFilePath(ID, processing = 'SQ10N', scale = 'lin', mount = TRUE))
    message('Getting signal for', basename(path(file)))
    gr <- prom %>% pull(promoter) %>% unique %>% GRanges()
    signal <- unlist(summary(file, gr))
    signal_tbl <- data_frame(promoter=paste(signal), x=signal$score)
    
    info <- rbeads:::ParseName(basename(path(file)))
    nam <- paste(info$Factor, info$Strain, ID, sep='_') 
    
    colnames(signal_tbl) <- c('promoter', nam)
    out <- left_join(prom, signal_tbl)
} 


data <- wb_prom %>% get_signal('RC009') %>% get_signal('RC010')
data %<>% get_signal('RC007') %>% get_signal('RC008')
data %<>% get_signal('AA673') %>% get_signal('AA467')
data %<>% get_signal('AA778') %>% get_signal('AA779')
data %<>% get_signal('AA780') %>% get_signal('AA781')

ja_cfp1_ex <- left_join(ja_cfp1, data)
ja_cfp1_ex %<>% mutate(upreg=padj<0.05 & log2FoldChange>0)
ja_cfp1_ex %<>% mutate(downreg=padj<0.05 & log2FoldChange<0)


ja_cfp1_ex$staus[ja_cfp1_ex$upreg] <- 'UP'
ja_cfp1_ex$staus[ja_cfp1_ex$downreg] <- 'DOWN'
ja_cfp1_ex$staus[is.na(ja_cfp1_ex$staus)] <- 'ns'

writexl::write_xlsx(ja_cfp1_ex, 'ja_cfp1_ex_hda1.xlsx')
slackr_upload('ja_cfp1_ex_hda1.xlsx', channels = '#figures')

gg0 <- ja_cfp1_ex %>% 
    filter(baseMean > 5, !is.na(staus)) %>% 
    select(staus, matches(".+_[A-Z]{2}[0-9]{3}")) %>% 
    gather(ex, value, starts_with('RC'), matches(".+[A-Z]{2}[0-9]{3}")) %>% 
        ggplot(aes(x=staus, y=value, fill=staus)) + scale_y_log10() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        theme(legend.position="bottom") +
        guides(fill=guide_legend(nrow=1,byrow=TRUE, title="Regulated: ")) 







gg1 <-  gg0 + geom_boxplot(notch = TRUE) + geom_violin(alpha=0.3)

 #+ geom_jitter(alpha=0.2) #
gg1 <- gg1 +  stat_summary(fun.data = function(x){return(c(y = min(x)-.2, label = length(x)))}, geom = "text")
gg2 <- gg1 + facet_wrap(~ex, nrow = 1, labeller = as_labeller(function(x) gsub('_', '\n', x)))
ggslackr(gg2, width = 16, height = 10)
