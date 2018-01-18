##### coding genes from Ensembl #####
txdb_fls <- 'ftp://ftp.ensembl.org/pub/release-90/gtf/caenorhabditis_elegans/Caenorhabditis_elegans.WBcel235.90.gtf.gz'
txdb_gtf <- import.gff(txdb_fls, format='gtf')
txdb <- GenomicFeatures::makeTxDbFromGRanges(txdb_gtf[txdb_gtf$gene_biotype == 'protein_coding'])
wb <- GenomicFeatures::genes(txdb)

anno_add(
    wb, Type = 'Genes', Version = 'Ensembl90', 
    Description = 'protein_coding genes from Ensembl, identified by WB id',
    Source = txdb_fls
)

##### coding promoters from Ensembl #####
prom <- promoters(wb, 500, 500)
anno_add(
    prom, Type = 'PromotersWB', Version = 'Ensembl90', 
    Description = 'promoters +/-500 TSS for protein_coding Ensembl genes',
    Source = txdb_fls
)

#####  Dfam 2.0 repeats #####
reptest_ce11 <- function() {
    data("repeatModel", package='JADBtools')
    require(rtracklayer)
    require(magrittr)
    
    chain <- import.chain('/Users/przemol/Downloads/ce10ToCe11.over.chain')
    
    rep <- unlist(repeatModel)
    rep$name <- rep$id
    rep$all_text <- NULL
    rep$description <- NULL
    rep$id <- NULL
    rep$synonym_and_previous_id <- NULL
    
    rep$classification_tags -> tg
    tg %<>% gsub('DNA Transposon', 'DNA_Transposon', .)
    tg %<>% gsub('Cut and Paste', 'Cut_and_Paste', .)
    tg %<>% gsub('Rolling Circle', 'Rolling_Circle', .)
    tg %<>% gsub('Gene', '', .)
    tg %<>% gsub('  ', ' ', .)
    
    tgt <- tg %>% strsplit(' ') %>% do.call(rbind, .)
    colnames(tgt) <- c('type', 'class', 'superfamily')
    tgt[tgt[,'class'] == 'Undefined','class'] <- 'ncRNA'
    
    elementMetadata(rep) <- cbind(elementMetadata(rep), as.data.frame(tgt))
    rep$classification_tags <- NULL
    
    rep_ce11 <- liftOver(rep, chain) %>% reduce(min.gapwidth=50) %>%  unlist
    elementMetadata(rep_ce11) <- elementMetadata(rep)
    
    rep_ce11_s <- sort(rep_ce11, ignore.strand=T)
    rep_ce11_s$score <- as.integer(rep_ce11_s$score)
    rep_ce11_s$score[rep_ce11_s$score>1000] <- 1000
    
    
    cmap <- c(Cut_and_Paste='#8A0000', Helitron='#FEC0CB', LTR='#2570B5', LINE='#6DADD6', SINE='#C6DBEF', ncRNA='#FFFFFF', Unknown='#BEBEBE', Satellite='#000000')
    rep_ce11_s$itemRgb <- cmap[as.character(rep_ce11_s$class)]
    return(rep_ce11_s)
    
    #export.bed(rep_ce11_s, '~/BIGfiles/rep_dfam20_ce11_class.sorted.bed')
    #system('/Users/przemol/miniconda2/bin/bedToBigBed ~/BIGfiles/rep_dfam20_ce11_class.sorted.bed ~/BIGfiles/ce11 ~/BIGfiles/rep_dfam20_ce11_class.bb')
    

    
}
dfam_reps_ce11 <- reptest_ce11()

anno_add(
    dfam_reps_ce11, Type = 'Repeats', Version = 'Dfam2.0', 
    Description = 'All repeats in Dfam 2.0',
    Source = 'http://www.dfam.org/web_download/Release/Dfam_2.0/'
)