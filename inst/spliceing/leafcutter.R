devtools::install_github("davidaknowles/leafcutter/leafcutter")
require(leafcutter)

download.file('ftp://ftp.ensembl.org/pub/release-90/gtf/caenorhabditis_elegans/Caenorhabditis_elegans.WBcel235.90.gtf.gz', 'Caenorhabditis_elegans.WBcel235.90.gtf.gz')
cmd_e <- '
../scripts/gtf_to_exons.R Caenorhabditis_elegans.WBcel235.90.gtf.gz exons.txt.gz
'
system(cmd_e)
ee <- read_delim('exons.txt.gz', delim = '\t')
write_delim(ee, 'ee.txt', delim='\t')


ids <- dir('/Users/przemol/BIGfiles/HCF/RNAseq/YA/bw/20degC') %>% strsplit('_') %>% do.call(rbind, .) %>% tbl_df %>% .$V3
ids2 <- dir('/Users/przemol/BIGfiles/HCF/RNAseq/YA/bw/15degC') %>% strsplit('_') %>% do.call(rbind, .) %>% tbl_df %>% .$V3
#ids <- c('rAM057',  'rAM058', 'rAM061',  'rAM062')
jadb_download_files(ids, 'STAR') 


nm <- dir(pattern = 'STAR.+bam$')
tab <- nm %>% strsplit('_') %>% do.call(rbind, .) %>% tbl_df %>% mutate(path=nm) %>% select(path, strain=V4) 


cmd <- '
source ~/.bash_profile
for bamfile in `ls *.bam`
do
echo Converting $bamfile to $bamfile.junc
sh ../scripts/bam2junc.sh $bamfile $bamfile.junc
echo $bamfile.junc >> test_juncfiles.txt
done
'
system(cmd)

test_strain <- 'set25met2'
write_delim(
    tab %>% dplyr::mutate(path=paste0(path, '.junc')) %>% 
         filter(strain=='N2'|strain==test) %>% dplyr::select(path), 
    'test.txt', delim = '\t', col_names = FALSE
)
cmd2 <- '
python ../clustering/leafcutter_cluster.py -j test.txt -m 10 -o clusters -l 10000
'
system(cmd2)
dd <- read_delim('clusters_perind_numers.counts.gz',  delim=' ', col_names = FALSE, skip = 1)

write_delim(
    tab %>% filter(strain=='N2'|strain==test_strain), 
    'test_diff_introns.txt', delim = ' ', col_names = FALSE
)
cmd3 <- '
../scripts/leafcutter_ds.R --exon_file=ee.txt -i 1 -g 1 -c 10 clusters_perind_numers.counts.gz test_diff_introns.txt
'
system(cmd3)

#system('../leafviz/gtf2leafcutter.pl -o ce11_90 Caenorhabditis_elegans.WBcel235.90.gtf.gz')
cmd_v <- sprintf('
../leafviz/prepare_results.R -f 0.9 -m %s %s_perind_numers.counts.gz %s_cluster_significance.txt %s_effect_sizes.txt ce11_90
', 'test_diff_introns.txt', 'clusters', 'leafcutter_ds', 'leafcutter_ds')
system(cmd_v)

out <- read_delim('leafcutter_ds_cluster_significance.txt', delim='\t')
out$cluster <- sub('chr', '', out$cluster)
write_delim(out , 'leafcutter_ds_cluster_significance.txt', delim='\t')


system('cd ../leafviz && pwd && ./run_leafviz.R ../test/leafviz.Rdata')


counts_file <- 'clusters_perind_numers.counts.gz'
cluster_significance_file <- 'leafcutter_ds_cluster_significance.txt'
effect.sizes.file <- 'leafcutter_ds_effect_sizes.txt'
annotation_code <- 'ce11_90'





out <- read_delim('leafcutter_ds_cluster_significance.txt', delim='\t')
out2 <- read_delim('leafcutter_ds_effect_sizes.txt', delim='\t')


dd <- read_delim('clusters_perind_numers.counts.gz',  delim=' ', col_names = FALSE, skip = 1)
strain <- strsplit(readLines('clusters_perind_numers.counts.gz', 1), ' ')[[1]] %>% strsplit('_') %>% do.call(rbind, .) %>% tbl_df %>% .[[4]] %>% factor


M <- dd %>% dplyr::select(-X1) %>% as.matrix()
rownames(M) <- dd$X1
colnames(M) <- strsplit(readLines('clusters_perind_numers.counts.gz', 1), ' ')[[1]] %>% strsplit('_') %>% do.call(rbind, .) %>% tbl_df %>% .[[3]]

map <- map_clusters_to_genes(get_intron_meta(rownames(M)), ee) %>% tbl_df %>% dplyr::select(cluster=clu, genes)
pmap <- get_intron_meta(rownames(M)) %>% tbl_df() %>% mutate(cluster=paste0(chr, ':', clu), pos=paste0(chr, ':', start, '-', end)) %>% dplyr::select(cluster,pos)

diff <- differential_splicing(
    M, 
    strain,  
    min_samples_per_intron =1, 
    min_samples_per_group =1, 
    min_coverage = 2
)

res <- cluster_results_table(diff) %>% tbl_df()
res %<>% tbl_df %>% arrange(p) %>% filter(p<0.05)
res %<>% left_join(map)
res %<>% left_join(pmap)
return(res)

test_strain='set25met2'
require(leafcutter)
require(dplyr)
require(magrittr)
get_ds <- function(test_strain) {
    message(test_strain)
    diff <- differential_splicing(
        M[,strain == 'N2' | strain == test_strain], 
        strain[strain == 'N2' | strain == test_strain],  
        min_samples_per_intron =1, 
        min_samples_per_group =1, 
        min_coverage = 2
    )

    res <- cluster_results_table(diff) %>% tbl_df()
    res %<>% tbl_df %>% arrange(p) %>% filter(p<0.05)
    res %<>% left_join(map)
    res %<>% left_join(pmap)
    return(res)
}
RES <- lapply(levels(strain)[-4], get_ds)
names(RES) <- levels(strain)[-4]

hpl2 <- res


## DEX
require(DEXSeq)
dxd <- DEXSeqDataSetFromSE( SEdex, design= ~ sample + exon + strain:exon )
dxd <- estimateSizeFactors( dxd )
dxd <- estimateDispersions( dxd )

plotDispEsts( dxd )
dxd <- testForDEU( dxd )
dxd <- estimateExonFoldChanges( dxd, fitExpToVar="strain")
dxr1 <- DEXSeqResults( dxd )
dxr1 %>% as.data.frame %>% tbl_df %>% arrange(padj) %>% filter(padj<0.05) -> res

txdb_gtf <- import.gff('Caenorhabditis_elegans.WBcel235.90.gtf.gz', format='gtf')
txdb <- makeTxDbFromGRanges(txdb_gtf[txdb_gtf$gene_biotype == 'protein_coding'])
elementMetadata(txdb_gtf) %>% as.data.frame %>% tbl_df %>% 
    filter(gene_biotype == 'protein_coding', type=='exon') %>% 
    mutate(seq_name=(gsub('^([A-Za-z0-9]+\\.[0-9+]).+', '\\1', exon_id))) %>% 
    dplyr::select(wb=gene_id, gene_name, seq_name) %>% unique -> ANNO

resa <- right_join(ANNO %>% rename(groupID=wb), res)

plotDEXSeq( dxr1, "WBGene00001996", legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2,fitExpToVar="strain" )

plotDEXSeq( dxr1, "WBGene00006876", FDR=0.0000001, legend=TRUE, fitExpToVar="strain", displayTranscripts=TRUE, expression=TRUE, splicing=TRUE)

require(DEXSeq)
BPPARAM = MultiCoreParam(workers=4)
dxd = estimateSizeFactors( dxd )
dxd = estimateDispersions( dxd, BPPARAM=BPPARAM)
dxd = testForDEU( dxd, BPPARAM=BPPARAM)
dxd = estimateExonFoldChanges(dxd, BPPARAM=BPPARAM)


