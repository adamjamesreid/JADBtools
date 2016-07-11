require(DESeq2)
re



nn <- load('/Volumes/raid0/_PROJECTS/repeats/Repeat_polyA_vs_longCAP_tags.Rdata')
dds <- DESeqDataSet(SEuniq, design = ~ stage )



#table <- assays(SEuniq)$counts
table <-fpkm(dds, robust = FALSE)

nam <- make.unique(rowRanges(SEuniq)$seqID)
polyA_AM_EE <- rowMeans(table[,c('rAM007', 'rAM008')])
polyA_FB_EE <- rowMeans(table[,c('rFB001', 'rFB013', 'rFB026', 'rFB027')])
LongCappedRNA_YD_EMB <- rowMeans(table[,c('rYD004', 'rYD008', 'rYD017')])
pos <- paste0(seqnames(rowRanges(SEuniq)), ':', start(rowRanges(SEuniq)), '-', end(rowRanges(SEuniq)))
tic <- read.csv('/Volumes/raid0/_PROJECTS/repeats/Supp_TableS2 - Figure 1 - Source Data 2.xls.csv')
t <- read.csv('/Volumes/raid0/_PROJECTS/repeats/Supp_TableS2 - Figure 1 - Source Data 2.xls.csv')
t <- t[!is.na(t$mode_position),]
tic <- GRanges(paste0('chr',t$chr), IRanges(t$mode_position, width=1), strand = t$strand)
tic1k <- resize(tic, 1000, fix='center')
haveTICplusMinus500 <- rowRanges(SEuniq)  %over% tic1k 

out <- data.frame(polyA_AM_EE, polyA_FB_EE, LongCappedRNA_YD_EMB, haveTICplusMinus500)

out_anno <- data.frame(name=rowRanges(SEuniq)$seqID, unique_name=nam, pos=pos, out)
write.csv(out_anno, 'repeats_expression_emb_N2_plusTIC_RPKM.csv')

out_classes <- sapply(unique(out_anno$name), function(x) {
    val <- colSums(out_anno[out_anno$name==x,4:7])
    len <- nrow(out_anno[out_anno$name==x,])
    c(val, len)
})

out_classes <- t(out_classes)
row.names(out_classes) <- unique(out_anno$name)

write.csv(out_anno, 'classes_repeats_expression_emb_N2_plusTIC_RPKM.csv')



p <- ggplot(out, aes(factor(haveTICplusMinus500), polyA_AM_EE))
p + geom_violin()