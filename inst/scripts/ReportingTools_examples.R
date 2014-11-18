system('which bwa', wait=TRUE)


system( paste(system('which tophat2', intern=TRUE), '--help') )


str( arrange(filter(r, padj < 0.05, padj ) ))

dds <- ans$dds
require("org.Ce.eg.db")
require(ReportingTools)
require(biomaRt)

des2Report <- HTMLReport(
    shortName = 'RNAseq_analysis_with_DESeq2',
    title = 'RNA-seq analysis of differential expression using DESeq2',
    reportDirectory = "inst/tmp"
)
publish(
    object=dds, 
    des2Report,  pvalueCutoff=0.005,
    reportDir="inst/tmp",
    factor = colData(dds)$strain,
    .modifyDF = list(
        function(df, ...){
            w <- rownames(df)
            ID <- paste0(
                '<a href="http://www.ensembl.org/Caenorhabditis_elegans/Gene/',
                'Summary?db=core;g=', w, '">', w,"</a>"
            )
            data(gnmodel)
            nam <-  elementMetadata(gnmodel[w,'geneName'])
            df <- cbind(ID, nam, df[, 2:ncol(df)])
            rownames(df) <- w
            return(df)
        }, modifyReportDF
    ), #object = object,
    make.plots = TRUE
)
browseURL(finish(des2Report))


add2 <- 



add.anns <- function(df, ...)
{
    mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset="celegans_gene_ensembl", host="www.ensembl.org")
    nm <- rownames(df)
    anns <- getBM(
        attributes = c("ensembl_gene_id", "external_gene_name"), 
        filters = "ensembl_gene_id", values = nm, mart = mart)
    anns <- anns[match(nm, anns[, 1]), ]
    colnames(anns) <- c("ID", "Gene Symbol")
    w <- anns$ID
    anns$ID <- paste0('<a href="http://www.ensembl.org/Caenorhabditis_elegans/Gene/Summary?db=core;g=', w, '">', w,"</a>")
    df <- cbind(anns, df[, 2:ncol(df)])
    rownames(df) <- nm
    df
}
        
        
        rdb2 <- makeTranscriptDbFromBiomart(, 
                                            filters=list( "biotype"="protein_coding" ), biomart="ENSEMBL_MART_ENSEMBL", 
                                            host='www.ensembl.org')


##Go

library(GOstats)
library(org.Ce.eg.db)

res <-  results(dds)
res <- res[res$padj < 0.05 & !is.na(res$padj), ]
selectedIDs <- rownames(res)
selectedIDs <- mappedLkeys(org.Ce.egENSEMBL2EG[selectedIDs[selectedIDs %in% keys(org.Ce.egENSEMBL2EG)]])
universeIDs <- names(gnmodel)
universeIDs <- mappedLkeys(org.Ce.egENSEMBL2EG[universeIDs[universeIDs %in% keys(org.Ce.egENSEMBL2EG)]])

goParams <- new("GOHyperGParams", 
    geneIds = selectedIDs, 
    universeGeneIds = universeIDs, 
    annotation ="org.Ce.eg.db" , 
    ontology = "MF", 
    pvalueCutoff = 0.01,
    conditional = TRUE, 
    testDirection = "over")
goResults <- hyperGTest(goParams)



library(ReportingTools)
goReport <- HTMLReport(shortName = 'go_analysis_rnaseq',
    title = "GO analysis of mockRnaSeqData",
	reportDirectory = "./reports")
publish(goResults, goReport, selectedIDs=selectedIDs, annotation.db="org.Ce.eg.db", 
	pvalueCutoff= 0.05)
browseURL(finish(goReport))

##edger 

library(edgeR)
d <- DGEList(counts = counts(dds), group = colData(dds)$strain)
d <- calcNormFactors(d)
d <- estimateCommonDisp(d)
edgeR.de <- exactTest(d)



d <- calcNormFactors(d)
design <- model.matrix(~colData(dds)$strain)
d <- estimateGLMCommonDisp(d, design)
d <- estimateGLMTrendedDisp(d, design)
d <- estimateGLMTagwiseDisp(d, design)
fit <- glmFit(d,design)
edgeR.lrt <- glmLRT(fit, coef=2)

deReport <- HTMLReport(shortName = 'RNAseq_analysis_with_edgeR',
    title = 'RNA-seq analysis of differential expression using edgeR (LRT)',
    reportDirectory = "./reports")
publish(
    edgeR.lrt, deReport, countTable=counts(dds),
    annotation.db = NULL,
    conditions=colData(dds)$strain, 
    .modifyDF = list(
        function(df, ...){
            w <- rownames(df)
            ID <- paste0(
                '<a href="http://www.ensembl.org/Caenorhabditis_elegans/',
                'Gene/Summary?db=core;g=', w, '">', w,"</a>"
            )
            data(gnmodel)
            nam <-  elementMetadata(gnmodel[w,'geneName'])
            df <- cbind(ID, nam, df[, 2:ncol(df)])
            rownames(df) <- w
            return(df)
        }, modifyReportDF
    ), 
	pvalueCutoff = .05, name="edgeR"
)
browseURL(finish(deReport))


indexPage <- HTMLReport(shortName = "indexRNASeq",
    title = "Analysis of mockRnaSeqData",
    reportDirectory = "./reports")
publish(Link(list(deReport, goReport), report = indexPage), indexPage)
browseURL( finish(indexPage) )

#PFAM
library(Category)
params <- new("KEGGHyperGParams", 
     geneIds= selectedIDs, 
 	universeGeneIds=universeIDs, 
 	annotation="org.Ce.eg.db",
 	pvalueCutoff= 0.01,
 	testDirection="over")
KEGGResults <- hyperGTest(params)

KEGGReport <- HTMLReport(shortName = 'pfam_analysis_rnaseq',
    title = "PFAM analysis of mockRnaSeqData",
	reportDirectory = "./reports")
publish(KEGGResults, KEGGReport, selectedIDs=selectedIDs, annotation.db="org.Ce.eg.db",pvalueCutoff= 0.05)
finish(PFAMReport)

PFAMResults <- hyperGTest(params)

PFAMReport <- HTMLReport(shortName = 'pfam_analysis_rnaseq',
    title = "PFAM analysis of mockRnaSeqData",
	reportDirectory = "./reports")
publish(PFAMResults, PFAMReport, selectedIDs=selectedIDs, annotation.db="org.Mm.eg",categorySize=5)
finish(PFAMReport)


## ref 

data(mockRnaSeqData)

library(DESeq2)
conditions <- c(rep("case",3), rep("control", 3))
mockRna.dse <- DESeqDataSetFromMatrix(countData = mockRnaSeqData,
                         colData = as.data.frame(conditions), design = ~ conditions)
colData(mockRna.dse)$conditions <- factor(colData(mockRna.dse)$conditions, levels=c("control", "case"))
## ## Get a DESeqDataSet object
mockRna.dse <- DESeq(mockRna.dse)


des2Report <- HTMLReport(shortName = 'RNAseq_analysis_with_DESeq2',
    title = 'RNA-seq analysis of differential expression using DESeq2',
     reportDirectory = "./reports")
 publish(mockRna.dse,des2Report,
         .modifyDF = list(add.anns, modifyReportDF), mart = mart,
     reportDir="./reports")
 browseURL( finish(des2Report) )




des2Report <- HTMLReport(shortName = 'RNAseq_analysis_with_DESeq2',
                         title = 'RNA-seq analysis of differential expression using DESeq2',
                         reportDirectory = "./reports")
publish(y,des2Report, pvalueCutoff = 1, make.plots = FALSE)
browseURL( finish(des2Report) )



des2Report <- HTMLReport(shortName = 'RNAseq_analysis_with_DESeq2',
                         title = 'RNA-seq analysis of differential expression using DESeq2',
                         reportDirectory = "./reports")
publish(z,des2Report, pvalueCutoff = 1, make.plots = FALSE)
browseURL( finish(des2Report) )
### ALL ####


library(ReportingTools)
data(mockRnaSeqData)


###################################################
### code chunk number 2: run_edgeR (eval = FALSE)
###################################################
library(edgeR)
conditions <- c(rep("case",3), rep("control", 3))
d <- DGEList(counts = mockRnaSeqData, group = conditions)
d <- calcNormFactors(d)
d <- estimateCommonDisp(d)
## Get an edgeR object
edgeR.de <- exactTest(d)


###################################################
### code chunk number 3: edgeR_report (eval = FALSE)
###################################################
library(lattice)
rep.theme <- reporting.theme()
## Change symbol colors in plots
rep.theme$superpose.symbol$col <- c("blue", "red")
rep.theme$superpose.symbol$fill <- c("blue", "red")
lattice.options(default.theme = rep.theme)
## Publish a report of the top 10 genes with p-values < 0.05 and log-fold change > 2
## In this case, the plots contain the counts from mockRnaSeqData, which are not normalized.
## The publish function does not normalize counts for the countTable argument to allow for
## flexibility in plotting various units (e.g. RPKM instead of counts).

deReport <- HTMLReport(shortName = 'RNAseq_analysis_with_edgeR',
    title = 'RNA-seq analysis of differential expression using edgeR',
    reportDirectory = "./reports")
publish(edgeR.de, deReport, countTable=mockRnaSeqData,
    conditions=conditions, annotation.db = 'org.Mm.eg', 
	pvalueCutoff = .05, lfc = 2, n = 10, name="edgeR")
finish(deReport)

## If you would like to plot normalized counts, run the following commands instead:
mockRnaSeqData.norm <- d$pseudo.counts
publish(edgeR.de, deReport, mockRnaSeqData.norm, 
       conditions, annotation.db = 'org.Mm.eg', 
      pvalueCutoff = .05, lfc = 2, n = 10)
finish(deReport)




> rld <- rlog(dds)
> z=plotPCA(rld, intgroup = 'strain')
> plotPCA(rld, intgroup = 'strain')
> z=plotPCA(rld, intgroup = 'strain')
> plots <- HTMLReport(shortName = 'RNAseq_analysis_with_edgeR',
                      +                          title = 'RNA-seq analysis of differential expression using edgeR (LRT) - plots',
                      +                          reportDirectory = "./reports")
> publish(z, plots)
> browseURL(finish(plots))
