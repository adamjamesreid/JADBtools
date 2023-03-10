#' Run full differential expression pipeline
#' 
#' @param ID the vector of ContactExpIDs
#' @param mnt_point give the mount point for filesystem, use URLs if NULL, defaults to NULL
#'   
#' @return stats string
#' 
#' @author Przemyslaw Stempor
#' 
#' @family RNAseq
#' @export
#' 
ids2diffExpr <- function(ID, mnt_point=NULL) {
    fls <- sapply(ID, getFilePath, processing='aligned', format='bam', scale='NA', url=is.null(mnt_point), eq=F)
    if(!is.null(mnt_point)) fls <- file.path(mnt_point, fls)
    
    message('Sumarizing the experiments')
    e <- summarizeBAMs(fls)
    
    message('Getting annotations')
    colData(e)$strain <- factor( unlist(getStrain(colnames(e))) )
    
    message('Running DEseq2')
    res <- doDiffExpr(e, design = ~ strain)
    
    return( res )
    
}

#' Get counts from bam files using embended gene meodel
#' 
#' File URLs can be obtained like this:
#' fls <- sapply(c('rAM001', 'rAM002', 'rAM003', 'rAM004'), getFilePath, processing='aligned', format='bam')
#' 
#' @param fls file names vector
#'   
#' @return eset
#' 
#' @author Przemyslaw Stempor
#' 
#' @family RNAseq
#' @export
#' 
summarizeBAMs <-function(fls, model = 'gnmodel', ignore.strand=TRUE) {
    bfl <- BamFileList(fls)
    expset_repeats <- summarizeOverlaps( get(data(list = model)), bfl, ignore.strand=ignore.strand )
    return(expset_repeats)
    
}

#' Get strain from IDs
#' 
#' @param ContactExpID vector of IDs
#'   
#' @return eset
#' 
#' @author Przemyslaw Stempor
#' 
#' @family RNAseq
#' @export
#' 
getStrain <- function( ContactExpID, EXTABLE='labrnaseq'){
    
    con <- dbConnect(dbDriver(DRIVER), group = GROUP, default.file='~/.my.cnf')
    PK <- dbGetQuery(con, sprintf("SHOW INDEX FROM %s WHERE Key_name = 'PRIMARY'", gsub('view$', '', EXTABLE) ))[['Column_name']]
    
    strain <- sapply(ContactExpID, function(x) dbGetQuery(
        con, sprintf('SELECT %s FROM %s WHERE %s = "%s"', 'Strain', 'labrnaseq', PK, x)
    ))
    
    dbDisconnect(con)
    return(strain)
}

#' Get stage from IDs
#' 
#' @param ContactExpID vector of IDs
#'   
#' @return eset
#' 
#' @author Przemyslaw Stempor
#' 
#' @family RNAseq
#' @export
#' 
getStage <- function(ContactExpID, EXTABLE='labrnaseq'){
    
    con <- dbConnect(dbDriver(DRIVER), group = GROUP, default.file='~/.my.cnf')
    PK <- dbGetQuery(con, sprintf("SHOW INDEX FROM %s WHERE Key_name = 'PRIMARY'", gsub('view$', '', EXTABLE) ))[['Column_name']]
    
    strain <- sapply(ContactExpID, function(x) dbGetQuery(
        con, sprintf('SELECT %s FROM %s WHERE %s = "%s"', 'Stage', 'labrnaseq', PK, x)
    ))
    
    dbDisconnect(con)
    return(strain)
}

#' Get data filed from IDs
#' 
#' @param ContactExpID vector of IDs
#' @param field name of the nannotation field
#'   
#' @return eset
#' 
#' @author Przemyslaw Stempor
#' 
#' @family RNAseq
#' @export
#' 
getDBdataField <- function(ContactExpID, field='ContactExpID', EXTABLE='labrnaseq'){
    
    con <- dbConnect(dbDriver(DRIVER), group = GROUP, default.file='~/.my.cnf')
    PK <- dbGetQuery(con, sprintf("SHOW INDEX FROM %s WHERE Key_name = 'PRIMARY'", gsub('view$', '', EXTABLE) ))[['Column_name']]
    
    out <- sapply(ContactExpID, function(x) dbGetQuery(
        con, sprintf('SELECT %s FROM %s WHERE %s = "%s"', field, EXTABLE, PK, x)
    ))
    
    dbDisconnect(con)
    return(out)
}


#' Get SummarizedEperiments from IDs
#' 
#' @param ContactExpID vector of IDs
#'   
#' @return SummarizedEperiment
#' 
#' @author Przemyslaw Stempor
#' 
#' @family RNAseq
#' @export
#' 
getSummarizedEperiment <- function( ContactExpIDs="rAM022" ){
    out <- sapply(ContactExpIDs, function(x) {
        message(x)
        get(load(url(getFilePath(ID=x, processing='TagCounts'))))
    })
    return( do.call(cbind, out) )
}

#' Get spikeIN summary from IDs
#' 
#' @param ContactExpID vector of IDs
#'   
#' @return SummarizedEperiment
#' 
#' @author Przemyslaw Stempor
#' 
#' @family RNAseq
#' @export
#' 
spikeINsummary <- function( ContactExpIDs=paste0("rFB0", 26:32) ){
    se <- getSummarizedEperiment(ContactExpIDs) 
    rownames(se)[rownames(se) == ''] <- 
        as.character(elementMetadata(se)$geneName[rownames(se) == ''])
    colnames(se) <- sapply(strsplit(colnames(se), '_'), '[[', 10)
    
    
    
    
    
    
    
}

#' Get SummarizedEperiments as Rdata from IDs
#' 
#' @param ContactExpID vector of IDs
#'   
#' @return SummarizedEperiment
#' 
#' @author Przemyslaw Stempor
#' 
#' @family RNAseq
#' @export
#' 
downloadSummarizedEperiment <- function( ContactExpIDs="rAM022", dir='SE' ){
    dir.create(SE)
    out <- sapply(ContactExpIDs, function(x) {
        message("Downloading: ", dir)
        download.file(getFilePath(ID=x, processing='TagCounts'))
    })
    return( do.call(cbind, out) )
}



#' Run DEseq2 and output results
#' 
#' @param e eset
#' @param design experiment design
#'   
#' @return list
#' 
#' @author Przemyslaw Stempor
#' 
#' @family RNAseq
#' @export
#' 
doDiffExpr <- function(e, design = ~ strain) {
    
    dds <- DESeqDataSet(e, design = design )
    dds <- DESeq(dds)
    
    res <- results(dds)
    out <- as.data.frame(res)
    colnames(out) <- elementMetadata(res)$description
    
    M <- assays(e)$counts
    RPKM <- as.data.frame(
        M / ((sum(width(gnmodel)) / 1e3) %*% t(colSums(M) / 1e6))
    )
    colnames(RPKM) <- paste0('RPKM_', colnames(RPKM))
    
    out <- cbind(
        elementMetadata(gnmodel), 
        data.frame(exon_width_kb=(sum(width(gnmodel)) / 1e3)),
        RPKM, 
        out
    )         

    outname <- elementMetadata(res)$description[[2]]
    outname <- gsub(' ', '_', gsub('log2 fold change \\(MAP\\): ', '', outname))
    outname <- paste(paste(colnames(e), collapse='-'), outname, 'RPKM_and_DESeq2.csv', sep='_')
    write.csv(out, 'diff.csv')
    
    return(list(dds=dds, out=out, outname=outname))
    
}
#' Get DEseq2 redult from DB
#' 
#' @param ID id
#' @param path mmt path, return URL if null
#'   
#' @return list
#' 
#' @author Przemyslaw Stempor
#' 
#' @family RNAseq
#' @export
#' 
getDEseq2res <- function(ID, path=NULL) {
    if(is.null(path)) {
        addr <- getFilePath(
            ID, format='Rdata', processing='DESeq2', 
            scale='.', url=TRUE, eq=FALSE
        ) 
        return(get(load(url(addr))))
    } else {
        message("Loading Rdata, using local path ", path)
        addr <- getFilePath(
            ID, format='Rdata', processing='DESeq2', 
            scale='.', url=FALSE, eq=FALSE
        ) 
        addr <- gsub('files', path, addr)
        return(get(load(addr)))
    }  
}

#' Get DEseq2 report as HTML
#' 
#' @param ID id
#'   
#' @return file colation
#' 
#' @author Przemyslaw Stempor
#' 
#' @family RNAseq
#' @export
#' 
getDEreport <- function(ID, path=NULL, location='DEreport', show=TRUE) {
    
    message('Getting expression results')
    dds <- getDEseq2res(ID, path=path)
    getDEreportForDESeqDataSet(dds, location=location, show=show)
}

#' Get DEseq2 report as HTML
#' 
#' @param dds DESeqDataSet form DESeq2 package
#' @param location sub-directory to create the report
#' @param show should R open the report in the browser
#'   
#' @return file colation
#' 
#' @author Przemyslaw Stempor
#' 
#' @family RNAseq
#' @export
#' 
getDEreportForDESeqDataSet <- function(dds, location='DEreport', show=TRUE) {
    
    final <- function(df, ...){
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
    }
    
    ## DESeq2 results
    message('Generating DESeq2 results')
    des2Report <- HTMLReport(
        shortName = 'RNAseq_analysis_with_DESeq2',
        title = 'Differential expression using DESeq2, FDR less 0.05',
        reportDirectory = location
    )
    publish(
        object=dds, 
        des2Report,  pvalueCutoff=0.05,
        factor = colData(dds)$strain,
        .modifyDF = list(final, modifyReportDF),
        make.plots = TRUE
    )
    finish(des2Report)
    
    ## DESeq2 GO results
    message('Generating DESeq2 GO')
    
    res <-  as.data.frame(results(dds))
    deseqCSV <- CSVFile('deseqCSV', 'Full DESeq set as CSV file.', reportDirectory = location)
    publish(res[order(res$padj),], deseqCSV)
    
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
    
    message('Generating DESeq2 GO report')
    goReport <- HTMLReport(shortName = 'go_analysis_rnaseq',
                           title = "GO analysis for DESeq2 (differential expression FDR less 0.05, p-hyper less 0.01)",
                           reportDirectory = location)
    publish(goResults, goReport, selectedIDs=selectedIDs, annotation.db="org.Ce.eg.db", 
            pvalueCutoff= 0.01, make.plots=FALSE)
    finish(goReport)
    
    ## edgeR
    message('Generating edgeR report')
    d <- DGEList(counts = counts(dds), group = colData(dds)$strain)
    
    d <- calcNormFactors(d)
    design <- model.matrix(~colData(dds)$strain)
    d <- estimateGLMCommonDisp(d, design)
    d <- estimateGLMTrendedDisp(d, design)
    d <- estimateGLMTagwiseDisp(d, design)
    fit <- glmFit(d,design)
    edgeR.lrt <- glmLRT(fit, coef=2)
    
    edgeReport <- HTMLReport(shortName = 'RNAseq_analysis_with_edgeR',
                             title = 'Differential expression using edgeR likelihood ratio tests (LRT), FDR less 0.05',
                             reportDirectory = location)
    publish(
        edgeR.lrt, edgeReport, countTable=counts(dds),
        annotation.db = NULL,
        conditions=colData(dds)$strain, 
        .modifyDF = list(final, modifyReportDF), 
        pvalueCutoff = .05, name="edgeR"
    )
    finish(edgeReport)
    
    ## edgeR GO results
    
    res <-  topTags(edgeR.lrt, n = Inf, adjust.method = 'BH', sort.by = 'p.value')$table
    edgerCSV <- CSVFile('edgerCSV', 'Full edgeR set as CSV file.', reportDirectory = location)
    publish(res, edgerCSV)
    
    res <- res[res$FDR < 0.05 & !is.na(res$FDR), ]
    selectedIDs <- rownames(res)
    selectedIDs <- mappedLkeys(org.Ce.egENSEMBL2EG[selectedIDs[selectedIDs %in% keys(org.Ce.egENSEMBL2EG)]])
    
    goParams <- new("GOHyperGParams", 
                    geneIds = selectedIDs, 
                    universeGeneIds = universeIDs, 
                    annotation ="org.Ce.eg.db" , 
                    ontology = "MF", 
                    pvalueCutoff = 0.01,
                    conditional = TRUE, 
                    testDirection = "over")
    goResults <- hyperGTest(goParams)
    
    message('Generating edgeR GO report')
    goReportEdgeR <- HTMLReport(shortName = 'go_analysis_rnaseq_edgeR',
                                title = "GO analysis for edgeR (differential expression FDR less 0.05, p-hyper less 0.01)",
                                reportDirectory = location)
    publish(goResults, goReportEdgeR, selectedIDs=selectedIDs, annotation.db="org.Ce.eg.db", 
            pvalueCutoff= 0.01, make.plots=FALSE)
    finish(goReportEdgeR)
    
    #Plors
    rld <- rlog(dds)
    pca_plot <- plotPCA(rld, intgroup = 'strain')
    pca_plot_report <- HTMLReport(
        shortName = 'pca_plot',
        title = 'RNA-seq analysis of differential expression using DESeq2: principal component analysis (PCA) plot',
        reportDirectory = location
    )
    publish(pca_plot, pca_plot_report)
    finish(pca_plot_report)
    
    #Index
    
    indexPage <- HTMLReport(shortName = "indexRNASeq",
                            title = sprintf("Analysis of RNA-seq for %s experiment(s).", paste(colnames(dds), collapse=', ')),
                            reportDirectory = location)
    publish(Link(
        list(des2Report, goReport, pca_plot_report, edgeReport, goReportEdgeR), 
        report = indexPage), indexPage
    )
    
    publish(Link(
        c("", "CSV: Full DESeq diffreential expression results.", "CSV: Full edgeR diffreential expression results."), 
        c("#", paste0(location, c("/deseqCSV.csv", "/edgerCSV.csv"))),
        report = indexPage
    ), indexPage)
    
    
    
    addr <- finish(indexPage) 
        
    if(show) browseURL(addr)
    return(addr)
}



