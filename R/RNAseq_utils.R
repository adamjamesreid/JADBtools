combineGGpca <- function(out, sep=TRUE) {


    a=plotPCA(out$dds, 'strain', ntop=Inf)
    b=plotPCA(out$dds, 'stage', ntop=Inf)
    
    if(sep) {
        
        require(gridExtra)
        grid.arrange(arrangeGrob(a, b))
        
    } else {
    
        m <- ncol(a$data); n <- m <- ncol(b$data);
        dat <- cbind(rbind(a$data[,1:3], b$data[,1:3]), g=c(rep('strain', m), rep('stage', n)))
        ggplot(dat, aes(PC1, PC2, color = group)) + facet_wrap(~g, nrow=1) + geom_point()
    }

}

getStageCmb <- function(out, name='test') {

    require(openxlsx)
    wb <- createWorkbook()
    
    
    diff <- list()
    cmb <- combn(levels(colData(out$dds)$stage), 2)
    
    
    for(i in 1:ncol(cmb)) {
    
        c1 <- cmb[1,i]
        c2 <- cmb[2,i]
        
        message(paste0(c1, '_vs_', c2))
        
        addWorksheet(wb, paste0(c1, '_vs_', c2))
        
        res <- results(out$dds, c("stage",c1,c2), addMLE = TRUE)
        to_wb <- as.data.frame(res)
        colnames(to_wb) <- elementMetadata(res)$description
        
        writeData( wb, paste0(c1, '_vs_', c2), to_wb  )
        write.csv( to_wb, paste0(name, '_', c1, '_vs_', c2, '.csv') )
    
      
        diff[[i]] <- to_wb  
        names(diff)[i] <- paste0(c1, '_vs_', c2)  
            
    }
    openxlsx::saveWorkbook(wb, paste0(name, ".xlsx"), overwrite = TRUE)

}

getStrainParwise <- function(out) {
    
    require(openxlsx)
    wb <- createWorkbook()
    
    for(stage in levels(colData(out$dds)$stage) ) {
        
        message(stage)
        
        s <- SE[,colData(SE)$stage == stage]
        o0 <- doDiffExpr(s)
        res <- results(o0$dds, addMLE = TRUE)
        to_wb <- as.data.frame(res)
        colnames(to_wb) <- elementMetadata(res)$description
        
        addWorksheet(wb, stage)
        writeData( wb, stage, to_wb  )
        write.csv( to_wb, paste0(stage, '_N2_vs_let418ts', '.csv') )
    }
    
    openxlsx::saveWorkbook(wb, "stage_N2_vs_let418ts.xlsx", overwrite = TRUE)

}

getStrainParwiseReport <- function(ddsTC, SE, n=100) {
    
    for(stage in levels(colData(ddsTC)$stage) ) {
        
        message(stage)
        s <- SE[,colData(SE)$stage == stage]
        o0 <- doDiffExpr(s)
        res <- results(o0$dds, addMLE = TRUE)
        subReport(res, ddsTC, paste0('Wald test, ', stage, ' - strain: N2 vs let418ts'),  location)
        
       
    }
    
}

plotSignal <- function(gene, ddsTC, res, basedir) {
    
    require('hwriter')
    
    IMG <- file.path('exprs_img',paste0(gene, ".png"))
    PDF <- file.path('exprs_img',paste0(gene, ".pdf"))
    
    if( !file.exists(file.path(basedir, PDF)) ) {
        data <- plotCounts(ddsTC, gene, 
                           intgroup=c("stage","strain"), returnData=TRUE)
        info <- sapply(mcols(ddsTC[gene,])[,1:2], as.character)
        #bm <- round( mcols(ddsTC[gene,])[,3], 2) 
        g <- ggplot(data, aes(x=stage, y=count, color=strain, group=strain)) + 
            geom_point() + scale_y_log10() 
        
        g1 <- g + stat_smooth(se=FALSE,method="loess", size = 2) 
        ggsave(file.path(basedir, IMG), g1, height=2, width=6, dpi=100)
        
        g2 <- g + stat_smooth(se=FALSE,method="loess", size = 1.1) +
            ggtitle(paste0(
                info[1],  ' (', gene, ')\n ', info[2] #, '\n FDR=', res[gene,]$padj, '; BM=', bm
            ))
        ggsave(file.path(basedir, PDF), g2, height=6, width=12)
        cat('.')
    }
    return(
        hwriteImage(IMG, link = PDF, table = FALSE, height=75, width=150)
    )
}

final <- function(df, html, object, factor, DataSet, ...){
    w <- rownames(df)
    ID <- paste0(
        '<a href="http://www.ensembl.org/Caenorhabditis_elegans/',
        'Gene/Summary?db=core;g=', w, '">', w,"</a>"
    )
    data(gnmodel)
    nam <-  elementMetadata(gnmodel[w,'geneName'])
    dir.create(file.path(html$reportDirectory, 'exprs_img'))
    Figures <- sapply(w, plotSignal, DataSet, object, html$reportDirectory)
    
    baseMean <- object[w,]$baseMean
    df <- cbind(ID, nam, baseMean, df[, 3:ncol(df)], Figures)
    rownames(df) <- w
    return(df)
}


final2 <- function(df, html, object, factor, DataSet, ...){
    w <- rownames(df)
    ID <- paste0(
        '<a href="http://www.ensembl.org/Caenorhabditis_elegans/',
        'Gene/Summary?db=core;g=', w, '">', w,"</a>"
    )
    data(gnmodel)
    nam <-  elementMetadata(gnmodel[w,'geneName'])
    dir.create(file.path(html$reportDirectory, 'exprs_img'))
    Figures <- sapply(w, plotSignal, DataSet, object, html$reportDirectory)
    
    baseMean <- object[w,]$baseMean
    df <- cbind(ID, nam, baseMean, df[, 2:ncol(df)], Figures)
    rownames(df) <- w
    return(df)
}


subReport <- function(res, ddsTC, title,  location, n=100) {
    DEreport <- HTMLReport(
        shortName = gsub("[^0-9a-zA-Z_]", '', gsub(' ', '_', title)),
        title = title,
        reportDirectory = location
    )
    
    publish(elementMetadata(res)[2,2], DEreport)
    publish(elementMetadata(res)[5,2], DEreport)
    
    publish(
        object=res[order(res$pvalue),], DataSet=ddsTC,
        DEreport,  pvalueCutoff=1,
        .modifyDF = list(final2),
        make.plots = FALSE, n=n
    )
    
    deseqCSV <- CSVFile(gsub("[^0-9a-zA-Z_]", '', gsub(' ', '_', title)), 'CSV', reportDirectory = location)
    publish(as.data.frame(res[order(res$padj),]), deseqCSV)
    
    publish(hwrite('Get as CSV', link = paste0(deseqCSV@shortName, '.csv')), DEreport) 
    finish(DEreport)
}












