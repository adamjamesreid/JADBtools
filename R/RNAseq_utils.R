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











