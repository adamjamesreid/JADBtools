ggplot() + 
    geom_density(aes(x=x), data=data.frame(x=rnbinom(1000, 2.4, 0.8)), fill="red", alpha=0.5) + 
    geom_density(aes(x=x), data=data.frame(x=rnbinom(1000, 3.5, 0.9)), fill="blue", alpha=0.5) +
    geom_density(aes(x=x), data=data.frame(x=rnbinom(1000, 2, 0.5)), colour="green", alpha=0, linetype=3)


require(lmtest)
require(MASS)

y <- c(rnbinom(1000, 2.4, 0.8), rnbinom(1000, 3.5, 0.9))
x <- factor(c(rep('g1', 1000), rep('g2', 1000)))
dat <- data.frame(x=x, y=y)

M <- glm.nb(y ~ x, data=dat)
anova(M)
waldtest(M)
lrtest(M)


y <- c(rnbinom(1000, 2.4, 0.8), rnbinom(1000, 2.4, 0.8))
x <- factor(c(rep('g1', 1000), rep('g2', 1000)))
dat <- data.frame(x=x, y=y)

M <- glm.nb(y ~ x, data=dat)
anova(M)
waldtest(M)
lrtest(M)



y <- c(rnbinom(10000, 2.2, 0.8), rnbinom(10000, 1, 0.5))
x <- factor(c(rep('g1', 10000), rep('g2', 10000)))
dat <- data.frame(x=x, y=y)

M <- glm.nb(y ~ x, data=dat)
anova(M)
waldtest(M)
lrtest(M)

require(aod)
mod <- negbin( y ~ x, ~ 1, dat)
rnd <- negbin( y ~ 1, ~ 1, dat)

aod::anova(mod, rnd)
wald.test(b = coef(mod), Sigma = vcov(mod), Terms=2)



#### SIMPLE math

## fit NB distr from a sample X using concentrated logLik
fitNB <- function(X) {
    n <- length(X)
    loglik.conc <- function(r) {
        prob <- n*r / (sum(X) + n*r)
        sum( lgamma(r + X) - lgamma(r) - lgamma(X + 1) +
                 r * log(prob) + X * log(1 - prob) ) 
    }
    ## find 'r' with an 1D optim...
    res <- optimize(f = loglik.conc, interval = c(0.001, 1000),
                    maximum = TRUE)
    r <- res$maximum[1]
    params <- c(size = r, prob = n*r / (sum(X) + n*r))
    attr(params, "logLik") <- res$objective[1]
    params
}
## compute score vector and info matrix at params 'psi' using closed forms
scoreAndInfo <- function(psi, X) {
    size <- psi[1]; prob <- psi[2]
    n <- length(X)
    U <- c(sum(digamma(size + X) - digamma(size) + log(prob)),  
           sum(size / prob - X / (1-prob) ))
    I <- matrix(c(- sum(trigamma(size + X) - trigamma(size)),  
                  -n / prob, -n / prob,  
                  sum( size / prob^2  + X / (1-prob)^2)),
                nrow = 2, ncol = 2)
    names(U) <- rownames(I) <- colnames(I) <- c("size", "prob")
    LM <-  as.numeric(t(U) %*% solve(I) %*% U)
    list(score = U, info = I, LM = LM)
}
## continuing on the question code a is for "all" 

c.dots <- rnbinom(1000, 2.4, 0.8)
w.dots <- rnbinom(1000, 3.5, 0.9)

c.fit <- fitNB(X = c.dots)
w.fit <- fitNB(X = w.dots)
a.fit <- fitNB(X = c(c.dots, w.dots))
## LR test and p.value
D.LR <- 2 * ( attr(c.fit, "logLik") + attr(w.fit, "logLik") ) -
    2 * attr(a.fit, "logLik") 
p.LR <- pchisq(D.LR, df = 2, lower.tail = FALSE)
## use restricted parameter estimate to compute the LM contributions
c.sI <- scoreAndInfo(psi = a.fit, X = c.dots) 
w.sI <- scoreAndInfo(psi = a.fit, X = w.dots) 
D.LM <- c.sI$LM + w.sI$LM 
p.LM <- pchisq(D.LM, df = 2, lower.tail = FALSE)


y <- c(c.dots, w.dots)
x <- factor(c(rep('g1', 1000), rep('g2', 1000)))
dat <- data.frame(x=x, y=y)

M <- glm.nb(y ~ x, data=dat)
anova(M)
waldtest(M)
lrtest(M)

require(aod)
mod <- negbin( y ~ x, ~ 1, dat)
rnd <- negbin( y ~ 1, ~ 1, dat)

aod::anova(mod, rnd)
wald.test(b = coef(mod), Sigma = vcov(mod), Terms=2)


### DEseq

dds <- DESeqDataSet(SE, design = ~ strain + stage + strain:stage )
dds0 <- DESeq(dds)

ddsTC <- DESeqDataSet(SE, design = ~ strain + stage + strain:stage )
ddsTC <- DESeq(ddsTC, test="LRT", reduced = ~ strain + stage )

resTC <- results(ddsTC)
resTC$symbol <- mcols(ddsTC)$geneName
top10 <- head(resTC[order(resTC$pvalue),],10)


gene <- rownames(top10[6,])


data <- plotCounts(ddsTC, gene, 
                   intgroup=c("stage","strain"), returnData=TRUE)
g <- ggplot(data, aes(x=stage, y=count, color=strain, group=strain)) + 
    geom_point() + stat_smooth(se=FALSE,method="loess") +  scale_y_log10() +
    ggtitle(gene)

ggsave("length-hist.pdf")



ddsTC0 <- DESeqDataSet(SE, design = ~ strain + stage )
ddsTC0 <- DESeq(ddsTC0, test="LRT", reduced = ~ 1 )

resTC0 <- results(ddsTC0)
resTC0$symbol <- mcols(ddsTC0)$geneName
top10 <- head(resTC0[order(resTC0$pvalue),],10)


gene <- rownames(top10[1,])


data <- plotCounts(ddsTC, gene, 
                   intgroup=c("stage","strain"), returnData=TRUE)
ggplot(data, aes(x=stage, y=count, color=strain, group=strain)) + 
    geom_point() + stat_smooth(se=FALSE,method="loess", size = 1.2) +  scale_y_log10() +
    ggtitle(gene)



message('Generating DESeq2 results')
des2Report <- HTMLReport(
    shortName = 'RNAseq_analysis_with_DESeq2',
    title = 'Differential expression using DESeq2, FDR less 0.05',
    reportDirectory = location
)
publish(
    object=top10, DataSet=ddsTC,
    des2Report,  pvalueCutoff=.1,
    .modifyDF = list(final, modifyReportDF),
    make.plots = TRUE,
    conditions=colData(dds)$stage
)
finish(des2Report)









message('Generating DESeq2 results')
des2Report <- HTMLReport(
    shortName = 'RNAseq_TS',
    title = "Time series analysis",
    reportDirectory = location
)
publish("LRT p-value: '~ strain + stage + strain:stage' vs '~ strain + stage'", des2Report)
publish(
    object=resTC[order(resTC$pvalue),], DataSet=ddsTC,
    des2Report,  pvalueCutoff=1,
    .modifyDF = list(final),
    make.plots = FALSE, n = Inf
)

deseqCSV <- CSVFile('RNAseq_TS', 'CSV', reportDirectory = location)
publish(as.data.frame(resTC[order(resTC$pvalue),]), deseqCSV)

publish(hwrite('Get as CSV', link = paste0(deseqCSV@shortName, '.csv')), des2Report) 
finish(des2Report)




res <- results(ddsTC, contrast=c('strain', 'N2', 'let418ts'), test="Wald")
title <- 'Wald test, all samples - strain: N2 vs let418ts'
location


######### report ##########
location='testReport2'




makeNewImages<-function(object,...){
    imagename <- c()
    for (i in 1:nrow(object)){
        imagename[i] <- paste0("plotNew", i, ".png")
        png(filename = paste0("", imagename[i]))
        plot(object$baseMean[i], ylab = "logfc", xlab = object$symbol[i],
        main = "New logfc Plot", col = "red", pch = 15, cex=3)
        dev.off()
    }
    object$Image <- hwriteImage(imagename, link = imagename, table = FALSE, height=150, width=150)
    return(object)
}
## DESeq2 results
message('Generating DESeq2 results')
des2Report <- HTMLReport(
    shortName = 'RNAseq_analysis_with_DESeq2',
    title = 'Differential expression using DESeq2, FDR less 0.05',
    reportDirectory = location
)
publish(
    object=ddsTC, 
    des2Report,  pvalueCutoff=0.000005,
    .modifyDF = list(final, modifyReportDF),
    make.plots = FALSE
)
finish(des2Report)


## DESeq2 GO results
message('Generating DESeq2 GO')

res <-  as.data.frame(results(ddsTC))
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






