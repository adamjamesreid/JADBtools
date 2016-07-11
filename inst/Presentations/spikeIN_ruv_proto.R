se <- getSummarizedEperiment(ContactExpIDs) 
rownames(se)[rownames(se) == ''] <- 
    as.character(elementMetadata(se)$geneName[rownames(se) == ''])
colnames(se) <- sapply(strsplit(colnames(se), '_'), '[[', 10)



tab <- read.table(system.file('spikeIN/cms_095046.txt', package = 'JADBtools'), sep = '\t', header = TRUE, row.names = 2)
attomole <- 602200
infoIN <- unlist(getDBdataField(paste0("rFB0", 26:32), field = 'SpikeIN'))

require(DESeq2)
require(magrittr)
dds <- DESeqDataSet(se, ~strain)
plotPCA( 
    DESeqTransform( se ), intgroup = c("strain", "stage"), 
    ntop=Inf, returnData = TRUE 
)  %>% ggplot(., aes(PC1, PC2, col=strain, label=name, shape=stage))  + 
    geom_text(vjust=1.2) + geom_point(size=3)




require(reshape2)
rpkm <- fpkm(dds, robust = FALSE)
count <- melt(
    rpkm[grep('ERCC', rownames(rpkm)),], 
    varnames = c('spike','smpl'), value.name = 'count'
)

tab <- tab[rownames(ercc_rphm),]
cont <- melt( 
    tab[,grepl('2', infoIN)+3]*attomole,
    varnames = c('mix', 'spike'),value.name = 'cont'
)
data <- cbind( count, cont )
data  %>% mutate(cont=log2(cont), count=log2(count))  %>% 
    filter(!is.infinite(count)) -> data_log


require(ggplot2)
ggplot(
    data, 
    aes(cont, count, colour = smpl)
) + geom_point() + scale_y_continuous(trans=log2_trans()) + 
    scale_x_continuous(trans=log2_trans())



plotSpikeIN <- function(dd, xpos = 23, ypos=-10) {
    library(proto)
    stat_smooth_func <- function (mapping = NULL, data = NULL, geom = "smooth", position = "identity",
                                  method = "auto", formula = y ~ x, se = TRUE, n = 80, fullrange = FALSE,
                                  level = 0.95, na.rm = FALSE, xpos=NULL, ypos=NULL, ...) {
        StatSmoothFunc$new(mapping = mapping, data = data, geom = geom, position = position,
                           method = method, formula = formula, se = se, n = n, fullrange = fullrange,
                           level = level, na.rm = na.rm, xpos=xpos, ypos=ypos, ...)
    }
    
    StatSmoothFunc <- proto(ggplot2:::Stat, {
        objname <- "smooth"
        
        calculate_groups <- function(., data, scales, method="auto", formula=y~x, ...) {
            rows <- daply(data, .(group), function(df) length(unique(df$x)))
            
            if (all(rows == 1) && length(rows) > 1) {
                message("geom_smooth: Only one unique x value each group.",
                        "Maybe you want aes(group = 1)?")
                return(data.frame())
            }
            
            # Figure out what type of smoothing to do: loess for small datasets,
            # gam with a cubic regression basis for large data
            # This is based on the size of the _largest_ group.
            if (identical(method, "auto")) {
                groups <- count(data, "group")
                
                if (max(groups$freq) < 1000) {
                    method <- "loess"
                    message('geom_smooth: method="auto" and size of largest group is <1000,',
                            ' so using loess.',
                            ' Use \'method = x\' to change the smoothing method.')
                } else {
                    method <- "gam"
                    formula <- y ~ s(x, bs = "cs")
                    message('geom_smooth: method="auto" and size of largest group is >=1000,',
                            ' so using gam with formula: y ~ s(x, bs = "cs").',
                            ' Use \'method = x\' to change the smoothing method.')
                }
            }
            if (identical(method, "gam")) try_require("mgcv")
            
            .super$calculate_groups(., data, scales, method = method, formula = formula, ...)
        }
        
        calculate <- function(., data, scales, method="auto", formula=y~x, se = TRUE, n=80, fullrange=FALSE, xseq = NULL, level=0.95, na.rm = FALSE, xpos=xpos, ypos=ypos, ...) {
            data <- remove_missing(data, na.rm, c("x", "y"), name="stat_smooth")
            if (length(unique(data$x)) < 2) {
                # Not enough data to perform fit
                return(data.frame())
            }
            
            if (is.null(data$weight)) data$weight <- 1
            
            if (is.null(xseq)) {
                if (is.integer(data$x)) {
                    if (fullrange) {
                        xseq <- scale_dimension(scales$x, c(0, 0))
                    } else {
                        xseq <- sort(unique(data$x))
                    }
                } else {
                    if (fullrange) {
                        range <- scale_dimension(scales$x, c(0, 0))
                    } else {
                        range <- range(data$x, na.rm=TRUE)
                    }
                    xseq <- seq(range[1], range[2], length=n)
                }
            }
            if (is.character(method)) method <- match.fun(method)
            
            method.special <- function(...)
                method(formula, data=data, weights=weight, ...)
            model <- safe.call(method.special, list(...), names(formals(method)))
            
            predictdf(model, xseq, se, level)
            m = model
            eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                             list(a = format(coef(m)[1], digits = 3), 
                                  b = format(coef(m)[2], digits = 3), 
                                  r2 = format(summary(m)$r.squared, digits = 3)))
            func_string = as.character(as.expression(eq))
            
            data.frame(
                x=if(is.null(xpos)) min(data$x)*0.9 else xpos, 
                y=if(is.null(ypos)) max(data$y)*0.9 else ypos, 
                label=func_string
            )
        }
        
        required_aes <- c("x", "y")
        default_geom <- function(.) GeomSmooth
    })
    p <-  ggplot(dd, aes(cont, count, colour = smpl)) 
    p <- p + geom_point() 
    p <- p + geom_smooth(method=lm) 
    p <- p + facet_wrap(~smpl, nrow=2, scales="free")
    p <- p + xlab("RNA spike-in concentartion [RNA particles/ul]")
    p <- p + ylab("RNA expression measurement [RPKM]")
    #p + stat_smooth_func(geom="text",method="lm",hjust=-0.15,parse=TRUE, colour='black', size=4, vjust=11.5) 
    P <- p + facet_wrap(~smpl, nrow=2) + stat_smooth_func(geom="text",method="lm",parse=TRUE, colour='black', size=4, xpos = xpos, ypos=ypos) 
    show(P)
}



library(plyr)
data_log  %>% filter(count > 0)  %>% ddply(.(smpl), function(df) {
    m <- lm(count ~ cont, data=df)
    data.frame(a = coef(m)[1], b = coef(m)[2], pval=anova(m)$'Pr(>F)'[[1]] )
}) -> coefs

str(coefs)

counts(dds, normalized = FALSE, replaced = FALSE)[grep('ERCC', rownames(dds)),]  %>% colSums -> spike_sums

library(dplyr)
rpkm_vec <- unlist(as.data.frame(rpkm))
smpl <- factor(sapply(strsplit(names(rpkm_vec), '_'), '[[', 4))
dat <- data.frame(x=log2(rpkm_vec), smpl=smpl)
p <- ggplot(dat, aes(smpl, x, col=smpl))
p <- p +  geom_boxplot(alpha=0.7)  + geom_violin(alpha=0.5)
p + ylab('Expression [log2(RPKM)]') + xlab('Experiment')

spikes <- grep("^ERCC", rownames(rpkm), value = TRUE)
infoIN  %>%  sapply(rep, spikes  %>% length)  %>% melt  %>% 
    dplyr::select(value)  %>% rename(Mix=value)  %>% 
    cbind(melt(ercc_rphm[spikes,]), subgroup=tab[spikes,]$subgroup, .)  %>% 
    ggplot(aes(x=Var1, y=log2(value), col=subgroup, shape=Mix)) + 
    geom_point(size=4) + scale_shape_manual(values = c(20, 2)) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1))

infoIN  %>%  sapply(rep, spikes  %>% length)  %>% melt  %>% 
    dplyr::select(value)  %>% rename(Mix=value)  %>% 
    cbind(melt(ercc_rphm[spikes,]), subgroup=tab[spikes,]$subgroup, ., rat=tab[spikes,]$log2.Mix.1.Mix.2.)  %>% 
    ggplot(aes(x=rat, y=log2(value), col=subgroup, shape=Mix)) + 
    geom_violin(alpha=0.5) + geom_boxplot(alpha=0.5) +
    geom_jitter(size=4) + scale_shape_manual(values = c(20, 2)) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1))

#RUV
zfGenes <- assays(se)$counts
zfGenes[spikes,infoIN=="Spike-in 2"] <- zfGenes[spikes,infoIN=="Spike-in 2"]*tab[spikes,]$expected.fold.change.ratio

filter <- apply(zfGenes, 1, function(x) length(x[x>5])>=2)
filtered <- zfGenes[filter,]
genes <- rownames(filtered)[grep("^WB", rownames(filtered))]
spikes <- rownames(filtered)[grep("^ERCC", rownames(filtered))]

#colnames(zfGenes) <- sapply(strsplit(colnames(zfGenes), '_'), '[[', 10)

heatmap.2(log2(zfGenes[spikes,]+1), ColSideColors = rainbow(2)[colData(se)$strain], dendrogram = 'column', trace = 'none')
text(-0.6, 0, 'Log2 space\nafter spike-in\nconcentrations\nnormalization:\ncluster by Stage',)

filter <- apply(zfGenes, 1, function(x) length(x[x>5])>=2)
filtered <- zfGenes[filter,]
genes <- rownames(filtered)[grep("^WB", rownames(filtered))]
spikes <- rownames(filtered)[grep("^ERCC", rownames(filtered))]

ruv <- RUVg(round(zfGenes), spikes, k=1)

set <- newSeqExpressionSet(as.matrix(filtered), phenoData = data.frame(colData(se)[,1:2], row.names=colnames(filtered)))
set <- newSeqExpressionSet(as.matrix(filtered), phenoData = data.frame(colData(se)[,2], row.names=colnames(filtered)))


library(RColorBrewer)
colors <- brewer.pal(3, "Set2")
#plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors[pData(set)[[1]]])
plotPCA(set, col=colors[pData(set)[[1]]], cex=1.2)

# set <- betweenLaneNormalization(set, which="upper")
# plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors[x], main='UC norm')
# plotPCA(set, col=colors[x], cex=1.2, main='UC norm')

set1 <- RUVg(set, spikes, k=1); #pData(set1)
#plotRLE(set1, outline=FALSE, ylim=c(-4, 4), col=colors[pData(set)[[1]]], main='RUVg with spikes')
plotPCA(set1, col=colors[pData(set)[[1]]], cex=1.2, main='RUVg with spikes')

design <- model.matrix(~strain , data=pData(set1))

y <- DGEList(counts=normCounts(set1), group=pData(set)[[1]])
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
fit <- glmFit(y, design)
lrtN <- glmLRT(fit, coef=2)
plotSmear(lrtN, smooth.scatter=TRUE, lowess=TRUE, normalize=FALSE, de.tags = spikes)

y <- DGEList(counts=counts(set), group=pData(set)[[1]])
y <- calcNormFactors(y, method="upperquartile")
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef=2)
plotSmear(lrt, smooth.scatter=TRUE, lowess=TRUE, normalize=FALSE, de.tags = spikes)







ggplot(melt(zfGenes[spikes,]), aes(x=Var1, y=log2(value), col=Var2)) + geom_point()
