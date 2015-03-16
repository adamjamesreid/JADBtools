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


gene <- rownames(top10[3,])
data <- plotCounts(ddsTC, gene, 
                   intgroup=c("stage","strain"), returnData=TRUE)
ggplot(data, aes(x=stage, y=count, color=strain, group=strain)) + 
    geom_point() + stat_smooth(se=FALSE,method="loess") +  scale_y_log10() +
    ggtitle(gene)

ggsave("length-hist.pdf")






