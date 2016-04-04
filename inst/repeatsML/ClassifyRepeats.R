fls <- dir()


data <- lapply(fls, function(f) {
    message(f)
    bwf <- BigWigFile(f)
    vec <- extarct_vector(bwf, which = unlist(repeatModel), size = 1L)
    return(vec)
})
m <- do.call(cbind, data)
M <- as.data.frame(m) 

M <- cbind(M, class=unlist(repeatModel)$id)
M <- M[rowMeans(M[,1:9]) >= 1, ]

n <- fls  %>% strsplit('_')  %>% sapply('[[', 1)
nn <- fls  %>% strsplit('_')  %>% sapply('[[', 7)
nam <- paste(n, nn, sep='_')
colnames(M)[1:9] <- nam

low <- names(table(M$class)[table(M$class) < 10])
M <- M[!(M$class %in% low),]
M$class <- droplevels(M$class)

smp_size <- floor(0.75 * nrow(M))
train_ind <- sample(seq_len(nrow(M)), size = smp_size)

train <- M[train_ind, ]
test <- M[-train_ind, ]

Z <- as.data.frame( t(scale(  t(M[,1:9])  )) )
Z$class <- M$class

train <-Z[train_ind, ]
test <- Z[-train_ind, ]


#Logistic gregression
train <- M[train_ind, ]
test <- M[-train_ind, ]
train$class <- factor(train$class == 'CELE1') 
test$class <- factor(test$class == 'CELE1')
model <- glm(class ~., family=binomial(link='logit'),data=train)
fitted.results <- predict(model,newdata=test, type='response')
fitted.results <- ifelse(fitted.results > 0.5,1,0)


library(ROCR)
p <- predict(model, newdata=test[,-10], type="response")
pr <- prediction(p, test$class)
prf <- performance(pr, measure = "tpr", x.measure = "fpr")
plot(prf, main="")
text(0.70, 0.25, paste0('AUC = ',round(performance(pr, measure = "auc")@y.values[[1]], 3)), cex=3)
abline(0, 1, col='gray')
title(main = 'CELE1 vs. ALL')
anova(model, test="Chisq")


require('e1071')
svm.model <- svm(class ~., data = train)
svm.pred <- predict(svm.model, test, type = "response")

pr <- prediction(as.numeric(svm.pred), as.numeric(test$class))
prf <- performance(pr, measure = "tpr", x.measure = "fpr")
plot(prf, main="")
text(0.70, 0.25, paste0('AUC = ',round(performance(pr, measure = "auc")@y.values[[1]], 3)), cex=3)
abline(0, 1, col='gray')


summary(model)


model <- glm(class ~., family=quasipoisson(link = "log"),data=train)
summary(model)


fitted.results <- predict(model,newdata=test, type = "terms")
fitted.results <- ifelse(fitted.results > 0.5,1,0)

misClasificError <- mean(fitted.results != test$Survived)
print(paste('Accuracy',1-misClasificError))

library(nnet)
mod <- multinom(class ~., train, MaxNWts = 10000)
mod.p <- predict(rp.iris, type = "class", newdata = test)

sum(mod.p == test$class)

library(rpart)

trainp <- train[ train$class  %in% predictable, ]
rp.iris <- rpart(class ~., train)
rpart.plot(rp.iris)

pred1 <- predict(rp.iris, type = "class")
pred1  %>% table  %>% .[.>0]  %>% names -> predictable

testp <- test[ test$class  %in% predictable, ]
pred1 <- predict(rp.iris, type = "class", newdata = testp)
pred1  %>% table  %>% .[.>0]  %>% names -> predictable

sum(testp$class == pred1)/nrow(testp)

cm = as.matrix(table(Actual = droplevels(testp$class), Predicted = droplevels(pred1)))

head(pred1)



library(nnet)
mod <- multinom(class ~., train)





rowMeans(M[,1:9]) < 1




unlist(repeatModel)$id  %>%  head



########### VIS #########


plot(prcomp(M[, 1:9])$rotation[, 1:2])
text(prcomp(M[, 1:9])$rotation[, 1:2], colnames(M)[1:9], cex=.5)

plot(prcomp(t(M[, 1:9]))$rotation[, 1:2], col=M$class)
plot(prcomp(M[, 1:9])$x[, 1:2], col=M$class)

########## H20 ############

library(h2o)
h2o.init(nthreads=-1)
D <- as.h2o(M)
h2o.summary(D)

data = h2o.splitFrame(D,ratios=c(.7,.15),destination_frames = c("train","test","valid"))
names(data) <- c('Train', 'Test', 'Valid')


response <- "class"
predictors <- setdiff(names(D), response)

m1 <- h2o.deeplearning(
    model_id="dl_model_first", 
    training_frame=data$Train, 
    validation_frame=data$Valid,   ## validation dataset: used for scoring and early stopping
    x=predictors,
    y=response,
    #activation="Rectifier",  ## default
    #hidden=c(200,200),       ## default: 2 hidden layers with 200 neurons each
    epochs=1,
    variable_importances=TRUE    ## not enabled by default
)
summary(m1)

m2 <- h2o.deeplearning(
    model_id="dl_model_faster", 
    training_frame=data$Train, 
    validation_frame=data$Valid,
    x=predictors,
    y=response,
    hidden=c(32,32,32),                  ## small network, runs faster
    epochs=1000000,                      ## hopefully converges earlier...
    score_validation_samples=10000,      ## sample the validation dataset (faster)
    stopping_rounds=2,
    stopping_metric="misclassification", ## could be "MSE","logloss","r2"
    stopping_tolerance=0.01
)



D = h2o.importFile(path = normalizePath("repeat_data.csv"))





nsfa <- "/Users/przemol/Documents/MATLAB/nsfa_vanilla_git/All_repeats_mean_signal_NSFA_unmasked_200iter.mat"
nr <- nr$finalsample
names(nr) <- attributes(nr)$dimnames[[1]]


ans <- sapply( levels(M$class), function(x) {
    colSums(nr$G[M$class==x, ])
})

sum((as.matrix(M[,-10]) - (nr$G%*%nr$X))  ^2) / length(as.matrix(M[,-10]))

