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

M$class <- droplevels(M$class)

smp_size <- floor(0.75 * nrow(M))
train_ind <- sample(seq_len(nrow(M)), size = smp_size)

train <- M[train_ind, ]
test <- M[-train_ind, ]

model <- glm(class ~., family=binomial(link='logit'),data=train)
summary(model)

fitted.results <- predict(model,newdata=test)
fitted.results <- ifelse(fitted.results > 0.5,1,0)

misClasificError <- mean(fitted.results != test$Survived)
print(paste('Accuracy',1-misClasificError))

library(nnet)
mod <- multinom(class ~., train)



library(nnet)
mod <- multinom(class ~., train)



rowMeans(M[,1:9]) < 1




unlist(repeatModel)$id  %>%  head