#Lets construct two well corelated vectos with random noise

b <- 1:100+10*rnorm(100)
a <- 1:100+10*rnorm(100)

#Plotting as scatterplot
plot(a,b)

#Correlation coefficient
cor(a,b)

plot(a,b)
cor(a,b)