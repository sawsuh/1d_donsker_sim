vals <- read.csv('res.csv', header=F)
x <- vals$V1
hist(x, breaks=50)
