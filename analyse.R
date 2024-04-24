library(ggplot2)
vals <- read.csv('res.csv', header=F)
x <- vals$V1

ggplot(vals, aes(x=V1, y=after_stat(density))) +
	geom_histogram(bins=250) +
	ylim(0,1.4)
