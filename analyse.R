library(ggplot2)
vals <- read.csv('res.csv', header=F)
x <- vals$V1

binwidth <- 0.01*4
args <- commandArgs(trailingOnly=T)
if (length(args) > 1) {
	if (args[2] %in% c("_FIG1", "_FIG2", "_FIG3")) {
		binwidth <- 0.02*4
	}
	else if (args[2] %in% c("_FIG4", "_FIG5")) {
		binwidth <- 0.05*4
	}
}

ggplot(vals, aes(x=V1, y=after_stat(density))) +
	geom_histogram(binwidth=binwidth) +
	xlab("X_T")
# +
	# ylim(0,1.4)
