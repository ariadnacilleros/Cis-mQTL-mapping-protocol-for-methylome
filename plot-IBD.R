#This script will plot the IBD and write a text file
#with the name of the samples with a PI-HAT > 0.18

arg <- commandArgs(trailingOnly = T)
data <- read.table(paste0(arg, '.genome'), header = T)

namepdf <- paste0(arg,'-hist.pdf')
pdf(namepdf)
hist(
  data$PI_HAT, ylim = c(0, 100), col = 'red', breaks = 100,
  xlab = 'Estimated mean pairwise IBD', main = '')
dev.off()

out <- data[data$PI_HAT > 0.18,]
nametable <- paste0(arg, '-fail-IBD-check.txt')
write.table(out, nametable, sep = '\t', quote = F, row.names = F)