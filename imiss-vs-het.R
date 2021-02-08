#Objective: Plot the missing call rate vs heterozygosity and write a text file with the names of the samples with an heterozygosity > +/- 4 x Standard Deviation 
#Version: 08-02-2021
#Contact: acilleros001@ikasle.ehu.eus

library('geneplotter')

arg <- commandArgs(trailingOnly = T)
imiss_file <- fread(arg[1], sep = '\t')
het_file <- fread(arg[2], sep = '\t')

imiss <- read.table(imiss_file, header = T) #CHANGE FILE'S NAME
imiss$logF_MISS <- log10(imiss[, 6])
het <-  read.table(het_file, header = T) #CHANGE FILE'S NAME
het$meanHet <- (het$N.NM. - het$O.HOM.) / het$N.NM.
colors  <- densCols(imiss$logF_MISS, het$meanHet)

pdf('imiss-vs-het.pdf')
plot(imiss$logF_MISS, het$meanHet, col = colors, xlim = c(-3, 0),
     ylim = c(0, 0.5), pch = 20, xlab = 'Proportion of missing genotypes',
     ylab = 'Heterozygosity rate', axes = F, cex = 0.5)
axis(2, at = c(0, 0.05, 0.10, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5), tick = T)
axis(1, at = c(-3, -2, -1, 0), labels = c(0.001, 0.01, 0.1, 1))
abline(h = mean(het$meanHet) - (4 * sd(het$meanHet)), col = 'red', lty = 2)
abline(h = mean(het$meanHet) + (4 * sd(het$meanHet)), col = 'red', lty = 2)
abline(v = -1.522879, col = 'red', lty = 2)
dev.off()

t1 <- mean(het$meanHet) + (4 * sd(het$meanHet))
t2 <- mean(het$meanHet) - (4 * sd(het$meanHet))

cat(sum(het$meanHet > t1), 'samples above 4*SD from the average heterozigosity.\n')
cat(sum(het$meanHet < t2), 'samples below 4*SD from the average heterozigosity.\n')

het <- het[het$meanHet > t1 | het$meanHet < t2, ]
output <- het[,c(1,2)]
write.table(output, 'filter-het.txt', col.names = F,
            row.names = F, quote = F, sep = '\t')
