# Create covariates text file
subs <- read.csv('example.csv', sep = ';')
subs$sex[subs$sex == 1] <- 2  #change female code for PLINK
subs$sex[subs$sex == 0] <- 1  #change maele code for PLINK
subs <- subs[, c(2,8)] #select sex and IID (sample names) columns 
tsubs <- t(subs) #transpose dataframe

write.table(tsubs, file = 'covariates.txt', col.names = F,
            row.names = F, quote = F, sep = '\t')