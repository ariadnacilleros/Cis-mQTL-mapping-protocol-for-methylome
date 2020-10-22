# Add sex to GIP & AST cohort
fam <- read.delim('example.fam', header = F) #CHANGE FILE'S NAME
subs <- read.csv('pheno_example_data.csv', sep = ';') #CHANGE FILE'S NAME
subs$sex[subs$sex == 1] <- 2
subs$sex[subs$sex == 0] <- 1
subs <- subs[, c(2,8)]

tmp <- merge(fam, subs, by.x = 'V2', by.y = 'plink_IID', all.x = TRUE, sort = F)
tmp$V5 <- ifelse(is.na(tmp$sex), tmp$V5, tmp$sex)
tmp$sex <- NULL
tmp <- tmp[,c(2,1,3:6)]

library(stringr)
tmp <- tmp[str_order(tmp$V2, numeric = T),]
tmp <- tmp[str_order(tmp$V1, numeric = T),]

write.table(tmp, file = 'example.fam', col.names = F, #CHANGE OUTPUT FILE'S NAME FOR THE DESIRED ONE
            row.names = F, quote = F, sep = '\t')