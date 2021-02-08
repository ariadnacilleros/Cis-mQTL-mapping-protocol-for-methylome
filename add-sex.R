# Objective of the script: Add the sex of the samples. 
# Last version: 08-02-2021
# Contact: acilleros001@ikasle.ehu.eus

#Step 1: Read the fam file from our genotype data
fam <- read.delim('whole_genome.fam', header = F) 
#Step 2: Read or upload the data were is contained information regarding the sex of the samples
subs <- read.csv('pheno_example_data.csv', sep = ';') #CHANGE FILE'S NAME WITH SEX INFO

#Step 3: Make sure that female is coded by 2 and male coded by 1. In our case we need to change {0,1} for {1,2}. 
subs$sex[subs$sex == 1] <- 2
subs$sex[subs$sex == 0] <- 1
subs <- subs[, c(2,8)] #subset the whole dataframe with only the column 2 (name of the samples --> plink_IID) and column 8 (sex of the samples)

#Step 4: Merge both data frames (fam file and sex) by the name of the samples (e.g., plink_IID column)
tmp <- merge(fam, subs, by.x = 'V2', by.y = 'plink_IID', all.x = TRUE, sort = F)

#Step 5: Add sex into the column 5 from fam file to the final dataframe
tmp$V5 <- ifelse(is.na(tmp$sex), tmp$V5, tmp$sex) 
tmp$sex <- NULL #remove sex column (now we have column 5 with this information)
tmp <- tmp[,c(2,1,3:6)] #select the columns according to the fam file format of PLINK (have a look at the documentation)

#Step 6: Sort dataframe by IID and then by FID
library(stringr)
tmp <- tmp[str_order(tmp$V2, numeric = T),]
tmp <- tmp[str_order(tmp$V1, numeric = T),]

#Step 7: obtain the new version of the fam file with the sex uploaded
write.table(tmp, file = 'whole_genome.fam', col.names = F, 
            row.names = F, quote = F, sep = '\t')
