#Objective: Add the sex of the samples to the fam file.
#Last version: 17-02-2021
#Contact: acilleros001@ikasle.ehu.eus

#load libraries
library(data.table)
library(dplyr)

#Step 1: Upload the fam file from whole_genome and the sample sheet (or any file)
#containing the sex of the samples
gip <- fread('inter/whole_genome.fam')
pheno <- fread('pheno-from-excels.csv', na.strings = '') 

#Step 2: Make sure that male is coded by 1 and female by 2
recode1 <- recode(pheno$sex, '0' = 'Male', '1' = 'Female')
pheno[, sex := recode1]
recode1 <- recode(pheno$sex, 'Male' = '1', 'Female' = '2')
pheno[, sex := recode1]

#Step 3: Assign unkown (0) to samples with missing sex
gip[, V5 := pheno[match(gip$V2, pheno$plink_IID), sex]]
gip[is.na(V5), V5 := 0]

#Step 4: Rewrite the fam file with the sex added
fwrite(gip, 'inter/whole_genome.fam', col.names = F, row.names = F, quote = F, sep = ' ')
