setwd("./covariates/")

#Step 1: Load file with the sex of each sample
sex_repor <- read.delim("met-samples-sex.txt", header=FALSE)

#Step 2: Load Principal Component Values
PC <- read.table("whole_genome_maf05_filt_samples_prunned.PCs.eigenvec", sep=" ")

#Step 3: Merge both files by sample name or IID
cov <- merge(x = sex_repor[,c(1,3)], y = PC[,c(2:7)],by.x="IID", by.y="V2")

#Step 4: Change column names by PC followed by a number
colnames(cov)[3:7] <- c("PC1","PC2", "PC3", "PC4", "PC5")

#Step 5: Write the final text file on covariates folder
write.table(x = t(cov), file = "covariates_sex_PC5.txt", quote = F, row.names = T, col.names = F, sep = "\t") 
