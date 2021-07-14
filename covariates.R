#Objective: analyse results and extract the necessary tables to send 
#Latest version: 14-02-21
#Contact: acilleros001@ikasle.ehu.eus

#set your working directory inside covariates folder
setwd("./covariates/")

#Step 1: Load file with the sex of each sample
sex_repor <- read.delim("met-samples-sex.txt", header=FALSE)

#Step 2: Load Principal Component Values
PC <- read.table("whole_genome_maf05_filt_samples_prunned.PCs.eigenvec", sep=" ")

#Step 3: Merge both files by sample name or IID
cov <- merge(x = sex_repor[,c(1,3)], y = PC[,c(2:7)],by.x="IID", by.y="V2")

#Step 4: Change column names by PC followed by a number
colnames(cov)[3:7] <- c("PC1","PC2", "PC3", "PC4", "PC5")

#Step 5: Load Planet values (omegas)
eda <- readRDS('./EPIC/cohort_analysisdate_Output/exp_idat_diff_samples.RDS')
planet <- eda$Omega
rm(eda)

#You will need to filter planet data.frame by the final samples set on ./final_list_basename.txt or ./final_list_samples.txt
# and change the sample names being the same as the genotype, the PCs' and the sex data.frame...

#Step 6: Merge Planet values
cov <- merge(x = cov, y = planet, by.x="IID", by.y="row.names")
colnames(cov)[1] <- "id"

#Step 5: Write the final text file inside covariates folder
write.table(x = t(cov), file = "covariates.txt", quote = F, row.names = T, col.names = F, sep = "\t") 
