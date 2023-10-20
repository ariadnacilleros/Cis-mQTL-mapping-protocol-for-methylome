#Objective: GWAS of the residualized methylation PCs with the final genotype
#Version: 19/10/2023
#Contact: ariadna.cilleros@ehu.eus

library(data.table)
library(dplyr)

#Step 1: Load fam file
setwd("./rnt_model")
fam <- fread("../whole_genome_definitive/whole_genome_maf05_filt_samples.fam", col.names = c("FID","IID","fID","mID","sex","pheno"))

#Step 2: Load mPCs
mpcs <- fread("./reisudalized_mPCs.txt")
colnames(mpcs)[1] <- "FID"

#Step 3: Merge fam file with mPCs
fam_modified <- merge(fam, mpcs, by= "FID", sort = F)

#Step 4: Creating one fam file per mPC:
for (i in colnames(mpcs)[-1]) {
  print(paste0("working on ", i))
  data <- fam_modified
  data$pheno <- subset(data,select = i) # assign as phenotype column an mPC
  data <- data[, c(1:6)] # select only PLINK fam file columns
  fwrite(x = data, file = paste0("./whole_genome_maf05_filt_samples_", i, ".fam"), # write fam file with phenotype being an mPC
              col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
}

#Step 5: Perform GWAS
for (i in colnames(mpcs)[-1]) {
  file <- paste("./whole_genome_maf05_filt_samples_", i, ".fam", sep="") # create string with name of the modified fam file
  print(file)
  output_file <- paste("./gwas_", i, sep="") # create string with name of the output file
  print(output_file)
  gwas_command <- paste(paste(paste("plink1.9 --bfile ../whole_genome_definitive/whole_genome_maf05_filt_samples --fam ", file, sep="")," --assoc --allow-no-sex --out ",sep=""), output_file, sep="") # write string with GWAS command
  print(gwas_command)
  system(gwas_command) # perform GWAS using the string with the command
}

#Step 6: Filter all the GWAS 
system("awk '{ if ( $9 < 1e-7) print $0 }' ./*.qassoc | wc -l")
# If you get any SNP with p-value < 1e-7, exclude from the analysis the mPC containing that association