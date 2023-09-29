#Objective: Obtain a BED file from Alexandra Binder's R package output and transform beta-values to inverse variance rank-transform 
#Version: 29-09-2023
#Contact: acilleros001@ikasle.ehu.eus

library(minfi)
library(stringr)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

#Step 1: Move your working directory to EPIC folder that we created previously. 
setwd("EPIC/")
betas <- readRDS("./EPIC/cohort_analysisdate_Output/outlier.rds") #load output from the QC of PACEanalysis
dim(betas)

#Make sure we don't have NAs
table(is.na(betas))

#Step 2: Read file with all the metadata from the samples
pheno <- read.delim("./Phenodata_basename_EPIC_BMI_all_names.txt", header = T, sep="\t") #adapt it

#Filter by the samples that pass the methylation QC
pheno <- pheno[pheno$Basename %in% colnames(betas), ]

#Step 3: Intersection between methylation and genotype samples
samples_geno <- read.table("../imputed-rsq09/sample_list_geno.txt", quote="\"", comment.char="", sep = "_")

#Obtain a final list of samples combining the final list of each type of data
final_list_samples <- intersect(samples_geno$V1, pheno$FID_IID)

#Filter pheno and betas data.frames
pheno <- pheno[pheno$FID_IID %in% final_list_samples, ]
betas <- betas[,pheno$Basename]

#Write a text file that allows you to filter the genotype data (PLINK) as in the Step 3 from the whole protocol 
write.table(x = cbind(final_list_samples,final_list_samples), file = "./final_list_samples.txt", sep="\t", row.names = F, col.names = F, quote = F)

#Write a text file with the FID_IID and basenames of the final samples dataset
write.table(x = pheno[,c("FID_IID","Basename")], file = "./final_list_basename.txt", sep="\t", row.names = F, col.names = F, quote = F)

#Step 4: Transform beta values
rnt=qnorm(t(apply(betas, 1, rank, ties.method = "average"))/ (ncol(betas)+1))

#Step 5: Annotate CpGs
#Get annotation
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
annot <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
  
#Merge annotation and beta data.frames by CpG IDs
annot_rnt <- merge(x=annot[,c(1,2,2,4)], y=rnt, by="row.names")
rm(annot)
rownames(annot_rnt) <- annot_rnt$Row.names #Make CpG names' as row.names of the data.frame 
annot_rnt <- annot_rnt[,-1]  #Remove column with CpG names'

#Step 6: Tune the data.frame according to a BED format.
#Step 6.1: Change "chr1" by "1" 
annot_rnt$chr <- str_remove(annot_rnt$chr, "chr")

#Step 6.2: Change column names
colnames(annot_rnt)[1:4] <- c("#Chr", "start", "end", "ID")

for (i in 5:372){ #change basename to FID_IID 
  colnames(annot_rnt)[i] <- pheno[pheno$Basename == colnames(annot_rnt)[i], "FID_IID"]
}

#Step 6.3: Remove CpGs on sexual Chr
annot_rnt_sex <- annot_rnt[annot_rnt$`#Chr`!=("Y"),]
annot_rnt_sex <- annot_rnt_sex[annot_rnt_sex$`#Chr`!=("X"),]

#Step 6.4: Define the cis-window (end = start + 1)
annot_rnt_sex$end <- annot_rnt_sex$start + 1

#Step 7: Obtain the final dataframe inside a text file following a BED structure format
write.table(x = annot_rnt_sex, file = "./methylome_BED_rnt.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
