#Objective: Obtain a BED file from the Alexandra Binder's R package output. 
#Version: 14-07-2021
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
samples_geno <- read.table("../imputed-rsq09/sample_list_geno.txt", quote="\"", comment.char="")

#Obtain a final list of samples combining the final list of each type of data
final_list_samples <- intersect(samples_geno$V1, pheno$FID_IID)

#Filter pheno and betas data.frames
pheno <- pheno[pheno$FID_IID %in% final_list_samples, ]
betas <- betas[,pheno$Basename]

#Write a text file that allows you to filter the genotype data (PLINK) as in the Step 3 from the whole protocol 
write.table(x = c(final_list_samples,final_list_samples), file = "./final_list_samples.txt", sep="\t", row.names = F, col.names = F)

#Write a text file with the FID_IID and basenames of the final samples dataset
write.table(x = pheno[,c("FID_IID","Basename")], file = "./final_list_basename.txt", sep="\t", row.names = F, col.names = F, quote = F)

#Step 4: Annotate CpGs
#Get annotation
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
annot <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
  
#Merge annotation and beta data.frames by CpG IDs
annot_beta <- merge(x=annot[,c(1,2,2,4)], y=betas, by="row.names")
rm(annot)
rownames(annot_beta) <- annot_beta$Row.names #Make CpG names' as row.names of the data.frame 
annot_beta <- annot_beta[,-1]  #Remove column with CpG names'

#Step 5: Tune the data.frame according to a BED format.
#Step 5.1: Change "chr1" by "1" 
annot_beta$chr <- str_remove(annot_beta$chr, "chr")

#Step 5.2: Change column names
colnames(annot_beta)[1:4] <- c("#Chr", "start", "end", "ID")

for (i in 5:372){ #change basename to FID_IID 
  colnames(annot_beta)[i] <- pheno[pheno$Basename == colnames(annot_beta)[i], "FID_IID"]
}

#Step 5.3: Remove CpGs on sexual Chr
annot_beta_sex <- annot_beta[annot_beta$`#Chr`!=("Y"),]
annot_beta_sex <- annot_beta_sex[annot_beta_sex$`#Chr`!=("X"),]

#Step 6: Obtain the final dataframe inside a text file following a BED structure format
write.table(x = annot_beta_sex, file = "./methylome_BED.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

#############################################################################################################################################################################
#Step 7: Calculate the variance of the probes
#Step 7.1: Method difference between 10 and 90 percentiles
Variation<-function(x) {quantile(x, c(0.9), na.rm=T)[[1]]-quantile(x, c(0.1), na.rm=T)[[1]]} #define function
#apply the function per row (cpg)
diff_percent <-as.numeric(lapply(1:nrow(annot_beta_sex), function(x) Variation(annot_beta_sex[x,5:ncol(annot_beta_sex)])))
names(diff_percent) <- annot_beta_sex$ID
df_diff_percent <- data.frame(diff_percent, row.names = names(diff_percent)) #create dataframe for the values

#Step 7.2: Method variances
variances <- as.vector(apply(annot_beta_sex[,5:ncol(annot_beta_sex)], 1, FUN = var)) #calculate varainces per row (cpg)
names(variances) <- annot_beta_sex$ID
df_variances <- data.frame(variances, row.names = names(variances)) #create dataframe for the values

#Step 7.3: Write variability stats into a text file 
df <- merge(df_diff_percent, df_variances, by ="row.names")
write.table(file = "EPIC/all_cpg_variances.txt", x = df, sep = "\t", quote = F, row.names = F, col.names = T, dec = ".", na = "NA")
