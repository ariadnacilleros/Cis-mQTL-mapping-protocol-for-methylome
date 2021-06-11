#Objective: Obtain a BED file from the ESet object of Alexandra Binder's R package. 
#Version: 08-02-2021
#Contact: acilleros001@ikasle.ehu.eus

library(minfi)
library(stringr)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

#Step 1: Move your working directory to EPIC folder that we created previously. 
setwd("EPIC/")
processedOut <- readRDS("./EPIC/cohort_analysisdate_Output/preproc_exp_idat.RDS") #load the ESet object from the QC of PACEanalysis

#Step 2: Convert pheno data into a data.frame
pheno <- data.frame(pData(processedOut$mset))

#At this point you need to have the same sample names between methylation and genotype, and removed the duplicates (if this is the case). 
#To remove the duplicates we could do the following subset by discarding the rows 168, 300 and 340, were in our case we have the duplicates: 
no_duplicates <- pheno[-c(168,300,340),]

#Step 3: Filter samples. 
samples_geno <- read.table("../imputed-rsq09/sample_list_geno.txt", quote="\"", comment.char="")
final_list_samples <- intersect(samples_geno$V1, no_duplicates$Sample_Name)
write.table(x = final_list_samples, file = "./final_list_samples.txt", sep="\t", row.names = F, col.names = F)
no_duplicates <- no_duplicates[no_duplicates$Sample_Name %in% final_list_samples, ]

#Step 4: Extract beta values and CpG annotation. 
beta <- processedOut$processedBetas
annot <- getAnnotation(processedOut$mset)

#Step 5: Merge annotation and beta data.frames by CpG IDs.
annot_beta <- merge(annot[,c(1,2,2)], beta, by=0, all=TRUE) #**

#Step 6: Order the data.frame with first the chr, second the start (pos), third the end (pos.1), 
#fourth the CpG id (Row names) and the beta values from the filtered samples of the pheno data data.frame.
df <- annot_beta[,c("chr","pos","pos.1","Row.names",no_duplicates$Basename)]

#Step 7: Tune the data.frame according to a BED format.
#Step 7.1: Change "chr1" by "1" 
chr_vector <- str_remove(df$chr, "chr")
df$chr <- chr_vector

#Step 7.2: Change column names. On most cases, the sample names of the beta data.frame 
#will correspond to the basename of the array, therefore at this point we changed them 
#by the correct sample names, which should be the same as the genotype.
colnames(df)[1:4] <- c("#Chr", "start", "end", "ID")
for (i in 5:ncol(df)){ 
  colnames(df)[i] <- no_duplicates[no_duplicates$Basename == colnames(df)[i], "Sample_Name"]
}

#Step 8: Filter CpGs. 
#Step 8.1: Filter CpGs located on sexual chromosomes.
sexual_probes <- df[df$`#Chr`==c("X","Y"),4]
df_filt_sex <- df[df$`#Chr`!=("X") & df$`#Chr`!=("Y"),]

#Step 8.2: Filter CpGs with SNPs (MAF < 5% in EUR)
SNP <- read.table(file = "1-s2.0-S221359601630071X-mmc1.txt", header = T, sep = "\t")
SNP_probes <- SNP[SNP$EUR_AF>0.05,1]
df_filt_sex_snp <- subset(df_filt_sex, !(ID %in% SNP_probes))

#Step 8.3: Filter cross-hybridizing CpGs
CpG_cross <- read.table(file = "1-s2.0-S221359601630071X-mmc2.txt", sep = "\t")
df_filt_sex_snp_cross <- subset(df_filt_sex_snp, !(ID %in% CpG_cross$V1))

#Look for cross-hybridizing non-CpG positions
nonCpG_cross <- read.table(file = "1-s2.0-S221359601630071X-mmc3.txt", sep = "\t")
df_filt_sex_snp_cross <- subset(df_filt_sex_snp_cross, !(ID %in% nonCpG_cross$V1))

#Step 9: Obtain the final dataframe inside a text file following a BED structure format
write.table(x = df_filt_sex_snp_cross, file = "./methylome_BED.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

#Step 10: Write Planet values (Omegas) on text file
omegas <- processedOut$Omega
#filter by samples
omegas_filt <- omegas[rownames(no_duplicates),]
#change sample names 
for (i in 1:370){ #change 370 for the total ammount of samples that you have at this point
  rownames(omegas_filt)[i] <- no_duplicates[no_duplicates$Basename == rownames(omegas_filt)[i], "Sample_Name"]
}
#write omegas 
write.table(x = omegas_filt, file = "EPIC/planet_values.txt", quote = FALSE, sep = "\t", row.names = T, col.names = TRUE)


#Step 11: Calculate the variance of the probes
#Step 11.1: Method difference between 10 and 90 percentiles
Variation<-function(x) {quantile(x, c(0.9), na.rm=T)[[1]]-quantile(x, c(0.1), na.rm=T)[[1]]} #define function
#apply the function per row (cpg)
diff_percent <-as.numeric(lapply(1:nrow(df_filt_sex_snp_cross), function(x) Variation(df_filt_sex_snp_cross[x,5:ncol(df_filt_sex_snp_cross)])))
names(diff_percent) <- df_filt_sex_snp_cross$ID
df_diff_percent <- data.frame(diff_percent, row.names = names(diff_percent)) #create dataframe for the values

#Step 11.2: Method variances
variances <- as.vector(apply(df_filt_sex_snp_cross[,5:ncol(df_filt_sex_snp_cross)], 1, FUN = var)) #calculate varainces per row (cpg)
names(variances) <- df_filt_sex_snp_cross$ID
df_variances <- data.frame(variances, row.names = names(variances)) #create dataframe for the values

#Step 10.3: Write variability stats into a text file 
df <- merge(df_diff_percent, df_variances, by ="row.names")
write.table(file = "EPIC/all_cpg_variances.txt", x = df, sep = "\t", quote = F, row.names = F, col.names = T, dec = ".", na = "NA")
