setwd("")#set your wd
#Steps to perform before: 
#1)Load .IDAT files
#2)Obtain betas 
#3)Quality control
#For these steps you can use minfi Package for example.
#Output must be a GenomicRatioSet object!

library(minfi)
library(stringr)

class(processedOut$mset)
#GenomicRatioSet: class holding M or/and Beta values together with associated genomic coordinates

#Change sample's name: 
#In our case we stablished a pattern for sample's name: FID_IID 
#FID = cohort (GIP, VAL, SAB...)
#IID = C_number (C_0428,...)
#Ex: GIP_C_0428

##A.Standarize samples between genotype and methylome
#A.1.Obtain the samples data
samples_data <- data.frame(pData(processedOut$mset))
samples_data[samples_data$IID == "", ]

#For loop to build sample's name
for (i in 1:387){ #our total number of samples is 387!

  num <- scan(text=samples_data$Sample_Name[i], sep="_", what="", quiet=TRUE)[4] 
  
  if (num == "D90"){
    num <- scan(text=samples_data$Sample_Name[i], sep="_", what="", quiet=TRUE)[3]
  }
  
  samples_data$IID[i] <- (paste("C",sprintf("%04d", as.numeric(num)), sep="_"))
  samples_data$FID[i] <- samples_data$Cohort[i]
}

#A.2.We also removed duplicates: 168 300 340
no_duplicates <- samples_data[-c(168,300,340),]

#A.3.Create a text file with the IID of our samples to filter genotype data
write.table(x = no_duplicates$IID, file = "IID_methylome.txt", quote = F, sep = "\t", row.names = F, col.names = F)


##B.Start building BED file
#B.1.Extract beta values
beta <- getBeta(processedOut$mset) 

#B.2.Obtain CpG's annotation
annot <- getAnnotation(processedOut$mset) 

#B.3.Merge both dataframes by CpG id
annot_beta <- merge(annot[,c(1,2,2)], beta, by=0, all=TRUE)  

#B.4.Order columns: chr, start, end, ID, sample1, sample2...
df <- annot_beta[,c(2,3,4,1,(5:391))] 

#B.5.Change chr name (e.g. chr1 --> 1)
chr_vector <- str_remove(df$chr, "chr") 
df$chr <- chr_vector

#B.6.Change samples names from betas dataframe according to previous samples_data
names_vector <- c()
for (i in 5:391){
  names_vector <- append(names_vector, paste(no_duplicates[colnames(df)[i],c(32)], no_duplicates[colnames(df)[i],c(30)], sep = "_"))
}

#B.7.Define colnames
colnames(df) <- c("#Chr", "start", "end", "ID", names_vector)


#Jump to step 2.3.Filter samples of genotype data(Wiki), and filter genotype data and obtain the name of the samples


##C.Filter samples and SNPs
#C.1.Filter by genotype samples
samples_imp_vcf <- read.table("samples_imp_vcf.txt", quote="\"", comment.char="")
df_filt_imp <- df[,c("#Chr", "start", "end", "ID",samples_imp_vcf$V1)]

#C.2.Filter low variable CpGs, in this case we choose according to the quartiles of the beta's variance
quantile(apply(df_filt_imp[,c(5:377)], 1, FUN = var))
df_filt_imp_var <- df_filt_imp[apply(df_filt_imp[,c(5:377)], 1, FUN = var) > 2.44e-04,] #use the variance of the first quantile (25%)

#C.3.Look for CpGs with SNPs (MAF < 5% in EUR)
SNP <- read.table(file = "McCarthy/1-s2.0-S221359601630071X-mmc1.txt", header = T, sep = "\t")
table(SNP$EUR_AF>0.05)
SNP_probes <- SNP[SNP$EUR_AF>0.05,1]
df_filt_imp_var_sex_SNP <- subset(df_filt_imp_var_sex, !(ID %in% SNP_probes))

#C.4.Look for cross-hybridizing CpGs
CpG_cross <- read.table(file = "McCarthy/1-s2.0-S221359601630071X-mmc2.txt", sep = "\t")
df_filt_imp_var_sex_SNP_CpG <- subset(df_filt_imp_var_sex_SNP, !(ID %in% CpG_cross$V1))

#C.5.Look for cross-hybridizing non-CpG positions
nonCpG_cross <- read.table(file = "McCarthy/1-s2.0-S221359601630071X-mmc3.txt", sep = "\t")
df_filt_imp_var_sex_SNP_CpG_cross <- subset(df_filt_imp_var_sex_SNP_CpG, !(ID %in% nonCpG_cross$V1))

#Take the total number of CpGs remaining to calculate the statistical power 

#WRITE BED FILE IN TEXT FILE
write.table(x = df_filt_imp_var_sex_SNP_CpG_cross, file = "/data/Genomics_projects/Placenta_FastQTL/EPIC/whole_genome_imp_var_bed.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)




