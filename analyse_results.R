#Objective: analyse results and extract the necessary tables to send 
#Latest version: 14-02-21
#Contact: acilleros001@ikasle.ehu.eus

#set your working directory to TensorQTL folder
setwd("TensorQTL/")

#Step 1: Read results text file
results <- read.table("cis_tensorQTL_maf05_PC5_sex_INMA.txt", header = T, sep = "\t")

#Step 2: Adjust by multiple-testing (Bonferroni)
results$Bonferroni <- p.adjust(results$pval_beta, method = "bonferroni")

#Number of significant mQTLs
dim(results[results$Bonferroni<0.05, ])
#68.361 sign

#Step 3: Obtain the variance of the significant CpGs from the all_cpg_variances text file (EPIC folder)
var <- read.table(file = "../EPIC/all_cpg_variances.txt", header = T, sep = "\t")

#Step 3.1: Write a text file with the variance of each CpG
df_cpgs <- var[var$Row.names %in% results$phenotype_id, ]
write.table(x = df_cpgs, file = "variance_table_cpgs_INMA.txt", 
            sep = "\t", dec = ".", row.names = F, col.names = T, quote = F)

#Step 4: Obtain the number of counts from the significant SNPs from the freqx filr (whole_genome_definitive folder)
counts <- read.table(file = "../whole_genome_definitive/whole_genome_maf05_filt_samples.frqx", 
                     header = T, sep="\t")

#Step 4.1: Merge counts and MAF information of the SNPs
#we took the columns with the SNP id (SNP), the minor allele (should be A1), 
#the major allele (should be A2) and the homozigous minor allele counts (C.HOM.A1.)
df_snps <- counts[counts$SNP %in% results$variant_id, c("SNP", "A1","A2","C.HOM.A1.") ]
df_snps <- merge(x = results[,c("variant_id", "maf")], y = df_snps, by.x = "variant_id", by.y = "SNP")
write.table(x = df_snps, file = "maf_counts_table_snps_INMA.txt", 
            sep = "\t", dec = ".", row.names = F, col.names = T, quote = F)

