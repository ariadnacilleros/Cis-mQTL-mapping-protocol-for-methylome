# Objective: Get methylation PCs
# Version: 19/10/2023
#Contact: ariadna.cilleros@ehu.eus

library(data.table)

#Step 1: Load RNT transform values
setwd("./rnt_model")
rnt <- fread("./methylome_BED_RNT.txt")

#Step2: Transform BED file into data.frame with columns as samples, and rows as CpGs
cpgs <- rnt$ID  # get vector with CpG's ID
rnt <- rnt[,-c(1:4)] # remove chr, start, end and ID columns
rnt <- as.data.frame(rnt) # convert to dataframe
rownames(rnt) <- cpgs # set rownames as CpG's ID

#Step 3: Load list of final samples with FID_IID and Basename (obtained by GSet_to_Bed_RNT.R, Step 5.1 from GitHub)
pheno <- fread("./final_list_basename.txt")

#Step 4: Load covariates files (obtained by covariates.R, Step 4 from GitHub)
cov <- fread("../covariates/covariates.txt")

#Step 5: Create a matrix with methylation values and known confounders (being columns the variables; methylation values-CpGs, sex, planet, genotype PCs; and rows the samples)
#Step 5.1: Covariates file
colnames_cov <- cov$id #get vector with covariates name
cov_t <- as.data.frame(t(cov)[-1,]) # transpose datatable excluding id column
cov_t[,1] <- as.numeric(cov_t[,1]) # make sure that sex column is numeric
colnames(cov_t) <- colnames_cov # assign column names as covariates names

#Step 5.2: RNT file
rnt_t <- t(rnt) #transpose datatable 
dim(rnt_t)
# columns as CpGs, rows as samples

#Step 6: Merge covariates and RNT data frames  
merged_df <- merge(x = rnt_t, y = cov_t, by ="row.names")
rownames(merged_df) <- merged_df$Row.names # assign as rownames sample FID_IID

#Step 7: Perform mulitple linear model (methylaion ~ sex + PC1 + PC2 + PC3 + PC4 + PC5 + Trophoblasts + Stromal + Hofbauer + Endothelial + nRBC + Syncytiotrophoblast)
#Step 7.1: Get dataframe with only methylation values
cpgs <- rownames(rnt) # get vector with CpG's ID
subset_cpgs <- merged_df[, cpgs] # subset dataframe by CpGs
rownames(subset_cpgs) <- rownames(merged_df) # set sample's name as rownames
(sapply(subset_cpgs[,c(1:6)], class)) # make sure methylation values are considered as numeric 

#Step 7.2: Apply as.numeric to all covariates
merged_df[colnames_cov] <- sapply(merged_df[colnames_cov],as.numeric)
(sapply(merged_df[,c(747488:747499)], class))# make sure that methylation values are numeric (you can do it with a subset instead of all the CpGs)

#Step 7.3: Perform the multiple linear model
fit <- lm(as.matrix(subset_cpgs) ~ sex + PC1 + PC2 + PC3 + PC4 + PC5 + Trophoblasts + Stromal + Hofbauer + Endothelial + nRBC + Syncytiotrophoblast, data = merged_df)

#Step 7.4: Get residuals from multiple linear model
residuals <- fit$residuals

#Step 8: Perform PCA on the residuals 
pca <- prcomp((na.omit(residuals)))
colnames(pca$x) <- paste("m", colnames(pca$x), sep="")

#Step 8.1: Get variance explained per PC 
prop<-pca$sdev[1:25]^2/sum(pca$sdev^2)*100
cumprop<-cumsum(pca$sdev[1:25]^2)/sum(pca$sdev^2)*100
pcn<-c("mPC1","mPC2","mPC3","mPC4","mPC5","mPC6","mPC7","mPC8","mPC9","mPC10","mPC11",
       "mPC12","mPC13","mPC14","mPC15","mPC16","mPC17","mPC18","mPC19","mPC20",
       "mPC21","mPC22","mPC23","mPC24","mPC25")
pc<-data.frame(pcn,prop,cumprop)

#Step 8.2: Select PCs explaining 80% of the accumulative variance 
pc <- pc[cumprop > 80, ]

#Step 9: Write text file with the mPCs (select a maximum of 20 mPCs)
pca6 <- pca$x[,c(1:20)]
write.table(x = pca6, file = "./residualized_mPCs.txt", col.names = T, row.names = T, sep="\t", quote = F)

