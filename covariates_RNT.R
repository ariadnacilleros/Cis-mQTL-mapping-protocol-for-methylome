#Objective: Add mPCs to covariates file with
#Version: 19/10/2023
#Contact: ariadna.cilleros@ehu.eus

setwd("./rnt_model")
library(data.table)

#Step 1: Load methyaltion PCs
mpcs <- fread("./reisudalized_mPCs.txt")
colnames(mpcs)[1] <- "FID_IID"

#Step 2: Load covariates file (obtained by covariates.R, Step 4 from GitHub) 
cov <- fread("../covariates/covariates.txt", header = T)
cov_t <- t(as.matrix(cov[,-c(1)])) # transpose datatable excluding id column
colnames(cov_t) <- cov$id # assign colnames as the names of the variables
cov_t <- as.data.frame(cov_t) # convert to dataframe
cov_t$FID_IID <- rownames(cov_t) # create a column with samples name
table(sapply(cov_t, class))  # make sure that the class of the different variables is correct (all of them must be numeric, sex could be numeric or character)

#Step 3: Merge covariates and mPCs dataframe 
merged <- merge(x = cov_t, y = mpcs, by = "FID_IID")
colnames(merged)[1] <- "id"

#Step 4: Write text file with the following covariates: sex, genotype PCs, plante cell type proportions and residualized mPCs
write.table(x = t(merged), file = "./covariates_RNT.txt", quote = F, col.names = F, row.names = T, sep="\t")
