setwd("") #set your wd where you have the unzipped files from the server (.vcf.gz and .info.gz)

for (i in c("X", 1:22)){ 

require(data.table) #Load data.table package 
  
# Reading data from generated files info 
  
# NOTE: Filenames are being split by i into the part before the chromosome number and after the chromosome number, the chromosome number will be added automatically by the script

dataInfo <- read.delim(file = (paste(getwd(),"/Chr/chr",i,".info.gz",sep="")), header = T, sep = "\t", dec = ".")

## Make dataframe of wanted columns and change the column names

dataInfoClean <- data.frame(dataInfo$SNP,dataInfo$REF.0.,dataInfo$ALT.1.,dataInfo$Rsq,dataInfo$Genotyped)
names(dataInfoClean)[names(dataInfoClean)=="dataInfo.SNP"] <- "SNP"
names(dataInfoClean)[names(dataInfoClean)=="dataInfo.Rsq"] <- "Rsq"

dataInfoClean$SNP <- as.character(dataInfoClean$SNP)
colnames(dataInfoClean) <- c("SNP","REF_ALLELE","ALT_ALLELE","Rsq","SOURCE")

# Filtering Files (Rsq > 0.9)

filtered_data <- dataInfoClean[dataInfoClean$Rsq > 0.9,]

#Discard multiallelic SNPs

single_allelic <-names(which(table(filtered_data$CHR_POS)==1))
SNP_keep <- filtered_data[filtered_data$CHR_POS %in% single_allelic,]

#Write SNP text file to keep 

write.table(SNP_keep$SNP,paste(getwd(), "/filtered_imputed_SNP_Chr",i,".txt",sep=""),quote=F,row.names=F,col.names=F)

}
