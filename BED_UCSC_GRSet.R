setwd("")#set your wd
load("") #load your .RData object with the GRSet (in case it is necessary)

require(minfi)
require(stringr)

class(processedOut$mset)
#GenomicRatioSet: class holding M or/and Beta values together with associated genomic coordinates

beta <- getBeta(processedOut$mset) #extract beta values from RGSet
annot <- getAnnotation(processedOut$mset) #extract annotation from RGSet

CpG <- rownames(beta) #extract CpG IDs

annot_beta <- merge(annot[,c(1,2,2)], beta, by=0, all=TRUE) #merge chr, start, end and beta values per sample, by CpG id

df <- annot_beta[,c(2,3,4,1,(5:391))] #order columns into BED UCSC format file
chr_vector <- str_remove(df$chr, "chr") #change chr name (e.g. chr1 --> 1)
df$chr <- chr_vector
colnames(df) <- c("#Chr", "start", "end", "ID", colnames(df)[5:391]) #define colnames

write.table(x = df, file = "INMA_preprocessed_bed.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

