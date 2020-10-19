sex <- read.csv("Final_imp.sexcheck", sep="")

df$FID <- paste(df$FID,df$FID, sep = "_")

df_t <- t(df)
rownames(df_t)[1] <-("id")

write.table(x = df_t, file = "COV_imp.txt", sep = "\t", quote = F, col.names = F)
