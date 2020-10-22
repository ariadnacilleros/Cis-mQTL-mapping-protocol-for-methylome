library(data.table)
arg <- commandArgs(trailingOnly = T)
print(arg)
indfile <- fread(arg[1], sep = '\t')


ind <- c(rbind(indfile[, IID1], indfile[,IID2]))
geno <- fread('example.imiss') #CHANGE FILE'S NAME WITH .IMISS FILE FROM STEP 1.1.1.
merged <- merge(as.data.frame(ind), geno, by.x = 'ind', by.y = 'IID', sort = F)


fwrite(merged, 'duplicate_samples_miss.txt', sep = "\t",
       quote = F, row.names = F)

# Anotating duplicates with higher n_miss for removal
loophelp <- 1:nrow(merged)
loophelp <- loophelp[loophelp %% 2 == 1]

iids <- NULL
fids <- NULL
for (i in loophelp){
  temprow <- merged[c(i, i+1), ]
  higher_n <- max(temprow$N_MISS)
  selectedrow <- temprow[temprow$N_MISS == higher_n, ]
  iids <- c(iids, selectedrow$ind)
  fids <- c(fids, selectedrow$FID)
}

outdf <- data.frame('FID' = fids, 'IID' = iids)
fwrite(outdf, arg[2], sep = '\t', col.names = F)
