#Objective: Convert Variant ID column from a PLINK's map file to "chr:position" format, using chromosome and position columns
#Version: 08-02-2021
#Contact: acilleros@ikasle.ehu.eus

# OVERRIDES .BIM FILES

arg <- commandArgs(trailingOnly = T)

cat('\nExecuting change-rsid.R script.')

# List of dataframes
myfiles <- lapply(arg, read.delim, header = F)

# This section converts Variant ID column from a PLINK's map file
# to "chr':position" format, using chromosome and position columns

library(data.table)
chrpos <- function(x){
  x <- as.data.table(x)
  x[, V2 := paste(V1, V4, sep = ':')]
}

# Apply conversion function to each dataframe in myfiles
res <- lapply(myfiles, chrpos) # result

cat('\n')

for (i in 1:length(res)) {
  fwrite(res[[i]], arg[i], quote = F, sep = '\t',
         row.names = F, col.names = F)
  cat('File', arg[i], 'has been overwritten.\n')
}

cat('\n')
