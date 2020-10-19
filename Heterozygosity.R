het <- read.table("Oxford_H.het", header=TRUE)
het$meanHet <- (het$N.NM. - het$O.HOM.)/het$N.NM.
SD <- sd(het$meanHet)
mean_Het <- mean(het$meanHet)
upper <- (SD*4)+mean_Het
lower <- mean_Het - (SD*4) 
data <- data.frame(row.names=which(upper < het$meanHet), which(upper < het$meanHet))
write.table(x = data, file = "4SD.txt", sep = "\t", row.names=TRUE, col.names=TRUE)
