#This script allows you to calculate the power of your data
#and remove SNPs with low MAF

#library(devtools);install_github("sterding/powerEQTL")  
library('powerEQTL')

# Set sample size
N <- c(100,200,300,400,500)
nn <- length(N)

# Set MAF
MAF <- seq(0.5,20,0.1)/100
nq <- length(MAF)

# number of SNPs tested (20 SNPs per gene x 600.000 CpGs in total), change the number of CpGs for yours
nSNP=(20 * 600.000)

# significant level (FP)
a=0.01

# obtain power
power_unbalanced <- array(numeric(nn*nq), dim=c(nn,nq))
for (i in 1:nn){
  for (j in 1:nq){
    # unbalanced
    result <- powerEQTL.ANOVA(MAF=MAF[j], 
                              typeI=a, 
                              nTests=nSNP, 
                              myntotal=N[i], 
                              mystddev=0.03116973, #change 0.03116973 values by the mean of the standard deviance per CpG (beta) 
                              deltaVec = c(0.03116973, 0.03116973),
                              verbose = F)
    power_unbalanced[i,j] <-result;
  }
}

# set up graph
xrange <- range(MAF*100)
yrange <- c(0:1)
colors <- rainbow(length(N))
pdf(file = "power.pdf")
plot(xrange, yrange, log='x', type="n", 
     xlab="MAF (%)",
     ylab="Power",
     main="Power Estimation for eQTL Studies, Sig=8.3e-7, nCpG=12.000.000 (one-way unbalanced ANOVA)")

abline(v=0, h=seq(0,1,.1), lty=2, col="grey89")
abline(h=0, v=c(1:10), lty=2,col="grey89")

# add power curves
for (i in 1:nn){
  lines(MAF*100, power_unbalanced[i,], type="l", lwd=4, col=colors[i])
}

legend("topleft", title="Sample size (n)", as.character(N),fill=colors, bty='n')

dev.off()
