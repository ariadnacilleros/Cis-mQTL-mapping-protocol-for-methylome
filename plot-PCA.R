#This filw will plot the Prinipal Components by different 
# pheno traits: sex, ethnicity, cohort.
# You will need to adapt the script considering your inputs. 

library(dplyr)
library(ggplot2)
library(data.table)

df1 <- fread('example.csv', sep = ';') #we need a file with first two columns FID and IID, followed by pheno data (sex,age...)

recode1 <- recode(df1$sex,'0' = 'Male', '1' = 'Female') #we make sure that Male = 0 and Female = 1 according to PLINK
df1[,sex:=recode1] 
recode2 <- recode(df1$ethnic_origin, '2' = 'Asiatic', '3' = 'Black', '4' = 'Arabe', #change codes for ethnic origin
                  '5' = 'Gipsy', '6' = 'American indios', '7' = 'Others',
                  '11' = 'Both parents white and Spanish',
                  '12' = 'Both parents white, not Spanish, but European',
                  '13' = 'Both parents white, but at least one not European')
df1[,ethnic_origin:=recode2]

pcs <- fread('PIHAT-clean.PCs.eigenvec') #upload PCs data
names(pcs) <- c('FID', 'IID', paste0('PC', 1:20)) #change colum names

p <- merge(df1, pcs, all.y = T) #merge both dataframes: pheno data + PCs

#Plot PCs and pheno variables 

pdf('PIHAT-clean-PC1vsPC2-cohort.pdf')
qplot(PC1, PC2, colour = cohort8_name,
      data = p, alpha = 0.2) + theme_light()
dev.off()

pdf('PIHAT-clean-PC3vsPC4-cohort.pdf')
qplot(PC3, PC4, colour = cohort8_name,
      data = p, alpha = 0.2) + theme_light()
dev.off()

pdf('PIHAT-clean-PC1vsPC2-ethn.pdf')
qplot(PC1, PC2, colour = ethnic_origin,
      data = p) + theme_light()
dev.off()

pdf('PIHAT-clean-PC3vsPC4-ethn.pdf')
qplot(PC3, PC4, colour = ethnic_origin,
      data = p) + theme_light()
dev.off()

pdf('PIHAT-clean-PC1vsPC2-sex.pdf')
qplot(PC1, PC2, colour = sex,
      data = p) + theme_light()
dev.off()

pdf('PIHAT-clean-PC3vsPC4-sex.pdf')
qplot(PC3, PC4, colour = sex,
      data = p) + theme_light()
dev.off()


