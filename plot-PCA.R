#Objective: Plot the Prinipal Components by different pheno traits: sex, ethnicity, cohort.
#Version: 08-02-2021
#Contact: acilleros001@ikasle.ehu.eus

library(dplyr)
library(ggplot2)
library(data.table)

#Step 1: Read a a file with first two columns FID and IID, followed by pheno data (sex,age...)
df1 <- fread('example.csv', sep = ';') 

#Step 2: Make sure that Male = 0 and Female = 1 according to PLINK
recode1 <- recode(df1$sex,'0' = 'Male', '1' = 'Female') 
df1[,sex:=recode1] 

#Step 3: It is not necessary, but we changed the ethnic origin numbers by explanatory names
recode2 <- recode(df1$ethnic_origin, '2' = 'Asiatic', '3' = 'Black', '4' = 'Arabe', 
                  '5' = 'Gipsy', '6' = 'American indios', '7' = 'Others',
                  '11' = 'Both parents white and Spanish',
                  '12' = 'Both parents white, not Spanish, but European',
                  '13' = 'Both parents white, but at least one not European')
df1[,ethnic_origin:=recode2]

#Step 4: Upload PCs data
pcs <- fread('clean-PIHAT-prunned.eigenvec') 
names(pcs) <- c('FID', 'IID', paste0('PC', 1:20)) #change colum names

#Step 5: Merge both dataframes (pheno data + PCs)
p <- merge(df1, pcs, all.y = T) 

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


