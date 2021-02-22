#Objective: Preprocess the methylation data.
#Last version: 22-02-2021
#Contact: acilleros001@ikasle.ehu.eus

#load library and set working directory
setwd("./EPIC")
library(PACEanalysis)

#Step 1: Load your phenodata (sample information file/sample sheet )file with the basenames 
pheno <- read.delim("./Phenodata_basename_EPIC_BMI_all_names.txt", header = T, sep="\t") #adapt it

#Step 2: Loading Samples
#Have a look at the reference manual for the package and adapt it to your data
idat <- loadingSamples(PhenoData = pheno, 
                       IDlink = "FID_IID", 
                       SEXvar = "sex", 
                       FemaleInd = "female", 
                       MaleInd = "male", 
                       ETHNICvar = "ethnic_origin2c", 
                       BWTvar = "peso", 
                       GESTvar = "sges", 
                       HEADCIRCUMvar = "pc", 
                       BATCHvar = "Sample_Plate", 
                       SamplePlacement = NULL,
                       BIRTHLENGTHvar = "talla", 
                       destinationfolder = "./EPIC", 
                       IDATdir = "/data/RAW_DATA/EPIC", 
                       cohort = "INMA", 
                       analysisdate = "03022021", 
                       savelog = T)
#the funtion will create a folder named: cohort_analysisdate_Output
#inside we will store the rds objects
saveRDS(idat, "./EPIC/cohort_analysisdate_Output/idat.RDS")

#CLOSE THE R SESSION AND START A NEW ONE

#Step 3: Exploratory analysis
idat <- readRDS("./EPIC/cohort_analysisdate_Output/idat.RDS")

#you can set any covariates of your interest
globalvarexplore <- c('sex', 'Cohort', 'ethnic_origin2c',
                      'Sample_Group', 'Sentrix_ID')

#Have a look at the reference manual for the package and adapt it to your data
exp_idat <- ExploratoryDataAnalysis(RGset = idat, 
                                    cohort = "INMA", 
                                    destinationfolder = "./EPIC/PACEanalysis/cohort_analysisdate_Output",
                                    analysisdate =  "03022021", 
                                    globalvarexplore = globalvarexplore)

#the funtion will create a folder named EDA/ inside cohort_analysisdate_Output
saveRDS(exp_idat, "./EPIC/cohort_analysisdate_Output/exp_idat_diff_samples.RDS")

#CLOSE THE R SESSION AND START A NEW ONE

#Step 4: Perform Sample Pre-processing
idat <- readRDS("./EPIC/cohort_analysisdate_Output/idat.RDS")
exp_idat <- readRDS("./EPIC/cohort_analysisdate_Output/exp_idat_diff_samples.RDS")

#Have a look at the reference manual for the package and adapt it to your data
preproc <- preprocessingofData(RGset = idat,
                               SamplestoRemove = exp_idat$SamplestoRemove,
                               ProbestoRemove = exp_idat$ProbestoRemove,
                               DetectionPvals = exp_idat$DetectionPval,
                               compositeCellType = 'Placenta',
                               KchooseManual = NULL,
                               cohort = "INMA",
                               analysisdate = "03022021",
                               destinationfolder = "./EPIC/PACEanalysis/cohort_analysisdate_Output")

saveRDS(preproc, "./EPIC/cohort_analysisdate_Output/preproc_exp_idat.RDS")