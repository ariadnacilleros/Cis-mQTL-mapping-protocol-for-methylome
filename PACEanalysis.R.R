#Objective: Preprocess the methylation data.
#Last version: 14-07-2021
#Contact: acilleros001@ikasle.ehu.eus

#load library and set working directory
setwd("./EPIC")
library(PACEanalysis)

#Have a look at the reference manual for the package and adapt it to your data

# ==== LOADING DATA ====

#Step 1: Load your phenodata (sample information file/sample sheet )file with the basenames 
pheno <- read.delim("./Phenodata_basename_EPIC_BMI_all_names.txt", header = T, sep="\t") #adapt it

#Step 2: Loading Samples
cohort <- 'INMA' 
destf <- './EPIC' 
analysisdate <- '03022021' # YearMonthDay

idat <- loadingSamples(PhenoData = pheno,
                       IDlink = 'FID_IID',
                       BWTvar = 'peso',
                       BATCHvar = 'Sample_Plate',
                       SEXvar = 'sex',
                       FemaleInd = 'female',
                       MaleInd = 'male',
                       ETHNICvar = 'ethnic_origin2c',
                       GESTvar = 'sges',
                       BIRTHLENGTHvar = 'talla',
                       HEADCIRCUMvar = 'pc',
                       IDATdir = '/data/RAW_DATA/EPIC',
                       destinationfolder = destf,
                       cohort = cohort,
                       analysisdate = analysisdate)

saveRDS(idat, "./EPIC/cohort_analysisdate_Output/idat.RDS")

#CLOSE THE R SESSION AND START A NEW ONE

# ==== EXPLORATORY DATA ANALYSIS ====

#Step 3: Exploratory analysis
idat <- readRDS("./EPIC/cohort_analysisdate_Output/idat.RDS")

#exploratory variables
globalvarexplore <- c('sex', 'Cohort', 'ethnic_origin2c',
                      'Sample_Group', 'Sentrix_ID')

#Have a look at the reference manual for the package and adapt it to your data
exp_idat <- ExploratoryDataAnalysis(RGset = idat,
                                    globalvarexplore = globalvarexplore,
                                    DetectionPvalMethod = 'SeSAMe',
                                    DetectionPvalCutoff = 0.05,
                                    minNbeads = 3,
                                    FilterZeroIntensities = TRUE,
                                    destinationfolder = destf,
                                    savelog = TRUE,
                                    cohort = cohort,
                                    analysisdate = analysisdate)

saveRDS(exp_idat, "./EPIC/cohort_analysisdate_Output/exp_idat_diff_samples.RDS")

#CLOSE THE R SESSION AND START A NEW ONE

# ==== PREPROCESSING ====

# -- Samples to remove --
samplestr <- read.csv('EPIC/cohort_analysisdate_Output/EDA/INMA_03022021_Recommended_Samples_to_Remove.csv')


# Duplicated samples
clustered <- c('203740910098_R08C01', '203748260039_R07C01', '203748260158_R01C01',
               '203751250089_R08C01', '203740920065_R08C01')

# Sex inconsistencies
sexinc <- samplestr[samplestr$Sex_Wrong == 'Yes', 'Basename']

# Too many failed probes
failprobes <- samplestr[samplestr$TooManyFailedProbes == 'Yes', 'Basename']

# Contamination
contamination <- samplestr[samplestr$Meanlog2oddsContamination > -1, 'Basename']

# Concatenating all samples to remove into a single vector
str <- unique(c(clustered, sexinc, failprobes, contamination))

#Step 4: Perform Sample Pre-processing
idat <- readRDS("./EPIC/cohort_analysisdate_Output/idat.RDS")

preproc <- preprocessingofData(RGset = idata,
                               SamplestoRemove = str,
                               ProbestoRemove = eda$ProbestoRemove,
                               compositeCellType = 'Placenta',
                               KchooseManual = NULL,
                               cohort = cohort,
                               analysisdate = analysisdate,
                               destinationfolder = destf)

saveRDS(preproc, "./EPIC/cohort_analysisdate_Output/preproc_exp_idat.RDS")

#CLOSE THE R SESSION AND START A NEW ONE

# ==== OUTLIERS ====


#Step 5: Treat outliers values
preproc <- readRDS("./EPIC/cohort_analysisdate_Output/preproc_exp_idat.RDS")

outlier <- outlierprocess(processedBetas = preproc$processedBetas,
                          quantilemethod = "EmpiricalBeta",
                          trimming = FALSE,
                          pct = 0.005,
                          destinationfolder = destf,
                          cohort = cohort,
                          analysisdate = analysisdate)

saveRDS(object = outlier, file = './EPIC/cohort_analysisdate_Output/outlier.rds')