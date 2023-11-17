# Cis-mQTL mapping protocol for placental methylome

In this GitHub repository, you will find the protocol elaborated by the Immunogenetics Research Lab (IRLab) from the University of the Basque Country (UPV/EHU), to map placental cis-mQTLs using the command-line program TensorQTL. On the one hand, all the commands and scripts used are available in the Readme but be careful, you will need to custom some of them, mainly the Rscripts executed from RStudio and not from the command line. **Also, during the pipeline, we will be creating new directories for the outputs, but you should always work from outside of them! Have a look at the [diagram below](https://github.com/ariadnacilleros/Cis-mQTL-mapping-protocol-for-methylome/blob/main/directory_diagram_v2.jpg).** 

![Image of the working directory](https://github.com/ariadnacilleros/Cis-mQTL-mapping-protocol-for-methylome/blob/main/directory_diagram_v2.jpg)

On the other hand, in the [Wiki](https://github.com/ariadnacilleros/Cis-mQTL-mapping-protocol-for-methylome/wiki), you will find each step explained in more detail. To complete the pipeline, you need to have installed: [R and RStudio](https://rstudio-education.github.io/hopr/starting.html), [Plink1.9](https://www.cog-genomics.org/plink/1.9/) and [Plink2](https://www.cog-genomics.org/plink/2.0/), [TensorQTL v1.0.5 and its dependencies](https://github.com/broadinstitute/tensorqtl), [bcftools](http://samtools.github.io/bcftools/bcftools.html) and [vcftools](http://vcftools.sourceforge.net/index.html).

**In case you already performed the HRC genotype imputation at the Michigan Server and the preprocessing of the methylome data using PACEanalysis R package by A. Binder with the arguments specified in Step 2.1., you can jump to Step 2.2.**

## Step 1. Genotype data quality control 

### Step 1.1. Quality control of genotype data

Create a directory for the bfiles: 

`mkdir inter`

Conversion of long format files to PLINK binary format file: 

`plink1.9 --file {your filename} --make-bed --out inter/whole_genotype`

**You should adapt this command depending on the [input file](https://www.cog-genomics.org/plink/1.9/input) that do you have. Take in consideration how do you want PLINK to read your sample names as FID and IID.**

To know if the sex of the samples is reported on the [fam file](https://www.cog-genomics.org/plink/1.9/formats#fam): 

`head inter/whole_genotype.fam`

If the fifth column of the file contains only 0, you will need to add the sex of the samples adapting the following script to the sample sheet that you may have with the sex of the samples reported: \
[add-sex.R](https://github.com/ariadnacilleros/Cis-mQTL-mapping-protocol-for-methylome/blob/main/add-sex.R) 

Calculate frequencies: 

`plink1.9 --bfile inter/whole_genotype --freq --out inter/whole_genotype`

Download and unzip [Haplotype Reference Consortium (HRC) site list](http://www.haplotype-reference-consortium.org/site): 
```
wget ftp://ngs.sanger.ac.uk/production/hrc/HRC.r1-1/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz
gunzip HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz
```

Execute [Will Rayner's pre-imputation check script](https://www.well.ox.ac.uk/~wrayner/tools/): 
```
perl HRC-1000G-check-bim.pl -b inter/whole_genotype.bim -f inter/whole_genotype.frq -r HRC.r1-1.GRCh37.wgs.mac5.sites.tab -h -v
```

Execute Will Rayner’s output bash script: 
```
plink1.9 --bfile inter/whole_genotype --exclude inter/Exclude-whole_genotype-HRC.txt --make-bed --out inter/TEMP1
plink1.9 --bfile inter/TEMP1 --update-map inter/Chromosome-whole_genotype-HRC.txt --update-chr --make-bed --out inter/TEMP2
plink1.9 --bfile inter/TEMP2 --update-map inter/Position-whole_genotype-HRC.txt --make-bed --out inter/TEMP3
plink1.9 --bfile inter/TEMP3 --flip inter/Strand-Flip-whole_genotype-HRC.txt --make-bed --out inter/TEMP4
plink1.9 --bfile inter/TEMP4 --a2-allele inter/Force-Allele1-whole_genotype-HRC.txt --make-bed --out inter/whole_genotype-updated
rm inter/TEMP*
```

Change rsIDs of SNPs to ‘chr:position’ format with [change-rsid.R](https://github.com/ariadnacilleros/Cis-eQTL-mapping-protocol-for-methylome/blob/main/change-rsid.R): 

`Rscript change-rsid.R inter/whole_genotype-updated.bim`

#### Step 1.1.1. Filter SNPs

Make a directory for Quality Control: 

`mkdir qc`

Calculate frequencies: 

`plink1.9 --bfile inter/whole_genotype-updated --freq --out inter/whole_genotype-updated`

Calculate missing call rate: 

`plink1.9 --bfile inter/whole_genotype-updated --missing --out inter/whole_genotype-updated`

Remove markers by MAF/geno (missing call rate)/HWE thresholds: 
```
plink1.9 --bfile inter/whole_genotype-updated --geno 0.05 --hwe 1e-06 --maf 0.01 --make-bed --out qc/wg-updated-marker
```


#### Step 1.1.2. Filter samples

Calculate heterozygosity: 

`plink1.9 --bfile qc/wg-updated-marker --het --out qc/wg-updated-marker`

Plot missing call rate vs heterozygosity and subset individuals with > ± 4 x standard deviation (SD) using [imiss-vs-het.R](https://github.com/ariadnacilleros/Cis-eQTL-mapping-protocol-for-methylome/blob/main/imiss-vs-het.R): 

`Rscript imiss-vs-het.R inter/whole_genotype-updated.imiss qc/wg-updated-marker.het`

Remove selected individuals (> ± 4 x SD): 

`plink1.9 --bfile qc/wg-updated-marker --remove filter-het.txt --make-bed --out qc/wg-updated-marker-rmhet` 

Calculate missing call rate: 

`plink1.9 --bfile qc/wg-updated-marker-rmhet --missing --out qc/wg-updated-marker-rmhet`

Remove individuals with 0.03 missing markers: 

`plink1.9 --bfile qc/wg-updated-marker-rmhet --mind 0.03 --make-bed --out qc/wg-updated-marker-rmhet-ind`

Obtain the genotype sex from the samples: 
```
plink1.9 --bfile qc/wg-updated-marker-rmhet-ind --check-sex --out qc/wg-updated-marker-rmhet-ind
```

On the previous step, once we have calculated the sex of the samples with the [--check-sex](https://www.cog-genomics.org/plink/1.9/basic_stats#check_sex) flag, we should see the consistency that our reported sex has with the new genotype sex. 
```
awk '{if ($3 != $4) print $1,"\t",$2,"\t",$3,"\t",$4}' qc/wg-updated-marker-rmhet-ind.sexcheck >> qc/sex_inconsistencies.txt
```

In case of observing sex inconcistencies, we recommend you to talk with the researchers who obtained the samples and prepared the databases before removig those samples, to see if there could be a mistake. Depending on each case, you might ant to reassign the sex of the sample in the .fam file with [--make-just-fam](https://www.cog-genomics.org/plink/2.0/data#make_just_pvar) flag or discard the sample using [--remove](https://www.cog-genomics.org/plink/1.9/filter#indiv) flag from PLINK.


Calculate relatedness by [Identity-by-descent (IBD)](https://www.cog-genomics.org/plink/1.9/ibd): 

`plink1.9 --bfile qc/wg-updated-marker-rmhet-ind --genome --make-bed --out qc/wg-updated-marker-rmhet-ind`

Plot IBD values and subset individuals with PI_HAT > 0.18 with [plot-IBD.R](https://github.com/ariadnacilleros/Cis-eQTL-mapping-protocol-for-methylome/blob/main/plot-IBD.R): 

`Rscript plot-IBD.R qc/wg-updated-marker-rmhet-ind`

Select one sample per pair (with lower genotyping freq.) to remove with [rm-pihat018.R](https://github.com/ariadnacilleros/Cis-eQTL-mapping-protocol-for-methylome/blob/main/rm-pihat018.R): 
```
Rscript rm-pihat018.R qc/wg-updated-marker-rmhet-ind-fail-IBD-check.txt qc/wg-updated-marker-rmhet.imiss qc/rmpihat018.txt
```

Remove one from each pair: 

`plink1.9 --bfile qc/wg-updated-marker-rmhet-ind --remove qc/rmpihat018.txt --make-bed --out qc/clean-PIHAT`

Calculate PCs:

`plink1.9 --bfile qc/clean-PIHAT --indep-pairwise 50 5 0.2 --out qc/clean-PIHAT-prunned`

Perform PCA:

`plink1.9 --bfile qc/clean-PIHAT --extract qc/clean-PIHAT-prunned.prune.in --pca --out qc/clean-PIHAT.PCs`

Plot PCs with individuals information (e.g. sex), the following script contains the main commands, but you should adapt it to your data:\
[plot-PCA.R](https://github.com/ariadnacilleros/Cis-eQTL-mapping-protocol-for-methylome/blob/main/plot-PCA.R)

#### Step 1.1.3. Prepare data for imputation

Make a directory for chromosome files:

`mkdir chr`

Obtain separate VCF files per chromosome:
```
for i in {1..23};
do
  # Split chromosomes
  plink1.9 --bfile qc/clean-PIHAT --make-bed --chr $i --out chr/clean-PIHAT-chr$i
  # Make VCF files per chromosome
  plink1.9 --bfile chr/clean-PIHAT-chr$i --recode vcf --out chr/chr$i
  # Sort and compress VCF files
  vcf-sort chr/chr${i}.vcf | bgzip -c > chr/chr$i.vcf.gz ;
done
```

### Step 1.2. Impute SNP genotypes in the [Michigan Imputation Server](https://imputationserver.sph.umich.edu/index.html#!)
-	Reference Panel: HRC r1. 2016 (GRCh37/hg19)
-	Array Build: GRCh37/hg19
-	rsq Filter: off
-	Phasing: Eagle v2.4 (phased output)
-	Population: EUR
-	Mode: Quality Control & Imputation 

We will also check AES 256 encryption. 

### Step 1.3. Post-imputation quality control
To download the imputed genotype data from the server, use the commands indicated in the [website](https://imputationserver.readthedocs.io/en/latest/getting-started/#download-results) (e.g.  curl -sL https://imputationserver.sph.umich.edu/get/...) from inside the following folder: 
```
mkdir chr_imp
cd chr_imp
```
Once you have donwloaded the files inside `chr_imp` folder, get out from it by: 

`cd ../`

Decompress downloaded folders with the corresponding password (you will receive it by email at the address you used for registration):
```
for i in {1..22}; do 7z x -p'PASSWORD' chr_imp/chr_${i}; done
7z x -p'PASSWORD' chr_imp/chr_X
```
Removing zips: 

`rm -r chr_imp/*.zip`

Create folder for filtered VCFs: 

`mkdir imputed-rsq09`

Filter SNPs by Rsq(R2) > 0.9: 
```
for i in {1..22};
do
bcftools filter -i 'R2>0.9' -o imputed-rsq09/chr${i}.vcf.gz -Oz chr_imp/chr${i}.dose.vcf.gz &
done
bcftools filter -i 'R2>0.9' -o imputed-rsq09/chrX.vcf.gz -Oz chr_imp/chrX.dose.vcf.gz
```

Make a list of the VCF files to be concatenated: 
```
for i in {1..22};
do
echo imputed-rsq09/chr${i}.vcf.gz >> imputed-rsq09/tmp-concat.txt
done
echo imputed-rsq09/chrX.vcf.gz >> imputed-rsq09/tmp-concat.txt
```

Concatenate VCF files: 
```
bcftools concat -f imputed-rsq09/tmp-concat.txt --threads 16 -o imputed-rsq09/chrALL.vcf.gz -Oz imputed-rsq09/tmp-concat.txt
```

Execute [Will Rayner’s post-imputation QC script](https://www.well.ox.ac.uk/~wrayner/tools/Post-Imputation.html). You can use the scripts you will find available in the repository, but it is preferable to have a look at the documentation at the previous link and install the latest version of the whole package to avoid errors: 
```
mkdir qc-vcfparse
mkdir qc-ic
perl vcfparse.pl -d chr_imp -o qc-vcfparse
perl ic.pl -d qc-vcfparse -r HRC.r1-1.GRCh37.wgs.mac5.sites.tab -o qc-ic
```

Obtain the final list of samples of the genotype file: 

`bcftools query -l imputed-rsq09/chrALL.vcf.gz >> imputed-rsq09/sample_list_geno.txt`



## Step 2. Methylome data quality control
### Step 2.1. Betas values quality control 
The first step to be done with the methylation data is the quality control and the preprocess of it. To do so, we will use the PACEanalysis R package from A.Binder that must be already installed. The functions that we will execute take a long time and the inputs and the outputs are heavy, therefore we strongly recommend you to consider these advices: 
- Execute the functions in an R session from the command line, **avoid using RStudio**. 
- Use a **different R sessions** per function to be executed, saving the R objects with `saveRDS` at the end of it, and loading them with `readRDS` in the next step. If not, your R session could be killed because of the object's size. 
- As much as possible, use **`screen` method** to execute the functions, in case of losing the connection with the server, the jobs will continue running. 

For the output of these steps, we will create a new directory: 

`mkdir EPIC` 

This script contains **the functions and the specific arguments that must be executed**, but, as always, you will need to adapt them to your data. \
[PACEanalysis.R](https://github.com/ariadnacilleros/Cis-mQTL-mapping-protocol-for-methylome/blob/main/PACEanalysis.R)


### Step 2.2. Prepare BED file for TensorQTL mapping
In the next R script, you will find the commands used by our group to obtain the final BED file for our data into a .txt, some of them can be used directly, but others will need an adaptation to your data or won't be needed. The main steps are: 
- Make sure the same samples are present in the methylation and genotype files.
- Make sure the same sample names are the same in methylation and genotype files. *Important: in the genotype file, sample name is the IID, not FID_IID!*
- Annotate the CpGs by chr, start and end using the Illumina's R package. **Important, to define cis-window: end = start + 1**
- Filter CpGs located on sexual chromosomes. 

The output of this step should be a text file with the CpGs in the rows and the chr, start, end, CpG ID and beta values per sample in the columns. Here you have an [example](https://github.com/ariadnacilleros/Cis-mQTL-mapping-protocol-for-methylome/blob/main/example_bed_file_format.txt).\
[GSet_to_BED.R](https://github.com/ariadnacilleros/Cis-mQTL-mapping-protocol-for-methylome/blob/main/GSet_to_BED.R)

Sort BED file: 
```
(head -n1 EPIC/methylome_BED.txt && sort -k1,1V -k2,2n -k3,3n <(tail -n+2 EPIC/methylome_BED.txt)) > EPIC/methylome_sorted.bed
```

Zip BED file: 

`bgzip EPIC/methylome_sorted.bed` 



## Step 3. Obtain final format for genotype data

Create folder for the final version of the genotype: 

`mkdir whole_genome_definitive`

Filter VCF by MAF > 5%, HWE p-value > 0.05 and by the final list of samples: 
```
plink1.9 --vcf  imputed-rsq09/chrALL.vcf.gz --maf 0.05 --hwe 0.05 --keep EPIC/final_list_samples.txt --keep-allele-order --make-bed --out whole_genome_definitive/whole_genome_maf05_filt_samples
```

*Sometimes --keep doesn't work and keeps all the samples except the ones listed in final_list_samples.txt, if this is the case, change --keep by --remove:* 
```
plink1.9 --vcf imputed-rsq09/chrALL.vcf.gz --maf 0.05 --hwe 0.05 --remove EPIC/final_list_samples.txt --keep-allele-order --make-bed --out whole_genome_definitive/whole_genome_maf05_filt_samples
```

Be careful, when PLINK converts a VCF to a binary PLINK file set, it subsets the name of the samples from the VCF into FID and IID on the binary plink file (.bim, .fam, .bed) by searching for a separator which by default is __ . In case of troubles with this, we provide links to instruction on how to [customize the way PLINK reads sample names](https://www.cog-genomics.org/plink/2.0/input#sample_id_convert) or [how to change the name of the samples once you already have the binary PLINK file set](https://www.cog-genomics.org/plink/1.9/data#update_indiv). In our case, we tend to use [--const-fid](https://www.cog-genomics.org/plink/1.9/input#double_id) flag, which allows you to set FID as 0 in all the samples and the IID as the whole sample name coming from the VCF, or [--double-id](), which causes both family and individual IDs to be set to the sample ID. Remember that to run TensorQTL, you should match the IID from PLINK with the sample name on the BED file from the methylome. 

An extra step that we will perform at this point is to calculate the homozygous and heterozygous counts, and get the linkage disequilibrium information for each SNP: 
```
plink1.9 --bfile whole_genome_definitive/whole_genome_maf05_filt_samples --freqx --keep-allele-order --out whole_genome_definitive/whole_genome_maf05_filt_samples
plink1.9 --bfile whole_genome_definitive/whole_genome_maf05_filt_samples --r2 --keep-allele-order --out whole_genome_definitive/whole_genome_maf05_filt_samples
```


## Step 4. Prepare the covariates file for TensorQTL mapping

In this analysis, the covariates that we are going to use are the sex of the individuals, the Planet values (cell type proportions), and the Principal Components of our genotype. The number of PCs is going to be defined by your data, but at least you should include the first five. In this protocol we have used five PCs as an example. Therefore, we had to perform a Principal Component Analysis (PCA) with PLINK considering the last version of the genotype: 

```
mkdir covariates

plink1.9 --bfile whole_genome_definitive/whole_genome_maf05_filt_samples --indep-pairwise 50 5 0.2 --out covariates/whole_genome_maf05_filt_samples_prunned

plink1.9 --bfile whole_genome_definitive/whole_genome_maf05_filt_samples --extract covariates/whole_genome_maf05_filt_samples_prunned.prune.in --pca --out covariates/whole_genome_maf05_filt_samples_prunned.PCs 

```

The format for the covariates should be a text file in which the first line corresponds to the IID of the sample, and the next rows the covariates as in this [example](https://github.com/ariadnacilleros/Cis-mQTL-mapping-protocol-for-methylome/blob/main/example_covariates_file.txt). In the following script are the main commands used by our group to obtain the text file: \
[covariates.R](https://github.com/ariadnacilleros/Cis-mQTL-mapping-protocol-for-methylome/blob/main/covariates.R)


## Step 5. Additional analysis with rank-based inverse normal transformed beta values

We will also run two additional models; these will use the rank-based inverse normal transformed beta values (RNT). As covariates, the first model will only consider the known confounders from the previous model (sex, genotype PCs and Planet cell type estimations), and the second model will additionally take in count methylation PCs (mPCs) calculated in the residualized methylation data, with the aim of correcting for extra technical variation. We have generated the code to perform these steps based on outputs that you have already obtained in the previous steps. We hope that this will ensure the simplicity of these extra analyses and their execution.

To perform these extra analyses, we will create a new folder: 

`mkdir rnt_model`

### Step 5.1. Obtention of the RNT values and the corresponding bed file

To get the RNT values from the methylation beta values, we have prepared a new R script based on [GSet_to_Bed.R](https://github.com/ariadnacilleros/Cis-mQTL-mapping-protocol-for-methylome/blob/main/GSet_to_BED.R) from [Step 2.2](https://github.com/ariadnacilleros/Cis-mQTL-mapping-protocol-for-methylome/tree/main#step-22-prepare-bed-file-for-tensorqtl-mapping). The new R script allows you to transform the methylation values and to get the preliminary BED file with the CpGs in the rows and the chr, start, end, CpG ID and RNT values per sample in the columns.\
[GSet_to_Bed_RNT.R](https://github.com/ariadnacilleros/Cis-mQTL-mapping-protocol-for-methylome/blob/main/GSet_to_BED_RNT.R)

Please notice that the transformation of the methylation values is in line 40 of the script, and the rest remains as it was. There is no need to calculate the variance of these values as we did with the betas.  

Once you have obtained the BED file, you have to sort it and zip it as we did previously:

```
(head -n1 rnt_model/methylome_BED_RNT.txt && sort -k1,1V -k2,2n -k3,3n <(tail -n+2 rnt_model/methylome_BED_RNT.txt)) > rnt_model/methylome_sorted_RNT.bed

bgzip rnt_model/methylome_sorted_RNT.bed
```

### Step 5.2. Get residualized mPCs 

To get the mPCs avoiding any correlation with the known confounders (sex, genotype PCs, and very especially, Planet cell type estimations) and therefore, multicollinearity, we will perform the PCA on the residuals from multiple linear regression. This regression will consider the transformed methylation values as the outcome, and the known confounders as predictors. Once we get the residuals, we will perform a PCA on them, and select all mPCs that cumulatively explain > 80% of the variance. **ATTENTION: Once we have filtered the mPCs by the variance, we will keep a maximum number of 20 mPCs.** All these steps are included in the following script, which should be easy to execute since all the inputs have been obtained in previous steps. **ATTENTION: For the multiple linear model you should use as many genetic PCs as you have considered at the end of [Step 4]().** \
[get_residualized_mPCs.R](https://github.com/ariadnacilleros/Cis-mQTL-mapping-protocol-for-methylome/blob/main/get_residualized_mPCs.R)

Make sure that all the variables from the multiple regression are detected as numeric, this will allow you to get a matrix of residuals (Nº samples x Nº CpGs) with continuous values. 

### Step 5.3. Perform a GWAS with the mPCs

We also want to avoid potential correlations of the mPCs with the genotype, to be sure of this, we will perform a GWAS with our final genotype ([Step 3](https://github.com/ariadnacilleros/Cis-mQTL-mapping-protocol-for-methylome/tree/main#step-3-obtain-final-format-for-genotype-data)) and each of the mPCs calculated in the previous step as phenotypes. We have prepared a script that will allow you to create a new PLINK fam file per mPC and perform the GWAS. \
[gwas_mPCs.R](https://github.com/ariadnacilleros/Cis-mQTL-mapping-protocol-for-methylome/blob/main/gwas_mPCs.R)

In the case that any SNP is associated with a particular mPC (p-value < 1e-7), we will exclude that mPC. 

### Step 5.4. Get covariate file

Finally, we will include the selected mPCs in a new covariate file. With this aim, we have developed a simple script, where we add the mPCs to the previous covariates.txt file obtained in [Step 4](https://github.com/ariadnacilleros/Cis-mQTL-mapping-protocol-for-methylome/tree/main#step-4-prepare-the-covariates-file-for-tensorqtl-mapping). \
[covariates_RNT.R](https://github.com/ariadnacilleros/Cis-mQTL-mapping-protocol-for-methylome/blob/main/covariates_RNT.R)

The output from this script should be the same table as covariates.txt, but containing the values from the selected mPCs. 


## Step 6. Mapping with [TensorQTL](https://github.com/broadinstitute/tensorqtl)

### Step 6.1. Compute mQTLs with the main model 

Create a folder for TensorQTL results: 

`mkdir tensorQTL`

Open python3 module: 

`python3`

Mapping cis-mQTLs using covariates with TensorQTL: 
```
#Load packages
import pandas as pd
import torch
import tensorqtl
from tensorqtl import genotypeio, cis, trans
print('PyTorch {}'.format(torch.__version__))
print('Pandas {}'.format(pd.__version__))

#Define paths to data
plink_prefix_path = 'whole_genome_definitive/whole_genome_maf05_filt_samples'
expression_bed = 'EPIC/methylome_sorted.bed.gz'
covariates_file = 'covariates/covariates.txt'
prefix = 'tensorQTL/maf05_hwe05_nominal_INMA_18022021_A' #For this variable, read below this code block
```
The `prefix` variable, should follow the pattern: 

`cis_tensorQTL_maf05_hwe05_NOMINAL_(cohort)_(ddmmyyyy)_(model)`

For example, if model A of the analysis had been performed by the INMA cohort on 18/02/21, the prefix variable should contain: 

`cis_tensorQTL_maf05_hwe05_NOMINAL_INMA_18022021_A`

```
#Load phenotypes and covariates
phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(expression_bed)
covariates_df = pd.read_csv(covariates_file, sep='\t', index_col=0).T

#PLINK reader for genotypes
pr = genotypeio.PlinkReader(plink_prefix_path)
genotype_df = pr.load_genotypes()
variant_df = pr.bim.set_index('snp')[['chrom', 'pos']]

#Sort phenotype sample names
phenotype_df = phenotype_df.reindex(sorted(phenotype_df.columns), axis=1)
genotype_df = genotype_df.reindex(sorted(genotype_df.columns), axis=1)
covariates_df = covariates_df.sort_index()

#Run TensorQTL
cis.map_nominal(genotype_df, variant_df,
                phenotype_df,
                phenotype_pos_df,
                prefix, covariates_df=covariates_df, window=500000)
 
 ```                
The results will be written in your working directory as parquet files (one per chromosome). 

### Step 6.2. Compute mQTLs with RNT values 

For the computation of the mQTLs for these additional models, we will also work on the root directory, and the prefix variable will specify "RNT_wo_mPCs" for the model adjusted by only the known confounders, and “RNT” for the model adjusted by both the known confounders and the residualized mPCs. For example, if the analyses had been performed in the INMA cohort on 18/02/21, the prefix variable should contain: \
`cis_tensorQTL_maf05_hwe05_NOMINAL_INMA_18022021_RNT_wo_mPCs` \
`cis_tensorQTL_maf05_hwe05_NOMINAL_INMA_18022021_RNT`

Here is the code with the input files and the command to run TensorQTL for the **"RNT" model**. *Please note that now, the covariate file and the bed file are as follows: covariates_RNT.txt and methylome_sorted_RNT.bed.gz.* 

```
#Load packages
import pandas as pd
import torch
import tensorqtl
from tensorqtl import genotypeio, cis, trans
print('PyTorch {}'.format(torch.__version__))
print('Pandas {}'.format(pd.__version__))

#Define paths to data
plink_prefix_path = 'whole_genome_definitive/whole_genome_maf05_filt_samples'
expression_bed = 'rnt_model/methylome_sorted_RNT.bed.gz'
covariates_file = 'rnt_model/covariates_RNT.txt'
prefix = 'tensorQTL/maf05_hwe05_nominal_INMA_18022021_RNT' 

#Load phenotypes and covariates
phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(expression_bed)
covariates_df = pd.read_csv(covariates_file, sep='\t', index_col=0).T

#PLINK reader for genotypes
pr = genotypeio.PlinkReader(plink_prefix_path)
genotype_df = pr.load_genotypes()
variant_df = pr.bim.set_index('snp')[['chrom', 'pos']]

#Sort phenotype sample names
phenotype_df = phenotype_df.reindex(sorted(phenotype_df.columns), axis=1)
genotype_df = genotype_df.reindex(sorted(genotype_df.columns), axis=1)
covariates_df = covariates_df.sort_index()

#Run TensorQTL
cis.map_nominal(genotype_df, variant_df,
                phenotype_df,
                phenotype_pos_df,
                prefix, covariates_df=covariates_df, window=500000)
```

Finally, here is the code with the input files and the command to run TensorQTL for the **"RNT_wo_RNT" model**. *Please note that now, the covariate file and the bed file are as follows: covariates.txt and methylome_sorted_RNT.bed.gz.* 

```
#Load packages
import pandas as pd
import torch
import tensorqtl
from tensorqtl import genotypeio, cis, trans
print('PyTorch {}'.format(torch.__version__))
print('Pandas {}'.format(pd.__version__))

#Define paths to data
plink_prefix_path = 'whole_genome_definitive/whole_genome_maf05_filt_samples'
expression_bed = 'rnt_model/methylome_sorted_RNT.bed.gz'
covariates_file = 'covariates/covariates.txt'
prefix = 'tensorQTL/maf05_hwe05_nominal_INMA_18022021_RNT_wo_mPCs' 

#Load phenotypes and covariates
phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(expression_bed)
covariates_df = pd.read_csv(covariates_file, sep='\t', index_col=0).T

#PLINK reader for genotypes
pr = genotypeio.PlinkReader(plink_prefix_path)
genotype_df = pr.load_genotypes()
variant_df = pr.bim.set_index('snp')[['chrom', 'pos']]

#Sort phenotype sample names
phenotype_df = phenotype_df.reindex(sorted(phenotype_df.columns), axis=1)
genotype_df = genotype_df.reindex(sorted(genotype_df.columns), axis=1)
covariates_df = covariates_df.sort_index()

#Run TensorQTL
cis.map_nominal(genotype_df, variant_df,
                phenotype_df,
                phenotype_pos_df,
                prefix, covariates_df=covariates_df, window=500000)
```


**Again, the results will be written in the tensorQTL folder as parquet files (one per chromosome).**


## Step 7. Send the results

Finally, you have to send us the following files: 
- TensorQTL results (`cis_tensorQTL_maf05_hwe05_NOMINAL_(cohort)_(ddmmaaaa)_(model).cis_qtl_pairs.chr{1:22}.parquet`) 
- CpGs variability information (`all_cpg_variances.txt`)
- SNPs MAF and counts information (`whole_genome_maf05_filt_samples.frqx`)
- SNPs LD information (`whole_genome_maf05_filt_samples.ld`)

### Step 7.1. Send the results of the RNT analysis

For the RNT analyses, you only have to send the TensorQTL results of this part of the README: `cis_tensorQTL_maf05_hwe05_NOMINAL_(cohort)_(ddmmaaaa)_RNT.cis_qtl_pairs.chr{1:22}.parquet`
`cis_tensorQTL_maf05_hwe05_NOMINAL_(cohort)_(ddmmaaaa)_RNT_wo_mPCs.cis_qtl_pairs.chr{1:22}.parquet`

