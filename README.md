# Cis-mQTL mapping protocol for placental methylome

In this GitHub repository, you will find the protocol elaborated by the Immunogenetics Research Lab (IRLab) from the University of the Basque Country (UPV/EHU), to map placental cis-mQTLs using the command-line program FastQTL. On the one hand, all the commands and scripts used are available in the Readme but be careful, you will need to custom some of them, mainly the Rscripts executed from RStudio and not from the command line. **Also, during the pipeline, we will be creating new directories for the outputs, but you should always work from outside of them! Have a look at the [diagram below](https://github.com/ariadnacilleros/Cis-mQTL-mapping-protocol-for-methylome/blob/main/directory_diagram.PNG).** 

![Image of the working directory](https://github.com/ariadnacilleros/Cis-mQTL-mapping-protocol-for-methylome/blob/main/directory_diagram.PNG)

On the other hand, in the Wiki, you will find each step explained in more detail. To complete the pipeline, you need to have installed: [R and RStudio](https://rstudio-education.github.io/hopr/starting.html), [Plink1.9](https://www.cog-genomics.org/plink/1.9/) and [Plink2](https://www.cog-genomics.org/plink/2.0/), [TensorQTL and its dependencies](https://github.com/broadinstitute/tensorqtl), [bcftools](http://samtools.github.io/bcftools/bcftools.html) and [vcftools](http://vcftools.sourceforge.net/index.html).

**In case that you already performed the genotype imputation by the Michigan Server and the preprocessing of the methylome by Alexandra Binder’s R package, you can jump to Step 2.2.**

## Step 1. Genotype data quality control 

### Step 1.1. Quality control of genotype data

Create a directory for the bfiles: \
`mkdir inter`

Conversion of long format files to PLINK binary format file: \
`plink1.9 --file {your filename} --make-bed --out inter/whole_genotype`

Change rsIDs of SNPs to ‘chr:position’ format with [change-rsid.R](https://github.com/ariadnacilleros/Cis-eQTL-mapping-protocol-for-methylome/blob/main/change-rsid.R): \
`Rscript change-rsid.R inter/whole_genotype.bim`

Calculate frequencies: \
`plink1.9 --bfile inter/whole_genotype --freq --out inter/whole_genotype`

Download and unzip [Haplotype Reference Consortium (HRC) site list](http://www.haplotype-reference-consortium.org/site): 
```
wget ftp://ngs.sanger.ac.uk/production/hrc/HRC.r1-1/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz
gunzip HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz
```

Execute [pre-imputation Will Rayner’s script](https://www.well.ox.ac.uk/~wrayner/tools/): 
```
perl HRC-1000G-check-bim.pl -b inter/whole_genome.bim -f inter/whole_genome.frq -r HRC.r1-1.GRCh37.wgs.mac5.sites.tab -h -v
```

Execute Will Rayner’s bash script output: 
```
plink1.9 --bfile inter/whole_genome --exclude inter/Exclude-whole_genome-HRC.txt --make-bed --out inter/TEMP1
plink1.9 --bfile inter/TEMP1 --update-map inter/Chromosome-whole_genome-HRC.txt --update-chr --make-bed --out inter/TEMP2
plink1.9 --bfile inter/TEMP2 --update-map inter/Position-whole_genome-HRC.txt --make-bed --out inter/TEMP3
plink1.9 --bfile inter/TEMP3 --flip inter/Strand-Flip-whole_genome-HRC.txt --make-bed --out inter/TEMP4
plink1.9 --bfile inter/TEMP4 --a2-allele inter/Force-Allele1-whole_genome-HRC.txt --make-bed --out inter/whole_genome-updated
rm inter/TEMP*
```

Make a directory for Quality Control: 
```
mkdir qc
```

Add missing sex like in the following R script, but adapt it to your inputs: \
[add-sex.R](https://github.com/ariadnacilleros/Cis-eQTL-mapping-protocol-for-methylome/blob/main/add-sex.R)

Merge sexual chromosomes: \
`plink1.9 --bfile inter/whole_genome-updated --merge-x no-fail --make-bed --out qc/wg-updated-mxy` 

#### Step 1.1.1. Filter SNPs

Calculate frequencies: \
`plink1.9 --bfile qc/wg-updated-mxy --freq --out qc/wg-updated-mxy`

Calculate missing call rate: \
`plink1.9 --bfile qc/wg-updated-mxy --missing --out qc/wg-updated-mxy`

Remove markers by MAF/geno (missing call rate)/HWE thresholds: \
`plink2 --bfile qc/wg-updated-mxy --geno 0.05 --hwe 1e-06 --maf 0.01 --make-bed --out qc/wg-updated-marker`

#### Step 1.1.2. Filter samples

Calculate heterozygosity: \
`plink1.9 --bfile qc/wg-updated-marker --het`

Plot missing call rate vs heterozygosity and subset individuals with > ± 4 x standard deviation (SD) using [imiss-vs-het.R](https://github.com/ariadnacilleros/Cis-eQTL-mapping-protocol-for-methylome/blob/main/imiss-vs-het.R): \
`Rscript imiss-vs-het.R qc/wg-updated-mxy.imiss qc/wg-updated-marker.het`

Remove selected individuals (> ± 4 x SD): \
`plink1.9 --bfile qc/wg-updated-marker --remove filter-het.txt --make-bed --out qc/wg-updated-rmhet` 

Remove individuals with 0.03 missing markers: \
`plink1.9 --bfile qc/wg-updated-rmhet --mind 0.03 --make-bed --out qc/wg-updated-ind`

Check sex concordance: \
`plink1.9 --bfile qc/wg-updated-ind --check-sex --out qc/wg-updated-ind`

Calculate relatedness by [Identity-by-descent (IBD)](https://www.cog-genomics.org/plink/1.9/ibd): \
`plink1.9 --bfile qc/wg-updated-ind --genome --make-bed --out qc/wg-updated-IBD`

Plot IBD values and subset individuals with PI_HAT > 0.18 with [plot-IBD.R](https://github.com/ariadnacilleros/Cis-eQTL-mapping-protocol-for-methylome/blob/main/plot-IBD.R): \
`Rscript plot-IBD.R qc/wg-updated-IBD`

Remove individuals PI_HAT > 0.18 w/less genotype:
-	Open the file with related individual pairs (qc/wg-updated-IBD-fail-IBD-check.txt)
-	Make a list of all samples involved (PIHAT018.txt) 

Select one sample per pair (with lower genotyping freq.) to remove with [rm-pihat018.R](https://github.com/ariadnacilleros/Cis-eQTL-mapping-protocol-for-methylome/blob/main/rm-pihat018.R): \
`Rscript rm-pihat018.R qc/wg-updated-IBD-fail-IBD-check.txt qc/wg-updated-mxy.imiss qc/rmpihat018.txt`

Remove one from each pair: \
`plink1.9 --bfile qc/wg-updated-IBD --remove qc/rmpihat018.txt --make-bed --out qc/clean-PIHAT`

Calculate PCs:\
`plink1.9 --bfile qc/clean-PIHAT --indep-pairwise 50 5 0.2 --out qc/clean-PIHAT-prunned`

Perform PCA:\
`plink1.9 --bfile qc/clean-PIHAT --extract qc/clean-PIHAT-prunned.prune.in --pca --out qc/clean-PIHAT.PCs`

Plot PCs with individuals information (e.g. sex), the following script contains the main commands, but you should adapt to your data:\
[plot-PCA.R](https://github.com/ariadnacilleros/Cis-eQTL-mapping-protocol-for-methylome/blob/main/plot-PCA.R)

#### Step 1.1.3. Prepare data for imputation

Make directory for chromosome files:  
```
mkdir chr
```

Obtain VCF file per chromosome:
```
for i in {1..23};
do
  # Split chromosomes
  plink1.9 --bfile qc/clean-PIHAT --reference-allele inter/Force-Allele1-whole_genome-HRC.txt  \
  --make-bed --chr $i --out chr/clean-PIHAT-chr$i
  # Make VCF files per chromosome
  plink1.9 --bfile chr/clean-PIHAT-chr$i --recode vcf --out chr/chr$i
  # Sort and compress VCF files
  vcf-sort chr/chr${i}.vcf | bgzip -c > chr/chr$i.vcf.gz ;
done
```

### Step 1.2. Impute genotype in [Michigan Imputation Server](https://imputationserver.sph.umich.edu/index.html#!)
-	Reference Panel: HRC r1. 2016 (GRCh37/hg19)
-	Array Build: GRCh37/hg19
-	rsq Filter: off
-	Phasing: Eagle v2.4 (phased output)
-	Population: EUR
-	Mode: Quality Control & Imputation 

We will also check AES 256 encryption. 

### Step 1.3. Post-imputation quality control
To download the imputed genotype data from the server, use the commands indicated in the [website](https://imputationserver.readthedocs.io/en/latest/getting-started/#download-results) (e.g.  curl -sL https://imputationserver.sph.umich.edu/get/...) inside the following folder: 
```
mkdir chr_imp
cd chr_imp
```
Once you donwloaded the files inside `chr_imp` folder, get out from it by: \
`cd ../`

Decompress downloaded folders with the corresponding password (email):
```
for i in {1..22}; do unzip -P 'PASSWORD' chr_imp/chr_${i}; done
unzip -P 'PASSWORD' chr_imp/chr_X
```
Decompress info files: 
```
for i in {1..22}; do bgzip -d chr_imp/chr${i}.info.gz; done
bgzip -d chr_imp/chrX.info.gz
```
Execute [Will Rayner’s post-imputation QC](https://www.well.ox.ac.uk/~wrayner/tools/Post-Imputation.html) which scripts you will find available on the repository, but it is preferable to have a look at the documentation of the previous link and install the whole package to avoid errors: 
```
mkdir qc-vcfparse
mkdir qc-ic
perl vcfparse.pl -d chr_imp -o qc-vcfparse
perl ic.pl -d qc-vcfparse -r HRC.r1-1.GRCh37.wgs.mac5.sites.tab -o qc-ic
```
Filter SNPs by Rsq from chromosome 1 to 22: 
```
mkdir qc-results
for i in {1..22};
do
wc -l chr_imp/chr${i}.info
awk '{if($7 > 0.9) {print $1}}' chr_imp/chr${i}.info | tail -n +2 > qc-results/filt${i}.snps
wc -l qc-results/filt${i}.snps
vcftools --gzvcf  chr_imp/chr${i}.dose.vcf.gz --snps qc-results/filt${i}.snps --recode --out qc-results/chr-filtered-${i} &
done
```
Filter SNPs by Rsq from chromosome X:
```
i=X
wc -l chr_imp/chr${i}.info
awk '{if($7 > 0.9) {print $1}}' chr_imp/chr${i}.info | tail -n +2 > qc-results/filt${i}.snps
wc -l qc-results/filt${i}.snps
vcftools --gzvcf chr_imp/chr${i}.dose.vcf.gz --snps qc-results/filt${i}.snps --recode --out qc_results/chr-filtered-${i}
```
Record the name of all the VCFs file to concatenate: 
```
mkdir qc-results/tmp-concat
for i in {1..22};
do
echo qc-results/chr-filtered-${i}.recode.vcf >> qc-results/tmp-concat/tmp.concat
done
echo qc-results/chr-filtered-X.recode.vcf >> qc-results/tmp-concat/tmp.concat
```
Concatenate VCF files into one file: \
`bcftools concat -f qc-results/tmp-concat/tmp.concat -o qc-results/concat-allchr.vcf`

Obtain the final list of samples of the genotype: \
`bcftools query -l qc-results/concat-allchr.vcf >> qc-results/sample_list_geno.txt`

## Step 2. Methylome data quality control
### Step 2.1. Betas quality control 
Alexandra Binder’s R package to be defined. 

### Step 2.2. Prepare BED file for FastQTL mapping
In the next R script, you will find the commands used to obtain the final text file filtered with all the annotation data of the CpGs and the samples. As you will see the following script contains the main commands used by our group to obtain a bed file for our data in a text file, using the data contained in an ExpressionSet R object which is the output from Alexandran Binder's R package. Some of the commands can be used directly, but others will need an adaptation to your data or won't be needed. The main steps are: 
- Have the same sample names between methylation and genotype. *Important: on the genotype, we understand as sample names the IID, not the FID_IID!*
- Remove duplicates (in case of necessity). 
- Have the same samples between methylation and genotype. 
- Filter CpGs with cross-hybridizing potential, with SNPs (European MAF < 5%) and the ones located in sexual chromosomes. *Important: here are the commands to download the CpG lists*
```
wget https://ars.els-cdn.com/content/image/1-s2.0-S221359601630071X-mmc1.txt
wget https://ars.els-cdn.com/content/image/1-s2.0-S221359601630071X-mmc2.txt
wget https://ars.els-cdn.com/content/image/1-s2.0-S221359601630071X-mmc3.txt
```
- Annotate the CpGs by chr, start and end. 
**In the last lines of the script, you will find the code to obtain the variability information of the CpGs that will be subset and sent to us once the mapping has been done.**

For the output of these steps, we will create a new directory: \
`mkdir EPIC` 

The output of this step should be a text file with the CpGs on the rows and the chr, start, end, CpG ID and beta values per sample on the columns. Here you have an [example](https://github.com/ariadnacilleros/Cis-mQTL-mapping-protocol-for-methylome/blob/main/example_bed_file_format.txt).\
[GSet_to_BED.R](https://github.com/ariadnacilleros/Cis-mQTL-mapping-protocol-for-methylome/blob/main/GSet_to_BED.R)

Sort BED file: 
```
(head -n1 EPIC/methylome_BED.txt && sort -k1,1V -k2,2n -k3,3n <(tail -n+2 EPIC/methylome_BED.txt)) > EPIC/methylome_sorted.bed
```

Zipp BED file: \
`bgzip EPIC/methylome_sorted.bed` 

Index BED file: \
`tabix -p bed EPIC/methylome_sorted.bed.gz`

## Step 3. Obtain final format for genotype data

Create folder for the final version of the genotype: \
`mkdir whole_genome_definitive`

Filter VCF by the 5% MAF and by the final list of samples: 
```
plink2 --vcf qc-results/concat-allchr.vcf --maf 0.05 --keep EPIC/final_list_samples.txt --make-bed --out whole_genome_definitive/whole_genome_maf05_filt_samples
```

*Sometimes --keep doesn't work and keeps all the samples except the ones listes in final_list_samples.txt, if this is the case, change --keep by --remove:* 
```
plink2 --vcf qc-results/concat-allchr.vcf --maf 0.05 --remove EPIC/final_list_samples.txt --make-bed --out whole_genome_definitive/whole_genome_maf05_filt_samples
```

Be careful, when PLINK converts a VCF to a binary PLINK file set, it subsets the name of the samples from the VCF into FID and IID on the binary plink file (.bim, .fam, .bed) by searching a separator which by default is _ . In case of having troubles with it, we leave here a link to [costumize the read of the sample names by PLINK](https://www.cog-genomics.org/plink/2.0/input#sample_id_convert) or [how to change the name of the samples once you already have the binary PLINK file set](https://www.cog-genomics.org/plink/1.9/data#update_indiv). In our case, we always set [--const-fid](https://www.cog-genomics.org/plink/1.9/input#double_id) flag which allows you to set FID as 0 in all the samples and the IID as the whole sample name coming from the VCF. Remember that to run TensorQTL, you should match the IID from PLINK with the sample name on the BED file from the methylome. 

An extra step that we will perform at this point is to calculate the homozygous and heterozygous counts for each SNP: 
```
plink --bfile whole_genome_definitive/whole_genome_maf05_filt_samples --freqx --out whole_genome_definitive/whole_genome_maf05_filt_samples
```


## Step 4. Prepare covariates file for TensorQTL mapping
In this analysis, the covariates that we are going to use are the sex of the samples and the first five Principal Components of our genotype. Therefore, we will need to perform a Principal Component Analysis (PCA) with PLINK: 

```
mkdir covariates

plink --bfile whole_genome_definitive/whole_genome_maf05_filt_samples --indep-pairwise 50 5 0.2 --out covariates/whole_genome_maf05_filt_samples_prunned --double-id

plink --bfile whole_genome_definitive/whole_genome_maf05_filt_samples --extract covariates/whole_genome_maf05_filt_samples_prunned.prune.in --pca --out covariates/whole_genome_maf05_filt_samples_prunned.PCs --double-id 

```

The format file for the covariates should be a text file in which the first line corresponds to the IID of the sample, being the next rows the other covariates as in this [example](https://github.com/ariadnacilleros/Cis-mQTL-mapping-protocol-for-methylome/blob/main/example_covariates_file.txt). In the following script are the main commands used to obtain the text file: \
[covariates_sex_PC5.R](https://github.com/ariadnacilleros/Cis-mQTL-mapping-protocol-for-methylome/blob/main/covariates_sex_PC5.R)

An extra step that could be done is to compute the sex of the samples from the genotype by [--check-sex](https://www.cog-genomics.org/plink/1.9/basic_stats#check_sex) and compare if this one matches with your notes. It is not clear if TensorQTL takes into count the sex of the samples provided by the .fam file of the binary PLINK set, in case that you want to take it into count for your analysis, we recommend you to have it described in both places, covariates text file and .fam file. 

## Step 5. Mapping with [TensorQTL](https://github.com/broadinstitute/tensorqtl)

Change timestamps from index files: \
`touch whole_genome_definitive/whole_genome_maf05_filt_samples.bed.gz.tbi`

Create a folder for TensorQTL results: \
`mkdir tensorQTL`

Open python3 module: \
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
covariates_file = 'covariates/covariates_sex_PC5.txt'

#Load phenotypes and covariates
phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(expression_bed)
covariates_df = pd.read_csv(covariates_file, sep='\t', index_col=0).T

#PLINK reader for genotypes
pr = genotypeio.PlinkReader(plink_prefix_path)
genotype_df = pr.load_genotypes()
variant_df = pr.bim.set_index('snp')[['chrom', 'pos']]

#Sort phenotype sample names
phenotype_df = phenotype_df.reindex(sorted(phenotype_df.columns), axis=1)

#Run TensorQTL
cis_df = cis.map_cis(genotype_df, variant_df, phenotype_df, phenotype_pos_df, covariates_df=covariates_df, seed=123456789)

#Write TensorQTL results
cis_df.to_csv('tensorQTL/cis_tensorQTL_maf05_PC5_sex.txt', header=True, index=True, sep='\t')
```
