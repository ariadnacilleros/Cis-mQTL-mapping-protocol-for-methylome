# Cis-mQTL mapping protocol for placental methylome

In this GitHub repository, you will find the protocol elaborated by the Immunogenetics Research Lab (IRLab) from the University of the Basque Country (UPV/EHU), to map placental cis-mQTLs using the command-line program FastQTL. On the one hand, all the commands and scripts used are available in the Readme, but you will need to custom them depending on your data. On the other hand, in the Wiki you will find each step explained in more detail. To complete the pipeline, you need to have installed: [R and RStudio](https://rstudio-education.github.io/hopr/starting.html), [Plink1.9](https://www.cog-genomics.org/plink/1.9/) and [Plink2](https://www.cog-genomics.org/plink/2.0/), [FastQTL](http://fastqtl.sourceforge.net/), [bcftools](http://samtools.github.io/bcftools/bcftools.html) and [vcftools](http://vcftools.sourceforge.net/index.html).
**In case that you already performed the genotype imputation by Michigan Server and the preprocessment of the methylome by Alexandra Binder’s R package, you can jump to Step 2.2.**

<p align="center">
<img src="https://github.com/ariadnacilleros/Cis-eQTL-mapping-protocol-for-methylome/blob/main/Imp_scheme_vertical.jpg" width="500">
</p>

## Step 1. Genotype data quality control 

### Step 1.1. Quality control of genotype data

Conversion of long format files to PLINK binary format file: \
`plink1.9 --file {filename} --make-bed --out {filename}`

Change rsIDs of SNPs to ‘chr:position’ format with [change-rsid.R](https://github.com/ariadnacilleros/Cis-eQTL-mapping-protocol-for-methylome/blob/main/change-rsid.R): \
`Rscript change-rsid.R {filename}.bim`

Calculate frequencies: \
`plink1.9 --bfile {filename} --freq --out {filename}`

Execute [pre-imputation Will Rayner’s script](https://www.well.ox.ac.uk/~wrayner/tools/): \
`perl HRC-1000G-check-bim.pl -b {filename}.bim -f {filename}.frq -r HRC.r1-1.GRCh37.wgs.mac5.sites.tab -h -v`

Execute Will Rayner’s bash script output: \
`bash Run-plink.sh`

Add missing sex like in the following R script, but adapt it to your inputs: \
[add-sex.R](https://github.com/ariadnacilleros/Cis-eQTL-mapping-protocol-for-methylome/blob/main/add-sex.R)

Merge sexual chromosomes: \
`plink1.9 --bfile {filename} --merge-x no-fail --make-bed --out {filename}-mxy` 

#### Step 1.1.1. Filter SNPs

Calculate frequencies: \
`plink1.9 --bfile {filename}-mxy --freq --out {filename}-mxy`

Calculate missing call rate: \
`plink1.9 --bfile {filename}-mxy --missing --out {filename}-mxy`

Remove markers by MAF/geno (missing call rate)/HWE thresholds: \
`plink2 --bfile {filename}-mxy --geno 0.05 --hwe 1e-06 --maf 0.01 --make-bed --out {filename}-marker`

#### Step 1.1.2. Filter samples

Calculate heterozygosity: \
`plink1.9 --bfile {filename}-marker --het`

Plot missing call rate vs heterozygosity and subset individuals with > ± 4 x standard deviation (SD) using [imiss-vs-het.R](https://github.com/ariadnacilleros/Cis-eQTL-mapping-protocol-for-methylome/blob/main/imiss-vs-het.R): \
`Rscript imiss-vs-het.R {filename}.imiss {filename}.het`

Remove selected individuals (> ± 4 x SD): \
`plink1.9 --bfile {filename}-marker --remove filter-het.txt --make-bed --out {filename}-rmhet` 

Remove individuals with 0.03 missing markers: \
`plink1.9 --bfile {filename}-rmhet --mind 0.03 --make-bed --out {filename}-ind`

Check sex concordance: \
`plink1.9 --bfile {filename}-ind --check-sex --out {filename}-ind`

Calculate relatedness by [Identity-by-descent (IBD)](https://www.cog-genomics.org/plink/1.9/ibd): \
`plink1.9 --bfile {filename}-ind --genome --make-bed --out {filename}-IBD`

Plot IBD values and subset individuals with PI_HAT > 0.18 with [plot-IBD.R](https://github.com/ariadnacilleros/Cis-eQTL-mapping-protocol-for-methylome/blob/main/plot-IBD.R): \
`Rscript plot-IBD.R {filename}-IBD`

Remove individuals PI_HAT > 0.18 w/less genotype:
-	Open the file with related individual pairs (...fail-IBD-check.txt)
-	Make a list of all samples involved (PIHAT018.txt) 

Select one sample per pair (with lower genotyping freq.) to remove with [rm-pihat018.R](https://github.com/ariadnacilleros/Cis-eQTL-mapping-protocol-for-methylome/blob/main/rm-pihat018.R): \
`Rscript rm-pihat018.R {filename}-fail-IBD-check.txt {filename}.imiss rmpihat018.txt`

Remove one from each pair: \
`plink1.9 --bfile {filename}-IBD --remove rmpihat018.txt --make-bed --out clean-PIHAT`

Calculate PCs:\
`plink1.9 --bfile clean-PIHAT --indep-pairwise 50 5 0.2 --out clean-PIHAT-prunned`

Perform PCA:\
`plink1.9 --bfile clean-PIHAT --extract {filename}-clean-prunned.prune.in --pca --out clean-PIHAT.PCs`

Plot PCs with individuals information (e.g. sex):\
[plot-PCA.R](https://github.com/ariadnacilleros/Cis-eQTL-mapping-protocol-for-methylome/blob/main/plot-PCA.R)

#### Step 1.1.3. Prepare data for imputation
```
for i in {1..23};
do
  # Split chromosomes
  plink1.9 --bfile clean-PIHAT --reference-allele Force-Allele1-{filename}.txt  \
  --make-bed --chr $i --out clean-PIHAT-chr$i
  # Make VCF files per chromosome
  plink1.9 --bfile clean-PIHAT-chr$i --recode vcf --out chr$i
  # Sort and compress VCF files
  vcf-sort chr${i}.vcf | bgzip -c > chr$i.vcf.gz ;
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
To download the imputed genotype data from the server, use the commands indicated in the website (e.g.  curl -sL https://imputationserver.sph.umich.edu/get/...). \
Decompress downloaded folders with the corresponding password (email):
```
for i in {1..22}; do unzip -P 'PASSWORD' chr_${i}; done
unzip -P 'PASSWORD' chr_X
```
Decompress info files: 
```
for i in {1..22}; do bgzip -d chr${i}.info.gz; done
bgzip -d chrX.info.gz
```
Execute [Will Rayner’s post-imputation QC](https://www.well.ox.ac.uk/~wrayner/tools/Post-Imputation.html): 
```
perl vcfparse.pl -d {directory path with imputed VCF files} -o {output directory name}
perl ic.pl -d {output directory from vcfparse} -r HRC.r1-1.GRCh37.wgs.mac5.sites.tab -o {output directory name}
```
Filter SNPs by Rsq from chromosome 1 to 22: 
```
for i in {1..22};
do
wc -l chr${i}.info
awk '{if($7 > 0.9) {print $1}}' chr${i}.info | tail -n +2 >filt${i}.snps
wc -l filt${i}.snps
vcftools --gzvcf  chr${i}.dose.vcf.gz --snps filt${i}.snps --recode --out chr-filtered-${i} &
done
```
Filter SNPs by Rsq from chromosome X:
```
i=X
wc -l chr${i}.info
awk '{if($7 > 0.9) {print $1}}' chr${i}.info | tail -n +2 > filt${i}.snps
wc -l filt${i}.snps
vcftools --gzvcf chr${i}.dose.vcf.gz --snps filt${i}.snps --recode --out chr-filtered-${i}
```
Concatenate VCF files into one file: 
```
for i in {1..22};
do
echo chr-filtered-${i}.recode.vcf >> tmp.concat
done
echo chr-filtered-X.recode.vcf >> tmp.concat
```
Convert final file into VCF: \
`bcftools concat -f tmp.concat -o concat-allchr.vcf`

### Step 1.4. Prepare covariates file for FastQTL mapping
In this case our covariate for this analysis is the sex of the samples, and as we said, we took it from an excel sample sheet. Into the following R script you will see the commands that we used, but depending in your input, they could change, but the output text file should have the same format. \
[covariates.R](https://github.com/ariadnacilleros/Cis-eQTL-mapping-protocol-for-methylome/blob/main/covariates.R)

## Step 2. Methylome data quality control
### Step 2.1. Betas quality control 
Alexandra Binder’s R package to be defined. 

### Step 2.2. Prepare BED file for FastQTL mapping
In the next R script, you will find the commands used to obtain the final text file filtered with all the annotation data of the CpGs and the samples: \
[BED_UCSC_GRSet.R](https://github.com/ariadnacilleros/Cis-eQTL-mapping-protocol-for-methylome/blob/main/BED_UCSC_GRSet.R)

Sort BED file: 
```
(head -n1 methylome_var_bed.txt && sort -k1,1V -k2,2n -k3,3n <(tail -n+2 methylome_var_bed.txt)) > methylome_var_sorted.bed
```

Zipp BED file: \
`bgzip methylome_var_sorted.bed` 

Index BED file: \
`tabix -p bed methylome_var_sorted.bed.gz`

***Jump to 2.4. Calculate statistical power to filter genotype data***

### Step 2.3. Filter samples of genotype data
Filter by methylation samples: \
`bcftools view -S IID_methylome.txt --force-samples concat-allchr.vcf -o concat-allchr-metfilt.vcf`

Obtain sample names for methylome filtering: \
`bcftools query -l concat-allchr-metfilt.vcf  > samples_imp_vcf.txt`

***Return to BED_UCSC_GRSet.R***

### Step 2.4. Calculate statistical power to filter genotype data
Calculate the statistical power of your data and decide the MAF to filter genotype data. You will need to change some parameters depending on your statistics: \
[power.R](https://github.com/ariadnacilleros/Cis-eQTL-mapping-protocol-for-methylome/blob/main/power.R)

Filter VCF by the MAF selected (power.pdf plot). For example, in our case we choose 8% of MAF: \
`plink2 --vcf concat-allchr-metfilt.vcf --maf 0.08 --recode vcf --out whole_genome_imp_maf`

Compress and index VCF: 
```
bgzip whole_genome_imp_maf.vcf	
tabix -p vcf whole_genome_imp_maf.vcf.gz
```

## Step 3. Mapping [FastQTL](http://fastqtl.sourceforge.net/)
Change timestamps from index files: 
```
touch methylome_var_sorted.bed.gz.tbi
touch whole_genome_imp_maf.vcf.gz.tbi
```
Mapping cis-mQTLs with FastQTL:
```
for j in $(seq 1 1000); do fastQTL --vcf whole_genome_imp_maf.vcf.gz --bed methylome_var_sorted.bed.gz --cov COV.txt --permute 1000 10000 --seed 123456789 --out permutations.imp.${j}.txt.gz --chunk $j 1000;done
```
Zipp all chunks into one file: \
`zcat permutations.imp.*.txt.gz | gzip -c > permutations.all.chunks.txt.gz`

## Extra step: 
In the case of having an error in one chunk due to low variable CpGs and/or SNPs like the following: 
```
  * Number of variants in cis = 28499
  * Best correlation = 0.5440
  * Number of permutations = 53 / 1000
gsl: beta.c:44: ERROR: domain error
Default GSL error handler invoked.
Aborted
```

Create file and write inside the name of the problematic CpGs or SNPs: \
`nano file.exc`

Run FastQTL for chunk 807 excluding CpGs: 
```
fastQTL --vcf  whole_genome_imp_maf.vcf.gz  --bed methylome_var_sorted.bed.gz --cov COV.txt --permute 1000 10000 --seed 123456789 --out permutations.imp.807.txt.gz --chunk 807 1000 --exclude-phenotypes file.exc
```
Run FastQTL for chunk 807 excluding SNPs: 
```
fastQTL --vcf  whole_genome_imp_maf.vcf.gz  --bed methylome_var_sorted.bed.gz --cov COV.txt --permute 1000 10000 --seed 123456789 --out permutations.imp.807.txt.gz --chunk 807 1000 --exclude-sites file.exc
```
