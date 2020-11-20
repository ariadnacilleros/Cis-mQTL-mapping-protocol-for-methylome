# Cis-eQTL mapping protocol for methylome

(Introduction)

## Step 1. Genotype data quality control 

### Step 1.1. Pre-imputation quality control

Conversion long format files to PLINK binary format file: \
`plink1.9 --file {filename} --make-bed --out {filename}`

Change rsIDs from SNPs to ‘chr:position’ format with [change-rsid.R](https://github.com/ariadnacilleros/Cis-eQTL-mapping-protocol-for-methylome/blob/main/change-rsid.R): \
`Rscript change-rsid.R {filename}.bim`

Calculate frequencies: \
`plink1.9 --bfile {filename} --freq --out {filename}`

Execute Will Rayner’s script: \
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
`plink1.9 --bfile {file}-marker --het`

Plot missing call rate vs heterozygosity and subset individuals with > ± 4 x standard deviation (SD) using [imiss-vs-het.R](https://github.com/ariadnacilleros/Cis-eQTL-mapping-protocol-for-methylome/blob/main/imiss-vs-het.R): \
`Rscript imiss-vs-het.R {filename}.imiss {filename}.het`

Remove selected individuals (> ± 4 x SD): \
`plink1.9 --bfile {filename}-marker --remove filter-het.txt --make-bed --out {filename}-rmhet` 

Remove individuals with 0.03 missing markers: \
`plink1.9 --bfile {filename}-rmhet --mind 0.03 --make-bed --out {filename}-ind`

Check sex concordance: \
`plink1.9 --bfile {filename}-ind --check-sex --out {filename}-ind`

Split mothers and rest:
```
plink1.9 --bfile {filename}-ind –keep {text file with mother IDs} --make-bed –out {file}-mothers
plink1.9 --bfile {filename}-ind --remove {text file with mother IDs} --make-bed --out {file}-nomothers
```

Calculate relatedness by IBD: 
```
plink1.9 --bfile {file}-mothers --genome --make-bed --out {file}-mothers-IBD
plink1.9 --bfile {file}-nomothers --genome --make-bed --out {file}-nomothers-IBD
```

Plot IBD values and subset individuals with PI_HAT > 0.18 with [plot-IBD.R](https://github.com/ariadnacilleros/Cis-eQTL-mapping-protocol-for-methylome/blob/main/plot-IBD.R):
```
Rscript plot-IBD.R {file}-mothers-IBD
Rscript plot-IBD.R {file}-nomothers-IBD
```
Remove individuals PI_HAT > 0.18 w/less genotype:
-	Open the file with related individual pairs (...fail-IBD-check.txt)
-	Make a list of all samples involved (PIHAT018.txt) 

Select one sample per pair (with lower genotyping freq.) to remove with [rm-pihat018.R](https://github.com/ariadnacilleros/Cis-eQTL-mapping-protocol-for-methylome/blob/main/rm-pihat018.R): 
```
Rscript rm-pihat018.R {file}-fail-IBD-check.txt {filename}.imiss rmpihat018-mothers.txt
Rscript rm-pihat018.R {file}-fail-IBD-check.txt {filename}.imiss rmpihat018-nomothers.txt
```
Remove one from each pair:
```
plink1.9 --bfile {file}-mothers-IBD --remove rmpihat018-mothers.txt --make-bed --out {file}-mothers-rmpihat
plink1.9 --bfile {file}-nomothers-IBD --remove rmpihat018-nomothers.txt --make-bed --out {file}-nomothers-rmpihat
```
Merge both files again: 
```
echo {file}-mothers-rmpihat > mergelist2
echo {file}-nomothers-rmpihat >> mergelist2
plink1.9 --merge-list mergelist2 –out merged-cohorts-clean
```
Calculate PCs:\
`plink1.9 --bfile merged-cohorts-clean --indep-pairwise 50 5 0.2 --out merged-cohorts-clean-prunned`

Perform PCA:\
`plink1.9 --bfile merged-cohorts-clean --extract {filename}-clean-prunned.prune.in --pca --out merged-cohorts-clean.PCs`

Plot PCs with individuals information (e.g. sex):\
[plot-PCA.R](https://github.com/ariadnacilleros/Cis-eQTL-mapping-protocol-for-methylome/blob/main/plot-PCA.R)

#### Step 1.1.3. Prepare data for imputation
```
for i in {1..23};
do
  # Split chromosomes
  plink1.9 --bfile merged-cohorts-clean --reference-allele Force-Allele1-{filename}.txt  \
  --make-bed --chr $i --out merged-cohorts-clean-chr$i
  # Make VCF files per chromosome
  plink1.9 --bfile merged-cohorts-clean-chr$i --recode vcf --out chr$i
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
To download the imputed genotype data from the server, use the commands indicated in the website from it (e.g.  curl -sL https://imputationserver.sph.umich.edu/get/...). \
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
Execute Will Rayner’s post-imputation QC: 
```
perl vcfparse.pl -d {directory path with imputed VCF files} -o {output directory name}
perl ic.pl -d {output directory from vcfparse} -r HRC.r1-1.GRCh37.wgs.mac5.sites.tab -o {output directory name}
```
Filter by Rsq chromosomes from 1 to 22: 
```
for i in {1..22};
do
wc -l chr${i}.info
awk '{if($7 > 0.9) {print $1}}' chr${i}.info | tail -n +2 >filt${i}.snps
wc -l filt${i}.snps
vcftools --gzvcf  chr${i}.dose.vcf.gz --snps filt${i}.snps --recode --out chr-filtered-${i} &
done
```
Filter by Rsq chromosome X:
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
In the next R script, you will find the commands used for us to obtain the final text file filtered with all the annotation data of the CpGs and the samples, to use them as examples or guidelines: \
`BED_UCSC_GRSet.R` 

Sort BED file: 
```
(head -n1 whole_genome_imp_var_bed.txt && sort -k1,1V -k2,2n -k3,3n <(tail -n+2 whole_genome_imp_var_bed.txt)) > whole_genome_var_sorted.bed
```

Zipp BED file: \
`bgzip whole_genome_var_sorted.bed` 

Index BED file: \
`tabix -p bed whole_genome_var_sorted.bed.gz`

***Jump to 2.4. Calculate statistical power to filter genotype data***

### Step 2.3. Filter samples of genotype data
Filter by methylation samples: \
`bcftools view -S IID_methylome.txt --force-samples concat-allchr.vcf -o concat-allchr-metfilt.vcf`

Obtain samples names for methylome filtering: \
`bcftools query -l concat-allchr-metfilt.vcf  > samples_imp_vcf.txt`

***Return to BED_UCSC_GRSet.R***

### Step 2.4. Calculate statistical power to filter genotype data
Calculate the statistical power of your data and decide the MAF to filter genotype, but you will need to change some parameters depending on your statistics: \
[power.R](https://github.com/ariadnacilleros/Cis-eQTL-mapping-protocol-for-methylome/blob/main/power.R)

Filter VCF by the MAF selected by you from the power.pdf plot. For example, in our case we choose 8% of MAF: \
`plink2 --vcf concat-allchr-metfilt.vcf --maf 0.08 --recode vcf --out whole_genome_imp_maf`

Compress and index VCF: 
```
bgzip Geno_imp/whole_genome_imp_maf.vcf	
tabix -p vcf Geno_imp/whole_genome_imp_maf.vcf.gz
```

## Step 3. Mapping [FastQTL](http://fastqtl.sourceforge.net/)
Change timestamps from index files: 
```
touch whole_genome_var_sorted.bed.gz.tbi
touch whole_genome_imp_maf.vcf.gz.tbi
```
Mapping cis-mQTLs with FastQTL:
```
for j in $(seq 1 1000); do fastQTL --vcf Geno_imp/whole_genome_imp_maf.vcf.gz --bed EPIC/whole_genome_var_sorted.bed.gz --cov Covariates/COV.txt --permute 1000 10000 --seed 123456789 --out permutations.imp.${j}.txt.gz --chunk $j 1000;done
```
## Extra step: 
Create file and write inside the name of the CpGs or SNPs: \
`nano file.exc`

Run FastQTL for chunk 807 excluding CpGs: 
```
fastQTL --vcf  whole_genome_imp_maf.vcf.gz  --bed whole_genome_var_sorted.bed.gz --cov COV.txt --permute 1000 10000 --seed 123456789 --out permutations.imp.807.txt.gz --chunk 807 1000 --exclude-phenotypes file.exc
```
Run FastQTL for chunk 807 excluding SNPs: 
```
fastQTL --vcf  whole_genome_imp_maf.vcf.gz  --bed whole_genome_var_sorted.bed.gz --cov COV.txt --permute 1000 10000 --seed 123456789 --out permutations.imp.807.txt.gz --chunk 807 1000 --exclude-sites file.exc
```
