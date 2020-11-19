# Cis-eQTL mapping protocol for methylome

(Introduction)

## Step 1. Genotype data quality control 

### Step 1.1. Pre-imputation quality control

Conversion long format files to PLINK binary format file: 
`plink1.9 --file {filename} --make-bed --out {filename}`

Change rsIDs from SNPs to ‘chr:position’ format: 
`Rscript change-rsid.R`

Calculate frequencies: 
`plink1.9 --bfile {filename} --freq --out {filename}`

Execute Will Rayner’s script:   
`perl HRC-1000G-check-bim.pl -b {filename}.bim -f {filename}.frq -r HRC.r1-1.GRCh37.wgs.mac5.sites.tab -h -v`

Execute Will Rayner’s bash script output: 
`bash Run-plink.sh`

Add missing sex using excel sheet:
`Rscript add-sex.R`

Merge sexual chromosomes: 
`plink1.9 --bfile {filename} --merge-x no-fail --make-bed --out {filename}-mxy`



#### Step 1.1.1. Filter SNPs

Calculate frequencies: 
`plink1.9 --bfile {filename}-mxy --freq --out {filename}-mxy`

Calculate missing call rate: 
`plink1.9 --bfile {filename}-mxy --missing --out {filename}-mxy`

Remove markers by MAF/geno (missing call rate)/HWE thresholds: 
`plink2 --bfile {filename}-mxy --geno 0.05 --hwe 1e-06 --maf 0.01 --make-bed --out {filename}-marker`


