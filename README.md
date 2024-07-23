# Local-adaptation-in-sugar-maple

Workflow for detection of adaptive genetic variation associated with climate related synthetic variables in 24 natural populations of sugar maple (Acer saccharum Marshall) (N=641 individuals). 

## Data
In the SRA (BioProject PRJNA1130992) you will found access to 641 individual fastqs to use as input for this workflow.

## Bioinformatic steps from raw FASTQs to set of raw SNPs 
Execute script fastq2snps.sh

## Genotyping error
Calculate genotyping error for replicated individuals in raw SNPs filtered for four values different minor allele frequency (MAF).
Run genotyping_error.sh : bash genotyping_error.sh path/to/input.vcf path/to/replicates.txt path/to/output_directory
Replicates.txt should have the same format of the provided example.
Output of the script is a count of diffences in genotypes for each pair of replicates that you can use to calculate genotyping error rate per pairs of replicates par MAF values (diff_count/total number of SNPs)

## Filtering steps toward clean data set of SNPs
### First round of filters and pruning
Filters SNPs according to depth, missing data and MAF then pruns SNPs for linkage desequilibirum (LD)
Run filters_pruning.sh : bash path/to/input.vcf path/to/output_directory
### Create VCF of filtered and pruned SNPs 
Last step of filters_pruning.sh produces a bed file containing informations on genetoyprs of individuals for SNPs remaining after pruning. SNPs are listed into the associated .bim file. 
Next steps are filters on the VCF file of pruned positions. To create this VCF we need to extract pruned positions of VCF outputed by filters_pruning.sh.
To do so, create a "pruned_positions_data.txt" listing pruned SNPs in two tab-separated columns, no header : CHROM  POS. Make sure that CHROM are numbered/designated as in the VCF as plink renamed them during pruning.
Then use VCFTOOLS to extract pruned positions : 
vcftools --vcf filtered.vcf --positions pruned_positions_data.txt --out filtered.pruned --recode --remove-filtered-all
### Second round of filtering on pruned SNPs
The second round of filtering consists in : i) removing replicated individuals that were used for calculating genotyping error rates, ii) removing SNPs suspected of paralogy,  iii) removing individuals suspected to be clones.
  i) We are keeping one replicate out of two, the one with less missing data. So get missing data from filtered.pruned.vcf created in precedent step.  
  vcftools --vcf filtered.pruned.vcf --missing-indv --out filtered.pruned
  From filtered.pruned.imiss identify replicates to remove and list them in remove.txt file so that we can remove them from VCF. 
  vcftools --vcf filtered.pruned.vcf --remove remove.txt --out filtered.pruned.norep --recode --remove-filtered-all
  ii) Checking for paralogy necessitates GT info of individuals for SNPs :
  vcftools --vcf filtered.pruned.norep.vcf --extract-FORMAT-info GT --out filtered.pruned.norep
  From filtered.pruned.norep.GT.FORMAT file check for SNPs with unrealistic heterozygotes rate : per SNPs, count heterozygotes (o/1) and divide it per total number of individuals.
  List SNPs (CHROM\tPOS) with Het rate >= 85% in snps2remove.txt
  vcftools --vcf filtered.pruned.norep.vcf --exclude-positions snps2remove.txt --out filtered.pruned.norep.nohet --recode --remove-filtered-all
  iii) Clones will be excessively related. Compute relatedness values from filtered.pruned.norep.vcf to identify putative clones.
  vcftools --vcf filtered.pruned.norep.nohet.recode.vcf --relatedness2
  Check in out.relatedness2 for pairs of distinct individuals that have relatedness_phi value > 0.4. Make a list of putative clones and remove it from VCF
  vcftools --vcf filtered.pruned.norep.nohet.recode.vcf --remove clones.txt --out filtered.pruned.norep.nohet.noclones --recode --remove-filtered-all

## Output final datasets (filtered.pruned.norep.nohet.noclones.vcf) from  to serve as input for next analyses
### For all individuals : we need bed file for input to Admixture and .raw format for input in R adegenet
bash final_data.sh path/to/input.vcf path/to/output_directory
### For pairs of populations : we need bed files for input in R pcadapt
You first need to list indivduals in pairs of populations in txt files and store them in a folder that will serve as output directory in the next script.
bash pairs.sh path/to/input.vcf path/to/output_directory

## Admixture
### run cross validation for choosing K
cd admixture
for K in {1..23}; do admixture --cv ../path/to/filtered.pruned.norep.nohet.noclones.bed $K | tee log${K}.out; done
### visualize CV errors
grep -h CV log*.out

## VCF statistics
### output stats of final vcf
IN_VCF=path/to/filtered.pruned.norep.nohet.noclones.vcf
for i in {"depth","site-mean-depth"}; do name=`basename $IN_VCF .recode.vcf`; vcftools --vcf $IN_VCF --$i --out $name; done 

## Genetic statistics (He, Ho, Fis), PCA, admixture composition plots, Fst calculations and Mantel test, DAPC, pcadapt and RDA
Follow Rmarkdown

## Candidate lists
for i in `cat candidates/list`; do bedtools intersect -a ~/Data/Acer/Genome/Aesc.1_0/Aesc.1_0.gff.gz -b candidates/$i.bed -wo -F 0.1 > candidates/$i.gff.txt; done
