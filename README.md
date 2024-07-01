# Local-adaptation-in-sugar-maple

Workflow for detection of adaptive genetic variation associated with climate related synthetic variables in 24 sugar maple natural populations. 

# Data
In the SRA (project) you will found access to 641 individual fastqs to use as input for this workflow.

# Mapping 
for f in reads/*R1.fq.gz; do
    sample=`basename $f _R1.fq.gz`
    echo "Trimming $sample"
    fastp -i reads/${sample}_R1.fq.gz -I reads/${sample}_R2.fq.gz \
          -o trimmed/${sample}_trimmed_R1.fq.gz \
          -O trimmed/${sample}_trimmed_R2.fq.gz \
          -w 24 -p -c -l 20 -M 30 -r --detect_adapter_for_pe \
          -h trimmed/${sample}_trimming.html \
          -j trimmed/${sample}_trimming.json -V
    bwa mem -t 24 -R "@RG\tID:${sample}\tSM:${sample}" \
            reference/Aesc.1_0.fa.gz \
            trimmed/${sample}_trimmed_R1.fq.gz \
            trimmed/${sample}_trimmed_R2.fq.gz |
            samtools sort -n -O sam - | samtools fixmate -m -O bam - - |
            samtools sort -O bam - | samtools view -h -b -f 3 - > \
            mapping/${sample}.mapped.sorted.concordant.bam
    bamtools index -in mapping/${sample}.mapped.sorted.concordant.bam
done

# Statistics on bams
for f in mapping/*.bam;
        do sample=`basename $f .bam`
        echo "stats on $sample "
        samtools flagstats $f > bamstats/$sample.stats
done
cd bamstats/
echo -e "sample\ttotal\tproperly_paired\tmapped\tduplicates\tsupp" > bamstats.txt
for f in *.stats; do
        sample=`basename $f .mapped.sorted.concordant.stats`
        total=`awk '{print $1}' $f | sed '1q;d'`
        paired=`awk '{print $1}' $f | sed '12q;d'`
        mapped=`awk '{print $1}' $f | sed '7q;d'`
        duplicates=`awk '{print $1}' $f | sed '5q;d'`
        supp=`awk '{print $1}' $f | sed '4q;d'`
        echo -e "$sample\t$total\t$paired\t$mapped\t$duplicates\t$supp" >> bamstats.txt
done

# Variant calling
ref_map.pl --samples mapping/ --popmap popmap.txt -o stacks -X populations:"-p 5 -r 0.85 --fstats --write-random-snp --vcf"

# Genotyping error
vcftools --vcf stacks/populations.snps.vcf --maf 0.01 --out genotyping_error/populations.snps.maf001 --recode --remove-filtered-all
vcftools --vcf stacks/populations.snps.vcf --maf 0.05 --out genotyping_error/populations.snps.maf005 --recode --remove-filtered-all
vcftools --vcf stacks/populations.snps.vcf --maf 0.005 --out genotyping_error/populations.snps.maf0005 --recode --remove-filtered-all
vcftools --vcf stacks/populations.snps.vcf --maf 0.001 --out genotyping_error/populations.snps.maf0001 --recode --remove-filtered-all

cd genotyping_error
split -l 2 replicates.txt

for i in x*; do vcftools --vcf populations.snps.maf001.recode.vcf --keep $i --out maf001/$i.populations.snps.001 --recode --remove-filtered-all; done
for i in x*; do vcftools --vcf populations.snps.maf005.recode.vcf --keep $i --out maf005/$i.populations.snps.005 --recode --remove-filtered-all; done
for i in x*; do vcftools --vcf npopulations.snps.maf0001.recode.vcf --keep $i --out maf0001/$i.populations.snps.0001 --recode --remove-filtered-all; done
for i in x*; do vcftools --vcf populations.snps.maf0005.recode.vcf --keep $i --out maf0005/$i.populations.snps.0005 --recode --remove-filtered-all; done

for i in maf001/*.populations.snps.001.recode.vcf; do vcftools --vcf $i --extract-FORMAT-info GT --out $i; done
for i in maf005/*.populations.snps.005.recode.vcf; do vcftools --vcf $i --extract-FORMAT-info GT --out $i; done
for i in maf0001/*.populations.snps.0001.recode.vcf; do vcftools --vcf $i --extract-FORMAT-info GT --out $i; done
for i in maf0005/*.populations.snps.0005.recode.vcf; do vcftools --vcf $i --extract-FORMAT-info GT --out $i; done

for i in maf001/*.populations.snps.001.recode.vcf.GT.FORMAT; do awk '{ if ($3 == $4) { print "same"; } else { print "different"; } }' $i > $i.diff; done
for i in maf005/*.populations.snps.005.recode.vcf.GT.FORMAT; do awk '{ if ($3 == $4) { print "same"; } else { print "different"; } }' $i > $i.diff; done
for i in maf0001/*.populations.snps.0001.recode.vcf.GT.FORMAT; do awk '{ if ($3 == $4) { print "same"; } else { print "different"; } }' $i > $i.diff; done
for i in maf0005/*.populations.snps.0005.recode.vcf.GT.FORMAT; do awk '{ if ($3 == $4) { print "same"; } else { print "different"; } }' $i > $i.diff; done

for i in maf001/*.populations.snps.001.recode.vcf.GT.FORMAT.diff; do grep -c "different" $i >> random_snps.diff_genotypes.001.txt; done
for i in maf005/*.populations.snps.005.recode.vcf.GT.FORMAT.diff; do grep -c "different" $i >> random_snps.diff_genotypes.005.txt; done
for i in maf0001/*.populations.snps.0001.recode.vcf.GT.FORMAT.diff; do grep -c "different" $i >> random_snps.diff_genotypes.0001.txt; done
for i in maf0005/*.populations.snps.0005.recode.vcf.GT.FORMAT.diff; do grep -c "different" $i >> random_snps.diff_genotypes.0005.txt; done

# filtering
vcftools --vcf stacks/populations.snps.vcf --min-meanDP 20 --out filters/populations.snps.DP20 --recode --remove-filtered-all
vcftools --vcf filters/populations.snps.DP20.recode.vcf --max-meanDP 500 --out filters/populations.snps.DP20-500 --recode --remove-filtered-all
vcftools --vcf filters/populations.snps.DP20-500.recode.vcf --max-missing 0.85 --out filters/populations.snps.DP20-500.miss15 --recode --remove-filtered-all
vcftools --vcf filters/populations.snps.DP20-500.miss15.recode.vcf --missing-indv --out filters/populations.snps.DP20-500.miss15
awk -F"\t" '$5>0.15' filters/populations.snps.DP20-500.miss15.imiss | cut -f 1 | sed '1d' > filters/exclude.txt
vcftools --vcf filters/populations.snps.DP20-500.miss15.recode.vcf --remove filters/exclude.txt --out filters/populations.snps.DP20-500.miss15-15 --recode --remove-filtered-all
vcftools --vcf filters/populations.snps.DP20-500.miss15-15.recode.vcf --maf 0.001 --out filters/populations.snps.DP20-500.miss15-15.0001 --recode --remove-filtered-all

# pruning
## making chrom map file for plink conversion (because plink will change name of chromosoms, I choose that it uses those that I want)
bcftools query -f '%CHROM\n' filters/populations.snps.DP20-500.miss15-15.0001.recode.vcf > chrom_maf0001.txt
grep acsa chrom_maf0001.txt | uniq > chrom_maf0001_list

a=0
for i in `cat chrom_maf0001_list`; do a=$((a+1)); echo -e "$i\t$a" >> chrom_maf0001_map; done

## Convert to PLINK raw format
vcftools --vcf filters/populations.snps.DP20-500.miss15-15.0001.recode.vcf --out pruning/populations.snps.20500_1515_0001 --plink --chrom-map chrom_maf0001_map

## prune linked snps
plink --file pruning/populations.snps.20500_1515_0001  --indep-pairwise 50 10 0.1 --allow-extra-chr --out plink-prune
plink --file pruning/populations.snps.20500_1515_0001  --extract plink-prune.prune.in --make-bed --out pruning/populations.snps.20500_1515_0001.prun --allow-extra-chr

## extract pruned positions of initial vcf(take data.prun.bim and in EXCEL make correspondance between chrom numbers with chrom map and extract right chrom name with positions and convert it to pruned_positions_data.txt)
vcftools --vcf filters/populations.snps.DP20-500.miss15-15.0001.recode.vcf --positions pruning/pruned_positions_populations.snps.20500_1515_0001.txt  --out pruning/populations.snps.20500_1515_0001.prun  --recode --remove-filtered-all

## get imiss data to indentify replicates to keep
vcftools --vcf pruning/populations.snps.20500_1515_0001.prun.recode.vcf --missing-indv --out pruning/noRIMO01.random_snps.20500_1515_0001.prun

## from imiss created after_pruning choose which replicate to exclude (with more missingness), create a list remove.txt
vcftools --vcf pruning/populations.snps.20500_1515_0001.prun.recode.vcf --remove remove.txt --out final_data/npopulations.snps.20500_1515_0001.prun.norep --recode --remove-filtered-all

## from final vcf create bed files and raw file
## making chrom map file for plink conversion (because plink will change name of chromosoms, I choose that it uses those that I want)
## You can reuse chrom_map used for pruning
#rm -v chrom*
#bcftools query -f '%CHROM\n' final_data/noRIMO01.20500_1515_001.prun.norep.recode.vcf > chrom_001.txt
#grep acsa chrom_001.txt | uniq > chrom_001.list
#a=0
#for i in `cat chrom_001.list`; do a=$((a+1)); echo -e "$i\t$a" >> chrom_001.map; done

## Convert to PLINK format and make BED
vcftools --vcf final_data/noRIMO01.random_snps.20500_1515_0001.prun.norep.recode.vcf --out final_data/noRIMO01.random_snps.20500_1515_0001.prun.norep --plink --chrom-map chrom_maf0001_map
plink --file final_data/noRIMO01.random_snps.20500_1515_0001.prun.norep --make-bed --out final_data/noRIMO01.random_snps.20500_1515_0001.prun.norep --allow-extra-chr

## convert plink format to raw format
plink --bfile final_data/noRIMO01.random_snps.20500_1515_0001.prun.norep --recode A --out final_data/noRIMO01.random_snps.20500_1515_0001.prun.norep --allow-extra-chr

## output stats of final_data
vcftools --vcf final_data/noRIMO01.random_snps.20500_1515_0001.prun.norep.recode.vcf --depth --out vcf-stats/noRIMO01.random_snps.20500_1515_0001.prun.norep
vcftools --vcf final_data/noRIMO01.random_snps.20500_1515_0001.prun.norep.recode.vcf --site-mean-depth --out vcf-stats/noRIMO01.random_snps.20500_1515_0001.prun.norep
vcftools --vcf final_data/noRIMO01.random_snps.20500_1515_0001.prun.norep.recode.vcf --het --out vcf-stats/noRIMO01.random_snps.20500_1515_0001.prun.norep
vcftools --vcf final_data/noRIMO01.random_snps.20500_1515_0001.prun.norep.recode.vcf --relatedness --out vcf-stats/noRIMO01.random_snps.20500_1515_0001.prun.norep
vcftools --vcf final_data/noRIMO01.random_snps.20500_1515_0001.prun.norep.recode.vcf --window-pi 10000 --out vcf-stats/noRIMO01.random_snps.20500_1515_0001.prun.norep
vcftools --vcf final_data/noRIMO01.random_snps.20500_1515_0001.prun.norep.recode.vcf --TajimaD 10000 --out vcf-stats/noRIMO01.random_snps.20500_1515_0001.prun.norep

## run cross validation for choosing K
INDIR=~/Analyses/Alix/ERS-localadapt/article/final_data/
cd admixture
for K in {1..23}; do admixture --cv $INDIR/noRIMO01.random_snps.20500_1515_0001.prun.norep.bed $K | tee log${K}.out; done

## visualize CV errors
grep -h CV log*.out

OUT_DIR=pairs

## Subset Final VCF for each pair of populations (make lists before)
for i in $OUT_DIR/*.txt; do name=`basename $i .txt`; vcftools --vcf final_data/noRIMO01.random_snps.20500_1515_0001.prun.norep.recode.vcf --keep $i --out $OUT_DIR/$name.noRIMO01.random_snps.20500_1515_0001.prun.norep --recode --remove-filtered-all; done

## from final vcf create bed files and raw file
## making chrom map file for plink conversion (because plink will change name of chromosoms, I choose that it uses those that I want)

rm $OUT_DIR/chrom_*
for i in $OUT_DIR/*.vcf; do name=`basename $i .recode.vcf`;bcftools query -f '%CHROM\n' $i | grep acsa | uniq > $OUT_DIR/chrom_$name.list; done
for f in $OUT_DIR/chrom_*.list; do name=`basename $f .list`; a=0; for i in `cat $f`; do a=$((a+1)); echo -e "$i\t$a" >> $OUT_DIR/$name.map; done; done

## Convert to PLINK format and make BED
for i in $OUT_DIR/*.vcf; do name=`basename $i .recode.vcf`; vcftools --vcf $i --out $OUT_DIR/$name --plink --chrom-map $OUT_DIR/chrom_$name.map; plink --file $OUT_DIR/$name --make-bed --out $OUT_DIR/$name --allow-extra-chr; done

for i in `cat candidates/list`; do bedtools intersect -a ~/Data/Acer/Genome/Aesc.1_0/Aesc.1_0.gff.gz -b candidates/$i.bed -wo -F 0.1 > candidates/$i.gff.txt; done

