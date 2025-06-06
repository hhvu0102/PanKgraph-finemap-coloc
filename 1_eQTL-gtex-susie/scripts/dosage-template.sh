#!/bin/bash

## fetch variants in the region and intersect UKBB and FUSION vcfs
for i in $t1_vcf; do tabix $$i $ukbb_fetch | awk '{if (($$0 !~ /^#/ && $$0 !~ /^chr/)) print "chr"$$0; else print $$0}' ; done | sort | uniq > ${susie_locus}.ukbb.genotypes
zcat $t1_vcf1 | head -10000 | awk '{if (($$0 ~ /^#/)) print $$0}' > ${susie_locus}.ukbb.header
cat ${susie_locus}.ukbb.header ${susie_locus}.ukbb.genotypes | bgzip -c > ${susie_locus}.ukbb.vcf.gz; tabix ${susie_locus}.ukbb.vcf.gz
rm ${susie_locus}.ukbb.genotypes ${susie_locus}.ukbb.header

## fetch UKBB dosages 
zcat ${susie_locus}.ukbb.vcf.gz | head -10000 | awk -F'\t' '{if (($$0 ~/^#CHROM/)) print $$0}' OFS='\t' | sed -e 's:#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT:ID:g' > ${susie_locus}.ukbb-header.txt 
bcftools query -f "%ID-%REF-%ALT[\t%DS]\n" ${susie_locus}.ukbb.vcf.gz | cat ${susie_locus}.ukbb-header.txt - > ${susie_locus}.ukbb-dosages.tsv 

## bgzip to save space
module load Bioinformatics
module load Bioinformatics  gcc/10.3.0-k2osx5y
module load samtools/1.13-fwwss5n

bgzip -@ 2 ${susie_locus}.ukbb-dosages.tsv

## cleanup
rm -rf ${susie_locus}.ukbb-header.txt ${susie_locus}.ukbb.vcf.gz*

