#!/bin/bash

## fetch variants in the region and intersect UKBB and FUSION vcfs
for i in $t1_vcf; do tabix $$i $ukbb_fetch | awk '{if (($$0 !~ /^#/ && $$0 !~ /^chr/)) print "chr"$$0; else print $$0}' ; done | sort | uniq > ${susie_locus}.ukbb.genotypes
zcat $t1_vcf1 | head -10000 | awk '{if (($$0 ~ /^#/)) print $$0}' > ${susie_locus}.ukbb.header
cat ${susie_locus}.ukbb.header ${susie_locus}.ukbb.genotypes | bgzip -c > ${susie_locus}.ukbb.vcf.gz; tabix ${susie_locus}.ukbb.vcf.gz
rm ${susie_locus}.ukbb.genotypes ${susie_locus}.ukbb.header

## fetch UKBB dosages 
zcat ${susie_locus}.ukbb.vcf.gz | head -10000 | awk -F'\t' '{if (($$0 ~/^#CHROM/)) print $$0}' OFS='\t' | sed -e 's:#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT:ID:g' > ${susie_locus}.ukbb-header.txt 
bcftools query -f "%ID-%REF-%ALT[\t%DS]\n" ${susie_locus}.ukbb.vcf.gz | cat ${susie_locus}.ukbb-header.txt - > ${susie_locus}.ukbb-dosages.tsv 

## fetch GWAS variants 
tabix -h $t1_summary $fetch > ${susie_locus}.gwas.tsv ;

## align GWAS alleles with UKBB reference and have consistent rsids
/scratch/scjp_root/scjp99/vthihong/2_PanKBase/colocGWAS_T1D/1_t1d-susie/scripts/align-gwas-refpanel-alleles.py ${susie_locus}.ukbb.vcf.gz ${susie_locus}.gwas.tsv ${susie_locus}.gwas.tsv

## cleanup
rm -rf ${susie_locus}.ukbb-header.txt ${susie_locus}.ukbb.vcf.gz*

