#!/usr/bin/env Rscript

library("argparse")
library("dplyr")
library("tidyr")

parser <- ArgumentParser(description= 'Translate GWAS summary stats from hg19 to hg38.\n
                                       Some part of this script is hard-coding col names, so check carefully.')
parser$add_argument('--ukbb_hg38', '-u', help= 'VCF file on hg38 containing SNPs. Example: "/scratch/scjp_root/scjp99/vthihong/2_PanKBase/colocGWAS_T1D/2_t1d-eQTL-coloc/results/hg38/ukbb-rsid/ukbb-rsid.vcf.gz')
parser$add_argument('--summStat', '-t', help= 'Text file containing GWAS signals and summary stats in hg19. Example: "/scratch/scjp_root/scjp99/vthihong/2_PanKBase/colocGWAS_T1D/2_t1d-eQTL-coloc/results/hg38/temp-ukbb/chr{chrom}_{analysis}.hg19.bed"')
parser$add_argument('--bed_hg38', '-b', help= 'Bed file containing GWAS signal locations in hg38, obtained using LiftOver. Example: /scratch/scjp_root/scjp99/vthihong/2_PanKBase/colocGWAS_T1D/2_t1d-eQTL-coloc/results/hg38/temp-ukbb/new.chr{chrom}_{analysis}.hg38')
parser$add_argument('--output', '-o', help= 'Output file.')

xargs <- parser$parse_args()

library(data.table)
library(dplyr)

complements <- list("A" = "T",
                    "T" = "A",
                    "G" = "C",
                    "C" = "G")

revcomp <- function(x) {
  out <- sapply(strsplit(x, "")[[1]], function(i) complements[[i]])
  return(paste(out, collapse = ""))
}

ukbb_hg38 <- fread(xargs$ukbb_hg38, skip = "CHROM")

summStat <- fread(xargs$summStat, header = F)
tmp <- which(colSums(is.na(summStat))<nrow(summStat))
summStat <- as.data.frame(summStat)
summStat <- summStat[, names(tmp)]
print(head(summStat))

if (length(grep("Exons", xargs$summStat)) > 0) {
	colnames(summStat) <- c("GeneName", "Strand", "GencodeLevel", "GeneType", "GeneID", "ChrPheno", "StartPheno", "EndPheno", "BestExonID", "NumExons", "NumVariantCis", "DistanceWithBest", "SNPid", "A1", "A2", "MAF", "SNPchr", "StartSNP", "EndSNP", "Nominal_Pval", "Slope", "EmpiricalAdjustedPval", "BetaAdjustedPval", "eQTLnum", "OrderByDistanceTSS", "OrderBySlope", "PvalueOrder", "Probability")
} else {
	colnames(summStat) <- c("GeneName", "Strand", "GencodeLevel", "GeneType", "GeneID", "ChrPheno", "StartPheno", "EndPheno", "NumSNPs", "DistanceWithBest", "SNPid", "A1", "A2", "MAF", "SNPchr", "StartSNP", "EndSNP", "Nominal_Pval", "Slope", "EmpiricalAdjustedPval", "BetaAdjustedPval", "DiscoveryOrder")
}

#colnames Exons: GeneName	Strand	GencodeLevel	GeneType	GeneID	ChrPheno	StartPheno	EndPheno	BestExonID	NumExons	NumVariantCis	DistanceWithBest	SNPid	A1	A2	MAF	SNPchr	StartSNP	EndSNP	Nominal_Pval	Slope	EmpiricalAdjustedPval	BetaAdjustedPval	eQTLnum	OrderByDistanceTSS	OrderBySlope	PvalueOrder	Probability
#colnames Genes: GeneName	Strand	GencodeLevel	GeneType	GeneID	ChrPheno	StartPheno	EndPheno	NumSNPs	DistanceWithBest	SNPid	A1	A2	MAF	chrSNP	StartSNP	EndSNP	Nominal_Pval	Slope	EmpiricalAdjustedPval	BetaAdjustedPval	DiscoveryOrder

summStat$SNPchr <- as.character(summStat$SNPchr)
summStat$SNPchr <- paste0("chr", as.character(summStat$SNPchr))

bed_hg38 <- read.table(xargs$bed_hg38, header = F)
colnames(bed_hg38) <- c("chr", "start", "end", "id")
bed_hg38 <- separate(bed_hg38, "id", c("snp", "GeneID"), sep = "_E")
bed_hg38$GeneID <- paste0("E", bed_hg38$GeneID)
bed_hg38 <- distinct(bed_hg38)
test <- inner_join(as.data.frame(ukbb_hg38[,1:5]), bed_hg38, by = c("#CHROM" = "chr", "POS" = "end")) #1:5 here means we will exclude QUAL, FILTER and INFO columns of ukbb file
test2 <- inner_join(test, summStat, by = c("#CHROM" = "SNPchr", "snp" = "SNPid", "GeneID" = "GeneID"))
#needed info - something like the following, depending on what GWAS file gave:
##snp_chrom	snp_start	snp_end	SNP	REF	ALT	minor_allele	maf	n_complete_samples	AC_hg19	ytx_hg19
##beta	se	tstat	pval	SNP_hg19	multiply

### check matching between ref/alt alleles in different genome built
# REF/ATL columns are ref/alt in reference genome, in this case hg38
# effect_allele/other_allele/effect_allele_frequency are in the original summ stat file, in this case hg19
#revEA <- unlist(lapply(test2$effect_allele, revcomp)) #obtain reverse compliments of the alleles
#revOA <- unlist(lapply(test2$other_allele, revcomp))

#if alt allele = reverse compliment of effect alleles & alt allele = reverse compliment of effect alleles, sign(beta) remains the same
test2$multiply <- ifelse(test2$REF == test2$A1 & test2$ALT == test2$A2, -1,
                         ifelse(test2$ALT == test2$A1 & test2$REF == test2$A2, 1, NA))
                                #ifelse(test2$ALT == revEA & test2$REF == revOA, 1, NA))) 
test2 <- test2[!is.na(test2$multiply),] #on lines where none of the cases in ifelse condition holds, drop these lines
test2$Slope <- test2$Slope * test2$multiply
if (length(grep("Exons", xargs$summStat)) > 0) {
  test2 <- test2[, c("#CHROM", "start", "POS", "ID", "REF", "ALT", "A1", "A2",
                     "GeneID", "GeneName", "BestExonID", "Nominal_Pval", "Slope")]
  colnames(test2) <- c("snp_chrom", "snp_start", "snp_end", "SNP",  "REF",  "ALT",  "effect_allele", "other_allele",
                       "GeneID", "GeneName", "BestExonID", "Nominal_Pval", "Slope")
} else {
  test2 <- test2[, c("#CHROM", "start", "POS", "ID", "REF", "ALT", "A1", "A2",
                     "GeneID", "GeneName", "Nominal_Pval", "Slope")]
  colnames(test2) <- c("snp_chrom", "snp_start", "snp_end", "SNP",  "REF",  "ALT",  "effect_allele", "other_allele",
                       "GeneID", "GeneName", "Nominal_Pval", "Slope")
}

write.table(test2, xargs$output, col.names = F, row.names = F, quote = F, sep = "\t")
