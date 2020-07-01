#!/bin/bash

# Download a vcf file containing variant information. In this case, common variants from dbSNP151. Make sure that the coordinates are for the 
# genome same assembly as the coordinates reported in the microarray annotation (This case hg19, aka GRCh37).
wget https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/00-common_all.vcf.gz
gunzip 00-common_all.vcf.gz

bedtools intersect -wa -wb -a /media/ben/DATA/Ben/3_genetic_artefacts/annotation/dbSNP_v151/00-common_all.vcf -b CpG_sites.bed > CpG_sites.vcf
bedtools intersect -wa -wb -a /media/ben/DATA/Ben/3_genetic_artefacts/annotation/dbSNP_v151/00-common_all.vcf -b SBE_sites_typeI.bed > SBE_sites_typeI.vcf
bedtools intersect -wa -wb -a /media/ben/DATA/Ben/3_genetic_artefacts/annotation/dbSNP_v151/00-common_all.vcf -b probe_sites_notCpG.bed > probe_sites_notCpG.vcf