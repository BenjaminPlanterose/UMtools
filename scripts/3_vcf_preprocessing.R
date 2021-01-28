library(data.table)

extract.maf <- function(metadata)
{
  partB = unlist(lapply(strsplit(metadata, split= "CAF="), function(x) x[2]))
  CAF = unlist(lapply(strsplit(partB, split= ";"), function(x) x[1]))
  MAF = unlist(lapply(strsplit(CAF, split = ","), function(x) min(suppressWarnings(sort(as.numeric(x), decreasing = T)[1:2]))))
  return(MAF)
}

extract.type <- function(metadata)
{
  partB = unlist(lapply(strsplit(metadata, split= "VC="), function(x) x[2]))
  type = unlist(lapply(strsplit(partB, split= ";"), function(x) x[1]))
  return(type)
}

###################################### CpG ####################################

setwd("/home/ben/Documents/Git/UMtools_dat/dat_450K/vcf/raw/")
setwd("/home/ben/Documents/Git/UMtools_dat/dat_EPIC/vcf/raw/")

#setwd("/media/ben/DATA/Ben/3_genetic_artefacts/annotation/UMtools_dat/dat_450K/vcf/raw/")

CpG_SNP = fread("CpG_sites.vcf"); dim(CpG_SNP) # 61733    12
head(CpG_SNP)

metadata = as.vector(CpG_SNP$V8); MAF = extract.maf(metadata)
CpG_SNP = CpG_SNP[MAF > 0.01,] # Filter > 1 % variants
dim(CpG_SNP) # 17842    12
metadata = as.vector(CpG_SNP$V8); type = extract.type(metadata)
table(type)
# DIV   SNV
# 1118 16724



indel = CpG_SNP[type == "DIV",]; dim(CpG_SNP) # 16724   12
CpG_SNP = CpG_SNP[type == "SNV",]; dim(indel) # 1118  12

setwd("/home/ben/Documents/Git/UMtools_dat/dat_450K/vcf/proc/")
setwd("/home/ben/Documents/Git/UMtools_dat/dat_EPIC/vcf/proc/")
fwrite(CpG_SNP, "CpG_SNP.vcf", quote = F, sep = "\t", col.names = F)
fwrite(indel, "CpG_indel.vcf", quote = F, sep = "\t", col.names = F)

###################################### SBE ####################################

setwd("/home/ben/Documents/Git/UMtools_dat/dat_450K/vcf/raw/")
setwd("/home/ben/Documents/Git/UMtools_dat/dat_EPIC/vcf/raw/")

SBE_SNP = fread("SBE_sites_typeI.vcf"); dim(SBE_SNP) # 2197   12
metadata = as.vector(SBE_SNP$V8); MAF = extract.maf(metadata)
SBE_SNP = SBE_SNP[MAF > 0.01,] # Filter > 1 % variants
dim(SBE_SNP) # 855  12
metadata = as.vector(SBE_SNP$V8); type = extract.type(metadata)
table(type)
# DIV SNV
# 293 562

setwd("/home/ben/Documents/Git/UMtools_dat/dat_450K/vcf/proc/")
setwd("/home/ben/Documents/Git/UMtools_dat/dat_EPIC/vcf/proc/")

fwrite(SBE_SNP[type == "SNV",], "SBE_SNP.vcf", quote = F, sep = "\t", col.names = F)
fwrite(SBE_SNP[type == "DIV",], "SBE_indel.vcf", quote = F, sep = "\t", col.names = F)

###################################### Probes ####################################


setwd("/home/ben/Documents/Git/UMtools_dat/dat_450K/vcf/raw/")
setwd("/home/ben/Documents/Git/UMtools_dat/dat_EPIC/vcf/raw/")

probe_SNP = fread("probe_sites_notCpG.vcf"); dim(probe_SNP) # 319133     12
metadata = as.vector(probe_SNP$V8); MAF = extract.maf(metadata)
probe_SNP = probe_SNP[MAF > 0.01,] # Filter > 1 % variants
dim(probe_SNP) # 114411     12
metadata = as.vector(probe_SNP$V8); type = extract.type(metadata)
table(type)
# DIV    SNV
# 10683 103728

setwd("/home/ben/Documents/Git/UMtools_dat/dat_450K/vcf/proc/")
setwd("/home/ben/Documents/Git/UMtools_dat/dat_EPIC/vcf/proc/")

fwrite(probe_SNP[type == "SNV",], "probe_SNP.vcf", quote = F, sep = "\t", col.names = F)
fwrite(probe_SNP[type == "DIV",], "probe_indel.vcf", quote = F, sep = "\t", col.names = F)

