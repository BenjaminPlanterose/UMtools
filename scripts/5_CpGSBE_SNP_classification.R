library(data.table)

extract.maf <- function(metadata)
{
  partB = unlist(lapply(strsplit(metadata, split= "CAF="), function(x) x[2]))
  CAF = unlist(lapply(strsplit(partB, split= ";"), function(x) x[1]))
  MAF = unlist(lapply(strsplit(CAF, split = ","), function(x) min(suppressWarnings(sort(as.numeric(x), decreasing = T)[1:2]))))
  return(MAF)
}

extract.caf <- function(metadata)
{
  partB = unlist(lapply(strsplit(metadata, split= "CAF="), function(x) x[2]))
  CAF = unlist(lapply(strsplit(partB, split= ";"), function(x) x[1]))
  return(CAF)
}

prep_alternate_allele <- function(CpG_SNP, Thr)
{
  # Extract metadata from sites with multiple minor alleles
  minor = CpG_SNP$m
  l_m = sapply(strsplit(minor, split = ","), length)
  where = which(l_m >= 2)
  metadata = as.vector(CpG_SNP$metadata)[where]
  CAF = extract.caf(metadata)

  # Eliminate alleles if below a threshold
  eliminate = lapply(strsplit(CAF, split =  ","), function(x) which(x == "." | as.numeric(x) < Thr)-1)
  length0 = which(sapply(eliminate, length) == 0)
  splitted = strsplit(minor[where], split =  ",")
  minor_cor = minor
  minor_cor[where] = unlist(lapply(1:length(splitted), function(x) paste((splitted[[x]])[-eliminate[[x]]], collapse = ",")))
  minor_cor[where][length0] = minor[where][length0]
  CpG_SNP$m = minor_cor
  return(CpG_SNP)
}

#### 450K
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
data(Manifest); head(Manifest)
data(Locations); head(Locations)
annotation <- cbind(Locations, Manifest); rm(Manifest, Locations); gc()
sum(!(rownames(annotation) == annotation$Name)) # 0
dim(annotation) # 485512     11
setwd("/home/ben/Documents/Git/UMtools_dat/dat_450K/vcf/proc/")
####


#### EPIC
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
data(Manifest); head(Manifest)
data(Locations); head(Locations)
annotation <- cbind(Locations, Manifest); rm(Manifest, Locations); gc()
sum(!(rownames(annotation) == annotation$Name)) # 0
dim(annotation) # 865859     11
setwd("/home/ben/Documents/Git/UMtools_dat/dat_EPIC/vcf/proc/")
####

# CpG-SNP
CpG_SNP = fread("CpG_SNP.vcf")
CpG_SNP = as.data.frame(CpG_SNP)
colnames(CpG_SNP) = c("chr", "start", "rs", "M", "m", "v6", "v7", "metadata", "chr2", "start2", "end2", "cg")
CpG_SNP$type = annotation[CpG_SNP$cg,"Type"]
CpG_SNP$strand = annotation[CpG_SNP$cg,"strand"]
CpG_SNP$real_pos = annotation[CpG_SNP$cg, "pos"]
metadata = as.vector(CpG_SNP$metadata); MAF = extract.maf(metadata)
CpG_SNP = CpG_SNP[MAF > 0.05,] # Filter > 5 % variants

# SBE-SNP
SBE_SNP = fread("SBE_SNP.vcf")
SBE_SNP = as.data.frame(SBE_SNP)
colnames(SBE_SNP) = c("chr", "start", "rs", "M", "m", "v6", "v7", "metadata", "chr2", "start2", "end2", "cg")
SBE_SNP$strand = annotation[SBE_SNP$cg,"strand"]
SBE_SNP$real_pos = annotation[SBE_SNP$cg, "pos"]
SBE_SNP$type = annotation[SBE_SNP$cg,"Type"]
metadata = as.vector(SBE_SNP$metadata); MAF = extract.maf(metadata)
SBE_SNP = SBE_SNP[MAF > 0.05,] # Filter > 5 % variants


###################################### Biallelic Vs triallelic CpG/SNP ####################################

CpG_SNP <- prep_alternate_allele(CpG_SNP, 0.05) # Eliminates alternate alleles < Thr
where = sapply(strsplit(CpG_SNP$m, split = ","), length) > 1
CpG_SNP_Tri = CpG_SNP[where, ]
CpG_SNP_Bi = CpG_SNP[!where, ]
rownames(CpG_SNP_Bi) = rownames(CpG_SNP_Tri) = NULL


table(CpG_SNP_Tri$M, CpG_SNP_Tri$m) # 450K
# A,C A,T C,T G,T
# C   0   4   0   4
# G   6   4   2   0

table(CpG_SNP_Bi$M, CpG_SNP_Bi$m) # 450K
#        A    C    G    T
# C    316    0  329 2570
# G   3366  437    0  407


setwd("/home/ben/Documents/Git/UMtools/data/")
save(CpG_SNP_Tri, file = "triallelic_CpG_SNP_450K.RData")
save(CpG_SNP_Tri, file = "triallelic_CpG_SNP_EPIC.RData")


###################################### SBE ####################################

SBE_SNP <- prep_alternate_allele(SBE_SNP, 0.05) # Eliminates alternate alleles < Thr
where = sapply(strsplit(SBE_SNP$m, split = ","), length) > 1; sum(where) # 0

table(SBE_SNP$M, SBE_SNP$m)
#    A  C  G  T
# A  0  8 35  9
# C 16  0 23 57
# G 56 28  0 20
# T  8 28  7  0


###################################### CpG/SBE-SNPs ####################################

# II
#   (+)
#       SNP1
#           C -> T/A          SNP is interpreted as unmethylated allele
#           C -> G            SNP is interpreted as methylated allele
#       SNP2
#           G -> A/C/T        3'-overhang: probe failure
#   (-)
#       SNP1
#           C -> A/G/T        3'-overhang: probe failure
#       SNP2
#           G -> T/A          SNP is interpreted as unmethylated allele
#           G -> C            SNP is interpreted as methylated allele

# I
#   (+)
#       SNP0
#           T <-> C <-> A     U/M Informative detection in the same channel
#           A/T/C <-> G       U/M Signal goes to the wrong channel (oob)
#       SNP1
#           C -> T            SNP is interpreted as unmethylated allele
#           C -> A/G          3’-overhang: probe failure for either U or M
#       SNP2
#           G -> A/C/T        3'-overhang: probe failure for either U or M
#   (-)
#       SNP1
#           C -> A/G/T        3'-overhang: probe failure for either U or M
#       SNP2
#           G -> A            SNP is interpreted as unmethylated allele
#           G -> C/T          3'-overhang: probe failure for either U or M
#       SNP3
#           T <-> C <-> A     U/M Informative detection in the same channel
#           A/T/C <-> G       U/M Signal goes to the wrong channel (oob)
#


names = c("II_plus_SNP1_C_TA", "II_plus_SNP1_C_G", "II_plus_SNP2_G_ACT",
          "II_minus_SNP1_C_AGT", "II_minus_SNP2_G_TA", "II_minus_SNP2_G_C",
          "I_plus_SNP0_A_C_T", "I_plus_SNP0_ACT_G", "I_plus_SNP1_C_T", "I_plus_SNP1_C_AG", "I_plus_SNP2_G_ACT",
          "I_minus_SNP1_C_AGT", "I_minus_SNP2_G_A", "I_minus_SNP2_G_CT", "I_minus_SNP3_A_G_T", "I_minus_SNP3_AGT_C")

classification = as.list(rep(NA, length(names)))
names(classification) = names

classification[["II_plus_SNP1_C_TA"]] = CpG_SNP_Bi[CpG_SNP_Bi$type == "II" &
                                  CpG_SNP_Bi$strand == "+" &
                                  CpG_SNP_Bi$M == "C" &
                                  CpG_SNP_Bi$m %in% c("T", "A"),]

classification[["II_plus_SNP1_C_G"]] = CpG_SNP_Bi[CpG_SNP_Bi$type == "II" &
                                    CpG_SNP_Bi$strand == "+" &
                                    CpG_SNP_Bi$M == "C" &
                                    CpG_SNP_Bi$m %in% c("G"),]

classification[["II_plus_SNP2_G_ACT"]] = CpG_SNP_Bi[CpG_SNP_Bi$type == "II" &
                                                    CpG_SNP_Bi$strand == "+" &
                                                    CpG_SNP_Bi$M == "G" &
                                                    CpG_SNP_Bi$m %in% c("A", "C", "T"),]

classification[["II_minus_SNP1_C_AGT"]] = CpG_SNP_Bi[CpG_SNP_Bi$type == "II" &
                                                      CpG_SNP_Bi$strand == "-" &
                                                      CpG_SNP_Bi$M == "C" &
                                                      CpG_SNP_Bi$m %in% c("A", "G", "T"),]

classification[["II_minus_SNP2_G_TA"]] = CpG_SNP_Bi[CpG_SNP_Bi$type == "II" &
                                                       CpG_SNP_Bi$strand == "-" &
                                                       CpG_SNP_Bi$M == "G" &
                                                       CpG_SNP_Bi$m %in% c("A", "T"),]

classification[["II_minus_SNP2_G_C"]] = CpG_SNP_Bi[CpG_SNP_Bi$type == "II" &
                                                      CpG_SNP_Bi$strand == "-" &
                                                      CpG_SNP_Bi$M == "G" &
                                                      CpG_SNP_Bi$m %in% c("C"),]

classification[["I_plus_SNP1_C_T"]] = CpG_SNP_Bi[CpG_SNP_Bi$type == "I" &
                                                     CpG_SNP_Bi$strand == "+" &
                                                     CpG_SNP_Bi$M == "C" &
                                                     CpG_SNP_Bi$m %in% c("T"),]

classification[["I_plus_SNP1_C_AG"]] = CpG_SNP_Bi[CpG_SNP_Bi$type == "I" &
                                                   CpG_SNP_Bi$strand == "+" &
                                                   CpG_SNP_Bi$M == "C" &
                                                   CpG_SNP_Bi$m %in% c("A", "G"),]

classification[["I_plus_SNP2_G_ACT"]] = CpG_SNP_Bi[CpG_SNP_Bi$type == "I" &
                                                    CpG_SNP_Bi$strand == "+" &
                                                    CpG_SNP_Bi$M == "G" &
                                                    CpG_SNP_Bi$m %in% c("A", "C", "T"),]

classification[["I_minus_SNP1_C_AGT"]] = CpG_SNP_Bi[CpG_SNP_Bi$type == "I" &
                                                     CpG_SNP_Bi$strand == "-" &
                                                     CpG_SNP_Bi$M == "C" &
                                                     CpG_SNP_Bi$m %in% c("A", "G", "T"),]

classification[["I_minus_SNP2_G_A"]] = CpG_SNP_Bi[CpG_SNP_Bi$type == "I" &
                                                      CpG_SNP_Bi$strand == "-" &
                                                      CpG_SNP_Bi$M == "G" &
                                                      CpG_SNP_Bi$m %in% c("A"),]

classification[["I_minus_SNP2_G_CT"]] = CpG_SNP_Bi[CpG_SNP_Bi$type == "I" &
                                                    CpG_SNP_Bi$strand == "-" &
                                                    CpG_SNP_Bi$M == "G" &
                                                    CpG_SNP_Bi$m %in% c("C", "T"),]


classification[["I_plus_SNP0_A_C_T"]] = SBE_SNP[SBE_SNP$type == "I" &
                                                SBE_SNP$strand == "+" &
                                                SBE_SNP$M %in% c("C", "T", "A") & SBE_SNP$m %in% c("C", "T", "A"),]

classification[["I_plus_SNP0_ACT_G"]] = SBE_SNP[SBE_SNP$type == "I" &
                                                  SBE_SNP$strand == "+" &
                                                  (SBE_SNP$M %in% c("C", "T", "A", "G") & SBE_SNP$m %in% c("G") |
                                                  SBE_SNP$M %in% c("G") & SBE_SNP$m %in% c("C", "T", "A", "G")),]

classification[["I_minus_SNP3_A_G_T"]] = SBE_SNP[SBE_SNP$type == "I" &
                                                  SBE_SNP$strand == "-" &
                                                  SBE_SNP$M %in% c("G", "T", "A") & SBE_SNP$m %in% c("G", "T", "A"),]

classification[["I_minus_SNP3_AGT_C"]] = SBE_SNP[SBE_SNP$type == "I" &
                                                  SBE_SNP$strand == "-" &
                                                  (SBE_SNP$M %in% c("G", "T", "A") & SBE_SNP$m %in% c("C") |
                                                  SBE_SNP$M %in% c("C") & SBE_SNP$m %in% c("G", "T", "A")),]


sapply(classification, nrow) # 450K
# II_plus_SNP1_C_TA    II_plus_SNP1_C_G  II_plus_SNP2_G_ACT II_minus_SNP1_C_AGT  II_minus_SNP2_G_TA
# 1131                 111                1266                1221                1916
# II_minus_SNP2_G_C   I_plus_SNP0_A_C_T   I_plus_SNP0_ACT_G     I_plus_SNP1_C_T    I_plus_SNP1_C_AG
# 188                  49                  97                 283                  98
# I_plus_SNP2_G_ACT  I_minus_SNP1_C_AGT    I_minus_SNP2_G_A   I_minus_SNP2_G_CT  I_minus_SNP3_A_G_T
# 364                 371                 355                 121                  62
# I_minus_SNP3_AGT_C
# 89

sum(sapply(classification, nrow)) # 7722; 450K
nrow(SBE_SNP) + nrow(CpG_SNP_Bi) # 7722; 450K


setwd("/home/ben/Documents/Git/UMtools/data/")
save(classification, file = "classification_CpG_SNP_450K.RData")
save(classification, file = "classification_CpG_SNP_EPIC.RData")
