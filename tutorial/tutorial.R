#### install-uninstall
setwd("/home/ben/Documents/Git/")
library(devtools); library(roxygen2); install("UMtools")
remove.packages("UMtools")
####


#### test datasets
data(annot_450K)
data(annot_EPIC)
data(classification_CpG_SNP_450K)
data(classification_CpG_SNP_EPIC)
data(CR_probes)
data(triallelic_CpG_SNP_450K)
data(triallelic_CpG_SNP_EPIC)
data(training_set)
####

#### check documentation
help(annot_450K)
help(annot_EPIC)
help(classification_CpG_SNP_450K)
help(classification_CpG_SNP_EPIC)
help(CR_probes)
help(triallelic_CpG_SNP_450K)
help(triallelic_CpG_SNP_EPIC)
help(training_set)
help(bGMM)
help(compute_BC_CV)
help(compute_CV)
help(GR_to_UM)
help(density_jitter_plot)
help(export_bigmat)
help(import_bigmat)
help(Kcall_CpG)
help(par_EW_Kcalling)
help(train_k_caller)
help(UM_plot)
help(Visualize_cometh)
####

#### Update documentation
library(roxygen2); library(devtools)
setwd("/home/ben/Documents/Git/UMtools/")
document()
####


########################################## 1. Installation ##########################################

# To install R-packages from github, you need devtools:
library(devtools)

# Dependencies from Bioconductor
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
BiocManager::install('minfi')
BiocManager::install('GEOquery')
# Dependencies from CRAN
devtools::install_version("modes", "0.7.0")
install.packages("scales")
install.packages("EMCluster")
install.packages("dbscan")
install.packages("RColorBrewer")
# Dependencies from Github
install_github("dphansti/Sushi")
install_github("BenjaminPlanterose/UMtools")

########################################## 2. Downloading example ##########################################

# Run in bash
# wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE104nnn/GSE104812/suppl/GSE104812_RAW.tar --wait=10 --limit-rate=50K
# tar -xvf GSE104812_RAW.tar
# find . -type f ! -name '*.idat.gz' -delete
# gunzip *.gz

########################################## 3. Peaking into an IDAT ##########################################

# Illuminaio is a dependency of minfi
library(illuminaio)
setwd("/media/ben/DATA/Ben/1_evCpGs/data/aging_children/GSE104812_RAW/") # Change this route to fit your system
example = illuminaio::readIDAT("GSM2808239_sample1_Grn.idat")
names(example)
# [1] "fileSize"  "versionNumber" "nFields"   "fields"  "nSNPsRead" "Quants" "MidBlock"
# [8] "RedGreen"  "Barcode"       "ChipType"  "RunInfo" "Unknowns"

head(example$RunInfo, 3)
#             RunTime          BlockType
# [1,] "1/12/2016 2:37:05 AM" "Decoding"
# [2,] "4/8/2016 6:03:57 PM"  "Scan"
# [3,] "4/8/2016 6:03:57 PM"  "Register"

head(example$Quants, 3)
#           Mean  SD  NBeads
# 10600313  284  137     13
# 10600322 9405 1363     14
# 10600328 3538  439     11

########################################## 4. Extracting fluorescence intensity matrices ##########################################

library(UMtools)
setwd("/media/ben/DATA/Ben/1_evCpGs/data/aging_children/GSE104812_RAW/") # Change this route to fit your system
rgSet = read.metharray.exp(getwd(), extended = TRUE)
TypeI.Red <- getProbeInfo(rgSet, type = "I-Red")
TypeI.Green <- getProbeInfo(rgSet, type = "I-Green")
TypeII <- getProbeInfo(rgSet, type = "II")
ctrls <- getProbeInfo(rgSet, type = "Control")
SnpI <- getProbeInfo(rgSet, type = "SnpI")
SnpII <- getProbeInfo(rgSet, type = "SnpII")
known_probes = c(SnpI$AddressA, SnpI$AddressB, SnpII$AddressA, ctrls$Address, TypeI.Red$AddressA,
                 TypeI.Red$AddressB, TypeI.Green$AddressA, TypeI.Green$AddressB, TypeII$AddressA)
length(known_probes) # 621926
all = rownames(rgSet); length(known_probes) # 622399
orphan = all[!(all %in% known_probes)]; length(orphan) # 473

head(TypeI.Green[, 1:3], n = 2)
#          Name    AddressA    AddressB
#   <character> <character> <character>
# 1  cg02004872    25785404    58629399
# 2  cg02050847    43656343    73683470

head(TypeII[, 1:2], n = 2)
#          Name    AddressA
#   <character> <character>
# 1  cg00035864    31729416
# 2  cg00061679    28780415

# Getting the probe annotation
library(IlluminaHumanMethylation450kanno.ilmn12.hg19) # 450K
annotation <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
# library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19) # 850K
# annotation <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

colnames(annotation)
# [1]  "chr"                      "pos"                      "strand"                   "Name"
# [5]  "AddressA"                 "AddressB"                 "ProbeSeqA"                "ProbeSeqB"
# [9]  "Type"                     "NextBase"                 "Color"                    "Probe_rs"
# [13] "Probe_maf"                "CpG_rs"                   "CpG_maf"                  "SBE_rs"
# [17] "SBE_maf"                  "Islands_Name"             "Relation_to_Island"       "Forward_Sequence"
# [21] "SourceSeq"                "Random_Loci"              "Methyl27_Loci"            "UCSC_RefGene_Name"
# [25] "UCSC_RefGene_Accession"   "UCSC_RefGene_Group"       "Phantom"                  "DMR"
# [29] "Enhancer"                 "HMM_Island"               "Regulatory_Feature_Name"  "Regulatory_Feature_Group"
# [33] "DHS"

Grn = assay(rgSet, "Green")       # Green mean across beads
Red = assay(rgSet, "Red")         # Red mean across beads
GrnSD = assay(rgSet, "GreenSD")   # Green SD across beads
RedSD = assay(rgSet, "RedSD")     # Red SD across beads
nBeads = assay(rgSet, "NBeads")   # Number of Beads across probes

M_U = GR_to_UM(Red = Red, Grn = Grn, rgSet = rgSet, what = "Mean")
M_U_sd = GR_to_UM(Red = RedSD, Grn = GrnSD, rgSet = rgSet, what = "SD")
nBeads_cg = GR_to_UM(nBeads = nBeads, rgSet = rgSet, what = "NBeads")

# PEAK
M_U$M[1:3, 1:3]; M_U$U[1:3, 1:3]
M_U_sd$M[1:3, 1:3]; M_U_sd$U[1:3, 1:3]
nBeads_cg[1:3, 1:3]


########################################## 5. Compute beta- and M-values ##########################################

offset = 100 # For numerical stability at low fluorescence intensities
beta_value = M_U$M/(M_U$M + M_U$U + offset)

offset = 1 # For numerical stability at low fluorescence intensities
M_value = log2((M_U$M + offset)/(M_U$U + offset))

rm(Grn, Red, GrnSD, RedSD, nBeads, nBeads_cg, M_value); gc()

########################################## 6. Quickly importing/exporting large matrices with data.table ##########################################

setwd("/media/ben/DATA/Ben/3_genetic_artefacts/R-packages/test/") # Change this route to fit your system
export_bigmat(M_U$M, "M.txt", nThread = 4)
list.files()
M = import_bigmat("2021-01-28_M.txt", nThread = 4)

########################################## 7. Quickly importing GEO phenotypes with GEOquery ##########################################

library(GEOquery)
setwd('/media/ben/DATA/Ben/1_evCpGs/data/aging_children/GSE104812_RAW/')
pheno_object <- getGEO('GSE104812', destdir=".", getGPL = FALSE)
pheno <- pheno_object[[1]]
pheno <- phenoData(pheno)
pheno <- pData(pheno)
pheno = data.frame(GEO_ID = as.character(rownames(pheno)),
                   sex = as.factor(pheno$`gender:ch1`),
                   age = as.numeric(pheno$`age (y):ch1`))
IDAT_IDs = sapply(strsplit(colnames(rgSet), split = "_"),function(x) x[1])
pheno <- pheno[match(pheno$GEO_ID, IDAT_IDs),] # Make sure samples in pheno are in the same order as in IDATs

########################################## 8. UMtools in action ##########################################

# 450K annotation
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
annotation <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# Y-chromosome targeting probe
density_jitter_plot(beta_value, "cg00050873", pheno$sex)
annotation["cg00050873", c("chr", "pos")] # chrY   9363356
UM_plot(M = M_U$M, U = M_U$U, CpG = "cg00050873", sex = pheno$sex)

# X-inactivation
UM_plot(M = M_U$M, U = M_U$U, CpG = "cg00026186", sex = pheno$sex)
annotation["cg00026186", c("chr", "pos")] # chrX  48367230

# X-inactivation escape
UM_plot(M = M_U$M, U = M_U$U, CpG = "cg04927982", sex = pheno$sex)
annotation["cg04927982", c("chr", "pos")] # chrX  53254653

# X-hypermethylation
UM_plot(M = M_U$M, U = M_U$U, CpG = "cg02973417", sex = pheno$sex)
annotation["cg02973417", c("chr", "pos")] # chrX  153220895


### Bivariate Gaussian Mixture Models (bGMMs)
set.seed(2); res = bGMM(M_U$M, M_U$U, "cg03398919", K = 2)
set.seed(3); res = bGMM(M_U$M, M_U$U, "cg00814218", K = 3)
set.seed(2); res = bGMM(M_U$M, M_U$U, "cg27024127", K = 4)
set.seed(4); res = bGMM(M_U$M, M_U$U, "cg23186955", K = 5)


### Quantifying epigenome-wide ambivalency in probe failure
CV = compute_CV(M_SD = M_U_sd$M, U_SD = M_U_sd$U, M = M_U$M, U = M_U_sd$U, alpha = 100)
BC_CV = compute_BC_CV(CV = CV)

density_jitter_plot(CV, "cg00050873", pheno$sex)
BC_CV["cg00050873"]
# cg00050873
#   1.128555
annotation["cg00050873", c("chr", "pos")] # chrY   9363356


### K-calling
Kcall_CpG("cg15771735", M_U$M, M_U$U, minPts = 5, eps = 0.1)
Kcall_CpG("cg03398919", M_U$M, M_U$U, minPts = 5, eps = 0.1)
Kcall_CpG("cg00814218", M_U$M, M_U$U, minPts = 5, eps = 0.1)
Kcall_CpG("cg27024127", M_U$M, M_U$U, minPts = 5, eps = 0.1)

chrY = rownames(annotation)[annotation$chr == "chrY"]



start_time <- Sys.time()
K_vec = par_EW_Kcalling(M_U$M[chrY,], M_U$U[chrY,], minPts = 5, eps = 0.1, nThread = 10)
end_time <- Sys.time()
end_time - start_time # Consumes 6 GB, 17.46199 secs

table(K_vec)
# K_vec
#  1   2
# 35 381

# K = 1 are cross-reactive
UM_plot(M = M_U$M, U = M_U$U, CpG = "cg02494853", sex = pheno$sex)
annotation["cg02494853", c("chr", "pos")] # chrY   4868397


# If aiming to employ epigenome-wide, it is very important to adjust {minPts, eps} to the
# sample size of the data employed.

# setwd("/home/ben/Documents/Git/UMtools/data/")
# load("training_set.Rdata")
data("training_set")

# Training annotated with dataset of 426 EUR MZ twin pairs. Training set may not be correctly annotated
# for other ancestries and other sample sizes.
Kcall_CpG(sample(training_set$k_1, 1), M_U$M, M_U$U, minPts = 5, eps = 0.1)
Kcall_CpG(sample(training_set$k_2, 1), M_U$M, M_U$U, minPts = 5, eps = 0.1)
Kcall_CpG(sample(training_set$k_3, 1), M_U$M, M_U$U, minPts = 5, eps = 0.1)
Kcall_CpG(sample(training_set$k_4, 1), M_U$M, M_U$U, minPts = 5, eps = 0.1)

# Lower maf variants are wrongly annotated in this dataset. Sample size not big enough.
train_k_caller(M_U$M, M_U$U, training_set, 3, 0.07, nThread = 10) # 0.7948261

### Comethylation plots

# meQTL
UM_plot(M = M_U$M, U = M_U$U, CpG = "cg14911689", sex = NULL)
annotation["cg14911689", c("chr", "pos", "UCSC_RefGene_Name")] # chr12    739980    NINJ2
annotation <- annotation[order(annotation$chr, annotation$pos),]
pos <- which(rownames(annotation) == "cg14911689")
UM_plot(M = M_U$M, U = M_U$U, CpG = rownames(annotation)[pos-1], sex = NULL)
UM_plot(M = M_U$M, U = M_U$U, CpG = rownames(annotation)[pos+1], sex = NULL)

res = Visualize_cometh(annotation = annotation, CpG = 'cg14911689', distance = 1000,
                       L_bound = 3, R_bound = 2, beta_mat = beta_value, max_y = 5)


# A genetic artefact
UM_plot(M = M_U$M, U = M_U$U, CpG = "cg11495604", sex = NULL)
res = Visualize_cometh(annotation = annotation, CpG = 'cg11495604', distance = 1000,
                       L_bound = 0, R_bound = 2, beta_mat = beta_value, max_y = 5)
annotation["cg11495604", c("chr", "pos", "UCSC_RefGene_Name")] # chr20  62053198 KCNQ2;KCNQ2;KCNQ2;KCNQ2


# HLA locus
res = Visualize_cometh(annotation = annotation, CpG = 'cg00211215', distance = 200,
                       L_bound = 3, R_bound = 0, beta_mat = beta_value, max_y = 5)
annotation["cg00211215", c("chr", "pos", "UCSC_RefGene_Name")] # chr6  32552246   HLA-DRB1
