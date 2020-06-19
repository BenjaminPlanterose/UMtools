# UMtools
## An R-package for analysing Illumina DNA Methylation microarrays at the fluorescence intensity level



#### Benjamin Planterose Jiménez, Manfred Kayser, Athina Vidaki

### Department of Genetic Identification, Erasmus MC University Medical Centre Rotterdam, The Netherlands

## License
[MIT](https://choosealicense.com/licenses/mit/)


## Tested on

    Ubuntu 18.04.4 LTS (bionic), R version 3.6.3 (2020-02-29) -- "Holding the Windsock"
    Contact b.planterosejimenez@erasmusmc.nl for any issues arising while running UMtools.
    
## Installation of Dependencies 

```r
# From Bioconductor
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install('minfi')
BiocManager::install('GEOquery')

# From CRAN
install.packages(modes)
install.packages(scales)
install.packages(EMCluster)
install.packages(dbscan)

# From Github
library("devtools")
install_github("dphansti/Sushi")
```

## Installation of UMtools:
```r
library(devtools)
devtools::install_github("BenjaminPlanterose/UMtools")
```

## About this tutorial
    
Do not attempt to perform this tutorial without at least 8GB of RAM. Working with fluorescence intensities
involves several large matrices

## 0) A word on the Beadchip microarray technology and the IDAT format: 

The microarray itself consists of a silica substrate with uniformly interspaced microwells.
Hundreds of thousands of copies of a specific oligonucleotide lie on the surface of silica beads.
During the manufacture of the chip, all 622,399 type of beads are pooled together and applied on
the microarray. Subsequently, beads automatically self-assemble on the microarray's microwell. As a result, both
the order and the number copies for a given bead type are random, hence requiring the decoding of the microarray.

Although targeting 485,512 cytosines, the 450K technology employs 622,399 probe oligonucleotides. There are
several reasons for this. To begin with, 450K combines three types of probes concerning detection:

* Type I (n = 135,476 x 2) - two bead types per cytosine 

  * Type-I Green (n = 46,289 x 2) - Quantification is informative only in the Green channel
  
  * Type-I Red (n = 89,187 x 2) - Quantification is informative only in the Red channel
  
* Type-II (n = 350,036): one bead type per cytosine, quantification is informative in both channels

In addition, there are a set of quality control probes (n = 848):
Staining (n = 4), extension (n = 4), hybridization (n = 3), target removal (n = 2), bisulfite conversion I (n = 12),
bisulfite conversion II (n = 4), specificity I (n = 12), specificity II (n = 3), non-polymorphic (n = 4),
negative (n = 613), restoration (n = 1) and normalization (n = 186).

Also, there are sample mix-up probes (n = 65):
SnpI (n = 25 x 2) and SnpII (n = 40)
Finally, there 473 orphan probes with placed on the array for unknown purposes

    473 + 25*2 + 40 + 848 + 46289 * 2 + 89187 * 2 + 350036 = 622399



## 1) Introduction to the IDAT format

The .IDAT extension (Intensity Data) is Illumina's proprietary format for storage of the fluorescence scanners'
raw output across several genome-wide platform. The IDAT format is an encrypted and non-human readable.
It was not until the birth the R-package illuminaio, that IDAT files could only be read via vendor's software.
We first need to download some example IDAT files from the GEO database. You may perform: 


```bash
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE104nnn/GSE104812/suppl/GSE104812_RAW.tar --wait=10 --limit-rate=50K
tar -xvf GSE104812_RAW.tar
find . -type f ! -name '*.idat.gz' -delete
gunzip *.gz
```

Or manually download at https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE104812

Back to R, we firstly load libraries

```r
library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)
library(GEOquery)
library(scales)
library(modes)
library(EMCluster)
library(dbscan)
library(Sushi)
```

To peak into an IDAT file, in this case of the green channel of random sample, we can ran:

```r
setwd("/media/ben/DATA/Ben/1_evCpGs/data/aging_children/GSE104812_RAW/")
example = illuminaio::readIDAT("GSM2808239_sample1_Grn.idat")
names(example)
# [1] "fileSize"      "versionNumber" "nFields"       "fields"        "nSNPsRead"     "Quants"        "MidBlock"
# [8] "RedGreen"      "Barcode"       "ChipType"      "RunInfo"       "Unknowns"
```

Information about the run is stored at:

```r
head(example$RunInfo)
#             RunTime          BlockType
# [1,] "1/12/2016 2:37:05 AM" "Decoding"
# [2,] "4/8/2016 6:03:57 PM"  "Scan"
# [3,] "4/8/2016 6:03:57 PM"  "Register"
```

But most importantly, number of beads, mean and standard deviation of the fluorescence intensity.

```r
head(example$Quants)
#           Mean  SD  NBeads
# 10600313  284  137     13
# 10600322 9405 1363     14
# 10600328 3538  439     11
```

## 2) Extracting raw intensities with the minfi package

The minfi library is a massive library. It contains utils to extract raw information. However, the use of
S4 object oriented language can make it hard for users to identify the right functions. In this tutorial
we have put them all together to ease the painstacking journey through minfi's dense documentation.

To read all IDAT files in a directory, we use the read.metharray.exp function. If additionally, we intend
to read additional information such as the number of beads or the standard deviation of the fluorescence
intensity channels, we will need to set the extended argument to TRUE.

```r
setwd("/media/ben/DATA/Ben/1_evCpGs/data/aging_children/GSE104812_RAW/")
rgSet = read.metharray.exp(getwd(), extended = TRUE)
IDAT_IDs = sapply(strsplit(colnames(rgSet), split = "_"),function(x) x[1])
```

From the resulting RGChannelSetExtended class object, it is possible to extract information for all probe types

```r
annotation <- getAnnotation(rgSet)
TypeI.Red <- getProbeInfo(rgSet, type = "I-Red")
TypeI.Green <- getProbeInfo(rgSet, type = "I-Green")
TypeII <- getProbeInfo(rgSet, type = "II")
ctrls <- getProbeInfo(rgSet, type = "Control")
SnpI <- getProbeInfo(rgSet, type = "SnpI")
SnpII <- getProbeInfo(rgSet, type = "SnpII")
```

Intringuingly, a set of 473 probes (orphan probes) have been placed on the 450K microarray for unknown purposes

```r
known_probes = c(SnpI$AddressA, SnpI$AddressB, SnpII$AddressA, ctrls$Address, TypeI.Red$AddressA, 
                 TypeI.Red$AddressB, TypeI.Green$AddressA, TypeI.Green$AddressB, TypeII$AddressA)
                 
length(known_probes)                  # 621926
all = rownames(rgSet); length(probes) # 622399
orphan = all[!(all %in% known_probes)]; length(missing) # 473
```

In any case, to extract Green/Red fluorescence mean/standard deviation or number of beads from an RGChannelSetExtended object, we employ the function assay:

```r
Grn = assay(rgSet, "Green")       # Green mean across beads
Red = assay(rgSet, "Red")         # Red mean across beads
GrnSD = assay(rgSet, "GreenSD")   # Green SD across beads
RedSD = assay(rgSet, "RedSD")     # Red SD across beads
nBeads = assay(rgSet, "NBeads")   # Number of Beads across probes
```

To convert from probes to CpG sites, we made it easier with the wrapper GR_to_UM (which internally employs unexported minfi:::.preprocessRaw function), return a list containing fluorescence intensities assigned to unmethylated and methylated epiallele.

```r
M_U = GR_to_UM(Red, Grn, rgSet)
M_U_sd = GR_to_UM(RedSD, GrnSD, rgSet)
```
    
To convert nBeads from probes to CpGs, a criteria for type-I probes is required. In beads_GR_to_UM, the minimum number of beads between address-A and -B is selected to represent a CpG targetted by each pair of type-I probes.

```r
nBeads_cg = beads_GR_to_UM(nBeads, rgSet)
```

Finally, to compute a matrix of raw beta-values, one can simply execute:

```r
beta_value = M_U$M/(M_U$M + M_U$U + 100)
```

However, raw beta-values should be avoided as these are affected by within and between array batch effects. Normalisation techniques are thus required, for which a wide variety of R-packages already exist. Here we name the most popular: minfi,  wateRmelon, ENmix, lumi, methylumi, ChAMP, meffil, preprocessCore and EWAStools. UMtools focuses on the analysis of methylation data employ U/M intensity signals rather than the beta-value.


## 3) Quickly importing/exporting with data.table


```r
setwd("/media/ben/DATA/Ben/3_genetic_artefacts/R-packages/test/")
export_bigmat(M_U$M, "M.txt", nThread = 4)
M = import_bigmat("2020-06-19_M.txt", nThread = 4)
```


## 4) UM tools

To start employing some of the functions in UM tools, we will need to extract the phenotypic information from GEO. GEOquery allows to parse from GEO in a minimum number of lines.

```r
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
```

### 4.1) U/M-plots

```r
UM_plot(M = M_U$M, U = M_U$U, CpG = "cg00050873", sex = pheno$sex)
UM_plot(M = M_U$M, U = M_U$U, CpG = "cg00026186", sex = pheno$sex)
```

### 4.2) CV jitter plots

```r
CV = compute_cv(M_U_sd$M, M_U_sd$U, M_U$M, M_U_sd$U)
```


```r
density_jitter_plot(CV, "cg00050873", pheno$sex)
density_jitter_plot(CV, "cg00214611", pheno$sex)
density_jitter_plot(CV, "cg02839557", pheno$sex)
density_jitter_plot(CV, "cg05544622", pheno$sex)
density_jitter_plot(beta_value, "cg00050873", pheno$sex)
```

### 4.3) Bivariate Gaussian Mixture Models (bGMMs)

```r
set.seed(1); bGMM(M_U$M, M_U$U, "cg13293246", 1) # K = 1
set.seed(2); bGMM(M_U$M, M_U$U, "cg03398919", 2) # K = 2
set.seed(3); bGMM(M_U$M, M_U$U, "cg00814218", 3) # K = 3
set.seed(2); bGMM(M_U$M, M_U$U, "cg27024127", 4) # K = 4
set.seed(6); bGMM(M_U$M, M_U$U, "cg23186955", 5) # K = 5
```

### 4.4) CV and BC(CV)
Compute CV per CpG and per sample


```r
BC_CV = compute_BC_CV(CV)
density_jitter_plot(CV, which.max(BC_CV), pheno$sex)
```

### 4.5.1) K-calling with visual output

```r
Kcall_CpG("cg15771735", M_U$M, M_U$U, minPts = 5, reach = seq(0.99, 1.01, 0.01)) # K = 1
Kcall_CpG("cg03398919", M_U$M, M_U$U, minPts = 5, reach = seq(0.99, 1.01, 0.01)) # K = 2
Kcall_CpG("cg00814218", M_U$M, M_U$U, minPts = 5, reach = seq(0.99, 1.01, 0.01)) # K = 3
Kcall_CpG("cg27024127", M_U$M, M_U$U, minPts = 5, reach = seq(0.99, 1.01, 0.01)) # K = 4
```

### 4.5.2) K-calling epigenome-wide

```r
chrY = rownames(annotation)[annotation$chr == "chrY"]
K_vec = par_EW_Kcalling(M_U$M[chrY,], M_U$U[chrY,], minPts = 5, reach = seq(0.99, 1.01, 0.01), R = 2)
```

### 4.6) Comethylation plots




### References and Supporting Information
B. Planterose *et al* (**2018**). Universal epigenetic dissimilarity: DNA methylation variation that is equivalent between monozygotic co-twins and unrelated individuals.






