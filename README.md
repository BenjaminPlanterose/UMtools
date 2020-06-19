# UMtools
## An R-package for analysing Illumina DNA Methylation microarrays at the fluorescence intensity level


#### Benjamin Planterose Jiménez, Manfred Kayser, Athina Vidaki

### Department of Genetic Identification, Erasmus MC University Medical Centre Rotterdam, The Netherlands

## License
[MIT](https://choosealicense.com/licenses/mit/)


## What is UMtools?


## Tested on

    Ubuntu 18.04.4 LTS (bionic), R version 3.6.3 (2020-02-29) -- "Holding the Windsock"
    Contact b.planterosejimenez@erasmusmc.nl for any issues arising while running UMtools.
    
## Installation 

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
install_github("BenjaminPlanterose/UMtools")
```


# UMtools tutorial
    
Do not attempt to perform this tutorial without at least 8GB of RAM. Working with fluorescence intensities
involves the use of several large matrices. To avoid any issues, we recommend to always monitor the consumption of RAM via htop:

```bash
htop
```


## 0) A word on the Beadchip microarray technology and the probes included in the 450K: 

The microarray itself consists of a silica substrate with uniformly interspaced microwells.
Hundreds of thousands of copies of a specific oligonucleotide lie on the surface of silica beads.
During the manufacture of the chip, a total of 622,399 types of beads are pooled together and depeosited on
the microarray. Subsequently, beads automatically self-assemble on the microarray's microwell. As a result, both
the order and the number copies for a given bead type are random, hence requiring the decoding of the microarray.

Although employing 622,399 probe oligonucleotides, the 450K technology targes 485,512 cytosines. We register here the count of probes. To begin with, 450K combines three types of probes concerning detection:

* Type I (n = 135,476 x 2) - two bead types per cytosine 

  * Type-I Green (n = 46,289 x 2) - Quantification is informative only in the Green channel
  
  * Type-I Red (n = 89,187 x 2) - Quantification is informative only in the Red channel
  
* Type-II (n = 350,036): one bead type per cytosine, quantification is informative in both channels

In addition, there are a wide range of quality control probes (n = 848):

* Staining (n = 4)

* extension (n = 4)

* hybridization (n = 3)

* target removal (n = 2)

* bisulfite conversion I and II (n = 12 and 4)

* specificity I and II (n = 12 and 3) 

* non-polymorphic (n = 4)

* negative control (n = 613)

* restoration (n = 1)

* normalization (n = 186).

As well, there are probes for targetting SNPs rather than CpGs for assessing sample mix-up (n = 65):

* SnpI (n = 25 x 2)

* SnpII (n = 40)

Finally, there 473 orphan probes, placed on the array for unknown purposes. In total, that makes:

    473 + 25*2 + 40 + 848 + 46,289 * 2 + 89,187 * 2 + 350,036 = 622,399 probes



## 1) Peaking into an IDAT file

The .IDAT extension (Intensity Data) corresponds to Illumina's proprietary format for storage of microarray scanners'
raw fluorescence output among several genome-wide platform. The IDAT format is encrypted and non-human readable. Before the R-package illuminaio was developed, no alternatives to vendor's software existed in order to read IDAT files.

We first need to download some example files from the GEO database. Throughout this tutorial, we will be using the data from Shi *et al* (GEO_ID = GSE104812), consisting of whole blood DNA methylation data from healthy children. To download the data, you may ran the following bash commands:

```bash
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE104nnn/GSE104812/suppl/GSE104812_RAW.tar --wait=10 --limit-rate=50K
tar -xvf GSE104812_RAW.tar
find . -type f ! -name '*.idat.gz' -delete
gunzip *.gz
```
Or if prefered, go to https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE104812 and download TAR (of IDAT) via http.

Back to R, we firstly load illuminaio (a dependancy of minfi that is automatically downloaded with it):

```r
library(illuminaio)
```

To peak into an IDAT file, in this case of the green channel of random sample, we can ran:

```r
setwd("/media/ben/DATA/Ben/1_evCpGs/data/aging_children/GSE104812_RAW/")
example = illuminaio::readIDAT("GSM2808239_sample1_Grn.idat")
```

An IDAT contains the following fields:

```r
names(example)
# [1] "fileSize"  "versionNumber" "nFields"   "fields"  "nSNPsRead" "Quants" "MidBlock"
# [8] "RedGreen"  "Barcode"       "ChipType"  "RunInfo" "Unknowns"
```

For example, information about the run is stored at:

```r
head(example$RunInfo)
#             RunTime          BlockType
# [1,] "1/12/2016 2:37:05 AM" "Decoding"
# [2,] "4/8/2016 6:03:57 PM"  "Scan"
# [3,] "4/8/2016 6:03:57 PM"  "Register"
```

But most importantly, it contains number of beads, mean and standard deviation of the fluorescence intensity.

```r
head(example$Quants)
#           Mean  SD  NBeads
# 10600313  284  137     13
# 10600322 9405 1363     14
# 10600328 3538  439     11
```

The information about the standard deviation of fluorescence intensities is highly valuable and has rarely made it to the literature.


## 2) Extracting raw intensities

The minfi R-package is a massive library that has set the standards of quality in epigenomics. 
However, its extensive use of S4-object oriented language can make it hard for users to find and repurpose functions.
In this tutorial, we have gathered together handy pieces of code for analysis of raw fluorescence intensities to ease the painstacking journey through minfi's dense documentation.

We firstly load all required libraries:

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

To read all IDAT files in a directory, we use the read.metharray.exp function. If additionally, we intend
to read additional information such as the number of beads or the standard deviation of the fluorescence
intensity channels, we will need to set the extended argument to TRUE.

```r
setwd("/media/ben/DATA/Ben/1_evCpGs/data/aging_children/GSE104812_RAW/")
rgSet = read.metharray.exp(getwd(), extended = TRUE)
```

From the resulting RGChannelSetExtended class object, it is possible to extract the annotation for all probe types

```r
TypeI.Red <- getProbeInfo(rgSet, type = "I-Red")
TypeI.Green <- getProbeInfo(rgSet, type = "I-Green")
TypeII <- getProbeInfo(rgSet, type = "II")
ctrls <- getProbeInfo(rgSet, type = "Control")
SnpI <- getProbeInfo(rgSet, type = "SnpI")
SnpII <- getProbeInfo(rgSet, type = "SnpII")
```

For example:

```r
# Type I probes have two addresses corresponding to methylated and unmethylated probes
head(TypeI.Green[, 1:3], n = 2)
#          Name    AddressA    AddressB
#   <character> <character> <character>
# 1  cg02004872    25785404    58629399
# 2  cg02050847    43656343    73683470

# Type II probes have only one address. Methylated and unmethylated corresponds 
# to the same probe in two different fluorescent channels
head(TypeII[, 1:3], n = 2)
#          Name    AddressA
#   <character> <character>
# 1  cg00035864    31729416
# 2  cg00061679    28780415
```

For are more thorough annotation of type-I and II probes, we can execute:

```r
annotation <- getAnnotation(rgSet)
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
```


As a side note, a set of 473 probes (orphan probes) have been placed on the 450K microarray for unknown purposes

```r
known_probes = c(SnpI$AddressA, SnpI$AddressB, SnpII$AddressA, ctrls$Address, TypeI.Red$AddressA, 
                 TypeI.Red$AddressB, TypeI.Green$AddressA, TypeI.Green$AddressB, TypeII$AddressA)
                 
length(known_probes)                  # 621926
all = rownames(rgSet); length(known_probes) # 622399
orphan = all[!(all %in% known_probes)]; length(orphan) # 473
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

Internally, it performs the following:

| Probe Type    |  Methylated       |     Unmethylated  |
|:-------------:|:-----------------:|:-----------------:|
| Type-II       | Green (addressA)  |  Red (addressA)   |
| Type-I Green  | Green (addressB)  | Green (addressA)  |
| Type-I Red    | Red (addressB)    |  Red (addressA)   |


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






