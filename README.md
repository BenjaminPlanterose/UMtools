# UMtools
## An R-package for analysing Illumina DNA Methylation microarrays at the fluorescence intensity level

#### Benjamin Planterose Jiménez, Manfred Kayser, Athina Vidaki
#### Department of Genetic Identification, Erasmus MC University Medical Centre Rotterdam, The Netherlands



## Why UMtools?

A great range of R-packages have already been developed to analyze data from Illumina's DNA methylation microarray platforms such as minfi, wateRmelon, ENmix, ChAMP, lumi, methylumi, meffil, EWAStools, etc.
Where does UMtools fit in this ecosystem? UMtools focuses on the low-level analysis of Illumina DNA methylation microarray data, at the level of fluorescence intensities. 

We believe that we can harvest much more from the IDAT file, Illumina's propietary format. For example, the standard deviation across beads has rarely been mentioned in the 
literature and could be used in applications studying to the technical noise of DNA methylation microarray platforms.

We have additionally included new tools such as CVlogT, BC(CVlogT), K-caller or comethylation plots that we developed for the verification of genetic artefacts in the 450K array 
(B. Planterose *et al* **2021**) but that could well be employed for other applications.

Finally, as most libraries tend to hide the initial steps of analysis (unexported functions, lack of documentation), we have rescued code (especially from the minfi 
R-package) and wrapped it to ease working at this level.

## Tested on

    Ubuntu 18.04.4 LTS (bionic)
    R version 3.6.3 (2020-02-29) -- "Holding the Windsock"
    
For issues arising while running UMtools, either [report an issue](https://github.com/BenjaminPlanterose/UMtools/issues) or simply contact b.planterosejimenez@erasmusmc.nl
    
## Installation 

To install an R-package from Github, the library devtools is required:

```r
install.packages("devtools")
library(devtools)
```

To install UMtools dependencies:

```r
# From Bioconductor
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install('minfi')

# From CRAN
install.packages(parallel)
install.packages(scales)
install.packages(EMCluster)
install.packages(dbscan)
install.packages(RColorBrewer)

# From Github
devtools::install_version("modes", "0.7.0")
devtools::install_github("dphansti/Sushi")
devtools::install_github("BenjaminPlanterose/UMtools")
```

To install other R-packages required for this tutorial:

```r
BiocManager::install('GEOquery')
BiocManager::install('IlluminaHumanMethylation450kanno.ilmn12.hg19')
BiocManager::install('IlluminaHumanMethylation450kmanifest')
```


# About the tutorial
    
With this tutorial we aim to increase the technical accesibility of the technology. We cover topics for which information is hard to find such as the Beadchip technology, 
control probes, the content of the IDAT file, as well as how to extract all of this information from an example dataset at GEO. In addition, we provide a summary of all 
the tools available at UMtools. 

However, make sure to have at least 8GB of RAM available. Working with fluorescence intensities involves large matrices, meaning heavy RAM usage. To avoid any issues, we recommend to always monitor resources via htop if working from a Linux machine:

```bash
sudo apt-get install htop
htop
```

Also, when deleting large objects in R, you may call the garbage collector to quickly repurpose RAM:

```r
help(gc)
gc()
```

All R commands are idented on the Wiki. But if prefered, you may find them all in one script at the following link
[tutorial.R](https://github.com/BenjaminPlanterose/UMtools/tree/master/tutorial/tutorial.R)


## A word on the Beadchip microarray technology and the probes on the 450K

Initially developed for genotyping arrays, a BeadChip microarray consists of a silicon substrate with regularly interspaced micro-wells. Silicon micro-beads randomly self-assemble in the wells. 
Each bead is covered by hundreds of thousands of copies of the same 50-nucleotide long probe. 

During the manufacture of the chip, a total of 622,399 types of beads, each bearing a different probe, are pooled together and deposited on
the microarray. Subsequently, beads automatically self-assemble on the microarray's microwell.

Detection of DNA methylation is done by coupling single-nucleotide variant detection with bisulfite conversion, by which unmethylated cytosines are converted to uraciles while leaving methylated cytosines unchanged. Two approaches are simultaneously performed on the microarray:

  * Infinium type-II: A single probe targets both epialleles by performing single base extension at CpG site positions +1 or +2 (depending on which strand is targetted). Thus, the addition of the nucleotide is dependent of the methylation status. As a result, Type-II probes are informative in both channels.
  
  * Infinium type-I: Relies on two probes, each targetting either the unmethylated converted CpG or the methylated unconverted CpG. Single base extension occurs one nucleotide before or after the CpG site at positions 0 or +3 (depending on which strand is targetted). Thus, the addition of the nucleotide is independent of the methylation status, and hence, both probes are registed in the same fluorescence channel but at two different beadtypes.
  
The addition of type-I probes together with the inclusion of control bead types that are not informative for CpG methylation explains why 485,512 cytosines are targeted by 622,399 different bead types. We have compiled a thorough count all the probes included in the 450K:

* Type I (n = 135,476 x 2) - two bead types per cytosine.

  * Type-I Green (n = 46,289 x 2) - Quantification is informative only in the Green channel.
  
  * Type-I Red (n = 89,187 x 2) - Quantification is informative only in the Red channel.
  
* Type-II (n = 350,036): one bead type per cytosine, quantification is informative in both channels.

* Control probes (n = 848).

  * Staining (n = 4).

  * extension (n = 4).

  * hybridization (n = 3).

  * target removal (n = 2).

  * bisulfite conversion I and II (n = 12 and 4).

  * specificity I and II (n = 12 and 3).

  * non-polymorphic (n = 4).

  * negative control (n = 613).

  * restoration (n = 1).

  * normalization (n = 186).

* SNP-targetting probes (n = 65) - included for assessing sample mix-up.

  * SnpI (n = 25 x 2) - Two bead types per SNP.

  * SnpII (n = 40) - One bead type per cytosine.

* orphan probes (n = 473) - placed on the array for unknown purposes. 


In total, that makes:

    473 + 25*2 + 40 + 848 + 46,289 * 2 + 89,187 * 2 + 350,036 = 622,399 probes


Also, during the manufacture, beads are pooled in equal ratio and then deposited on the array. 
As a result, the exact number of copies of each beadType is not controlled (varies from chip to chip). To assign the correspondence between microwells and bead types, decoding is required. 
This is done during manufacture via consecutive hybridizations with other sets of probes that target the address, a 23 nucleotide-long oligonucleotide handle (address) that links the bead 
to the probe in the form of a DMAP file (Nakabayashi 2017).


## The 450K protocol

The protocol to process samples goes as follows:

  * Genomic DNA extraction.
  
  * Bisulfite conversion - Cytosines<sup>Unmethylated</sup> become Uraciles and Cytosines<sup>Methylated</sup> remain Cytosines.
  
  * Whole-genome amplification - Uraciles become Thymines.
  
  * Enzymatic DNA fragmentation.
  
  * Hybridization to the microarray + Washing.
  
  * Infinium assay - incubation with a DNA polymerase and dideoxynucleotides-triphosphate (ddNTPs): ddATP and ddTTP labelled with biotin, ddCTP and ddGTP labelled with dinitrophenol (DNP). Upon single-based extension (SBE), elongation cannot continue due to the dideoxy nature of the incorporated nucleotide. Staining is carried out by incubating with fluorophore-labelled acceptors: Red-fluorescing Cy5-labelled anti-DNP (targetting ddA/T) and Green-fluorescing Cy3-labelled streptavidin (targeting ddC/G).
  
  * Fluorescence scanning of the microarray in the Green and Red channels with iScan/HiScan confocal laser microarray scanner. Fluorescence intensity information is stored as two IDAT files, one per fluorescence channel. 
  


## Peaking into an IDAT file

The .IDAT extension (Intensity Data) corresponds to Illumina's proprietary format for storage of microarray scanners'
raw fluorescence output among several genome-wide platform. The IDAT format is encrypted and non-human readable. Before the R-package illuminaio was developed, no alternatives to vendor's software existed in order to read IDAT files.

To unveil the content of the IDAT file, we first need to download some examples from the GEO database. Throughout this tutorial, we will be using the data from Shi *et al* (GEO_ID = GSE104812), consisting of whole blood DNA methylation data from healthy children. To download the data, you may ran the following bash commands:

```bash
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE104nnn/GSE104812/suppl/GSE104812_RAW.tar --wait=10 --limit-rate=50K
tar -xvf GSE104812_RAW.tar
find . -type f ! -name '*.idat.gz' -delete
gunzip *.gz
```
Or if prefered, it is also possible to go to the following url and manually download the TAR (of IDAT) via http https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE104812

Back to R, to peak into an IDAT file, in this case of the green channel of random sample, we can run:

```r
library(illuminaio)
#setwd("~/foo/") # Change this route to fit your system
example = illuminaio::readIDAT("GSM2808239_sample1_Grn.idat")
```

The output of readIDAT contains the following fields:

```r
names(example)
# [1] "fileSize"  "versionNumber" "nFields"   "fields"  "nSNPsRead" "Quants" "MidBlock"
# [8] "RedGreen"  "Barcode"       "ChipType"  "RunInfo" "Unknowns"
```

For example, deep information about the processing of the sample is stored at:

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


## Extracting fluorescence intensity matrices

For this stage, we will be mainly using functions from the minfi package and in some circumstances, functions from UMtools that encapsulate unexported code from the minfi's R-package.

We firstly load UMtools:

```r
library(UMtools)

```

To read all IDAT files in a directory, we use the *minfi::read.metharray.exp* function. If additionally, we intend
to read additional information such as the number of beads or the standard deviation of the fluorescence
intensity channels, we will need to set the *extended* argument to TRUE.

```r
setwd("~/foo/") # Change this route to fit your system
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

As a side note, a set of 473 probes (orphan probes) have been placed on the 450K microarray for unknown purposes

```r
known_probes = c(SnpI$AddressA, SnpI$AddressB, SnpII$AddressA, ctrls$Address, TypeI.Red$AddressA, 
                 TypeI.Red$AddressB, TypeI.Green$AddressA, TypeI.Green$AddressB, TypeII$AddressA)
                 
length(known_probes)                  # 621926
all = rownames(rgSet); length(known_probes) # 622399
orphan = all[!(all %in% known_probes)]; length(orphan) # 473
```

Diving into the annotation, Type I probes have two probes or addresses (A and B) corresponding to unmethylated and methylated probes, respectively:
```r
head(TypeI.Green[, 1:3], n = 2)
#          Name    AddressA    AddressB
#   <character> <character> <character>
# 1  cg02004872    25785404    58629399
# 2  cg02050847    43656343    73683470
```

Type II probes have only one probe or address (A) which covers both methylated and unmethylated intensities:
```r
head(TypeII[, 1:3], n = 2)
#          Name    AddressA
#   <character> <character>
# 1  cg00035864    31729416
# 2  cg00061679    28780415
```

For a biological annotation of the CpGs targetted by type-I and II probes, we can execute:

```r
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
```

In any case, to extract Green/Red fluorescence mean/standard deviation or number of beads from an RGChannelSetExtended object, we employ the function assay:

```r
Grn = assay(rgSet, "Green")       # Green mean across beads
Red = assay(rgSet, "Red")         # Red mean across beads
GrnSD = assay(rgSet, "GreenSD")   # Green SD across beads
RedSD = assay(rgSet, "RedSD")     # Red SD across beads
nBeads = assay(rgSet, "NBeads")   # Number of Beads across probes
```

To convert from probes to CpG sites, we made it easier with the wrapper GR_to_UM (which in turn employs the unexported function *minfi:::.preprocessRaw*), returning a list that contains methylated and unmethylated fluorescence intensities.
To convert nBeads from probes to CpGs, a criteria for type-I probes is required. In GR_to_UM, we select the smallest number of beads between addressA and addressB to represent a CpG targetted by type-I probe pairs.

```r
M_U = GR_to_UM(Red = Red, Grn = Grn, rgSet = rgSet, what = "Mean")
M_U_sd = GR_to_UM(Red = RedSD, Grn = GrnSD, rgSet = rgSet, what = "SD")
nBeads_cg = GR_to_UM(nBeads = nBeads, rgSet = rgSet, what = "NBeads")
```

Internally, it assigns the methylated and unmethylated intensity values as followed:

| Probe Type    |  Methylated       |     Unmethylated  |
|:-------------:|:-----------------:|:-----------------:|
| Type-II       | Green (addressA)  |  Red (addressA)   |
| Type-I Green  | Green (addressB)  | Green (addressA)  |
| Type-I Red    | Red (addressB)    |  Red (addressA)   |


Finally, to obtain a matrix of raw beta-values or M-values, one can simply compute:

```r
# Beta-values
offset = 100 # For numerical stability at low fluorescence intensities
beta_value = M_U$M/(M_U$M + M_U$U + offset)

# M-values
offset = 1 # For numerical stability at low fluorescence intensities
M_value = log2((M_U$M + offset)/(M_U$U + offset))
```

A wide range of tools are available to preprocess raw signals such as batch effect correction, background correction, colour channel balance, type-I and II correction or normalisation. In this tutorial we will however focus on
raw signals. Bare in mind that applications such as EWAS benefit enourmously from further preprocessing but that is discussed elsewhere; for more information, check the documentations for minfi, wateRmelon, ENmix, lumi, methylumi, ChAMP, meffil, preprocessCore and EWAStools. 





To clean up the workspace, we can perform:
```r
rm(Grn, Red, GrnSD, RedSD, nBeads, nBeads_cg, M_value); gc()
```

## Quickly importing/exporting large matrices with data.table

During the analysis of DNA methylation microarray data, reading IDATs can be time-consuming.  To reduce waiting time, we wrapped data.table's ultra-fast and RAM-efficient routines to conveniently import and export epigenomics matrices (import_bigmat and export_bigmat, respectively):

```r
setwd("~/foo/") # Change this route to fit your system
export_bigmat(M_U$M, "M.txt", nThread = 4)
list.files() # "2021-01-28_M.txt"
M = import_bigmat("2021-01-28_M.txt", nThread = 4)
```

## Importing phenotypes with GEOquery

In addition, GEOquery allows to parse phenotypes from the GEO database in minimum number of lines of code:

```r
library(GEOquery)
setwd('~/foo/')
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

## UMtools in action

UMtools was initially developed to characterise the influence of genetic artefacts on the fluorescence intensity signals of Illumina's DNA methylation microarrays. But its use for exploring the behaviour of probes can be further extended to any probe.

Moving from univariate methylation value to the bivariate the U/M plane not only offers a gain of resolution, but also the possibility to distinguish probe failure from intermediate methylation: when a probe fails, background fluorescence is acquired in both channels of roughly the same scale; as a result the ratio to the total intensity tends to 50 %. 

For example, the following probe targets the Y-chromosome and hence fails in females, giving rise to beta-values tending towards 0.5 (slightly skewed to the left by the established offset):
```r
density_jitter_plot(beta_value, "cg00050873", pheno$sex)
annotation["cg00050873", c("chr", "pos")] # chrY   9363356
```
![Alt text](img/img1.png?raw=true "cg00026186 U/M plot")

Unlike on the methylation scale, failed samples cluster are obvious in the U/M-plane:
```r
UM_plot(M = M_U$M, U = M_U$U, CpG = "cg00050873", sex = pheno$sex)
annotation["cg00050873", c("chr", "pos")] # chrY   9363356
```
![Alt text](img/img2.png?raw=true "cg00050873 U/M plot")


U/M plots are a highly interpretable way to visualize a wide variety of phenomena. For example, X-inactivation: due to the double gene dosage of Chromosome-X in females, one of the copies is randomly inactivated via large-scale targetted methylation. As a result, females are 50 % methylated while males are 0 or 100 % methylated:

```r
UM_plot(M = M_U$M, U = M_U$U, CpG = "cg00026186", sex = pheno$sex)
annotation["cg00026186", c("chr", "pos")] # chrX  48367230
```
![Alt text](img/img3.png?raw=true "cg00026186 U/M plot")


Also, some loci escape X-inactivation and as a result, are unmethylated in both males and females. However, given the double-copy of X-chromosomes in females, the unmethylated intensity is higher on average than in males
```r
UM_plot(M = M_U$M, U = M_U$U, CpG = "cg04927982", sex = pheno$sex)
annotation["cg04927982", c("chr", "pos")] # chrX  53254653
```
![Alt text](img/img4.png?raw=true "cg04927982 U/M plot")


Finally, some loci are hypermethylated in both females and males:
```r
UM_plot(M = M_U$M, U = M_U$U, CpG = "cg02973417", sex = pheno$sex)
annotation["cg02973417", c("chr", "pos")] # chrX  153220895
```
![Alt text](img/img5.png?raw=true "cg02973417 U/M plot")





### Bivariate Gaussian Mixture Models (bGMMs)

The formation of clusters in the U/M plane is a consistent feature of a genetic artefacts and meQTLs. 
To assign cluster identity to samples when the number of expected clusters is known (for example, after a U/M plot), we developed bivariate Gaussian mixture models (bGMMs).
In such approach, we wrapped the routines from the EMCluster. Here are some examples:

```r
set.seed(2); bGMM(M_U$M, M_U$U, "cg03398919", K = 2)
```
![Alt text](img/img6.png?raw=true "cg03398919 bGMM")


```r
set.seed(3); bGMM(M_U$M, M_U$U, "cg00814218", K = 3)
```
![Alt text](img/img7.png?raw=true "cg00814218 bGMM")


```r
set.seed(2); bGMM(M_U$M, M_U$U, "cg27024127", K = 4)
```
![Alt text](img/img8.png?raw=true "cg27024127 bGMM")


```r
set.seed(6); bGMM(M_U$M, M_U$U, "cg23186955", K = 5)
```
![Alt text](img/img9.png?raw=true "cg23186955 bGMM")


Knowing the assignation of samples to clusters can be used to estimate minor allele frequencies (MAF).



### Quantifying epigenome-wide ambivalency in probe failure

For probes suffering from a genetic variant that causes probe failure, a duality in probe efficiency can be observed: for some individuals it fails, for others it does not. We first define, CV, as the coefficient of variation of the log of the total intensity, defined as:

<img src="https://render.githubusercontent.com/render/math?math=CV_{ln(U %2B M)} = \dfrac{\hat{\sigma}_{ln(U %2B M)}}{\hat{\mu}_{ln(U %2B M)}}">

CV is a measure of noise-to-signal ratio. CV is highly bimodal when a probe fails due to a genetic artefact (females display large noise-to-signal ratio than males for chrY-targeting probes).

```r
CV = compute_CV(M_SD = M_U_sd$M, U_SD = M_U_sd$U, M = M_U$M, U = M_U_sd$U, alpha = 100)
density_jitter_plot(CV, "cg00050873", pheno$sex)
```
![Alt text](img/img10.png?raw=true "cg00050873 jitter")


Bimodality can be quantified by a *bimodality coefficient* as a function of the sample skewness (<img src="https://render.githubusercontent.com/render/math?math=\gamma">) and kurtosis (<img src="https://render.githubusercontent.com/render/math?math=\kappa">):

<img src="https://render.githubusercontent.com/render/math?math=BC(CV) = \dfrac{\hat{\gamma}_{CV} %2B 1}{\hat{\kappa}_{CV} %2B \dfrac{3(n-1)^2}{(n-2)(n-3)}}">

BC(CV) is defined in the range [0,1] and BC(CV)> 5/9 can be used as evidence for multimodality.

```r
BC_CV = compute_BC_CV(CV = CV)
BC_CV["cg00050873"]
# cg00050873 
#   1.128555  # higher than 1, because sample estimators are used to computed kurtosis and skewness and the sample size employed is small.
annotation["cg00050873", c("chr", "pos")] # chrY   9363356
```


### K-calling

Cluster formation in the UM plane cannot be tested epigenome-wide by bGMM since we do not know the target number of clusters. Because of this, we developed a K-caller, a tool able to automatically detect how many clusters are formed in the U/M plane based on the density-based spatial clustering of applications with noise (dbscan) algorithm.

On the one hand, the function *Kcall_CpG* provides of a visual output. Please note that outliers are not included in any class and that that K-calling is performed in a transformed scale to reduce the ellipticity of clusters in the U/M plane.

```r
Kcall_CpG("cg15771735", M_U$M, M_U$U, minPts = 5, eps = 0.1)
# [1] 1
```
![Alt text](img/img11.png?raw=true "cg15771735")


```r
Kcall_CpG("cg03398919", M_U$M, M_U$U, minPts = 5, eps = 0.1)
# [1] 2
```
![Alt text](img/img12.png?raw=true "cg03398919")


```r
Kcall_CpG("cg00814218", M_U$M, M_U$U, minPts = 5, eps = 0.1)
# [1] 3
```
![Alt text](img/img13.png?raw=true "cg00814218")



```r
Kcall_CpG("cg27024127", M_U$M, M_U$U, minPts = 5, eps = 0.1)
# [1] 4
```
![Alt text](img/img14.png?raw=true "cg27024127")


On the other hand, function *par_EW_Kcalling* is the parallel-computing version for epigenome-wide K-calling. Running it on all chrY probes (n = 416), we observe the following:

```r
chrY = rownames(annotation)[annotation$chr == "chrY"]
K_vec = par_EW_Kcalling(M_U$M[chrY,], M_U$U[chrY,], minPts = 5, eps = 0.1, nThread = 10)
table(K_vec)
# K_vec
# 1   2
# 35 381
```

It was expected that given that probe targeting ChrY fail for females, that two clusters are always formed. The observed K = 1 probes consist of cross-reactive probes and a minor component of K-calling mistakes. For example, this supposedly Y-Chr CpG clearly displays a strong methylated signal in females (thus, cross-reactive):

```r
names(which(K_vec == 1))[2] # "cg02494853"
UM_plot(M = M_U$M, U = M_U$U, CpG = "cg02494853", sex = pheno$sex)
annotation["cg02494853", c("chr", "pos")] # chrY   4868397
```
![Alt text](img/img15.png?raw=true "cg02494853")

As a word of notice, parameters *eps* and *minPts* depend on the sample size and hence require dataset-specific tuning. To enable the training of these parameters,
we manually compiled a training set that includes probes forming from 1 to 4 clusters in the U/M plane on a dataset of 426 MZ twin pairs of the E-risk cohort (British, EUR ancestry). 
However, please note that the number of clusters observed may be different in other datasets depending on ancestry and sample size and hence, the training set needs to be 
manually curated. In this specific dataset, we have a much smaller sample size and a different ancestry so we expect strong differences.
We can confirm our suspicions by: 

```r
data("training_set")
Kcall_CpG(sample(training_set$k_1, 1), M_U$M, M_U$U, minPts = 5, eps = 0.1) # U/M plot of a random probe in the set K = 1
Kcall_CpG(sample(training_set$k_2, 1), M_U$M, M_U$U, minPts = 5, eps = 0.1) # U/M plot of a random probe in the set K = 2
Kcall_CpG(sample(training_set$k_3, 1), M_U$M, M_U$U, minPts = 5, eps = 0.1) # U/M plot of a random probe in the set K = 3
Kcall_CpG(sample(training_set$k_4, 1), M_U$M, M_U$U, minPts = 5, eps = 0.1) # U/M plot of a random probe in the set K = 4
```

Though the training set has not been curated on this dataset, for the sake of displaying how to use this UMtools functionality, we display how performance for a set of *eps* and *minPts* can be computed:

```r
train_k_caller(M_U$M, M_U$U, training_set, 3, 0.07, nThread = 10)
```

The obtained confusion matrix, macroPrecision, macroRecall and macro-F1 score aim to assess the performance of the multi-class classification task.


### Comethylation plots

Finally, the formation of clusters is not a strict signature for a genetic artefact. It is possible that nearby CpGs, not affecting the Infinium assay in any way, influence the methylation status of the region: these are the so-called methylation quantitative trait loci (meQTL).


This CpG is forming three clusters:
```r
UM_plot(M = M_U$M, U = M_U$U, CpG = "cg14911689", sex = NULL)
annotation["cg14911689", c("chr", "pos", "UCSC_RefGene_Name")] # chr12    739980    NINJ2
```
![Alt text](img/img16.png?raw=true "cg14911689")


However, we observe neighbour CpGs form similar shapes:

```r
annotation <- annotation[order(annotation$chr, annotation$pos),] # sort annotation by chromosome and position
pos <- which(rownames(annotation) == "cg14911689")
UM_plot(M = M_U$M, U = M_U$U, CpG = rownames(annotation)[pos-1], sex = NULL)
UM_plot(M = M_U$M, U = M_U$U, CpG = rownames(annotation)[pos+1], sex = NULL)
```

![Alt text](img/img17.png?raw=true "cg26371957")
![Alt text](img/img18.png?raw=true "cg26654770")


Moreover, we can examine if neighbouring CpGs methylation status is correlated via:

```r
res = Visualize_cometh(annotation = annotation, CpG = 'cg14911689', distance = 1000,
                       L_bound = 3, R_bound = 2, beta_mat = beta_value, max_y = 5)
```

![Alt text](img/img19.png?raw=true "cg01201512")


This is not the case for a genetic artefact such as:

```r
UM_plot(M = M_U$M, U = M_U$U, CpG = "cg11495604", sex = NULL)
res = Visualize_cometh(annotation = annotation, CpG = 'cg11495604', distance = 1000,
                       L_bound = 0, R_bound = 2, beta_mat = beta_value, max_y = 5)
annotation["cg11495604", c("chr", "pos", "UCSC_RefGene_Name")] # chr20  62053198 KCNQ2;KCNQ2;KCNQ2;KCNQ2
```
![Alt text](img/img20.png?raw=true "cg11495604")
![Alt text](img/img21.png?raw=true "cg11495604")




A final word of notice, two assumptions are implicit under this approach:

  * Probes are available at a distance close enough to display co-methylation. This is not always the case
  
  * Nearby probes are not affected by genetic artefacts.

Asumption 2 is great violated at HLA loci, known to be hotspots for genetic variability. As a result, nearby CpGs atefactually affected by underlying genetic variants display
artefactual comethylation due to linkage disequilibrium between underlying variant. As an example:

```r
res = Visualize_cometh(annotation = annotation, CpG = 'cg00211215', distance = 200,
                       L_bound = 3, R_bound = 0, beta_mat = beta_value, max_y = 5)
annotation["cg00211215", c("chr", "pos", "UCSC_RefGene_Name")] # chr6  32552246   HLA-DRB1
```
![Alt text](img/img22.png?raw=true "cg11495604")


### Annotations included in UMtools

UMtools also contains some helpful annotations: 

```r
data(annot_450K) # Genetic variants associated to Illumina Infinium HumanMethylation450 Beadchip probes
data(annot_EPIC) # Genetic variants associated to Illumina Infinium MethylationEPIC Beadchip probes
data(classification_CpG_SNP_450K) # Classification of CpG/SBE-SNPs in the Illumina Infinium HumanMethylation450 Beadchip microarray
data(classification_CpG_SNP_EPIC) # Classification of CpG/SBE-SNPs in the Illumina Infinium MethylationEPIC Beadchip microarray
data(CR_probes) # List of in silico-predicted cross-reactive probes in the Illumina Infinium HumanMethylation450 Beadchip microarray
data(triallelic_CpG_SNP_450K) # Tri-allelic SNPs associated to Illumina Infinium HumanMethylation450 Beadchip probes
data(triallelic_CpG_SNP_EPIC) # Tri-allelic SNPs associated to Illumina Infinium MethylationEPIC Beadchip probes
data(training_set) # K-caller (Training Set; EUR, Illumina Infinium HumanMethylation450 Beadchip microarray)
```


### References

B. Planterose *et al* (**2021**). Revisiting genetic artefacts on DNA methylation microarrays: implications for meQTL mapping. *Genome research*






