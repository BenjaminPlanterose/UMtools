# UMtools
## An R-package for analysing Illumina DNA Methylation microarrays at the fluorescence intensity level

#### Benjamin Planterose Jiménez, Manfred Kayser, Athina Vidaki
#### Department of Genetic Identification, Erasmus MC University Medical Centre Rotterdam, The Netherlands



## Why UMtools?

Several R-packages have been developed to analyze data from Illumina's DNA methylation microarray platforms. In pursuit for the appealing of a wide range of end users, these however tend to sacrifice modularisation for simplicity: standarised pipelines are conceived to be deployed as one rather than being composed by repurposable modules. 

As a result, the average user performs rudimentary examination of the raw data dissuaded by the scarcity of available tools targetting the beginning of the analysis. Very often, these tools already exist but are embedded within established pipelines, are unexported to the main R-package and lack written documentation. 

For all the above, UMtools was developed as modular R-package that focuses on the low-level analysis of Illumina DNA methylation microarray data. Instead of the methylation ratio or beta-value, UMtools analyses fluorescence intensity means and standard deviation across beads, stored at the very heart of the IDAT file, Illumina's propietary format.


## Tested on

    Ubuntu 18.04.4 LTS (bionic)
    R version 3.6.3 (2020-02-29) -- "Holding the Windsock"
    
    Contact b.planterosejimenez@erasmusmc.nl for any issues arising while running UMtools.
    
## Installation 

```r
# From Bioconductor
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install('minfi')
BiocManager::install('GEOquery')

# From CRAN
install.packages(parallel)
install.packages(modes)
install.packages(scales)
install.packages(EMCluster)
install.packages(dbscan)

# From Github
library("devtools")
install_github("dphansti/Sushi")
install_github("BenjaminPlanterose/UMtools")
```

# About the tutorial
    
With this tutorial we aim to increase the technical accesibility of the technology. We cover topics for which information is hard to find such as the Beadchip technology, control probes, the content of the IDAT file, as well as how to extract all of this information from an example dataset at GEO. In addition, we provide a summary of all the tools available at UMtools. 

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

The Beadchip technology is the basis of Illumina DNA methylation microarrays. On the one hand, it is a probe-based approach where hundreds of thousands of copies of a specific 50 nucleotide-long probes lie on the surface of silica beads. On the other hand, The microarray itself consists of a silica substrate with uniformly interspaced microwells.
During the manufacture of the chip, a total of 622,399 types of beads, each bearing a different probe, are pooled together and deposited on
the microarray. Subsequently, beads automatically self-assemble on the microarray's microwell. As a result, both
the order and the number copies for each bead type are random.
To assign the correspondance between microwells and bead types, decoding is required. This is done during manufacture via consecutive hybridizations with other sets of probes that target the address, a 22 nucleotide-long oligonucleotide handle that links the bead to the probe (stored as a DMAP file).

Detection of DNA methylation is done by coupling single-nucleotide variant detection with bisulfite conversion, by which unmethylated cytosines are converted to uraciles while leaving methylated cytosines unchanged. Two approaches are simultaneously performed on the microarray:

  * Infinium type-II: A single probe targets both epialleles by performing single base extension at CpG site positions +1 or +2 (depending on which strand is targetted). Thus, the addition of the nucleotide is dependent of the methylation status. As a result, Type-II probes are informative in both channels.
  
  * Infinium type-I: Relies on two probes, each targetting either the unmethylated converted CpG or the methylated unconverted CpG. Single base extension occurs one nucleotide before or after the CpG site at positions 0 or +3 (depending on which strand is targetted). Thus, the addition of the nucleotide is independent of the methylation status, and hence, both probes are registed in the same fluorescence channel but at two different beadtypes.
  
The addition of type-I probes together with the inclusion of control bead types that are not informative for CpG methylation explains why 485,512 cytosines are targetted by 622,399 different bead types. We have compiled a thorough count all the probes included in the 450K:

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


## The 450K protocol

As a summary of the previous section, probes (50 nt) are tethered on the surface of each bead type (622,399 in total), connected via the address (23 nt). A random number for each bead type self-assemble on random micro-wells of the microarray. Decoding is performed to assign a bead type to each microwell, via consecutive hybridizations. Concerning CpG methylation detection, for the Infinium type-I assay, two bead types are required while for Infinium type-II assay, only a single bead type is required.

With a working microarray with known mapping between bead types and positions (in the shape of a DMAP file), the preparation of samples goes as follows:

  * Genomic DNA extraction.
  
  * Bisulfite conversion - Cytosines<sup>Unmethylated</sup> become Uraciles and Cytosines<sup>Methylated</sup> remain Cytosines.
  
  * Whole-genome amplification - Uraciles become Thymines.
  
  * Enzymatic DNA fragmentation.
  
  * Hybridization to the microarray + Washing.
  
  * Infinium assay - incubation with a DNA polymerase and dideoxynucleotides-triphosphate (ddNTPs): ddATP and ddTTP labelled with biotin, ddCTP and ddGTP labelled with dinitrophenol (DNP). Upon single-based extension, elongation cannot continue due to the dideoxy nature of the incorporated nucleotide. Staining is carried out by incubating with fluorophore-labelled acceptors: Red-fluorescing Cy5-labelled anti-DNP (targetting ddA/T) and Green-fluorescing Cy3-labelled streptavidin (targetting ddC/G).
  
  * Fluorescence scanning of the microarray in the Green and Red channels with iScan/HiScan.
  
  * Fluorescence intensity information is stored as two IDAT files, one per fluorescence channel.


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

Back to R, to peak into an IDAT file, in this case of the green channel of random sample, we ran:

```r
library(illuminaio)
setwd("/media/ben/DATA/Ben/1_evCpGs/data/aging_children/GSE104812_RAW/") # Change this route to fit your system
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

We particularly highlight the standard deviation of fluorescence intensities, as it is highly valuable and has rarely made it to the literature.


## Extracting fluorescence intensity matrices

The Bioconductor-based **minfi** library is a huge R-package that has set the standards of quality in methylomics data analysis. Though it contains a vast amount of code that could find other applications, its deep encapsullation via an S4 object-oriented implementation can make it hard for users to find and repurpose low-level functions, especially unexported functions, which are not even on the documentation. As a result, minfi is designed to run as a standarized pipeline rather than an adaptable toolset.

Given that UMtools depends on some basic processing operations already established in the minfi library, we have compiled in this tutorial some handy functions, many of which unexported or not included on minfi's documentation. 

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
library(UMtools)
library(parallel)
```

To read all IDAT files in a directory, we use the *minfi::read.metharray.exp* function. If additionally, we intend
to read additional information such as the number of beads or the standard deviation of the fluorescence
intensity channels, we will need to set the *extended* argument to TRUE.

```r
setwd("/media/ben/DATA/Ben/1_evCpGs/data/aging_children/GSE104812_RAW/") # Change this route to fit your system
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

In any case, to extract Green/Red fluorescence mean/standard deviation or number of beads from an RGChannelSetExtended object, we employ the function assay:

```r
Grn = assay(rgSet, "Green")       # Green mean across beads
Red = assay(rgSet, "Red")         # Red mean across beads
GrnSD = assay(rgSet, "GreenSD")   # Green SD across beads
RedSD = assay(rgSet, "RedSD")     # Red SD across beads
nBeads = assay(rgSet, "NBeads")   # Number of Beads across probes
```

To convert from probes to CpG sites, we made it easier with the wrapper GR_to_UM (which in turn employs the unexported function *minfi:::.preprocessRaw*), returning a list that contains methylated and unmethylated fluorescence intensities.

```r
M_U = GR_to_UM(Red, Grn, rgSet)
M_U_sd = GR_to_UM(RedSD, GrnSD, rgSet)
```

Internally, it assigns the methylated and unmethylated intensity values as followed:

| Probe Type    |  Methylated       |     Unmethylated  |
|:-------------:|:-----------------:|:-----------------:|
| Type-II       | Green (addressA)  |  Red (addressA)   |
| Type-I Green  | Green (addressB)  | Green (addressA)  |
| Type-I Red    | Red (addressB)    |  Red (addressA)   |


To convert nBeads from probes to CpGs, a criteria for type-I probes is required. In beads_GR_to_UM, we select the smallest number of beads between addressA and addressB to represent a CpG targetted by type-I probe pairs.

```r
nBeads_cg = beads_GR_to_UM(nBeads, rgSet)
```

Finally, to obtain a matrix of raw beta-values, one can simply compute:

```r
offset = 100 # For numerical stability at low fluorescence intensities
beta_value = M_U$M/(M_U$M + M_U$U + offset)
```

Also quite popular, one can compute M-values as:
```r
offset = 1 # For numerical stability at low fluorescence intensities
M_value = log2((M_U$M + offset)/(M_U$U + offset))
```

However, raw beta and M-values should be avoided for further analysis as these display strong within and between array batch effects, especially on large datasets. Normalisation techniques are thus required, for which a wide variety of R-packages can be deployed. Here, we name the most popular: minfi, wateRmelon, ENmix, lumi, methylumi, ChAMP, meffil, preprocessCore and EWAStools. Unlike the names R-packages, UMtools does not focus on the analysis of the methylation values, but rather on the U/M intensity signals directly.


To clean up the workspace, we can perform:
```r
rm(Grn, Red, GrnSD, RedSD, nBeads, nBeads_cg, M_value); gc()
```

## Quickly importing/exporting large matrices with data.table

During the analysis of DNA methylation microarray data, to avoid having to read IDATs again and again, it is prefarable to read once and export/import the rest of the times. As R-base cannot cope with large files, we made use of data.table ultra-fast and RAM-efficient routines to build up wrappers to conveniently import and export epigenomics matrices (import_bigmat and export_bigmat, respectively):

```r
setwd("/media/ben/DATA/Ben/3_genetic_artefacts/R-packages/test/") # Change this route to fit your system
export_bigmat(M_U$M, "M.txt", nThread = 4)
M = import_bigmat("2020-06-19_M.txt", nThread = 4)
```

## Quickly importing phenotypes with GEOquery

In addition, GEOquery allows to parse phenotypes from the GEO database in minimum number of lines of code:

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

## UMtools in action

UMtools was initially developed to characterise the influence of genetic artefacts on the fluorescence intensity signals of Illumina's DNA methylation microarrays. But its use for exploring the behaviour of probes can be further extended to any probe.

Moving from univariate methylation value to the bivariate the U/M plane not only offers a gain of resolution, but also the possibility to distinguish probe failure from intermediate methylation: when a probe fails, background fluorescence is acquired in both channels of roughly the same scale; as a result the ratio to the total intensity tends to 50 %. For example, the following probe targets the Y-chromosome and hence fails in females, giving rise to beta-values tending towards 0.5 (slightly skewed to the left by the established offset):
```r
density_jitter_plot(beta_value, "cg00050873", pheno$sex)
annotation["cg00050873", c("chr", "pos")] # chrY   9363356
```
![Alt text](img/jitter_betaval.png?raw=true "cg00026186 U/M plot")

Unlike on the methylation scale, failed samples cluster at the origin of the UM-plane:
```r
UM_plot(M = M_U$M, U = M_U$U, CpG = "cg00050873", sex = pheno$sex)
annotation["cg00050873", c("chr", "pos")] # chrY   9363356
```
![Alt text](img/UM.png?raw=true "cg00050873 U/M plot")

UM plots wide variety of phenomena such as X-inactivation: due to the double gene dosage of Chromosome-X in females, one of the copies is randomly inactivated via large-scale targetted methylation. As a result, females are 50 % methylated while males are 0 or 100 % methylated:

```r
UM_plot(M = M_U$M, U = M_U$U, CpG = "cg00026186", sex = pheno$sex)
annotation["cg00026186", c("chr", "pos")] # chrX  48367230
```
![Alt text](img/x_inact.png?raw=true "cg00026186 U/M plot")


However, some loci escape X-inactivation and as a result, are unmethylated in both males and females. However, given the double-copy of X-chromosomes in females, the unmethylated intensity is higher on average than in males
```r
UM_plot(M = M_U$M, U = M_U$U, CpG = "cg04927982", sex = pheno$sex)
annotation["cg04927982", c("chr", "pos")] # chrX  53254653
```
![Alt text](img/escape.png?raw=true "cg04927982 U/M plot")


Some probes are cross-reactive, e.g. they hybridize at several loci in the genome. Although supposedly targetting an autosomal locus, this probe looks exactly like X-inactivation due to its cross-reactivity towards the chromosome X:
```r
UM_plot(M = M_U$M, U = M_U$U, CpG = "cg20926353", sex = pheno$sex)
annotation["cg20926353", c("chr", "pos")] # chr9  84303358
```
![Alt text](img/CR.png?raw=true "cg20926353 U/M plot")


In this case, the targetted locus escapes X-inactivation but on top, the probe is cross-reactive to the Y-chromosome (males get an extra M-signal from the Y-chromosome):
```r
UM_plot(M = M_U$M, U = M_U$U, CpG = "cg26738106", sex = pheno$sex)
annotation["cg26738106", c("chr", "pos")] # chrX   3265038
```
![Alt text](img/CR_2.png?raw=true "cg26738106 U/M plot")


Genetic artefacts such as SNPs or indels can cause probe failure when homozygous. In this case, a SNP causes a 3'-overhang, making it impossible to develop the single-base extension step:
```r
UM_plot(M = M_U$M, U = M_U$U, CpG = "cg03398919", sex = NULL)
annotation["cg03398919", c("chr", "pos")] # chr2 173118470
```
![Alt text](img/PF.png?raw=true "cg03398919 U/M plot")


But other times, the SNP is confused for the U/M epiallele. When the biological context has the opposite state (SNP = U in a methylated region or viceversa), it gives rise to 3 clusters (homozygous and heterozygous for both allele). Here, SNP = U:
```r
UM_plot(M = M_U$M, U = M_U$U, CpG = "cg00814218", sex = NULL)
annotation["cg00814218", c("chr", "pos")] # chr14  37445440
```
![Alt text](img/SNP_U.png?raw=true "cg00814218 U/M plot")


And here, SNP = M:
```r
UM_plot(M = M_U$M, U = M_U$U, CpG = "cg17004290", sex = NULL)
annotation["cg17004290", c("chr", "pos")] # chr4 108853384
```
![Alt text](img/SNP_M.png?raw=true "cg17004290 U/M plot")


Sometimes, two SNPs interact giving rise to a mixture between probe failure and SNP confused for one of the epialleles:
```r
UM_plot(M = M_U$M, U = M_U$U, CpG = "cg27024127", sex = NULL)
annotation["cg27024127", c("chr", "pos")] # chr8  27522576
```
![Alt text](img/2SNP.png?raw=true "cg27024127 U/M plot")



Finally, we show a probe affected by a SNP on top of being cross-reactive:
```r
UM_plot(M = M_U$M, U = M_U$U, CpG = "cg23186955", sex = pheno$sex)
annotation["cg23186955", c("chr", "pos")] # chr3  16420793
```
![Alt text](img/CR_SNP.png?raw=true "cg23186955 U/M plot")



### Bivariate Gaussian Mixture Models (bGMMs)

As we have seen in the previous section, the formation of clusters in the UM plane is a consistent feature of a genetic artefact. For this reason, we tested several clustering techniques in the UM plane. The most succesful approach in cases where the number of expected clusters is known, were bivariate Gaussian mixture models (bGMMs). We wrapped the routines from the EMCluster library for straighforward deployment on epigenomic data. Here some case examples from K = {2, 3, 4, 5}.

```r
set.seed(2); bGMM(M_U$M, M_U$U, "cg03398919", K = 2)
```
![Alt text](img/bGMM2.png?raw=true "cg03398919 bGMM")


```r
set.seed(3); bGMM(M_U$M, M_U$U, "cg00814218", K = 3)
```
![Alt text](img/bGMM3.png?raw=true "cg00814218 bGMM")


```r
set.seed(2); bGMM(M_U$M, M_U$U, "cg27024127", K = 4)
```
![Alt text](img/bGMM4.png?raw=true "cg27024127 bGMM")


```r
set.seed(6); bGMM(M_U$M, M_U$U, "cg23186955", K = 5)
```
![Alt text](img/bGMM5.png?raw=true "cg23186955 bGMM")




### Quantifying epigenome-wide ambivalency in probe failure

For probes suffering from a genetic variant that causes probe failure, a duality in probe efficiency can be observed: for some individuals it fails, for others it does not. We first define, CV, as the coefficient of variation of the log of the total intensity, defined as:

<img src="https://render.githubusercontent.com/render/math?math=CV_{ln(U %2B M)} = \dfrac{\hat{\sigma}_{ln(U %2B M)}}{\hat{\mu}_{ln(U %2B M)}}">

CV is a measure of noise-to-signal ratio and can be simply computed by *compute_cv*. CV is highly bimodal when a probe fails due to a genetic artefact.

```r
CV = compute_CV(M_U_sd$M, M_U_sd$U, M_U$M, M_U_sd$U)
density_jitter_plot(CV, "cg00050873", pheno$sex)
```
![Alt text](img/jitter_CV.png?raw=true "cg00050873 jitter")


Bimodality can be quantified by a *bimodality coefficient*:

<img src="https://render.githubusercontent.com/render/math?math=BC(CV) = \dfrac{\hat{\gamma}_{CV} %2B 1}{\hat{\kappa}_{CV} %2B \dfrac{3(n-1)^2}{(n-2)(n-3)}}">

BC(CV) can be computed for all CpGs with *compute_BC_CV*, rendering a good measure for ambivalency in probe failure.

```r
BC_CV = compute_BC_CV(CV)
BC_CV["cg00050873"]
# cg00050873 
#   1.128741  
annotation["cg00050873", c("chr", "pos")] # chrY   9363356
```


### K-calling

Cluster formation in the UM plane cannot be tested epigenome-wide by bGMM since we do not know the target number of clusters. Because of this, we developed a K-caller, a tool able to automatically detect how many clusters are formed in the U/M plane based on the density-based spatial clustering of applications with noise (dbscan) algorithm.

On the one hand, the function *Kcall_CpG* provides of a visual output. Please note that outliers are not included in any class and that the scale of the plot has been transformed to reduce the ellipticity of clusters in the U/M plane.

```r
Kcall_CpG("cg15771735", M_U$M, M_U$U, minPts = 5, reach = seq(0.99, 1.01, 0.01))
# [1] 1
```
![Alt text](img/K_call1.png?raw=true "cg15771735")


```r
Kcall_CpG("cg03398919", M_U$M, M_U$U, minPts = 5, reach = seq(0.99, 1.01, 0.01))
# [1] 2
```
![Alt text](img/K_call2.png?raw=true "cg03398919")


```r
Kcall_CpG("cg00814218", M_U$M, M_U$U, minPts = 5, reach = seq(0.99, 1.01, 0.01))
# [1] 3
```
![Alt text](img/K_call3.png?raw=true "cg00814218")



```r
Kcall_CpG("cg27024127", M_U$M, M_U$U, minPts = 5, reach = seq(0.99, 1.01, 0.01))
# [1] 4
```
![Alt text](img/K_call4.png?raw=true "cg27024127")


On the other hand, function *par_EW_Kcalling* is the parallel-computing version for epigenome-wide K-calling. Running on all chrY probes (n = 416), we observe the following:

```r
chrY = rownames(annotation)[annotation$chr == "chrY"]
K_vec = par_EW_Kcalling(M_U$M[chrY,], M_U$U[chrY,], minPts = 5, reach = seq(0.99, 1.01, 0.01), R = 2)
table(K_vec)
# K_vec
# 1   2
# 38 378
```

It was expected that given that probe targeting ChrY fail for females, that two clusters are always formed. The observed K = 1 probes consist of cross-reactive probes and a minor component of K-calling mistakes. For example, this supposedly Y-Chr CpG clearly displays a strong methylated signal in females (thus, cross-reactive):

```r
names(which(K_vec == 1))[2] # "cg02494853"
UM_plot(M = M_U$M, U = M_U$U, CpG = "cg02494853", sex = pheno$sex)
annotation["cg02494853", c("chr", "pos")] # chrY   4868397
```
![Alt text](img/y_cr.png?raw=true "cg02494853")


### Comethylation plots

Finally, the formation of clusters is not a strict signature for a genetic artefact. It is possible that nearby CpGs, not affecting the Infinium assay in any way, influence the methylation status of the region: these are the so-called methylation quantitative trait loci (meQTL).


This CpG is forming three clusters:
```r
UM_plot(M = M_U$M, U = M_U$U, CpG = "cg14911689", sex = NULL)
annotation["cg14911689", c("chr", "pos", "UCSC_RefGene_Name")] # chr12    739980    NINJ2
```
![Alt text](img/mQTL.png?raw=true "cg14911689")


However, we observe neighbour CpGs form similar shapes:

```r
annotation <- annotation[order(annotation$chr, annotation$pos),]
pos <- which(rownames(annotation) == "cg14911689")
UM_plot(M = M_U$M, U = M_U$U, CpG = rownames(annotation)[pos-1], sex = NULL)
UM_plot(M = M_U$M, U = M_U$U, CpG = rownames(annotation)[pos+1], sex = NULL)
```

![Alt text](img/mQTL2.png?raw=true "cg26371957")
![Alt text](img/mQTL3.png?raw=true "cg26654770")


Moreover, we can examine if neighbouring CpGs methylation status is correlated via:

```r
res = Visualize_cometh(annotation = annotation, CpG = 'cg14911689', distance = 1000,
                       L_bound = 3, R_bound = 2, beta_mat = beta_value,
                       cgHeightLabel = -1, deltaposHeightLabel = -0.4, chrHeightLabel = -2,
                       max_y = 5)
```

![Alt text](img/cometh_mQTL.png?raw=true "cg01201512")


This is not the case for a genetic artefact such as:

```r
UM_plot(M = M_U$M, U = M_U$U, CpG = "cg11495604", sex = NULL)
res = Visualize_cometh(annotation = annotation, CpG = 'cg11495604', distance = 1000,
                       L_bound = 0, R_bound = 2, beta_mat = beta_value,
                       cgHeightLabel = -1, deltaposHeightLabel = -0.4, chrHeightLabel = -2,
                       max_y = 5)
annotation["cg11495604", c("chr", "pos", "UCSC_RefGene_Name")] # chr20  62053198 KCNQ2;KCNQ2;KCNQ2;KCNQ2
```
![Alt text](img/PF2.png?raw=true "cg11495604")

![Alt text](img/cometh_mQTL2.png?raw=true "cg11495604")


Only with the exception of fields of genetic variants such as in HLA loci
```r
res = Visualize_cometh(annotation = annotation, CpG = 'cg00211215', distance = 200,
                       L_bound = 3, R_bound = 0, beta_mat = beta_value,
                       cgHeightLabel = -1, deltaposHeightLabel = -0.4, chrHeightLabel = -2,
                       max_y = 5)
annotation["cg00211215", c("chr", "pos", "UCSC_RefGene_Name")] # chr6  32552246   HLA-DRB1
```
![Alt text](img/cometh_mQTL3.png?raw=true "cg11495604")






### References

B. Planterose *et al* (**2018**). Universal epigenetic dissimilarity: DNA methylation variation that is equivalent between monozygotic co-twins and unrelated individuals.






