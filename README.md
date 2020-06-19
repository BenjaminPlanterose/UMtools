# UMtools
## An R-package for analysing Illumina DNA Methylation microarrays at the fluorescence intensity level



#### Benjamin Planterose Jiménez, Manfred Kayser, Athina Vidaki

### Department of Genetic Identification, Erasmus MC University Medical Centre Rotterdam, The Netherlands


## Tested on:

    Ubuntu 18.04.4 LTS (bionic), R version 3.6.3 (2020-02-29) -- "Holding the Windsock"
    Contact b.planterosejimenez@erasmusmc.nl for any issues arising while running UMtools.
    
## Dependencies 

    library(minfi)
    library(modes)
    library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
    library(IlluminaHumanMethylation450kmanifest)
    library(GEOquery)
    library(scales)
    library(EMCluster)
    library(dbscan)
    library(Sushi)
    library(preprocessCore)

## About this tutorial
    
Do not attempt to perform this tutorial without at least 8GB of RAM. Working with fluorescence intensities
involves several large matrices

# 0) A word on the Beadchip microarray technology and the IDAT format: 

The microarray itself consists of a silica substrate with uniformly interspaced microwells.
Hundreds of thousands of copies of a specific oligonucleotide lie on the surface of silica beads.
During the manufacture of the chip, all 622,399 type of beads are pooled together and applied on
the microarray. Subsequently, beads automatically self-assemble on the microarray's microwell. As a result, both
the order and the number copies for a given bead type are random, hence requiring the decoding of the microarray.

Although targeting 485,512 cytosines, the 450K technology employs 622,399 probe oligonucleotides. There are
several reasons for this. To begin with, 450K combines three types of probes concerning detection:
Type-I Green (n = 46,289 x 2): two bead types per cytosine, quantification is informative only in the Green channel
Type-I Red (n = 89,187 x 2): two bead types per cytosine, quantification is informative only in the Red channel
Type-II (n = 350,036): one bead type per cytosine, quantification is informative in both channels

In addition, there are a set of quality control probes (n = 848):
Staining (n = 4), extension (n = 4), hybridization (n = 3), target removal (n = 2), bisulfite conversion I (n = 12),
bisulfite conversion II (n = 4), specificity I (n = 12), specificity II (n = 3), non-polymorphic (n = 4),
negative (n = 613), restoration (n = 1) and normalization (n = 186).

Also, there are sample mix-up probes (n = 65):
SnpI (n = 25 x 2) and SnpII (n = 40)
Finally, there 473 orphan probes with placed on the array for unknown purposes

    473 + 25*2 + 40 + 848 + 46289 * 2 + 89187 * 2 + 350036 = 622399



# 1) Introduction to the IDAT format

The .IDAT extension (Intensity Data) is Illumina's proprietary format for storage of the fluorescence scanners'
raw output across several genome-wide platform. The IDAT format is an encrypted and non-human readable.
It was not until the birth the R-package illuminaio, that IDAT files could only be read via vendor's software.
We first need to download some example IDAT files from the GEO database. You may perform: 

    wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE104nnn/GSE104812/suppl/GSE104812_RAW.tar --wait=10 --limit-rate=50K
    tar -xvf GSE104812_RAW.tar
    find . -type f ! -name '*.idat.gz' -delete
    gunzip *.gz

Or manually download at https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE104812

To peak into an IDAT file, in this case of the green channel of random sample, we can ran:

    setwd("/media/ben/DATA/Ben/1_evCpGs/data/aging_children/GSE104812_RAW/")
    example = illuminaio::readIDAT("GSM2808239_sample1_Grn.idat")
    names(example)
    # [1] "fileSize"      "versionNumber" "nFields"       "fields"        "nSNPsRead"     "Quants"        "MidBlock"
    # [8] "RedGreen"      "Barcode"       "ChipType"      "RunInfo"       "Unknowns"

Information about the run is stored at:

    head(example$RunInfo)
    #             RunTime          BlockType
    # [1,] "1/12/2016 2:37:05 AM" "Decoding"
    # [2,] "4/8/2016 6:03:57 PM"  "Scan"
    # [3,] "4/8/2016 6:03:57 PM"  "Register"
    
But most importantly, number of beads, mean and standard deviation of the fluorescence intensity.

    head(example$Quants)
      #           Mean  SD  NBeads
      # 10600313  284  137     13
      # 10600322 9405 1363     14
      # 10600328 3538  439     11


# 2) Extracting raw intensities with minfi

The minfi library is a massive library. It contains utils to extract raw information. However, the use of
S4 object oriented language can make it hard for users to identify the right functions. In this tutorial
we have put them all together to ease the painstacking journey through minfi's documentation.

To read all IDAT files in a directory, we use the read.metharray.exp function. If additionally, we intend
to read additional information such as the number of beads or the standard deviation of the fluorescence
intensity channels, we will need to set the extended argument to TRUE.

    setwd("/media/ben/DATA/Ben/1_evCpGs/data/aging_children/GSE104812_RAW/")
    rgSet = read.metharray.exp(getwd(), extended = TRUE)
    IDAT_IDs = sapply(strsplit(colnames(rgSet), split = "_"),function(x) x[1])

From the resulting RGChannelSetExtended class object, it is possible to extract information of all types of probes

    annotation <- getAnnotation(rgSet)
    TypeI.Red <- getProbeInfo(rgSet, type = "I-Red")
    TypeI.Green <- getProbeInfo(rgSet, type = "I-Green")
    TypeII <- getProbeInfo(rgSet, type = "II")
    ctrls <- getProbeInfo(rgSet, type = "Control")
    SnpI <- getProbeInfo(rgSet, type = "SnpI")
    SnpII <- getProbeInfo(rgSet, type = "SnpII")

Intringuingly, a set of 473 probes (orphan probes) have been placed on the microarray for unknown purposes

    known_probes = c(SnpI$AddressA, SnpI$AddressB, SnpII$AddressA, ctrls$Address, TypeI.Red$AddressA, TypeI.Red$AddressB,
                     TypeI.Green$AddressA, TypeI.Green$AddressB, TypeII$AddressA); length(known_probes) # 621926
    all = rownames(rgSet); length(probes) # 622399
    orphan = all[!(all %in% known_probes)]; length(missing) # 473


In any case, to extract Green/Red fluorescence mean/standard deviation or number of beads from an RGChannelSetExtended object, we employ the function assay:

    Grn = assay(rgSet, "Green")       # Green mean across beads
    Red = assay(rgSet, "Red")         # Red mean across beads
    GrnSD = assay(rgSet, "GreenSD")   # Green SD across beads
    RedSD = assay(rgSet, "RedSD")     # Red SD across beads
    nBeads = assay(rgSet, "NBeads")   # Number of Beads across probes

To convert from probes to CpG sites, we included wrappers GR_to_UM and beads_GR_to_UM.

    M_U = GR_to_UM(Red, Grn, rgSet)
    M_U_sd = GR_to_UM(RedSD, GrnSD, rgSet)
    nBeads_cg = beads_GR_to_UM(nBeads, rgSet)

GR_to_UM employs unexported minfi:::.preprocessRaw function.


# 3) UM tools
From this point on, we start relying on functions from UMtools. 








### References and Supporting Information
B. Planterose *et al* (**2018**). Universal epigenetic dissimilarity: DNA methylation variation that is equivalent between monozygotic co-twins and unrelated individuals.






