# UMtools
## An R-package for analysing Illumina DNA Methylation microarrays at the fluorescence intensity level

#### Benjamin Planterose Jiménez, Manfred Kayser, Athina Vidaki

### Department of Genetic Identification 
#### Erasmus MC University Medical Centre Rotterdam, The Netherlands


## Requirements

    Operating system: tested on Ubuntu 18.04.4 LTS (bionic)
    R: tested on R version 3.6.3 (2020-02-29) -- "Holding the Windsock"

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

### The content of an IDAT file

The .IDAT extension (Intensity Data) is Illumina's proprietary format for storage of the fluorescence scanners'
raw output across several genome-wide platform. The IDAT format is an encrypted and non-human readable.
It was not until the birth the R-package illuminaio, that IDAT files could only be read via vendor's software.

Here we read the green channel of random sample

    setwd("/media/ben/DATA/Ben/1_evCpGs/data/aging_children/GSE104812_RAW/")
    example = illuminaio::readIDAT("GSM2808239_sample1_Grn.idat")
    names(example)
    # [1] "fileSize"      "versionNumber" "nFields"       "fields"        "nSNPsRead"     "Quants"        "MidBlock"
    # [8] "RedGreen"      "Barcode"       "ChipType"      "RunInfo"       "Unknowns"

    head(example$RunInfo)
    #             RunTime          BlockType
    # [1,] "1/12/2016 2:37:05 AM" "Decoding"
    # [2,] "4/8/2016 6:03:57 PM"  "Scan"
    # [3,] "4/8/2016 6:03:57 PM"  "Register"

    head(example$Quants)
      #           Mean  SD  NBeads
      # 10600313  284  137     13
      # 10600322 9405 1363     14
      # 10600328 3538  439     11

Mean: mean Grn intensity across beads, SD Grn intensity across beads and number of beads



