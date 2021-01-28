#' Classification of CpG/SBE-SNPs in the Illumina Infinium MethylationEPIC Beadchip microarray
#'
#'
#'
#' @docType data
#'
#' @usage data(classification_CpG_SNP_EPIC)
#'
#' @format An object of class list
#'
#' @keywords datasets
#'
#' @references Revisiting Genetic artefacts on DNA methylation microarrays. Genome Research
#'
#' @description  Includes a list containing processed vcf files in 16 categories:
#' \itemize{
#'   \item II_plus_SNP1_C_TA: Infinium type-II probe, targeting (+)-strand, SNP at position 1, Alleles: C ⇔ A/T
#'   \item II_plus_SNP1_C_G: Infinium type-II probe, targeting (+)-strand, SNP at position 1, Alleles: C ⇔ G
#'   \item II_plus_SNP2_G_ACT: Infinium type-II probe, targeting (+)-strand, SNP at position 2, Alleles: G ⇔ A/C/T
#'   \item II_minus_SNP1_C_AGT: Infinium type-II probe, targeting (-)-strand, SNP at position 1, Alleles: C ⇔ A/G/T
#'   \item II_minus_SNP2_G_TA: Infinium type-II probe, targeting (-)-strand, SNP at position 2, Alleles: G ⇔ A/T
#'   \item II_minus_SNP2_G_C: Infinium type-II probe, targeting (-)-strand, SNP at position 2, Alleles: G ⇔ C
#'   \item I_plus_SNP0_A_C_T: Infinium type-I probe, targeting (+)-strand, SNP at position 0, Alleles: A ⇔ C ⇔ T
#'   \item I_plus_SNP0_ACT_G: Infinium type-I probe, targeting (+)-strand, SNP at position 0, Alleles: A/C/T ⇔ G
#'   \item I_plus_SNP1_C_T: Infinium type-I probe, targeting (+)-strand, SNP at position 1, Alleles: C ⇔ T
#'   \item I_plus_SNP1_C_AG: Infinium type-I probe, targeting (+)-strand, SNP at position 1, Alleles: C ⇔ A/G
#'   \item I_plus_SNP2_G_ACT: Infinium type-I probe, targeting (+)-strand, SNP at position 2, Alleles: G ⇔ A/C/T
#'   \item I_minus_SNP1_C_AGT: Infinium type-I probe, targeting (-)-strand, SNP at position 1, Alleles: A ⇔ A/G/T
#'   \item I_minus_SNP2_G_A: Infinium type-I probe, targeting (-)-strand, SNP at position 2, Alleles: G ⇔ A
#'   \item I_minus_SNP2_G_CT: Infinium type-I probe, targeting (-)-strand, SNP at position 2, Alleles: G ⇔ C/T
#'   \item I_minus_SNP3_A_G_T: Infinium type-I probe, targeting (-)-strand, SNP at position 3, Alleles: A ⇔ G ⇔ T
#'   \item I_minus_SNP3_AGT_C: Infinium type-I probe, targeting (-)-strand, SNP at position 3, Alleles: A/G/T ⇔ C
#' }
#'
#' Categories are mutually exclusive. Triallelic SNPs have been excluded at this point. For more information
#' concerning the logic of this classification, please check the associated paper (see References).
#' Column names correspond to the following:
#' \itemize{
#'   \item SNP_chr: SNP chromosome number
#'   \item SNP_pos: SNP position (hg19/GRCh37)
#'   \item rs: SNP identifier in dbSNP
#'   \item REF: Reference allele
#'   \item ALT: Alternate allele
#'   \item QUAL: Quality (not available on dbSNP, hence, simply .)
#'   \item FILTER: Filter (not available on dbSNP, hence, simply .)
#'   \item INFO: Contains metadata of the variant; for more information, consult .vcf formating rules (https://www.ncbi.nlm.nih.gov/snp/docs/products/vcf/redesign/).
#'   Here, in summary, we highlight: GENEINFO (closest gene to the SNP), CAF (Common Allele frequency) and VC (variant class: either SNV (single-nucleotide variant) or DIV (deletion and insertion variants))
#'   \item CpG_chr: CpG chromosome number
#'   \item CpG_start: probe region start (hg19/GRCh37), input into bedtools (consider that bedtools employed intervals from the form (a,b]
#'   \item CpG_end: probe region end (hg19/GRCh37), input into bedtools (consider that bedtools employed intervals from the form (a,b]
#'   \item cg: CpG identifier (see IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
#' }
#'
#' @examples
#' data(classification_CpG_SNP_EPIC)
"classification_CpG_SNP_EPIC"
