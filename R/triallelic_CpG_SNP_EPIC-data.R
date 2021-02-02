#' Tri-allelic SNPs associated to Illumina Infinium MethylationEPIC Beadchip probes
#' @docType data
#' @usage data(triallelic_CpG_SNP_EPIC)
#' @format An object of class data.frame
#' @keywords datasets
#' @references Revisiting Genetic artefacts on DNA methylation microarrays. Genome Research
#' @description Includes a matrix that was obtained by intersecting 850K probe coordinates with common variants in dbSNP151
#' file 00-common_all.vcf. We considered tri-allelic those variants for which three alleles display allelic frequencies larger than 0.05.
#' Scripts employed can be found at UMtools/scripts/.
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
#' @examples
#' data(triallelic_CpG_SNP_EPIC)
"triallelic_CpG_SNP_EPIC"
