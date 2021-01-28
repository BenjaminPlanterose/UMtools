#' Genetic variants associated to Illumina Infinium HumanMethylation450 Beadchip probes
#'
#'
#'
#' @docType data
#'
#' @usage data(annot_450K)
#'
#' @format An object of class list
#'
#' @keywords datasets
#'
#' @references Revisiting Genetic artefacts on DNA methylation microarrays. Genome Research
#'
#' @description Includes a list containing processed vcf files in 6 categories:
#' \itemize{
#'   \item CpG_SNP
#'   \item SBE_SNP
#'   \item probe_SNP
#'   \item CpG_indel
#'   \item SBE_indel
#'   \item probe_indel
#' }
#'
#' These files were obtained by intersecting 450K probe coordinates with common variants in dbSNP151
#' file 00-common_all.vcf. Additionally, we segregated SNPs from indels and kept only variants with
#' MAF > 0.01. All scripts employed can be found at UMtools/scripts/.
#'
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
#'   \item cg: CpG identifier (see IlluminaHumanMethylation450kanno.ilmn12.hg19)
#' }
#'
#'
#' @examples
#' data(annot_450K)
"annot_450K"
