#' Gene data
#'
#' A dataset containing Gene data
#'
#' @docType data
#'
#' @usage data(Genes.df.example)
#'
#' @format Dataframe, one row per SNP, with 6 columns
#' \describe {
#'  \item{Gene}{Gene symbol/name for which the Coordinate data refers to. Data type: character 
#'              Note: gene symbol/name must match entries in eQTL.df to ensure proper matching.}
#'  \item{CHR}{Chromosome the gene is on. Data type: integer 
#'             Note: do not include a "car" prefix, and sex chromosomes should be coded numerically.}
#'  \item{Start}{Chromosomal coordinate of start position (in basepairs) to use for gene. Data type: integer 
#'               Note: this should be the smaller of the two values between Start and Stop.}
#'  \item{Stop}{Chromosomal coordinate of end position (in basepairs) to use for gene. Data type: integer 
#'              Note: this should be the larger of the two values between Start and Stop.}
#'  \item{build}{The genome build (either "hg19" or "hg38") for the coordinate data. Data type: character 
#'               The default Genes.df dataframe contains entries for both genome builds for each gene, 
#'               and the script will select the appropriate entry based on the specified gbuild 
#'               (default is hg19)).}
#' }
#'
"Genes.df.example"
