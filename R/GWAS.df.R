#' GWAS data
#' 
#' A dataset containing GWAS data
#'
#' @docType data
#'
#' @usage data(GWAS.df.example)
#' 
#'
#' @format Dataframe, one row per SNP, with 6 columns (PLINK --logistic/--linear output format)
#' \describe{
#'   \item{CHR}{Chromosome for SNP (sex chromosomes coded numerically). Data type: integer}
#'   \item{POS}{Chromosomal position for each SNP, in base pairs. Data type: integer}
#'   \item{SNP}{Variant ID (such as dbSNP ID "rs..."). Data type: character
#'              Note: Must be the same naming scheme as used in eQTL.df to ensure proper matching.}
#'   \item{P}{p-value for the SNP from GWAS analysis. Data type: numeric}
#'   \item{BETA}{beta for the SNP from GWAS analysis. Data type: numeric}
#'   \item{P}{OPTIONAL Name of the phenotype for which the GWAS data refers. Data type: character
#'            This column is optional and is useful if your GWAS.df contains data for multiple phenotypes, 
#'            such as one might obtain from a PheWAS. If GWAS.df does not contain a PHE column, eQTpLot 
#'            will assume all the supplied GWAS data is for a single phenotype, with a name to be 
#'            specified with the trait argument. }
#' }
#'
"GWAS.df.example"
