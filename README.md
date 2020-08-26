# eQTpLot
## Visualization of Colocalization Between eQTL and GWAS Data

eQTpLot is an R package developed for the visualization of colocalization between eQTL and GWAS data. It requires a number of R packages (biomaRt , dplyr, GenomicRanges, ggnewscale, ggplot2, ggpubr, gridExtra, Gviz, patchwork) and takes as input two data frames (one of GWAS data, and the other of eQTL data), with the user specifying the name of the gene to be analyzed, the GWAS trait to be analyzed (useful if the GWAS data contains information on multiple associations, as one might obtain from a PheWAS), and the tissue type to use for the eQTL analysis (useful if eQTL data is available on multiple tissue types. A PanTissue analysis can be specified as well, combining data across tissue types for each variant). Additional parameters may be specified, including the p-value thresholds for GWAS or eQTL significance, the genomic range to be displayed, axis/layout modifications for the resultant graphs, etc. This data is then used to generate and output a series of plots visualizing colocalization, correlation, and enrichment between eQTL and GWAS signals for a given gene-trait pair.


## Installation

eQTpLot can be install using `devtools`, either directly from GitHub,

`devtools::install_github(RitchieLab/eQTpLot)`

or by downloading the repository to your computer, unzipping, and installing the `eQTpLot` folder.

`devtools::install("eQTpLot")`

*Note: For issues installing dependencies, try running the folling code prior to onstallation.

`Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS=TRUE)`

## Input files
eQTpLot requires two data files, read into R as data frames and passed to the function as the arguments `GWAS.df` and `eQTL.df`. Optionally, a third data frame can be passed to the function as `Genes.df` (discussed more below). A description of the required data is contained here:

#### GWAS.df

`GWAS.df` takes as input a data frame of standard GWAS data that must contain the following 6 columns:

  `CHR` Chromosome for SNP (X coded numerically as 23)  
  `POS` Chromosomal position for each SNP, in base pairs  
  `SNP` Variant ID (such as dbSNP ID "rs...". *Note: Must be the same naming scheme as used in `eQTL.df` to ensure proper matching)*  
  `pvalue_GWAS` p-value for the SNP from GWAS analysis  
  `beta_GWAS` beta for the SNP from GWAS analysis  
  `Trait` name of trait for which GWAS was run. This column must be present, with a value for every row, even if your `GWAS.df` contains data for only one trait. *Note: `GWAS.df` can contain multiple traits (such as one might obtain from a PheWAS. Specifying the trait to be analyzed will filter only for GWAS entries for this trait*
  
**Example:**
    CHR	POS		SNP		beta_GWAS	pvalue_GWAS	Trait
    14		88500045	rs17123978	0.000479064	0.227412	eGFR
    14		88500645	rs60623686	0.000118531	0.256264	eGFR

  

#### eQTL.df
`eQTL.df` takes as input a data frame of standard eQTL data (as might be downloaded from the GTEx Portal) that must contain the following 5 columns:

`SNP` Variant ID (such as dbSNP ID "rs...". *Note: Must be the same naming scheme as used in the `GWAS.df` to ensure proper matching)*  
`Gene` Gene symbol/name to which the eQTL expression data refers *Note: gene symbol/name must match entries in Genes.df to ensure proper matching*  
`pvalue_eQTL` p-value for the SNP from eQTL analysis (such as one might download from the GTEx Portal)  
`NES_eQTL` NES (normalized effect size) for the SNP from eQTL analysis (Per GTEx: defined as the slope of the linear regression, and is computed as the effect of the alternative allele (ALT) relative to the reference allele (REF) in the human genome reference (i.e., the eQTL effect allele is the ALT allele). NES are computed in a normalized space where magnitude has no direct biological interpretation.)  
`Tissue` tissue type in which the eQTL p-value/beta were obtained *Note: `eQTL.df` can contain multiple tissue types. Specifying the tissue type to be analyzed will filter only for eQTL entires for this tissue type. Alternatively, setting tissue type to "all" (or leaving out the tissue type argument) will automatically pick the smallest eQTL p-value for each SNP across all tissues for a PanTissue analysis*
  
**Example:**
    Gene	SNP		pvalue_eQTL	NES_eQTL	Tissue
    PTPN21	rs147470573	0.0479617	-0.156095	Adipose_Subcutaneous
    PTPN21	rs60688436	0.0479617	-0.156095	Adipose_Subcutaneous



#### Genes.df
`Genes.df` takes an input a data frame of chromosomal locations for the genes to be analyzed, and, if supplied, must contain the following 5 columns. Note: eQTpLot automatically loads a default Genes.df dataframe containing information for thousands of genes (identified by gene symbol) in both genomic builds hg19 and hg38, but you may wish to specify our own Genes.df dataframe if your gene of interest is not included in the default dataframe, or if your eQTL data uses a different gene naming scheme (for example, Gencode ID instead of gene symbol)

`Gene` Gene symbol/name for which the Coordinate data (below) refers to *Note: gene symbol/name must match entries in eQTL.df to ensure proper matching)  
`CHR` Chromosome the gene is on (X coded numerically as 23).  
`Start` Chromosomal coordinate of start position (in basepairs) to use for gene *Note: this should be the smaller of the two values between Start and Stop*  
`Stop` Chromosomal coordinate of end position (in basepairs) to use for gene *Note: this should be the larger of the two values between Start and Stop*  
`Build` The genome build (either hg19 or hg38) for the location data -- the default `Genes.df` dataframe contains entries for both genome builds for each gene, and the script will select the appropriate entry based on the specified gbuild (default is hg19)).

**Example:**
    Gene		CHR	Start		Stop		Build
    EML5		14	89078491	89259096	hg19
    KCNK10	14	88646451	88793256	hg19
    


## Function arguments
To run `eQTpLot`, a number of arguments must be specified. A number of optional arguments are available as well to customize and adjust the resultant plots.

*Required Arguments*

The following arguments must be specified to run eQTpLot:   
`GWAS.df` The name of GWAS.df, in quotes, defined as above  
`eQTL.df` The name of eQTL.df, in quotes, defined as above  
`gene` The name/symbol of gene to analyze, in quotes (must be present in both `Genes.df` and `eQTL.df`)  
`trait` The name of trait to analyze, in quotes (must be present in `GWAS.df`)  
`sigpvalue_GWAS` The GWAS p-value significance threshold to use (this value will be used for a horizontal line in plot A, and to define GWAS significant/non-significant SNPs for plot C)


*Optional Arguments*

The following arguments have default settings, which may be overridden to customize the resulting eQTpLot graphs:  
`Genes.df` The name of `Genes.df`, in quotes, defined as above. The default `Genes.df` contains chromosomal coordinates for most genes (identified by gene symbol) for both genome builds hg19 and hg38.  
`sigpvalue_eQTL`  The eQTL p-value significance threshold to use (eQTL data with a p-value larger than this threshold will be excluded from the analysis). The default value is 0.05.  
`tissue` The tissue, in quotes, from which to derive the eQTL information (the specified tissue must be present in `eQTL.df`). The default setting is "all" to run a Pan-Tissue analysis.  
`range` The range, in kB, to extend the analysis window on either side of the gene boundry of the gene of interest. The default value is 200 kb.  
`gbuild` The genome build used for SNP cooridnates specified in GWAS.df, in quotes. This information will be used to obtain the appropriate gene coordinates from.  
`Gene.df`, and will be used to fetch the genomic information for panel B. The default setting is `"hg19"` but can be changed to `"hg38"` if needed.   
`NESeQTLRange` the maximum and minimum limits (in the format c(0,2), for example (without quotes)) to display for the NES value in  `eQTL.df`. The default setting will adjust the size limits automatically for your data, whereas specifying the limits can keep them consistent between plots.   
`ylima` Used to manually adjust the y axis limit in plot A, if needed.  
`ylimd` Used to manually adjust the y axis limit in plot D, if needed.  
`xlimd` Used to manually adjust the x axis limit in plot D, if needed.  
`genometrackheight` used to set the height of the genome track panel (B), with default setting of 1.5. Gene-dense regions may require more plotting space, whereas gene-poor regions may look better with less plotting space   
`res` resolution of the output plot image (default is 300 dpi)  
`hgt` height of the output plot image (default is 12 inches)  
`wi` width of the output plot image (default is 14 inches)  
`getplot` default is TRUE. If set to false, script will not dsiplay the generated plot it in the viewer  
`saveplot` default is FALSE. If set to true, script will save the generated plot with the name `"gene.trait.tissue.eQTL.png"` using the supplied variables


## Notes on analysis:
A PanTissue analysis can be specified, combining data across all tissue types for each variant. This is the default analysis, and can be overridden by specifying the tissue to use with the tissue argument.

Not that, for all analyses, variants are divided into two groups – congruous (those with the same direction of effect on gene expression and the GWAS trait (e.g., both negative)) in blue, and incongruous (those with opposite directions of effect on gene expression and the GWAS trait (e.g., one negative, one positive)), in red. The division between congruous and incongruous variants provides a more nuanced view of the relationship between gene expression level and GWAS associations – a variant associated with increased expression of a candidate gene and increased risk for a given GWAS trait would seem to be operating through different mechanisms that a variant that is similarly associated with increased expression of the same candidate gene, but a decreased risk for the same GWAS trait. 

## Examples:
Using the supplied example data frames GWAS.df.example and eQTL.df.example, we can generate a PanTissue eQTL analysis for the gene SPATA7 and the trait eGFR, using a GWAS significance threshold of 5e-8, and an eQTL significance threshold of 5e-4, as follows:

    library(eQTpLot)
    data(GWAS.df.example)
    data(eQTL.df.example)
    eQTpLot(GWAS.df = GWAS.df.example, eQTL.df = eQTL.df.example, 
    gene = "SPATA7", trait = "eGFR", sigpvalue_GWAS = 5e-8, 
    sigpvalue_eQTL = 5e-4)

Which generates the following plot:

PLOT1

We can carry out the same analysis, confining out analysis to only kidney specific eQTLs as follows:

     eQTpLot(GWAS.df = GWAS.df.example, eQTL.df = eQTL.df.example, 
     gene = "SPATA7", trait = "eGFR", sigpvalue_GWAS = 5e-8, 
     sigpvalue_eQTL = 5e-4, tissue = "Artery_Coronary")

Which generates the following plot:

PLOT2

Switching our analysis to a different gene, PTPN21, we can compare results between our SPATA7 analysis and PTPN21:
   
    eQTpLot(GWAS.df = GWAS.df.example, eQTL.df = eQTL.df.example, 
    gene = "PTPN21", trait = "eGFR", sigpvalue_GWAS = 5e-8, 
    sigpvalue_eQTL = 5e-4)

Which generates the following plot:

PLOT3

Multiple additional modifications to the plots can be specified, as noted above.
