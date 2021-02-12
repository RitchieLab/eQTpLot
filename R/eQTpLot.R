#' An eQTpLot function
#'
#' This function creates an eQTpLot composite
#' @param GWAS.df Dataframe, one row per SNP, as one might obtain from a PLINK association analysis, with columns :
#'     CHR          Chromosome for SNP (X coded numerically as 23)
#'     BP           Chromosomal position for each SNP, in base pairs
#'     SNP          Variant ID (such ws dbSNP ID "rs...". Note: Must be the same naming scheme as  used in the eQTL.df to ensure proper matching)
#'     P            P value for the SNP from GWAS analysis
#'     BETA         Beta for SNP from GWAS analysis
#'     PHE          (Optional) Name of the phenotype for which the GWAS was run. This column is optional, and is useful if your GWAS.df contains data for multiple phenotypes,
#'                  such as one might obtain from a PheWAS). The "trait" parameter is subsequently used to filter in the GWAS.df entries for only this phenotype.
#'                  If GWAS.df does not contain a "PHE" column, eQTpLot will assume all the supplied GWAS data is for the phenotype specified by the "trait" parameter
#' @param eQTL.df Dataframe, one row per SNP, as one might obtain from the GTEx Portal, with columns:
#'     SNP.Id       Variant ID (such ws dbSNP ID "rs...". Note: Must be the same naming scheme as used in the GWAS.df to ensure proper matching)
#'     Gene.Symbol  Gene symbol/name to which the eQTL expression data refers (Note: gene symbol/name must match entries in Genes.df to ensure proper matching)
#'     P.value      pvalue for the SNP from eQTL analysis (such as one might download from the GTEx Portal)
#'     NES          NES (normalized effect size) for the SNP from eQTL analysis
#'                  (Per GTEx: defined as the slope of the linear regression, and is computed as the effect of the alternative allele (ALT) relative to the reference allele (REF) in the human genome reference
#'                  (i.e., the eQTL effect allele is the ALT allele). NES are computed in a normalized space where magnitude has no direct biological interpretation.)
#'     Tissue       Tissue type to which the eQTL pvalue/beta refer
#'                  (Note: eQTL.df can contain multiple tissue types. Specifying the tissue type to be analyzed will filter only for eQTL entires for this tissue type.
#'                  Alternatively, setting tissue type to "all" (the default setting) will automatically pick the smallest eQTL pvalue for each SNP across all tissues for a PanTissue analysis)
#' @param Genes.df Dataframe, one row per gene, with the following columns
#'         (Note: eQTpLot automatically loads a default Genes.df database containing information for both genomic builds hg19 and hg38,
#'         but you may wish to specify our own Genes.df dataframe if your gene of interest is not included in the default dataframe, or if your eQTL data uses a different gene naming scheme
#'         (for example, Gencode ID instead of gene symbol)):
#'     Gene         Gene symbol/name for which the Coordinate data (below) refers to
#'                  (Note: gene symbol/name must match entries in eQTL.df to ensure proper matching)
#'     CHR          Chromosome the gene is on (X coded numerically as 23)
#'     Start        Chromosomal coordinate of start position (in basepairs) to use for gene (Note: this should be the smaller of the two values between Start and Stop)
#'     Stop         Chromosomal coordinate of end position (in basepairs) to use for gene (Note: this should be the larger of the two values between Start and Stop)
#'     Build        The genome build (either hg19 or hg38) for the location data -
#'                 the default Genes.df dataframe contains entries for both genome builds for each gene, and the script will select the appropriate entry based on the specified gbuild (default is hg19)).
#' @param LD.df Optional dataframe of SNP linkage data, one row per SNP pair, with columns as one might obtain from a PLINK linkage analysis using the PLINK --r2 option:
#'         (If no LD.df is supplied, eQTpLot will plot data without LD information)
#'     BP_A         Basepair position of the first SNP in the LD pair
#'     SNP_A        Variant name of the first SNP in the LD pair (Not: only SNPs that also appear in the GWAS.df SNP column will be used for LD analysis)
#'     BP_B         Basepair position of the second SNP in the LD pair
#'     SNP_B        Variant name of the second SNP in the LD pair (Note: only SNPs that also appear in the GWAS.df SNP column will be used for LD analysis)
#'     R2           Squared correlation measure of linkage between the two SNPs in the pair
#' @param gene name/symbol of gene to analyze, in quotes
#' @param trait name of GWAS phenotype to analyze, in quotes. If all the data in GWAS.df is for a single phenotype and no PHE column is present, this parameter will be used as the name for the analyzed phenotype. 
#'        If GWAS.df contains information on multiple phenotypes, as specified in the optional GWAS.df PHE column, this parameter will be used to filter in GWAS.df entries for only this phenotype.
#' @param sigpvalue_GWAS GWAS pvalue significance threshold to use (this value will be used for a horizontal line in plot A, and to define GWAS significant/non-significant SNPs for plot C)
#' @param sigpvalue_eQTLeQTL pvalue significance threshold to use (eQTL data with a pvalue larger than this threshold will be excluded from the analysis)
#' @param tissue tissue type, in quotes, for eQTL data to use (from eQTL.df, default setting is "all" for Pan-Tissue analysis)
#' @param range range, in kB, to extend analysis window on either side of gene boundry. Default is 200 kb
#' @param NESeQTLRange the maximum and minimum limits in the format c(min,max), to display for the NES value in eQTL.df.
#'          The default setting will adjust the size limits automatically for your data, whereas specifying the limits
#'          can keep them consistent between plots.
#' @param congruence if set to TRUE, variants with congruent and incongruent effects will be plottes separately. Default is FALSE.         
#' @param R2min R^2 values less than this value in the the LD.df (is supplied) will not be displayed. Default is 0.1
#' @param LDmin for the LDheatmap panel, only variants that are in LD (with R^2 > R2min) with at least this many other variants will be displayed. The default is 10
#' @param LDcolor the default color scheme for the LDheatmap is in the viridis color palate. This option can be set to "black" to plot the LDheatmap in grayscale
#' @param leadSNP if LD.df data is included, this parameter is used to specify the lead SNP to use for plotting LD information in the P-P plots. The SNP ID must be present in both the GWAS.df and LD.df dataframes.
#' @param ylima used to manually adjust the y axis limit in plot A, if needed
#' @param ylimd used to manually adjust the y axis limit for the P-P plot, if needed
#' @param xlimd used to manually adjust the x axis limit for the P-P plot, if needed
#' @param genometrackheight used to set the height of the genome track panel (B), with default setting of 1.5.
#'           Gene-dense regions may require more plotting space, whereas gene-poor regions may look better with less plotting space.
#' @param gbuild the genome build to use, in quotes, for fetching genomic information for panel B.
#'            Default is "hg19" but can specify "hg38" if needed. This build should match the genome build used for "CHR" and "BP" in the GWAS.df
#' @param res resolution of the output plot image (default is 300 dpi)
#' @param wi the width of the output plot image (the height is calculated from this to maintain the correct aspect ratio)
#' @param getplot default is TRUE. If set to false, script will not dsiplay the generated plot in the viewer
#' @param saveplot default is TRUE, script will save the generated plot in the working diretory with the name
#'            "gene.trait.tissue.eQTL.png" using the supplied variables
#' @importFrom magrittr "%>%"
#' @export
#' @examples
#' Saves to current directory
#' eQTpLot(Genes.df = Genes.df.example, GWAS.df = GWAS.df.example,
#'         eQTL.df = eQTL.df.example, gene = "SPATA7", trait = "eGFR",
#'         getplot=FALSE)
#' eQTpLot()

eQTpLot <- function(GWAS.df, eQTL.df, Genes.df, LD.df = TRUE, gene, trait,
                    sigpvalue_GWAS = 5e-8, sigpvalue_eQTL = 0.05,
                    tissue = "all", range = 200, NESeQTLRange = c(NA,NA), 
                    congruence = FALSE, R2min = 0.2, LDmin = 10, leadSNP = TRUE,
                    LDcolor = "color", ylima = NA, ylimd = NA, xlimd = NA,
                    genometrackheight = 1.5, gbuild = "hg19",
                    res = 300, wi = "wi",
                    getplot = TRUE, saveplot = TRUE){
  
  
  
  ########################
  ### Required Packages
  #require(ggnewscale)
  #require(patchwork)
  #require(Gviz)
  #require(GenomicRanges)
  #require(biomaRt)
  #require(ggpubr)
  #require(LDheatmap)
  #require(ggplotify)
  
  
  
  ########################
  ### Check Data
  print("Checking input data...")
  
  if(missing(Genes.df)){Genes.df <- eQTpLot:::genes}
  
  if(all(c("CHR", "BP", "SNP", "BETA", "P") %in% colnames(GWAS.df))==FALSE) {
    stop("The data supplied to GWAS.df must contain columns 'CHR', 'BP', 'SNP', 'BETA', and 'P'")
  }
  
  if("PHE" %in% colnames(GWAS.df)){
    print(paste(sep = "", 'PHE column found in GWAS.df. Analyzing data for phenotype ', trait))
  } else {print(paste(sep = "", 'PHE column not found in GWAS.df. Assuming all data in GWAS.df is for phenotype ', trait))
    GWAS.df$PHE <- trait
  }
  
  if((trait %in% GWAS.df$PHE) == FALSE) {
    stop('Sorry, the  phenotype ', paste(trait), ' does not exist in the PHE column of the GWAS.df dataframe. Phenotypes included in the data supplied to GWAS.df are:\n',
         paste("'",as.character(unique(GWAS.df$PHE)),"'",collapse=", ",sep=""))
  }
  
  if(all(c("SNP.Id", "Gene.Symbol", "P.Value", "NES", "Tissue") %in% colnames(eQTL.df))==FALSE) {
    stop("The data supplied to eQTL.df must contain columns 'SNP.Id', 'Gene.Symbol', 'P.Value', 'NES', and 'Tissue'")
  }
  
  if(all(c("Gene", "CHR", "Start", "Stop", "Build") %in% colnames(Genes.df))==FALSE) {
    stop("The data supplied to Genes.df must contain columns 'Gene', 'CHR', 'Start', 'Stop', 'Build'")
  }
  
  if(is.numeric(eQTL.df$P.Value) == FALSE | is.numeric(eQTL.df$NES) == FALSE) {
    stop('Sorry, the eQTL.df dataframe must contain only numeric data for P.Value and NES')
  }
  
  if(is.numeric(GWAS.df$P) == FALSE | is.numeric(GWAS.df$BETA) == FALSE | is.integer(GWAS.df$BP) == FALSE | is.integer(GWAS.df$CHR) == FALSE) {
    stop('Sorry, the GWAS.df dataframe must contain only numeric data for CHR, BP, P, and BETA (Note: chromosomes must be coded numerically)')
  }
  
  if(is.integer(Genes.df$CHR) == FALSE | is.integer(Genes.df$Start) == FALSE | is.integer(Genes.df$Stop) == FALSE) {
    stop('Sorry, the Genes.df dataframe must contain only integer values for CHR, Start, and Stop (Note: chromosomes must be coded nuemrically)')
  }
  
  if(dim(Genes.df[which(Genes.df$Gene == gene & Genes.df$Build == gbuild), ])[1] == 0) {
    stop('Sorry, there is no information for the gene ', paste(gene), ' in the Genes.df dataframe for the genomic build ', paste(gbuild), '\nConsider supplying your own Genes.df input file with information on this gene.')
  }
  
  if(all(tissue == "all") == FALSE & (all(tissue %in% eQTL.df$Tissue)) == FALSE){
    stop('Sorry, at least one of the specified tissues, ', paste(tissue, " ", sep = "," ), ' does not exist in the eQTL.df dataframe. Tissues included in the data supplied to eQTL.df are:\n',
         paste("'",as.character(unique(eQTL.df$Tissue)),"'",collapse=", ",sep=""))
  }
  
  if(length(tissue) >= 2){
    eQTL.df <- eQTL.df[which(eQTL.df$Tissue %in% tissue), ]
  }
  
  print(paste(sep = "", 'eQTL analysis will be completed for tissues ',  paste("'",as.character(unique(tissue)),"'",collapse=", ",sep=""), ' and for gene ', paste(gene)))
  
  if(LDcolor != "color" & LDcolor != "black"){
    stop('Sorry, the argument LDcolor must be set to either "color" or "black"')
  }
  
  
  if(isTRUE(LD.df) == FALSE){
    if(all(c("BP_A", "SNP_A", "BP_B", "SNP_B", "R2") %in% colnames(LD.df))==FALSE) {
    stop("The data supplied to LD.df must contain columns 'BP_A', 'SNP_A', 'BP_B', 'SNP_B', and 'R2'")
    }
  
    if((is.integer(LD.df$BP_A) == FALSE | is.integer(LD.df$BP_B) == FALSE)) {
    stop('Sorry, the LD.df dataframe must contain only integer values for BP_A and BP_B')
    }
  
    if(is.numeric(LD.df$R2) == FALSE){
    stop('Sorry, the LD.df dataframe must contain only numeric values for R2')
    }
    
    if(isTRUE(leadSNP) == FALSE & (leadSNP %in% LD.df$SNP_A) == FALSE & (leadSNP %in% LD.df$SNP_B) == FALSE){
    stop('Sorry, the specified leadSNP is not present in your LD.df')
    }
    
    if(isTRUE(leadSNP) == FALSE & (leadSNP %in% GWAS.df$SNP) == FALSE){
      stop('Sorry, the specified leadSNP is not present in your GWAS.df')
    }
    
    if(isTRUE(leadSNP) == FALSE & (leadSNP %in% eQTL.df$SNP.Id) == FALSE){
      stop('Sorry, the specified leadSNP is not present in your eQTL.df')
    }
  }
  

  ########################
  ### Compile Data
  print("Compiling GWAS and eQTL data...")
  
  ### Set genomic ranges for data selection
  rangebp <- range*1000
  startpos <- (min(Genes.df[which(Genes.df$Gene == gene & Genes.df$Build == gbuild), ] %>% dplyr::select(Start)) - rangebp)
  stoppos <- (max(Genes.df[which(Genes.df$Gene == gene & Genes.df$Build == gbuild), ] %>% dplyr::select(Stop)) + rangebp)
  chromosome <- Genes.df$CHR[Genes.df$Gene == gene & Genes.df$Build == gbuild]
  
  ### Subset GWAS.df for gene of interest, check to make sure data is present
  gwas.data <- GWAS.df[which(GWAS.df$CHR == chromosome & GWAS.df$PHE == trait & GWAS.df$BP >= startpos & GWAS.df$BP <= stoppos & !(is.na(GWAS.df$P)) & !(is.na(GWAS.df$BETA))), ]
  if(dim(gwas.data)[1] == 0) stop('Sorry, GWAS.df does not contain data for any SNPs in the ', paste(range), 'kb flanking the gene ', paste(gene), ' for the trait ', paste(trait))
  if(dim(gwas.data[which(gwas.data$P <= sigpvalue_GWAS), ])[1] == 0) {
    NoFisher<-TRUE; print(paste(sep = "", 'WARNING: GWAS.df does not contain any SNPs with p-value < sigpvalue_GWAS within the ',
                                range, 'kb flanking the gene ', gene, ' for the trait ',
                                trait, '. eQTL Enrcihment Plot statistics will not be calculated'))
  } else {NoFisher <- FALSE}
  
  ### Subset eQTL.df for tissue and gene of interest, check to make sure data is present
  eqtl.data <- eQTL.df[which(eQTL.df$Gene.Symbol == gene & eQTL.df$P.Value <= sigpvalue_eQTL & !(is.na(eQTL.df$NES))), ]
  if(dim(eqtl.data)[1] == 0) stop('Sorry, eQTL.df does not have any data for the gene ', paste(gene), ' meeting your sigpvalue_eQTL threshold')
  if(any(tissue == "all") == TRUE | (length(tissue) >= 2) == TRUE){PanTissue <- TRUE}
  if(any(tissue == "all") == FALSE & (length(tissue) >= 2) == FALSE){PanTissue <- FALSE}
  if(PanTissue == TRUE) {
    eqtl.data <- eqtl.data[which(!(is.na(eqtl.data$SNP.Id))), ] %>% dplyr::group_by(SNP.Id) %>% dplyr::slice(which.min(P.Value));
    eqtl.data$Tissue <- "PanTissue"
  } else {
    eqtl.data <- eqtl.data[which(eqtl.data$Tissue == tissue), ]
  }
  if(dim(eqtl.data)[1] == 0) stop('Sorry, there are no eQTLs for the tissue', paste(tissue), ' with a p-value < sigeQTL')
  eqtl.data <- dplyr::ungroup(eqtl.data)
  
  ### Join GWAS and eQTL data, check to make sure there is at least some overlap in SNPs between the two datasets
  gwas.data$SNP <- as.factor(gwas.data$SNP)
  eqtl.data$SNP.Id <- as.factor(eqtl.data$SNP.Id)
  combinedSNPS <- sort(union(levels(gwas.data$SNP), levels(eqtl.data$SNP.Id)))
  Combined.eQTL.GWAS.Data <- dplyr::left_join(dplyr::mutate(gwas.data, SNP=factor(SNP, levels=combinedSNPS)),
                                              dplyr::mutate(eqtl.data, SNP.Id=factor(SNP.Id, levels=combinedSNPS)) %>% dplyr::rename(SNP = SNP.Id))
  if(dim(Combined.eQTL.GWAS.Data)[1] == 0) {
    stop('Sorry, for the gene ', paste(gene), ' and the trait ', paste(trait), ' there is no overlap between the SNPs in your GWAS.df and eQTL.df')
  }
  
  ### Determine directions of effect and congruence
  Combined.eQTL.GWAS.Data$DirectionOfEffect_GWAS <- ifelse(Combined.eQTL.GWAS.Data$BETA < 0, "Negative", ifelse(Combined.eQTL.GWAS.Data$BETA > 0, "Positive", NA))
  Combined.eQTL.GWAS.Data$DirectionOfEffect_eQTL <- ifelse(Combined.eQTL.GWAS.Data$NES < 0, "DOWN", ifelse(Combined.eQTL.GWAS.Data$NES > 0, "UP", NA))
  Combined.eQTL.GWAS.Data$Congruence <- (Combined.eQTL.GWAS.Data$BETA*Combined.eQTL.GWAS.Data$NES)
  Combined.eQTL.GWAS.Data$Congruence <- ifelse(Combined.eQTL.GWAS.Data$Congruence < 0, "Incongruent", ifelse(Combined.eQTL.GWAS.Data$Congruence > 0, "Congruent", NA))
  
  if(congruence == FALSE) {
    Combined.eQTL.GWAS.Data$Congruence <- ifelse(Combined.eQTL.GWAS.Data$Congruence == "Incongruent", "Congruent", ifelse(Combined.eQTL.GWAS.Data$Congruence == "Congruent", "Congruent", NA))
  }
  
  ### Build final dataframe
  Combined.eQTL.GWAS.Data$NeglogeQTLpValue <- -(log10(Combined.eQTL.GWAS.Data$P.Value))
  Combined.eQTL.GWAS.Data$Neglog10pvalue_GWAS <- -(log10(Combined.eQTL.GWAS.Data$P))
  Combined.eQTL.GWAS.Data <- Combined.eQTL.GWAS.Data[which(!(is.na(Combined.eQTL.GWAS.Data$P))), ]
  Combined.eQTL.GWAS.Data$significance <- ifelse(Combined.eQTL.GWAS.Data$P >= sigpvalue_GWAS, "Non-significant", "Significant")
  Combined.eQTL.GWAS.Data$Congruence[is.na(Combined.eQTL.GWAS.Data$Congruence)] <- "Non-eQTL"
  Combined.eQTL.GWAS.Data <- Combined.eQTL.GWAS.Data %>% dplyr::mutate(Congruence=factor(Congruence, levels=c("Non-eQTL", "Congruent", "Incongruent"), ordered=TRUE))
  if(dim(Combined.eQTL.GWAS.Data[ which( !is.na(Combined.eQTL.GWAS.Data$NES) & Combined.eQTL.GWAS.Data$Congruence == "Congruent") , ])[1] == 0){
    Congruentdata <- FALSE} else {Congruentdata <- TRUE
    }
  if(dim(Combined.eQTL.GWAS.Data[ which( !is.na(Combined.eQTL.GWAS.Data$NES) & Combined.eQTL.GWAS.Data$Congruence == "Incongruent") , ])[1] == 0){
    Incongruentdata <- FALSE} else {Incongruentdata <- TRUE
    }
  
  ########################
  ### LD Calculations
  if(isTRUE(LD.df) == FALSE){
    print("Compiling LD information...")
    Combined.eQTL.GWAS.Data$pvaluemult <- Combined.eQTL.GWAS.Data$P.Value*Combined.eQTL.GWAS.Data$P
    LD.df <- LD.df %>% dplyr::filter(SNP_A %in% levels(Combined.eQTL.GWAS.Data$SNP))
    LD.df <- LD.df %>% dplyr::filter(SNP_B %in% levels(Combined.eQTL.GWAS.Data$SNP))
    
    ### Select lead SNP for congruent data
    if(Congruentdata == TRUE){
      mostsigsnp.cong <-  as.character(Combined.eQTL.GWAS.Data %>% dplyr::filter(!is.na(pvaluemult)) %>% dplyr::filter(Congruence == "Congruent") %>% dplyr::filter(pvaluemult == min(pvaluemult)) %>% dplyr::pull(SNP))
      if(isTRUE(leadSNP) == FALSE){
        if((as.character(Combined.eQTL.GWAS.Data%>% dplyr::filter(SNP == leadSNP) %>% dplyr::pull(Congruence)) == "Congruent") == TRUE){
        mostsigsnp.cong <- leadSNP
      }}
      Combined.eQTL.GWAS.Data <- dplyr::left_join(Combined.eQTL.GWAS.Data, LD.df %>% dplyr::filter(SNP_A == mostsigsnp.cong) %>% dplyr::select(c("SNP_B","R2")), by = c("SNP" = "SNP_B"))
      names(Combined.eQTL.GWAS.Data)[names(Combined.eQTL.GWAS.Data) == "R2"] <- "R2cong"
      names(Combined.eQTL.GWAS.Data)[names(Combined.eQTL.GWAS.Data) == "SNP_B"] <- "SNP_Bcong"
      Combined.eQTL.GWAS.Data.cong <- subset(Combined.eQTL.GWAS.Data, SNP == mostsigsnp.cong)
    }
    
    ### Select lead SNP for incongruent data    
    if(Incongruentdata == TRUE){
      mostsigsnp.incong <-  as.character(Combined.eQTL.GWAS.Data %>% dplyr::filter(!is.na(pvaluemult)) %>% dplyr::filter(Congruence == "Incongruent") %>% dplyr::filter(pvaluemult == min(pvaluemult)) %>% dplyr::pull(SNP))
      if(isTRUE(leadSNP) == FALSE){
        if((as.character(Combined.eQTL.GWAS.Data%>% dplyr::filter(SNP == leadSNP) %>% dplyr::pull(Congruence)) == "Incongruent") == TRUE){
        mostsigsnp.incong <- leadSNP
      }}
      Combined.eQTL.GWAS.Data <- dplyr::left_join(Combined.eQTL.GWAS.Data, LD.df %>% dplyr::filter(SNP_A == mostsigsnp.incong) %>% dplyr::select(c("SNP_B","R2")), by = c("SNP" = "SNP_B"))
      names(Combined.eQTL.GWAS.Data)[names(Combined.eQTL.GWAS.Data) == "R2"] <- "R2incong"
      names(Combined.eQTL.GWAS.Data)[names(Combined.eQTL.GWAS.Data) == "SNP_B"] <- "SNP_Bincong"
      Combined.eQTL.GWAS.Data.incong <- subset(Combined.eQTL.GWAS.Data, SNP == mostsigsnp.incong)
    }
    
    ### Filter and combine LD data
    LD.df$R2[LD.df$R2 <= R2min] = NA 
    LD.df.matrix <- tidyr::spread(LD.df[(!duplicated(LD.df[,c("SNP_B","SNP_A")])),c("SNP_A","SNP_B","R2")], SNP_A, R2)
    rownames(LD.df.matrix) <- LD.df.matrix$SNP_B
    LD.df.matrix$SNP_B <- NULL
    colnames(LD.df.matrix[, which(colSums(!is.na(LD.df.matrix)) > LDmin )]) -> SNPsWithLDData1
    rownames(LD.df.matrix[which(colSums(!is.na(LD.df.matrix)) > LDmin ), ]) -> SNPsWithLDData2
    union(SNPsWithLDData1, SNPsWithLDData2) -> SNPsWithLDData
    if(length(SNPsWithLDData) < 2){
      stop('Sorry, after filtering the LD.df data by the supplied R2min and LDmin thresholds, fewer than 2 SNPs remain that are also present in GWAS.df')
    }
    
    ### Generate square LD matrix    
    LD.df.matrix[rownames(LD.df.matrix) %in% SNPsWithLDData, colnames(LD.df.matrix) %in% SNPsWithLDData] -> LD.df.matrix
    LD.df.matrix$startpos <-  NA
    LD.df.matrix$stoppos <- NA
    dat2 <- data.frame(matrix(nrow = 2, ncol = ncol(LD.df.matrix)))
    rownames(dat2) <- c("startpos", "stoppos")
    colnames(dat2) <- colnames(LD.df.matrix)
    dplyr::bind_rows(LD.df.matrix, dat2) -> LD.df.matrix
    LD.df.matrix[,c("startpos", "stoppos")] <- NA
    matrix <- as.matrix(LD.df.matrix)
    un1 <- unique(sort(c(colnames(LD.df.matrix), rownames(LD.df.matrix))))
    matrix2 <- matrix(NA, length(un1), length(un1), dimnames = list(un1, un1))
    matrix2[row.names(matrix), colnames(matrix)] <- matrix
    matrix <- t(matrix2)
    LD.df.matrix <- dplyr::coalesce(as.data.frame(matrix), as.data.frame(matrix2))
    LD.df.matrix[is.na(LD.df.matrix)] = 0
    SNPPositions <- unique(dplyr::bind_rows(unique(LD.df[which(LD.df$SNP_A %in% colnames(LD.df.matrix)), c('SNP_A', 'BP_A')]), unique(LD.df[which(LD.df$SNP_B %in% colnames(LD.df.matrix)), c('SNP_B', 'BP_B')] %>% dplyr::rename(SNP_A = 1, BP_A = 2))))
    SNPPositions[order(SNPPositions$BP_A),] -> SNPPositions
    rownames(SNPPositions) <- SNPPositions$SNP_A
    SNPPositions[order(SNPPositions$BP_A),]$SNP_A -> SNPorder
    SNPPositions[order(SNPPositions$BP_A),]$BP_A -> positions
    LD.df.matrix[SNPorder, SNPorder] -> LD.df.matrix
  }
  
  ########################
  ### Generate main plot
  print("Generating main plot...")
  
  ### Set plot limits
  if(is.na(ylima) == TRUE){ylima <- (max(Combined.eQTL.GWAS.Data %>% dplyr::select (Neglog10pvalue_GWAS), na.rm = TRUE) + 1)}
  minpos <- min(Combined.eQTL.GWAS.Data$BP, na.rm = TRUE)
  maxpos <- max(Combined.eQTL.GWAS.Data$BP, na.rm = TRUE)
  
  ### Generate plot 1
  p1 <-
    ggplot2::ggplot (data=Combined.eQTL.GWAS.Data, aes(stroke = 0)) +
    ggplot2::coord_cartesian(xlim = c(minpos, maxpos), expand = FALSE) +
    ggplot2::geom_point(data = subset(Combined.eQTL.GWAS.Data, Congruence == "Non-eQTL"), shape = 15, color = "black", alpha = 0.2,
                        aes(x=BP, y=Neglog10pvalue_GWAS)) +
    ggplot2::xlab("") +
    ggplot2::ylab(paste("-log10(GWAS p-value)  association with ", trait)) +
    ggplot2::scale_y_continuous(limits=c(NA,ylima)) +
    ggplot2::ggtitle(paste("Association study of ", trait, ", variants colored by eQTL data for ", gene,
                           "\n(GWAS significance threshold ", sigpvalue_GWAS, ", eQTL significance threshold ",
                           sigpvalue_eQTL, ")", sep = "")) +
    ggplot2::scale_shape_manual("GWAS Direction of Effect", values=c("Negative" = 25, "Positive" = 24), na.value = 22) +
    ggplot2::guides(alpha = FALSE,
                    size = guide_legend("eQTL Normalized Effect Size",
                                        override.aes = list(shape = 24, color = "black", fill ="grey"),
                                        title.position = "top", order = 2),
                    shape = guide_legend(title.position = "top",
                                         order = 1,
                                         override.aes = list(size= 3,
                                                             fill = "grey"))) +
    ggplot2::theme(legend.direction = "horizontal", legend.key = element_rect(fill = NA, colour = NA, size = 0.25)) +
    ggplot2::geom_hline(yintercept=-log10(sigpvalue_GWAS), linetype="solid", color = "red", size=0.5) +
    ggplot2::theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
    ggplot2::theme(plot.margin = unit(c(0,1,-0.8,0), "cm"))
  
  ### Add lead SNP data, if LD.df is supplied
  if(isTRUE(LD.df) == FALSE & Congruentdata == TRUE ){
    p1 <- p1 + 
      ggplot2::geom_point(data=Combined.eQTL.GWAS.Data.cong, 
                          aes(x=BP, y=Neglog10pvalue_GWAS), 
                          shape = 1, size = 4, stroke = 3, color = "#000099")
  }
  
  if(isTRUE(LD.df) == FALSE & Incongruentdata == TRUE ){
    p1 <- p1 + 
      ggplot2::geom_point(data=Combined.eQTL.GWAS.Data.incong, 
                          aes(x=BP, y=Neglog10pvalue_GWAS),
                          shape = 1, size = 4, stroke = 3, color = "#990000")
  }
  
  ### If congruence is TRUE, add congruent and incongruent data
  if(congruence == TRUE){
    if(Congruentdata == TRUE){
      p1 <- p1 + 
        ggplot2::geom_point(data = subset(Combined.eQTL.GWAS.Data, Congruence == "Congruent"),
                            alpha = 1,
                            aes(x=BP, y=Neglog10pvalue_GWAS, fill=NeglogeQTLpValue, alpha = 1,
                                shape = DirectionOfEffect_GWAS, size = abs(NES))) +
        ggplot2::scale_fill_gradient("-log10(eQTL p-value), SNPs with\nCongruous directions of effect",
                                     low="#000099", high="#33FFFF",
                                     guide = guide_colorbar(title.position = "top"),
                                     limits=c(min(Combined.eQTL.GWAS.Data %>%
                                                    dplyr::select (NeglogeQTLpValue), na.rm = TRUE),max(Combined.eQTL.GWAS.Data %>%
                                                                                                          dplyr::select (NeglogeQTLpValue), na.rm = TRUE)))
    }
    if(Congruentdata == TRUE & Incongruentdata == TRUE){
      p1 <- p1 + ggnewscale::new_scale_fill()
    }
    if(Incongruentdata == TRUE){
      p1 <- p1 + 
        ggplot2::geom_point(data = subset(Combined.eQTL.GWAS.Data, Congruence == "Incongruent"),
                            alpha = 1,
                            aes(x=BP, y=Neglog10pvalue_GWAS, fill=NeglogeQTLpValue, alpha = 1,
                                shape = DirectionOfEffect_GWAS, size = abs(NES))) +
        ggplot2::scale_fill_gradient("-log10(eQTL p-value), SNPs with\nIncongruous directions of effect",
                                     low="#990000", high="#FFCC33",
                                     guide = guide_colorbar(title.position = "top"),
                                     limits=c(min(Combined.eQTL.GWAS.Data %>%
                                                    dplyr::select (NeglogeQTLpValue), na.rm = TRUE),max(Combined.eQTL.GWAS.Data %>%
                                                                                                          dplyr::select (NeglogeQTLpValue), na.rm = TRUE)))
    }
  }
  
  ### If congruence is FALSE, add all data
  if(Congruentdata == TRUE & Incongruentdata == FALSE & congruence != TRUE) {
    p1 <- p1 + 
      ggplot2::geom_point(data = subset(Combined.eQTL.GWAS.Data, Congruence == "Congruent"),
                          alpha = 1,
                          aes(x=BP, y=Neglog10pvalue_GWAS, fill=NeglogeQTLpValue, alpha = 1, shape = DirectionOfEffect_GWAS, size = abs(NES))) +
      ggplot2::scale_fill_viridis_c("-log10(eQTL p-value)",
                                    option = "C",
                                    guide = guide_colorbar(title.position = "top"),
                                    limits=c(min(Combined.eQTL.GWAS.Data %>%
                                                   dplyr::select (NeglogeQTLpValue), na.rm = TRUE), max(Combined.eQTL.GWAS.Data %>%
                                                                                                          dplyr::select (NeglogeQTLpValue), na.rm = TRUE)))
  }
  
  ### Add size scale
  p1 <- p1 + ggplot2::scale_size_continuous(limits = NESeQTLRange)
  
  ### Add lead SNP labels, if LD.df is supplied
  if(isTRUE(LD.df) == FALSE & Congruentdata == TRUE){
    p1 <- p1 + ggrepel::geom_text_repel(aes(x=BP, y=Neglog10pvalue_GWAS, label = SNP, fontface = "bold"),
                                        color = "black", size =5, hjust="inward", vjust="inward", data = Combined.eQTL.GWAS.Data.cong)
  }
  
  if(isTRUE(LD.df) == FALSE & Incongruentdata == TRUE){
    p1 <- p1 + ggrepel::geom_text_repel(aes(x=BP, y=Neglog10pvalue_GWAS, label = SNP, fontface = "bold"),
                                        color = "black", size =5, hjust="inward", vjust="inward", data = Combined.eQTL.GWAS.Data.incong)
  }
  
  
  
  ########################
  ### Generate Gene Tracks
  print("Generating gene tracks...")
  
  ### Generate Gene Track Plot
  if(gbuild == "hg19"){hostname <- "https://grch37.ensembl.org"}
  if(gbuild == "hg38"){hostname <- "https://apr2020.archive.ensembl.org"}
  
  bm <- biomaRt::useMart(host = hostname,
                         biomart = "ENSEMBL_MART_ENSEMBL",
                         dataset = "hsapiens_gene_ensembl")
  
  biomTrack <- Gviz::BiomartGeneRegionTrack(genome = gbuild,
                                            chromosome = median(Combined.eQTL.GWAS.Data$CHR, na.rm = TRUE),
                                            start = minpos,
                                            end = maxpos,
                                            filter = list("with_refseq_mrna"=TRUE),
                                            name = "ENSEMBL",
                                            background.panel="gray95",
                                            biomart = bm,
                                            margin = c(-3,-3))
  
  gtrack <- Gviz::GenomeAxisTrack(fontcolor="#000000", fontsize=14, margin = c(-3,-3))
  
  genetracks <- patchwork::wrap_elements(panel = (grid::grid.grabExpr(Gviz::plotTracks(list(biomTrack, gtrack),
                                                                                       collapseTranscripts = "meta",
                                                                                       transcriptAnnotation = "symbol",
                                                                                       chromosome = median(Combined.eQTL.GWAS.Data$CHR, na.rm = TRUE),
                                                                                       from = min(Combined.eQTL.GWAS.Data$BP, na.rm = TRUE),
                                                                                       to= max(Combined.eQTL.GWAS.Data$BP, na.rm = TRUE),
                                                                                       showTitle = FALSE,
                                                                                       labelPos = "below",
                                                                                       distFromAxis = 10,
                                                                                       innermargin = 0,
                                                                                       maxHeight = (genometrackheight*10),
                                                                                       minHeight = (genometrackheight*10),
                                                                                       sizes=c(genometrackheight,1),
                                                                                       margin = c(-3,-3)))))
  
  ########################  
  ### Generate LDHeatMap
  if(isTRUE(LD.df) == FALSE){
    if(length(SNPsWithLDData) > 1000 & interactive()){
    notrun <- askYesNo(default = TRUE,
                       msg = paste(sep = "", length(SNPsWithLDData), ' variants being used to generate LDHeatMap.\nUsing more than 1000 variants can take a long time to run.\nYou can increase the values for LDmin and R2min to use fewer variants.\nDo you want to continue with ', length(SNPsWithLDData),  ' variants?'))
  } else {
    notrun <- TRUE }
  
  if(notrun == FALSE){
    opt <- options(show.error.messages = FALSE)
    on.exit(options(opt))
    print("Stopping analysis")
    stop()
  } else {
    if(notrun == TRUE){
      print(paste(sep = "", 'Generating LDHeatMap with ', length(SNPsWithLDData), ' variants'))
    }
    
  if(LDcolor == "color"){
    colorscale <-c(viridisLite::viridis(30, option = "C", direction = -1), "white")
    } else {if(LDcolor == "black"){
      colorscale <- c("grey10", "grey20", "grey30", "grey40","grey50", "grey60", "grey70", "grey80", "grey90","grey100", "white")
    }
      }
    
    
  LDheatmap::LDheatmap(as.matrix(LD.df.matrix), genetic.distances = positions, color = colorscale, flip = TRUE, add.map= TRUE, title = "", geneMapLabelX = NA, geneMapLabelY = NA, newpage = FALSE) -> LDmap
  }
  ###Update heatmap graphics
    LDmap$LDheatmapGrob$children$heatMap$children$heatmap$vp$width <- unit(0.7, "snpc")
    LDmap$LDheatmapGrob$children$heatMap$children$heatmap$vp$height <- unit(0.7, "snpc")
    LDmap$LDheatmapGrob$children$geneMap$children$segments$vp$width <- unit(0.7, "snpc")
    LDmap$LDheatmapGrob$children$geneMap$children$segments$vp$height <- unit(0.7, "snpc")
    LDmap$LDheatmapGrob$children$geneMap$children$diagonal$vp$width <- unit(0.7, "snpc")
    LDmap$LDheatmapGrob$children$geneMap$children$diagonal$vp$height <- unit(0.7, "snpc")
    LDmap$LDheatmapGrob$children$heatMap$children$heatmap$vp$y <- unit(0.875, "npc")
    LDmap$LDheatmapGrob$children$geneMap$children$diagonal$vp$y <- unit(0.875, "npc")
    LDmap$LDheatmapGrob$children$geneMap$children$segments$vp$y <- unit(0.875, "npc")
    LDmap$LDheatmapGrob$children$Key$vp$y <- unit(0.6, "npc")
    LDmap$LDheatmapGrob$children$Key$vp$x <- unit(0.95, "npc")
    LDmap$LDheatmapGrob$children$Key$vp$width <- unit(0.15, "npc")
    LDmap$LDheatmapGrob$children$Key$vp$height <- unit(0.05, "npc")
    LDmap$flipVP$justification <- c('center', 'top')
}
  
  ########################
  ###Generate eQTL Enrichment Plot
  print("Generating eQTL enrichment plot...")
  
  ###Fisher's exact test of eQTL enrichment
  Combined.eQTL.GWAS.Data$iseQTL <- Combined.eQTL.GWAS.Data$Congruence
  Combined.eQTL.GWAS.Data$iseQTL <- ifelse(Combined.eQTL.GWAS.Data$iseQTL=="Non-eQTL", "Non-eQTL", "eQTL")
  if(nrow(as.table(table(Combined.eQTL.GWAS.Data$iseQTL, Combined.eQTL.GWAS.Data$significance))) < 2 | ncol(as.table(table(Combined.eQTL.GWAS.Data$iseQTL, Combined.eQTL.GWAS.Data$significance))) < 2){
    NoFisher <- TRUE; print ("Not enough data to compute enrichment significance for Plot C")
  }
  if(NoFisher==FALSE){fisher <- fisher.test(table(Combined.eQTL.GWAS.Data$iseQTL,
                                                  Combined.eQTL.GWAS.Data$significance))}
  if(NoFisher==FALSE){fpvalue <- fisher$p.value}
  
  ### Generate plot 2
  p2 <- ggplot2::ggplot(Combined.eQTL.GWAS.Data) +
    ggplot2::aes(x=significance, y = 1, fill = (Congruence)) +
    ggplot2::geom_bar(stat = 'identity', position = "fill") +
    ggplot2::ggtitle(paste("Enrichment of eQTLs among\nGWAS-significant SNPs")) +
    ggplot2::ylab("Proportion of SNPs that are eQTLs\n") +
    ggplot2::xlab(paste("\nGWAS significance (threshold p <", sigpvalue_GWAS,")"))+
    ggplot2::ylim(0,1.2) +
    if(NoFisher==FALSE){ggpubr::geom_signif(y_position=c(1.1,1.1),
                                            xmin=c("Non-significant"),
                                            xmax=c("Significant"),
                                            annotation=(paste("p =", formatC(fpvalue, format = "e", digits = 2))),
                                            tip_length=0.05)}
  
  if(Congruentdata == TRUE & Incongruentdata == TRUE){
    p2 <- p2 + ggplot2::scale_fill_manual(values = c("Congruent" = "#000099", "Incongruent" = "#990000", "Non-eQTL" = "#C0C0C0")) + 
      ggplot2::guides(fill = guide_legend("eQTL Direction of Effect"))}
  
  if(Congruentdata == TRUE & Incongruentdata == FALSE & congruence == TRUE){
    p2 <- p2 + ggplot2::scale_fill_manual(values = c("Congruent" = "#000099", "Non-eQTL" = "#C0C0C0")) +
      ggplot2::guides(fill = guide_legend("eQTL Direction of Effect"))}
  
  if(Congruentdata == TRUE & Incongruentdata == FALSE & congruence == FALSE){
    p2 <- p2 + ggplot2::scale_fill_manual(labels = c("Non-eQTL", "eQTL"), values = c("Congruent" = "#ffee00", "Non-eQTL" = "#360f70")) + 
      ggplot2::guides(fill = guide_legend(""))}
  
  if(Congruentdata == FALSE & Incongruentdata == TRUE){
    p2 <- p2 + ggplot2::scale_fill_manual(values = c("Incongruent" = "#990000", "Non-eQTL" = "#C0C0C0")) +
      ggplot2::guides(fill = guide_legend("eQTL Direction of Effect"))}
  
  
  
  ########################
  ###Generate P-P plot
  print("Generating P-P plot...")
  
  p3 <- ggplot2::ggplot(data=Combined.eQTL.GWAS.Data[ which( !is.na(Combined.eQTL.GWAS.Data$NES)) , ]) +
    ggplot2::guides(color = guide_legend("Direction of Effect")) +
    ggplot2::xlab("\n-log10(eQTL p-value)") +
    ggplot2::ylab(paste("-log10(GWAS p-value)\n")) +
    ggplot2::scale_y_continuous(limits=c(NA,ylimd))+
    ggplot2::scale_x_continuous(limits=c(NA,xlimd))+
    ggplot2::theme(legend.position = "right")
  
  ### Add data and annotations to plot 3 if LD info is supplied
  ### For congruent data
  if(isTRUE(LD.df) == FALSE & Congruentdata == TRUE & nrow(Combined.eQTL.GWAS.Data[ which( !is.na(Combined.eQTL.GWAS.Data$NES) & Combined.eQTL.GWAS.Data$Congruence == "Congruent") , ]) >=2){
    pearson.congruent <- cor.test(Combined.eQTL.GWAS.Data[ which( !is.na(Combined.eQTL.GWAS.Data$NES) & Combined.eQTL.GWAS.Data$Congruence == "Congruent") , ]$NeglogeQTLpValue,
                                  Combined.eQTL.GWAS.Data[ which( !is.na(Combined.eQTL.GWAS.Data$NES) & Combined.eQTL.GWAS.Data$Congruence == "Congruent") , ]$Neglog10pvalue_GWAS,
                                  method = "pearson")
    
    p3.1 <- p3 + 
      ggplot2::geom_point(data=Combined.eQTL.GWAS.Data[ which(Combined.eQTL.GWAS.Data$Congruence == "Congruent" & is.na(Combined.eQTL.GWAS.Data$R2cong) == TRUE), ],
                          aes(y=Neglog10pvalue_GWAS, x=NeglogeQTLpValue, fill=R2cong), shape = 21, size = 4)  + 
      ggplot2::geom_point(data=Combined.eQTL.GWAS.Data[ which(Combined.eQTL.GWAS.Data$Congruence == "Congruent" & is.na(Combined.eQTL.GWAS.Data$R2cong) == FALSE) , ],
                          aes(y=Neglog10pvalue_GWAS, x=NeglogeQTLpValue, fill=R2cong), shape = 21, size = 4)  +
      ggplot2::geom_smooth(data=Combined.eQTL.GWAS.Data[ which(Combined.eQTL.GWAS.Data$Congruence == "Congruent") , ],
                           aes(y=Neglog10pvalue_GWAS, x=NeglogeQTLpValue), color = "black", method = "lm", formula = (y ~ x)) 
    
    p3.1 <- p3.1 + 
      ggplot2::geom_point(data=Combined.eQTL.GWAS.Data.cong,
                          aes(y=Neglog10pvalue_GWAS, x=NeglogeQTLpValue), shape = 23, size = 4, fill = "#33FF33") +
      ggplot2::geom_text(x= -Inf, y= Inf, 
                         label = paste(sep = "", "Pearson correlation: ", round(pearson.congruent$estimate, 3), ",\np-value: ", 
                                       formatC(pearson.congruent$p.value, format = "e", digits = 2)), hjust = -0.05, vjust= 1.1, family = "mono", color = "black")
  } else {if(isTRUE(LD.df) == FALSE){
    print("Not enough data to generate P-P plot for Congruent eQTLs")
  }}
  
  ### For incongruent data
  if(isTRUE(LD.df) == FALSE & Incongruentdata == TRUE & nrow(Combined.eQTL.GWAS.Data[ which( !is.na(Combined.eQTL.GWAS.Data$NES) & Combined.eQTL.GWAS.Data$Congruence == "Incongruent") , ]) >=2){
    pearson.incongruent <- cor.test(Combined.eQTL.GWAS.Data[ which( !is.na(Combined.eQTL.GWAS.Data$NES) & Combined.eQTL.GWAS.Data$Congruence == "Incongruent") , ]$NeglogeQTLpValue,
                                    Combined.eQTL.GWAS.Data[ which( !is.na(Combined.eQTL.GWAS.Data$NES) & Combined.eQTL.GWAS.Data$Congruence == "Incongruent") , ]$Neglog10pvalue_GWAS,
                                    method = "pearson")
    
    p3.2 <- p3 + 
      ggplot2::geom_point(data=Combined.eQTL.GWAS.Data[ which(Combined.eQTL.GWAS.Data$Congruence == "Incongruent" & is.na(Combined.eQTL.GWAS.Data$R2cong) == TRUE) , ],
                          aes(y=Neglog10pvalue_GWAS, x=NeglogeQTLpValue, fill=R2incong), shape = 21, size = 4)  + 
      ggplot2::geom_point(data=Combined.eQTL.GWAS.Data[ which(Combined.eQTL.GWAS.Data$Congruence == "Incongruent" & is.na(Combined.eQTL.GWAS.Data$R2cong) == FALSE) , ],
                          aes(y=Neglog10pvalue_GWAS, x=NeglogeQTLpValue, fill=R2incong), shape = 21, size = 4)  +
      ggplot2::scale_fill_gradient(bquote(atop(R^2 ~ "with lead SNP", ~{.(mostsigsnp.incong)} ~ ", in green")), low="#990000", high="#FFCC33", na.value = "grey80") +
      ggplot2::geom_smooth(data=Combined.eQTL.GWAS.Data[ which(Combined.eQTL.GWAS.Data$Congruence == "Incongruent") , ],
                           aes(y=Neglog10pvalue_GWAS, x=NeglogeQTLpValue), color = "black", method = "lm", formula = (y ~ x)) 
    
    p3.2 <- p3.2 + 
      ggplot2::geom_point(data=Combined.eQTL.GWAS.Data.incong,
                          aes(y=Neglog10pvalue_GWAS, x=NeglogeQTLpValue), shape = 23, size = 4, fill = "#33FF33") +
      ggplot2::geom_text(x= -Inf, y= Inf, 
                         label = paste(sep = "", "Pearson correlation: ", round(pearson.incongruent$estimate, 3), ",\np-value: ", 
                                       formatC(pearson.incongruent$p.value, format = "e", digits = 2)), hjust = -0.05, vjust = 1.1, family = "mono", color = "black") +
      ggplot2::ggtitle(paste("P-P plot, GWAS p-value vs. eQTL p-value\nfor incongruent variants"))
  } else {if(isTRUE(LD.df) == FALSE & congruence == TRUE){
    print("Not enough data to generate P-P plot for Incongruent eQTLs")
  }}
  
  ### Add color scale to plot 3.1 if LD info is supplied
  if(isTRUE(LD.df) == FALSE & congruence == FALSE){p3.1 <- p3.1 +ggplot2::scale_fill_viridis_c(bquote(atop(R^2 ~ "with lead SNP", ~{.(mostsigsnp.cong)} ~ ", in green")), option = "C", na.value = "grey80") +
    ggplot2::ggtitle(paste("P-P plot,\nGWAS p-value vs. eQTL p-value"))
  } else {
    if(isTRUE(LD.df) == FALSE & congruence == TRUE){p3.1 <- p3.1 + ggplot2::scale_fill_gradient(bquote(atop(R^2 ~ "with lead SNP", ~{.(mostsigsnp.cong)} ~ ", in green")), low="#000099", high="#33FFFF", na.value = "grey80") + 
      ggplot2::ggtitle(paste("P-P plot, GWAS p-value vs. eQTL p-value\nfor congruent variants"))
    }}
  
  ### Add data and annotations to plot 3 if LD data is not supplied
  if(isTRUE(LD.df) == TRUE){
    
    ### Annotations for congruent data points
    if(Congruentdata == TRUE & nrow(Combined.eQTL.GWAS.Data[ which( !is.na(Combined.eQTL.GWAS.Data$NES) & Combined.eQTL.GWAS.Data$Congruence == "Congruent") , ]) >=2){
      pearson.congruent <- cor.test(Combined.eQTL.GWAS.Data[ which( !is.na(Combined.eQTL.GWAS.Data$NES) & Combined.eQTL.GWAS.Data$Congruence == "Congruent") , ]$NeglogeQTLpValue,
                                    Combined.eQTL.GWAS.Data[ which( !is.na(Combined.eQTL.GWAS.Data$NES) & Combined.eQTL.GWAS.Data$Congruence == "Congruent") , ]$Neglog10pvalue_GWAS,
                                    method = "pearson")
      p3 <- p3 + 
        ggplot2::geom_text(x= -Inf, y= Inf,
                           label = paste(sep = "", "Pearson correlation: ", 
                                         round(pearson.congruent$estimate, 3),
                                         ",\np-value: ", 
                                         formatC(pearson.congruent$p.value, format = "e", digits = 2)),
                           color="#000099", hjust = -0.05, vjust= 1.1, family = "mono")
      
      ### Regression line for congruent data points  
      if(congruence == TRUE){
        p3 <- p3 + ggplot2::geom_smooth(data=Combined.eQTL.GWAS.Data[ which( !is.na(Combined.eQTL.GWAS.Data$NES) & Combined.eQTL.GWAS.Data$Congruence == "Congruent") , ],
                                        aes(y=Neglog10pvalue_GWAS, x=NeglogeQTLpValue, color=Congruence), method = "lm", formula = (y ~ x))
      } else {
        if(congruence == FALSE)     
          p3 <- p3 + ggplot2::geom_smooth(data=Combined.eQTL.GWAS.Data[ which( !is.na(Combined.eQTL.GWAS.Data$NES) & Combined.eQTL.GWAS.Data$Congruence == "Congruent") , ],
                                          aes(y=Neglog10pvalue_GWAS, x=NeglogeQTLpValue, color=Congruence), method = "lm", formula = (y ~ x), show.legend = FALSE, color = "#ffee00")
      }
      
      ### Data for congruent data points
      if(congruence == TRUE){
        p3 <- p3 +
          ggplot2::geom_point(data=Combined.eQTL.GWAS.Data[ which( !is.na(Combined.eQTL.GWAS.Data$NES) & Combined.eQTL.GWAS.Data$Congruence == "Congruent") , ],
                              aes(y=Neglog10pvalue_GWAS, x=NeglogeQTLpValue, color=Congruence)) 
      } else {
        p3 <- p3 + 
          ggplot2::geom_point(data=Combined.eQTL.GWAS.Data[ which( !is.na(Combined.eQTL.GWAS.Data$NES) & Combined.eQTL.GWAS.Data$Congruence == "Congruent") , ],
                              aes(y=Neglog10pvalue_GWAS, x=NeglogeQTLpValue), color = "#360f70") 
      }
    } else {
      print("Not enough data to generate P-P plot for Congruent eQTLs")
    }
    
    ### Annotations and data for incongruent data points
    if(Incongruentdata == TRUE & nrow(Combined.eQTL.GWAS.Data[ which( !is.na(Combined.eQTL.GWAS.Data$NES) & Combined.eQTL.GWAS.Data$Congruence == "Incongruent") , ]) >=2){
      pearson.incongruent <- cor.test(Combined.eQTL.GWAS.Data[ which( !is.na(Combined.eQTL.GWAS.Data$NES) & Combined.eQTL.GWAS.Data$Congruence == "Incongruent") , ]$NeglogeQTLpValue,
                                      Combined.eQTL.GWAS.Data[ which( !is.na(Combined.eQTL.GWAS.Data$NES) & Combined.eQTL.GWAS.Data$Congruence == "Incongruent") , ]$Neglog10pvalue_GWAS,
                                      method = "pearson")
      p3 <- p3 + ggplot2::geom_point(data=Combined.eQTL.GWAS.Data[ which( !is.na(Combined.eQTL.GWAS.Data$NES) & Combined.eQTL.GWAS.Data$Congruence == "Incongruent") , ],
                                     aes(y=Neglog10pvalue_GWAS, x=NeglogeQTLpValue, color=Congruence)) +
        ggplot2::geom_smooth(data=Combined.eQTL.GWAS.Data[ which( !is.na(Combined.eQTL.GWAS.Data$NES) & Combined.eQTL.GWAS.Data$Congruence == "Incongruent") , ],
                             aes(y=Neglog10pvalue_GWAS, x=NeglogeQTLpValue, color=Congruence), method = "lm", formula = (y ~ x)) + 
        ggplot2::geom_text(x= -Inf, y= Inf,
                           label = paste(sep = "", "\n\nPearson correlation: ",
                                         round(pearson.incongruent$estimate, 3),
                                         ",\np-value: ",
                                         formatC(pearson.incongruent$p.value, format = "e", digits = 2)),
                           color="#990000", hjust = -0.05, vjust= 1.1, family = "mono")
    } else {if(congruence == "TRUE"){
      print("Not enough data to generate P-P Plot for Incongruent eQTLs")
    }}
    
    ### Add correct color scales 
    if(Congruentdata == TRUE & Incongruentdata == TRUE){p3 <- p3 + scale_color_manual(values = c("Congruent" = "#000099", "Incongruent" = "#990000")) }
    
    if(Congruentdata == TRUE & Incongruentdata == FALSE){p3 <- p3 + scale_color_manual(values = c("Congruent" = "#000099"))}
    
    if(Congruentdata == FALSE & Incongruentdata == TRUE){p3 <- p3 + scale_color_manual(values = c("Incongruent" = "#990000"))} 
    
    p3 <- p3 + ggplot2::ggtitle(paste("P-P plot,\nGWAS p-value vs. eQTL p-value"))
  }
  
  
  
  ########################
  ###Combine plots
  print("Merging and plotting...")
  if(PanTissue == TRUE & length(tissue) >= 2){tissue <- "MultiTissue"}
  if(PanTissue == TRUE & tissue == "all"){tissue <- "PanTissue"}
  if(PanTissue == TRUE) tissuetitle <- paste(", ", tissue) else tissuetitle <- paste(", in", tissue)

  
  if(isTRUE(LD.df) == FALSE & Incongruentdata == FALSE){
    p4 <- ((p1 + genetracks + patchwork::wrap_elements(ggplotify::as.ggplot(LDmap$LDheatmapGrob) +theme(plot.margin=unit(c(-0.5,0,-2,0),"npc")) +coord_fixed()) + patchwork::plot_layout(nrow = 3, byrow = FALSE, heights = c(2,(genometrackheight/5),NA))) |
            (p2 / p3.1 / patchwork::plot_spacer())) +
    patchwork::plot_layout(ncol = 2, widths = c(2.5,1)) +
    patchwork::plot_annotation(tag_levels = 'A', tag_suffix = ".") & theme(plot.tag = element_text(size = 18), text = element_text(size = 12), plot.tag.position = c(0, 0.995))
  }
  
  if(isTRUE(LD.df) == FALSE & Congruentdata == FALSE){
    p4 <- ((p1 + genetracks + patchwork::wrap_elements(ggplotify::as.ggplot(LDmap$LDheatmapGrob) +theme(plot.margin=unit(c(-0.5,0,-2,0),"npc")) +coord_fixed()) + patchwork::plot_layout(nrow = 3, byrow = FALSE, heights = c(2,(genometrackheight/5),NA))) |
             (p2 / p3.2 / patchwork::plot_spacer())) +
      patchwork::plot_layout(ncol = 2, widths = c(2.5,1)) +
      patchwork::plot_annotation(tag_levels = 'A', tag_suffix = ".") & theme(plot.tag = element_text(size = 18), text = element_text(size = 12), plot.tag.position = c(0, 0.995))
  }
  
  if(isTRUE(LD.df) == FALSE & Congruentdata == TRUE & Incongruentdata == TRUE){
    p4 <- ((p1 + genetracks + patchwork::wrap_elements(ggplotify::as.ggplot(LDmap$LDheatmapGrob) +theme(plot.margin=unit(c(-0.5,0,-2,0),"npc")) +coord_fixed()) + patchwork::plot_layout(nrow = 3, byrow = FALSE, heights = c(2,(genometrackheight/5),NA))) |
             (p2 / p3.1 / p3.2)) +
      patchwork::plot_layout(ncol = 2, widths = c(2.5,1)) +
      patchwork::plot_annotation(tag_levels = 'A', tag_suffix = ".") & theme(plot.tag = element_text(size = 18), text = element_text(size = 12), plot.tag.position = c(0, 0.995))
  }
  
  if(isTRUE(LD.df) == TRUE) {
    p4 <- (p1 + genetracks + patchwork::plot_spacer() + (p2 + p3 + patchwork::plot_layout(ncol = 2, widths = c(2,3))) + 
             patchwork::plot_layout(ncol=1, height = c(4,genometrackheight, 0.1, 2)) + 
             patchwork::plot_annotation(tag_levels = 'A', tag_suffix = ".") & theme(plot.tag = element_text(size = 18), text = element_text(size = 12)))
  }
  
  pfinal <- p4 + patchwork::plot_annotation(title = paste("    Association study (",trait,") and eQTL analysis for ", gene, tissuetitle, sep = ""), theme=theme(plot.title = element_text(size = 19)))
  
  
  
  ########################  
  ### Export final plot
  if(isTRUE(LD.df) == FALSE & wi == "wi"){
    wi <- 18.5
  }
  if(isTRUE(LD.df) == TRUE & wi == "wi"){
    wi <- 12
  }
  if(isTRUE(LD.df) == FALSE){
    hgt <- ((1.33246*(wi) - 0.01427*(wi)^2 - 7.739))
  }
  if(isTRUE(LD.df) == TRUE){
    hgt <- wi
  }
  if(congruence == TRUE){congruence <- "WithCongruenceData"} else {congruence <- "WithoutCongruenceData"}
  if(isTRUE(LD.df) == FALSE){LDinfo <- "WithLinkageData"} else {LDinfo <- "WithoutLinkageData"}
  if(saveplot == TRUE){ggsave(pfinal, filename=paste(gene, trait, tissue, congruence, LDinfo, "eQTpLot", "png", sep="."), dpi=res, units="in", height=hgt, width=wi)}
  if(getplot == TRUE){return(pfinal)}
  
}