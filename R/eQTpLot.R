#' An eQTpLot function
#'
#' This function creates an eQTpLot composite
#' @param GWAS.df Dataframe, one row per SNP, with columns
#'     CHR Chromosome for SNP (X coded numerical ly as 23)
#'     POS Chromosomal position for each SNP, in base pairs
#'     SNP Variant ID (such ws dbSNP ID "rs...". Note: Must be the same naming scheme as  used in the eQTL.df to ensure proper matching)
#'     pvalue_GWAS pvalue for the SNP from GWAS analysis
#'     beta_GWAS beta for SNP from GWAS analysis
#'     Trait name of trait for which GWAS was run. This column must be present, with a value for every row, even if your GWAS.df contains data for only one trait.
#'         (Note: GWAS.df can contain multiple traits (such as one might get from a PheWAS. Specifying the trait to be analyzed will filter only for GWAS entires for this trait)
#' @param eQTL.df Dataframe, one row per SNP, with columns
#'     SNP Variant ID (such ws dbSNP ID "rs...". Note: Must be the same naming scheme as  used in the GWAS.df to ensure proper matching)
#'     Gene Gene symbol/name to which the eQTL expression data refers (Note: gene symbol/name must match entries in Genes.df to ensure proper matching)
#'     pvalue_eQTL pvalue for the SNP from eQTL analysis (such as one might download from the GTEx Portal)
#'     NES_eQTL NES (normalized effect size) for the SNP from eQTL analysis
#'         (Per GTEx: defined as the slope of the linear regression, and is computed as the effect of the alternative allele (ALT) relative to the reference allele (REF) in the human genome reference
#'         (i.e., the eQTL effect allele is the ALT allele). NES are computed in a normalized space where magnitude has no direct biological interpretation.)
#'     Tissue tissue type to which the eQTL pvalue/beta refer
#'         (Note: eQTL.df can contain multiple tissue types. Specifying the tissue type to be analyzed will filter only for eQTL entires for this tissue type.
#'         Alternatively, setting tissue type to "all" (the default setting) will automatically pick the smallest eQTL pvalue for each SNP across all tissues for a PanTissue analysis)
#' @param Genes.df Dataframe, one row per gene, with the following columns
#'         (Note: eQTpLot automatically loads a default Genes.df database containing information for both genomic builds hg19 and hg38,
#'         but you may wish to specify our own Genes.df dataframe if your gene of interest is not included in the default dataframe, or if your eQTL data uses a different gene naming scheme
#'         (for example, Gencode ID instead of gene symbol)):
#'     Gene Gene symbol/name for which the Coordinate data (below) refers to
#'         (Note: gene symbol/name must match entries in eQTL.df to ensure proper matching)
#'     CHR Chromosome the gene is on (X coded numerically as 23)
#'     Start Chromosomal coordinate of start position (in basepairs) to use for gene (Note: this should be the smaller of the two values between Start and Stop)
#'     Stop Chromosomal coordinate of end position (in basepairs) to use for gene (Note: this should be the larger of the two values between Start and Stop)
#'     Build The genome build (either hg19 or hg38) for the location data -
#'         the default Genes.df dataframe contains entries for both genome builds for each gene, and the script will select the appropriate entry based on the specified gbuild (default is hg19)).
#' @param gene name/symbol of gene to analyze, in quotes
#' @param trait name of trait to analyze, in quotes (from GWAS.df)
#' @param sigpvalue_GWAS GWAS pvalue significance threshold to use (this value will be used for a horizontal line in plot A, and to define GWAS significant/non-significant SNPs for plot C)
#' @param sigpvalue_eQTLeQTL pvalue significance threshold to use (eQTL data with a pvalue larger than this threshold will be excluded from the analysis)
#' @param tissue tissue type, in quotes, for eQTL data to use (from eQTL.df, default setting is "all" for Pan-Tissue analysis)
#' @param range range, in kB, to extend analysis window on either side of gene boundry. Default is 200 kb
#' @param NESeQTLRange the maximum and minimum limits in the format c(min,max), to display for the NES value in eQTL.df.
#'          The default setting will adjust the size limits automatically for your data, whereas specifying the limits
#'          can keep them consistent between plots.
#' @param ylima used to manually adjust the y axis limit in plot A, if needed
#' @param ylimd used to manually adjust the y axis limit in plot D, if needed
#' @param xlimd used to manually adjust the x axis limit in plot D, if needed
#' @param genometrackheight used to set the height of the genome track panel (B), with default setting of 1.5.
#'           Gene-dense regions may require more plotting space, whereas gene-poor regions may look better with less plotting space.
#' @param gbuild the genome build to use, in quotes, for fetching genomic information for panel B.
#'            Default is "hg19" but can specify "hg38" if needed. This build should match the genome build used for "CHR" and "POS" in the GWAS.df
#' @param res resolution of the output plot image (default is 300 dpi)
#' @param hgt height of the output plot image (default is 12 inches)
#' @param wi width of the output plot image (default is 14 inches)
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

eQTpLot <- function(GWAS.df, eQTL.df, Genes.df, gene, trait,
                    sigpvalue_GWAS = 5e-8, sigpvalue_eQTL = 0.05,
                    tissue = "all", range = 200, NESeQTLRange = c(NA,NA),
                    ylima = NA, ylimd = NA, xlimd = NA,
                    genometrackheight = 1.5, gbuild = "hg19",
                    res = 300, hgt = 18, wi = 18,
                    getplot = TRUE, saveplot = TRUE){

  ########################
  ### Required Packages
  #require(ggnewscale)
  #require(patchwork)
  #require(Gviz)
  #require(GenomicRanges)
  #require(biomaRt)
  #require(ggpubr)

  ########################
  ### Check Data
  print("Checking input data...")

  if(missing(Genes.df)){Genes.df <- eQTpLot:::genes}

  if(all(c("CHR", "POS", "SNP", "beta_GWAS", "pvalue_GWAS","Trait") %in% colnames(GWAS.df))==FALSE) {
    stop("The data supplied to GWAS.df must contain columns 'CHR', 'POS', 'SNP', 'beta_GWAS', 'pvalue_GWAS', and 'Trait'")
  }

  if(all(c("SNP", "Gene", "pvalue_eQTL", "NES_eQTL", "Tissue") %in% colnames(eQTL.df))==FALSE) {
    stop("The data supplied to eQTL.df must contain columns 'SNP', 'Gene', 'pvalue_eQTL', 'NES_eQTL', and 'Tissue'")
  }

  if(all(c("Gene", "CHR", "Start", "Stop", "Build") %in% colnames(Genes.df))==FALSE) {
    stop("The data supplied to Genes.df must contain columns 'Gene', 'CHR', 'Start', 'Stop', 'Build'")
  }


  if(is.numeric(eQTL.df$pvalue_eQTL) == FALSE | is.numeric(eQTL.df$NES_eQTL) == FALSE) {
    stop('Sorry, the eQTL.df dataframe must contain only numeric data for pvalue_eQTL and NES_eQTL')
  }

  if(is.numeric(GWAS.df$pvalue_GWAS) == FALSE | is.numeric(GWAS.df$beta_GWAS) == FALSE | is.integer(GWAS.df$POS) == FALSE | is.integer(GWAS.df$CHR) == FALSE) {
    stop('Sorry, the GWAS.df dataframe must contain only numeric data for CHR, POS, pvalue_GWAS, and beta_GWAS (Note: chromosomes must be coded numerically)')
  }

  if(is.integer(Genes.df$CHR) == FALSE | is.integer(Genes.df$Start) == FALSE | is.integer(Genes.df$Stop) == FALSE) {
    stop('Sorry, the Genes.df dataframe must contain only integer values for CHR, Start, and Stop (Note: chromosomes must be coded nuemrically)')
  }

  if(dim(Genes.df[which(Genes.df$Gene == gene & Genes.df$Build == gbuild), ])[1] == 0) {
    stop('Sorry, there is no information for the gene ', paste(gene), ' in the Genes.df dataframe for the genomic build ', paste(gbuild))
  }

  if((tissue %in% eQTL.df$Tissue) == FALSE & (tissue == "all") == FALSE){
    stop('Sorry, the tissue ', paste(tissue), ' does not exist in the eQTL.df dataframe. Tissues included in the data supplied to eQTL.df are:\n',
         paste("'",as.character(unique(eQTL.df$Tissue)),"'",collapse=", ",sep=""))
  }
  if((trait %in% GWAS.df$Trait) == FALSE) {
    stop('Sorry, the  trait ', paste(trait), ' does not exist in the GWAS.df dataframe. Traits included in the data supplied to GWAS.df are:\n',
         paste("'",as.character(unique(GWAS.df$Trait)),"'",collapse=", ",sep=""))
  }

  ########################
  ### Compile Data
  print("Compiling data...")

  ### Set genomic ranges for data selection
  rangebp <- range*1000
  startpos <- (min(Genes.df[which(Genes.df$Gene == gene & Genes.df$Build == gbuild), ] %>% dplyr::select(Start)) - rangebp)
  stoppos <- (max(Genes.df[which(Genes.df$Gene == gene & Genes.df$Build == gbuild), ] %>% dplyr::select(Stop)) + rangebp)
  chromosome <- Genes.df$CHR[Genes.df$Gene == gene & Genes.df$Build == gbuild]

  ### Subset GWAS.df for gene of interest, check to make sure data is present
  gwas.data <- GWAS.df[which(GWAS.df$CHR == chromosome & GWAS.df$Trait == trait & GWAS.df$POS >= startpos & GWAS.df$POS <= stoppos & !(is.na(GWAS.df$pvalue_GWAS)) & !(is.na(GWAS.df$beta_GWAS))), ]
  if(dim(gwas.data)[1] == 0) stop('Sorry, GWAS.df does not contain data for any SNPs in the ', paste(range), 'kb flanking the gene ', paste(gene), ' for the trait ', paste(trait))
  if(dim(gwas.data[which(gwas.data$pvalue_GWAS <= sigpvalue_GWAS), ])[1] == 0) {
    NoFisher<-TRUE; print(paste(sep = "", 'WARNING: GWAS.df does not contain any SNPs with p-value < sigpvalue_GWAS within the ',
                                range, 'kb flanking the gene ', gene, ' for the trait ',
                                trait, '. eQTL Enrcihment Plot statistics will not be calculated'))
  } else {NoFisher <- FALSE}

  ### Subset eQTL.df for tissue and gene of interest, check to make sure data is present
  eqtl.data <- eQTL.df[which(eQTL.df$Gene == gene & eQTL.df$pvalue_eQTL <= sigpvalue_eQTL & !(is.na(eQTL.df$NES_eQTL))), ]
  if(dim(eqtl.data)[1] == 0) stop('Sorry, eQTL.df does not have any data for the gene ', paste(gene), ' meeting your sigpvalue_eQTL threshold')
  PanTissue <- ifelse(tissue == "all", TRUE, FALSE)
  if(PanTissue == TRUE) {
    eqtl.data <- eqtl.data[which(!(is.na(eqtl.data$SNP))), ] %>% dplyr::group_by(SNP) %>% dplyr::slice(which.min(pvalue_eQTL));
    eqtl.data$Tissue <- "PanTissue"
  } else {
    eqtl.data <- eqtl.data[which(eqtl.data$Tissue == tissue), ]
  }
  if(dim(eqtl.data)[1] == 0) stop('Sorry, there are no eQTLs for the tissue', paste(tissue), ' with a p-value < sigpvalue_eQTL')
  eqtl.data <- dplyr::ungroup(eqtl.data)

  ### Join GWAS and eQTL data, check to make sure there is at least some overlap in SNPs between the two datasets
  gwas.data$SNP <- as.factor(gwas.data$SNP)
  eqtl.data$SNP <- as.factor(eqtl.data$SNP)
  combinedSNPS <- sort(union(levels(gwas.data$SNP), levels(eqtl.data$SNP)))
  Combined.eQTL.GWAS.Data <- dplyr::left_join(dplyr::mutate(gwas.data, SNP=factor(SNP, levels=combinedSNPS)),
                                              dplyr::mutate(eqtl.data, SNP=factor(SNP, levels=combinedSNPS)))
  if(dim(Combined.eQTL.GWAS.Data)[1] == 0) {
    stop('Sorry, for the gene ', paste(gene), ' and the trait ', paste(trait), ' there is no overlap between the SNPs in your GWAS.df and eQTL.df')
  }

  ### Determine directions of effect and congruence
  Combined.eQTL.GWAS.Data$DirectionOfEffect_GWAS <- ifelse(Combined.eQTL.GWAS.Data$beta_GWAS < 0, "Negative", ifelse(Combined.eQTL.GWAS.Data$beta_GWAS > 0, "Positive", NA))
  Combined.eQTL.GWAS.Data$DirectionOfEffect_eQTL <- ifelse(Combined.eQTL.GWAS.Data$NES_eQTL < 0, "DOWN", ifelse(Combined.eQTL.GWAS.Data$NES_eQTL > 0, "UP", NA))
  Combined.eQTL.GWAS.Data$Congruence <- (Combined.eQTL.GWAS.Data$beta_GWAS*Combined.eQTL.GWAS.Data$NES_eQTL)
  Combined.eQTL.GWAS.Data$Congruence <- ifelse(Combined.eQTL.GWAS.Data$Congruence < 0, "Incongruent", ifelse(Combined.eQTL.GWAS.Data$Congruence > 0, "Congruent", NA))

  ### Build final dataframe
  Combined.eQTL.GWAS.Data$NeglogeQTLpValue <- -(log10(Combined.eQTL.GWAS.Data$pvalue_eQTL))
  Combined.eQTL.GWAS.Data$Neglog10pvalue_GWAS <- -(log10(Combined.eQTL.GWAS.Data$pvalue_GWAS))
  Combined.eQTL.GWAS.Data <- Combined.eQTL.GWAS.Data[which(!(is.na(Combined.eQTL.GWAS.Data$pvalue_GWAS))), ]
  Combined.eQTL.GWAS.Data$significance <- ifelse(Combined.eQTL.GWAS.Data$pvalue_GWAS >= sigpvalue_GWAS, "Non-significant", "Significant")
  Combined.eQTL.GWAS.Data$Congruence[is.na(Combined.eQTL.GWAS.Data$Congruence)] <- "Non-eQTL"
  Combined.eQTL.GWAS.Data <- Combined.eQTL.GWAS.Data %>% dplyr::mutate(Congruence=factor(Congruence, levels=c("Non-eQTL", "Congruent", "Incongruent"), ordered=TRUE))
  if(dim(Combined.eQTL.GWAS.Data[ which( !is.na(Combined.eQTL.GWAS.Data$NES_eQTL) & Combined.eQTL.GWAS.Data$Congruence == "Congruent") , ])[1] == 0){
    Congruentdata <- FALSE} else {Congruentdata <- TRUE
  }
  if(dim(Combined.eQTL.GWAS.Data[ which( !is.na(Combined.eQTL.GWAS.Data$NES_eQTL) & Combined.eQTL.GWAS.Data$Congruence == "Incongruent") , ])[1] == 0){
    Incongruentdata <- FALSE} else {Incongruentdata <- TRUE
  }


  ########################
  ### Generate main plot
  print("Generating Plot A...")

  ### Set graph limits
  if(is.na(ylima) == TRUE){ylima <- (max(Combined.eQTL.GWAS.Data %>% dplyr::select (Neglog10pvalue_GWAS), na.rm = TRUE) + 1)}
  minpos <- min(Combined.eQTL.GWAS.Data$POS, na.rm = TRUE)
  maxpos <- max(Combined.eQTL.GWAS.Data$POS, na.rm = TRUE)

  ### Generate plot 1
  p1 <-
    ggplot2::ggplot (data=Combined.eQTL.GWAS.Data, aes(stroke = 0)) +
      ggplot2::coord_cartesian(xlim = c(minpos, maxpos), expand = FALSE) +
      ggplot2::geom_point(data = subset(Combined.eQTL.GWAS.Data, Congruence == "Non-eQTL"), shape = 15, color = "black", alpha = 0.2,
                          aes(x=POS, y=Neglog10pvalue_GWAS)) +
      ggplot2::xlab("") +
      ggplot2::ylab(paste("-log10(GWAS p-value)  association with ", trait)) +
      ggplot2::scale_y_continuous(limits=c(NA,ylima)) +
      ggplot2::ggtitle(paste("Association study of ", trait, ", variants colored by eQTL data for ", gene,
                             " (GWAS significance threshold ", sigpvalue_GWAS, ", eQTL significance threshold ",
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

  if(Congruentdata == TRUE & Incongruentdata == FALSE){
    p1 <- p1 + ggplot2::geom_point(data = subset(Combined.eQTL.GWAS.Data, Congruence == "Congruent"),
                          alpha = 1,
                          aes(x=POS, y=Neglog10pvalue_GWAS, fill=NeglogeQTLpValue, alpha = 1, shape = DirectionOfEffect_GWAS, size = abs(NES_eQTL)))+
                          ggplot2::scale_fill_gradient("-log10(eQTL p-value), SNPs with\nCongruous directions of effect",
                                                       low="#000099", high="#33FFFF",
                                                       guide = guide_colorbar(title.position = "top"),
                                                       limits=c(min(Combined.eQTL.GWAS.Data %>%
                                                                      dplyr::select (NeglogeQTLpValue), na.rm = TRUE), max(Combined.eQTL.GWAS.Data %>%
                                                                                                                             dplyr::select (NeglogeQTLpValue), na.rm = TRUE)))
    }

  if(Congruentdata == FALSE & Incongruentdata == TRUE){
    p1 <- p1 + ggplot2::geom_point(data = subset(Combined.eQTL.GWAS.Data, Congruence == "Incongruent"),
                         alpha = 1,
                         aes(x=POS, y=Neglog10pvalue_GWAS, fill=NeglogeQTLpValue, alpha = 1,
                             shape = DirectionOfEffect_GWAS, size = abs(NES_eQTL))) +
                ggplot2::scale_fill_gradient("-log10(eQTL p-value), SNPS with\nIncongruous directions of effect",
                                             low="#990000", high="#FFCC33",
                                             guide = guide_colorbar(title.position = "top"),
                                             limits=c(min(Combined.eQTL.GWAS.Data %>%
                                                            dplyr::select (NeglogeQTLpValue), na.rm = TRUE),max(Combined.eQTL.GWAS.Data %>%
                                                                                                                  dplyr::select (NeglogeQTLpValue), na.rm = TRUE)))
    }

  if(Congruentdata == TRUE & Incongruentdata == TRUE){
    p1 <- p1 + ggplot2::geom_point(data = subset(Combined.eQTL.GWAS.Data, Congruence == "Congruent"),
                          alpha = 1,
                          aes(x=POS, y=Neglog10pvalue_GWAS, fill=NeglogeQTLpValue, alpha = 1,
                              shape = DirectionOfEffect_GWAS, size = abs(NES_eQTL))) +
               ggplot2::scale_fill_gradient("-log10(eQTL p-value), SNPs with\nCongruous directions of effect",
                                            low="#000099", high="#33FFFF",
                                            guide = guide_colorbar(title.position = "top"),
                                            limits=c(min(Combined.eQTL.GWAS.Data %>%
                                                           dplyr::select (NeglogeQTLpValue), na.rm = TRUE),max(Combined.eQTL.GWAS.Data %>%
                                                                                                                 dplyr::select (NeglogeQTLpValue), na.rm = TRUE))) +
               ggnewscale::new_scale_fill() +
               ggplot2::geom_point(data = subset(Combined.eQTL.GWAS.Data, Congruence == "Incongruent"),
                 alpha = 1,
                 aes(x=POS, y=Neglog10pvalue_GWAS, fill=NeglogeQTLpValue, alpha = 1,
                     shape = DirectionOfEffect_GWAS, size = abs(NES_eQTL))) +
               ggplot2::scale_fill_gradient("-log10(eQTL p-value), SNPS with\nIncongruous directions of effect",
                                            low="#990000", high="#FFCC33",
                                            guide = guide_colorbar(title.position = "top"),
                                            limits=c(min(Combined.eQTL.GWAS.Data %>%
                                                           dplyr::select (NeglogeQTLpValue), na.rm = TRUE),max(Combined.eQTL.GWAS.Data %>%
                                                                                                                 dplyr::select (NeglogeQTLpValue), na.rm = TRUE)))
    }

  p1 <- p1 + scale_size_continuous(limits = NESeQTLRange)

  ########################
  ### Generate Gene Tracks
  print("Generating Plot B...")

  ### Generate Gene Track Plot
  if(gbuild == "hg19"){hostname <- "grch37.ensembl.org"}
  if(gbuild == "hg38"){hostname <- "apr2020.archive.ensembl.org"}

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
                                                                      from = min(Combined.eQTL.GWAS.Data$POS, na.rm = TRUE),
                                                                      to= max(Combined.eQTL.GWAS.Data$POS, na.rm = TRUE),
                                                                      showTitle = FALSE,
                                                                      labelPos = "below",
                                                                      distFromAxis = 10,
                                                                      innermargin = 0,
                                                                      maxHeight = (genometrackheight*10),
                                                                      minHeight = (genometrackheight*10),
                                                                      sizes=c(genometrackheight,1),
                                                                      margin = c(-3,-3)))))


  ########################
  ###Generate eQTL Enrichment Plot
  print("Generating Plot C...")

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
          ggplot2::ggtitle(paste("Enrichment of eQTLs among GWAS-significant SNPs")) +
          ggplot2::guides(fill = guide_legend("eQTL Direction of Effect")) +
          ggplot2::ylab("Proportion of SNPs that are eQTLs\n") +
          ggplot2::xlab(paste("\nGWAS significance (threshold p <", sigpvalue_GWAS,")"))+
          ggplot2::ylim(0,1.2) +
    if(NoFisher==FALSE){ggpubr::geom_signif(y_position=c(1.1,1.1),
                                    xmin=c("Non-significant"),
                                    xmax=c("Significant"),
                                    annotation=(paste("p =", formatC(fpvalue, format = "e", digits = 2))),
                                    tip_length=0.05)}

  if(Congruentdata == TRUE & Incongruentdata == TRUE){
    p2 <- p2 + ggplot2::scale_fill_manual(values = c("Congruent" = "#000099", "Incongruent" = "#990000", "Non-eQTL" = "#C0C0C0"))}

  if(Congruentdata == TRUE & Incongruentdata == FALSE){
    p2 <- p2 + ggplot2::scale_fill_manual(values = c("Congruent" = "#000099", "Non-eQTL" = "#C0C0C0"))}

  if(Congruentdata == FALSE & Incongruentdata == TRUE){
    p2 <- p2 + ggplot2::scale_fill_manual(values = c("Incongruent" = "#990000", "Non-eQTL" = "#C0C0C0"))}



  ########################
  ###Generate P-P plot
  print("Generating Plot D...")

  ### Regression Line Equation Function
  lm_eq <- function(df){

    m <- lm(NeglogeQTLpValue ~ Neglog10pvalue_GWAS, df);

    eq <- substitute(y == a + b %.% x*","~~r^2~"="~r2,

                     list(a = format(unname(coef(m)[1]), digits = 2),

                          b = format(unname(coef(m)[2]), digits = 2),

                          r2 = format(summary(m)$r.squared, digits = 3)))

    as.character(as.expression(eq));
  }

  ### Generate plot 3
  p3 <- ggplot2::ggplot(data=Combined.eQTL.GWAS.Data[ which( !is.na(Combined.eQTL.GWAS.Data$NES_eQTL)) , ]) +
          ggplot2::ggtitle(paste("P-P plot, associaton of GWAS p-value with eQTL p-value")) +
          ggplot2::guides(color = guide_legend("Direction of Effect")) +
          ggplot2::xlab("\n-log10(eQTL p-value)") +
          ggplot2::ylab(paste("-log10(GWAS p-value)\n")) +
          ggplot2::scale_y_continuous(limits=c(NA,ylimd))+
          ggplot2::scale_x_continuous(limits=c(NA,xlimd))+
          ggplot2::theme(legend.position = "right")

  if(Congruentdata == TRUE & nrow(Combined.eQTL.GWAS.Data[ which( !is.na(Combined.eQTL.GWAS.Data$NES_eQTL) & Combined.eQTL.GWAS.Data$Congruence == "Congruent") , ]) >=2) {
    pearson.congruent <- cor.test(Combined.eQTL.GWAS.Data[ which( !is.na(Combined.eQTL.GWAS.Data$NES_eQTL) & Combined.eQTL.GWAS.Data$Congruence == "Congruent") , ]$NeglogeQTLpValue,
                                  Combined.eQTL.GWAS.Data[ which( !is.na(Combined.eQTL.GWAS.Data$NES_eQTL) & Combined.eQTL.GWAS.Data$Congruence == "Congruent") , ]$Neglog10pvalue_GWAS,
                                  method = "pearson")
    p3 <- p3 + ggplot2::geom_point(data=Combined.eQTL.GWAS.Data[ which( !is.na(Combined.eQTL.GWAS.Data$NES_eQTL) & Combined.eQTL.GWAS.Data$Congruence == "Congruent") , ],
                                   aes(y=Neglog10pvalue_GWAS, x=NeglogeQTLpValue, color=Congruence)) +
               ggplot2::geom_smooth(data=Combined.eQTL.GWAS.Data[ which( !is.na(Combined.eQTL.GWAS.Data$NES_eQTL) & Combined.eQTL.GWAS.Data$Congruence == "Congruent") , ],
                                    aes(y=Neglog10pvalue_GWAS, x=NeglogeQTLpValue, color=Congruence), method = "lm", formula = (y ~ x)) +
               ggplot2::geom_text(x= -Inf, y= Inf,
                                  label = lm_eq(Combined.eQTL.GWAS.Data[ which( !is.na(Combined.eQTL.GWAS.Data$NES_eQTL) & Combined.eQTL.GWAS.Data$Congruence == "Congruent") , ]),
                                  parse = TRUE, color="#000099", hjust =-0.05, vjust= 1.2, family = "mono") +
               ggplot2::geom_text(x= -Inf, y= Inf,
                                  label = paste(sep = "", "Pearson correlation: ",
                                                round(pearson.congruent$estimate, 3),
                                                ", p-value: ",
                                                formatC(pearson.congruent$p.value,
                                                        format = "e",
                                                        digits = 2)),
                                  color="#000099", hjust =-0.1, vjust= 3.8, family = "mono")
  } else {
      print("Not enough data to generate Plot D for Congruent eQTLs")
  }


  if(Incongruentdata == TRUE & nrow(Combined.eQTL.GWAS.Data[ which( !is.na(Combined.eQTL.GWAS.Data$NES_eQTL) & Combined.eQTL.GWAS.Data$Congruence == "Incongruent") , ]) >=2) {
    pearson.incongruent <- cor.test(Combined.eQTL.GWAS.Data[ which( !is.na(Combined.eQTL.GWAS.Data$NES_eQTL) & Combined.eQTL.GWAS.Data$Congruence == "Incongruent") , ]$NeglogeQTLpValue,
                                  Combined.eQTL.GWAS.Data[ which( !is.na(Combined.eQTL.GWAS.Data$NES_eQTL) & Combined.eQTL.GWAS.Data$Congruence == "Incongruent") , ]$Neglog10pvalue_GWAS,
                                  method = "pearson")
    p3 <- p3 + ggplot2::geom_point(data=Combined.eQTL.GWAS.Data[ which( !is.na(Combined.eQTL.GWAS.Data$NES_eQTL) & Combined.eQTL.GWAS.Data$Congruence == "Incongruent") , ],
                         aes(y=Neglog10pvalue_GWAS, x=NeglogeQTLpValue, color=Congruence)) +
               ggplot2::geom_smooth(data=Combined.eQTL.GWAS.Data[ which( !is.na(Combined.eQTL.GWAS.Data$NES_eQTL) & Combined.eQTL.GWAS.Data$Congruence == "Incongruent") , ],
                                    aes(y=Neglog10pvalue_GWAS, x=NeglogeQTLpValue, color=Congruence),
                                    method = "lm", formula = (y ~ x)) +
               ggplot2::geom_text(x= -Inf, y= Inf,
                                  label = lm_eq(Combined.eQTL.GWAS.Data[ which( !is.na(Combined.eQTL.GWAS.Data$NES_eQTL) & Combined.eQTL.GWAS.Data$Congruence == "Incongruent") , ]),
                                  parse = TRUE, color="#990000", hjust =-0.05, vjust= 4.0, family = "mono") +
                ggplot2::geom_text(x= -Inf, y= Inf,
                                   label = paste(sep = "", "Pearson correlation: ",
                                                 round(pearson.incongruent$estimate, 3),
                                                 ", p-value: ",
                                                 formatC(pearson.incongruent$p.value, format = "e", digits = 2)),
                                                         color="#990000", hjust =-0.1, vjust= 8.5, family = "mono")
  } else {
     print("Not enough data to generate Plot D for Incongruent eQTLs")
  }

  if(Congruentdata == TRUE & Incongruentdata == TRUE){p3 <- p3 + scale_color_manual(values = c("Congruent" = "#000099", "Incongruent" = "#990000")    )}

  if(Congruentdata == TRUE & Incongruentdata == FALSE){p3 <- p3 + scale_color_manual(values = c("Congruent" = "#000099"))}

  if(Congruentdata == FALSE & Incongruentdata == TRUE){p3 <- p3 + scale_color_manual(values = c("Incongruent" = "#990000"))}



  ########################
  ###Combine plots
  print("Merging and plotting...")

  if(PanTissue == TRUE) tissuetitle <-", PanTissue" else tissuetitle <- paste(", in", tissue)
  if(PanTissue == TRUE) tissue <- "PanTissue"

  p4 <- (p1 + genetracks + patchwork::plot_spacer() + (p2 + p3 + patchwork::plot_layout(ncol = 2, widths = c(2,3))) +
              patchwork::plot_layout(ncol=1, height = c(4,genometrackheight, 0.1, 2)) +
              patchwork::plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = 18), text = element_text(size = 12)))

  pfinal <- p4 + patchwork::plot_annotation(title = paste("    Association study (",trait,") and eQTL analysis for ", gene, tissuetitle, sep = ""), theme=theme(plot.title = element_text(size = 18)))

  if(saveplot == TRUE){ggsave(pfinal, filename=paste(gene, trait, tissue, "eQTL", "png", sep="."), dpi=res, units="in", height=hgt, width=wi)}
  if(getplot == TRUE){return(pfinal)}

}
