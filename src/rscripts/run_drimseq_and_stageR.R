library(DRIMSeq, stageR)


drimseq_initialise_from_se <- function(se,
                                       group_name = "group",
                                       assay_to_test = "counts",
                                       save_to_path = "") {
  current_assay_counts <- assay(se, assay_to_test)
  
  # prep dm object for DIU
  
  transcripts_counts <- data.frame(current_assay_counts)
  transcripts_counts <- round(transcripts_counts)
  
  transcripts_counts$feature_id <- rownames(transcripts_counts)
  transcripts_counts$gene_id <- rowData(se)$GENEID
  
  
  samples_df <- data.frame(colData(se)) %>% dplyr::select(c(name, !!sym(group_name)))#%>%dplyr::mutate(group_name = droplevels(sym(group_name)))
  
  samples_df <- droplevels(samples_df)
  print(colnames(samples_df))
  samples_df$test_group <- samples_df[, group_name]
  samples_df$sample_id <- samples_df$name
  
  print(colnames(samples_df))
  group_ <- samples_df[, "test_group"]
  
  design <- model.matrix( ~ 0 + group_)
  
  # print(samples_df)
  d <- DRIMSeq::dmDSdata(counts = transcripts_counts, samples = samples_df)
  
  d <- DRIMSeq::dmFilter(
    d,
    min_samps_feature_expr = 2,
    min_feature_expr = 5,
    min_samps_gene_expr = 2,
    min_gene_expr = 10
  )
  
  
  
  
  # ======================================
  # Perform all DRIMSeq analysis stages
  # =====================================
  ## Calculate precision
  d <- DRIMSeq::dmPrecision(d, design = design)
  #
  # ## Fit full model proportions
  d <- DRIMSeq::dmFit(d, design = design)
  
  proportions_df <- data.frame(proportions(d))
  
  print(paste(save_to_path, assay_to_test, "_proportions.csv", sep = ""))
  
  write.csv(proportions_df,
            paste(save_to_path, assay_to_test, "_proportions.csv", sep = ""))
  
  saveRDS(d,
          paste(save_to_path, assay_to_test, "_drimSeqFitObject.rds", sep = ""))
  print(paste(save_to_path, assay_to_test, "_drimSeqFitObject.rds", sep = ""))
  return(d)
}

#' Perform dmTest on dmFit object and runs stage wise analysis
#'
#' @param d DRIMSeq dmFit object
#' @param contrast_to_test contrast to test
#'
#' @returns
#' @export
#'
#' @examples
dmContrastTest <- function(dmFitObj, contrast_to_test) {
  d <- dmTest(dmFitObj, contrast = contrast_to_test)
  
  
  ## Assign gene-level pvalues to the screening stage
  pScreen <- results(d)$pvalue
  
  # results(d)$pvalue[!is.na(results(d)$pvalue)]
  names(pScreen) <- results(d)$gene_id
  
  pScreen
  length(pScreen)
  ### =================have to  ensure no nan values are included in the list:
  pScreen <- results(d)$pvalue[!is.na(results(d)$pvalue)]
  names(pScreen) <- results(d)$gene_id[!is.na(results(d)$pvalue)]
  pScreen
  length(pScreen)
  # ========================================================
  
  ## Assign transcript-level pvalues to the confirmation stage
  # pConfirmation <- matrix(results(d, level = "feature")$pvalue, ncol = 1)
  # rownames(pConfirmation) <- results(d, level = "feature")$feature_id
  
  # pConfirmation
  # dim(pConfirmation)
  
  ### =================have to  ensure no nan values are included in the list:
  pConfirmation <- matrix(results(d, level = "feature")[!is.na(results(d, level = "feature")$pvalue), "pvalue"], ncol = 1)
  rownames(pConfirmation) <- results(d, level = "feature")[!is.na(results(d, level = "feature")$pvalue), "feature_id"]
  
  
  pConfirmation <- pConfirmation[]
  dim(pConfirmation)
  
  ## Create the gene-transcript mapping
  tx2gene <- results(d, level = "feature")[!is.na(results(d, level = "feature")$pvalue), c("feature_id", "gene_id")]
  
  ### =============== have to filter out genes with only 1 TX for the following test to work properly: ================
  ## Create the gene-transcript mapping
  # tx2gene <- results(d, level = "feature")[!is.na(results(d, level = "feature")$pvalue), c("feature_id", "gene_id")]
  
  # Remove genes with only one TX
  tx2keep <- tx2gene %>%
    dplyr::group_by(gene_id) %>%
    dplyr::mutate(n = n()) %>%
    dplyr::filter(n != 1) %>%
    ungroup() %>%
    dplyr::select(-n)
  
  
  
  # Consequently keep only the TX's corresponding to the updated list of genes:
  # pConfirmation <- pConfirmation[tx2keep$feature_id,]
  
  
  # note I have to do this a bit differently due to different handling of rownames in matrix vs data.frame data types
  # pConfirmation<-  matrix(results(d, level = "feature")[!is.na(results(d, level = "feature")$pvalue),"pvalue"], ncol = 1)
  temp_pConfirm <- results(d, level = "feature")[!is.na(results(d, level = "feature")$pvalue), c("gene_id", "feature_id", "pvalue")]
  rownames(temp_pConfirm) <- temp_pConfirm$feature_id
  temp_pConfirm <- temp_pConfirm[tx2keep$feature_id, ]
  dim(temp_pConfirm)
  
  pConfirmation <- matrix(temp_pConfirm[, "pvalue"])
  rownames(pConfirmation) <- temp_pConfirm$feature_id
  
  # pConfirmation <- pConfirmation[tx2gene$feature_id,]
  
  dim(pConfirmation)
  
  # Similarly keep only relevant genes:
  pScreen <- pScreen[unique(tx2keep$gene_id)]
  ## Create the stageRTx object and perform the stage-wise analysis
  
  
  stageRObj <- stageRTx(
    pScreen = pScreen,
    pConfirmation = pConfirmation,
    pScreenAdjusted = FALSE,
    tx2gene = tx2gene
  )
  stageRObj <- stageWiseAdjustment(object = stageRObj,
                                   method = "dtu",
                                   alpha = 0.05)
  getSignificantGenes(stageRObj)
  getSignificantTx(stageRObj)
  padj <- getAdjustedPValues(stageRObj,
                             order = TRUE,
                             onlySignificantGenes = FALSE)
  head(padj)
  
  padj$gene_05label <- 0
  padj$gene_05label[padj$gene < 0.05] <- 1
  
  padj$TX_05label <- 0
  padj$TX_05label[padj$transcript < 0.05] <- 1
  #
  # print(
  #   plotProportions(d,gene_id=padj$geneID[1] ,group_variable = "ordered_category",plot_type = "barplot",group_colors= category_colors)
  # )
  
  
  padj <- padj %>% left_join(DRIMSeq::mean_expression(d), by = join_by("geneID" == "gene_id"))
  # padj$mean_expression <-mean_expression(d)[!is.na(results(d, level = "gene")$pvalue),]
  return(list(padj, d))
}





run_drimseq_with_stageWise <- function(d, contrast_to_test) {
  dmResultsStageR_adj <- dmContrastTest(d, contrast_to_test = contrast_to_test)
  
  
  summarised_results <- summaraise_dmsq_stageR_results(dmResultsStageR_adj[[1]], contrast_to_test, d)
  
  
  
  return(summarised_results)
  
}



#' Consolidate drimSeq proportions and dmTest results after stageR for a specific contrast into a single df
#'
#' @param dmResultsStageR_adj
#' @param contrastToTest
#'
#' @returns
#' @export
#'
#' @examples
summaraise_dmsq_stageR_results <- function(dmResultsStageR_adj,
                                           tested_contrast,
                                           dmFitObj) {
  # Pull proportion results from drimseq obj for the relevant groups
  
  
  contrast_factors <- tested_contrast
  
  groups_in_contrast <- names(contrast_factors[contrast_factors != 0])
  plus_group_name <- names(contrast_factors[contrast_factors == 1])
  minus_group_name <- names(contrast_factors[contrast_factors == -1])
  # NOTE: This is a temporary solution for only the simple contrast case in which contrast = group1 - group2
  # TODO: generalise this to work for more complex contrasts too
  
  
  
  # join all adjusted drimseq dtu results with proportions (keep the order by gene adj_pval):
  design_matrix <- design(dmFitObj)
  all_groups_proportion_df <- data.frame(proportions(dmFitObj))
  
  all_groups_proportion_df <- all_groups_proportion_df %>%
    mutate(across(where(is.numeric), ~ round(.x, 5)))
  
  
  all_drimSeqResults <- dmResultsStageR_adj %>% left_join(all_groups_proportion_df, by = join_by("txID" == "feature_id"))
  
  
  # sample names are taken directly from the design matrix, based on the relevant contrast:
  no_TXlabels <- all_groups_proportion_df %>% dplyr::select(-c("gene_id", "feature_id"))
  
  plus_group_oneSampleColname <- colnames(no_TXlabels[, design_matrix[, plus_group_name] == 1][1])
  minus_group_oneSampleColname <- colnames(no_TXlabels[, design_matrix[, minus_group_name] == 1][1])
  
  
  
  
  
  # Pull out one proportions of one sample from each relevant group (column names obtained above). Keep only those, and the updated DTU test results:
  relevent_dmsq_results <- all_drimSeqResults %>%
    dplyr::rename("plusGroupProp" = plus_group_oneSampleColname) %>%
    dplyr::rename("minusGroupProp" = minus_group_oneSampleColname) %>%
    dplyr::select(all_of(c(
      colnames(dmResultsStageR_adj),
      "plusGroupProp",
      "minusGroupProp"
    ))) %>%
    dplyr::mutate(propDiff = plusGroupProp - minusGroupProp) %>%
    dplyr::mutate(propSum = plusGroupProp + minusGroupProp) %>%
    dplyr::mutate(isNovel = if_else(str_detect(txID, "BambuTx"), "novelTX", "knownTX")) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(directionalDIUlabel = getDIUlabel(TX_05label, plusGroupProp, minusGroupProp)) # Define a function that determines a directional DTU label
  
  
  
  
  numberOfTXperGene <- relevent_dmsq_results %>%
    dplyr::group_by(geneID) %>%
    dplyr::mutate(numOfTXinGene = n()) %>%
    ggplot(aes(x = numOfTXinGene, fill = factor(directionalDIUlabel))) +
    geom_histogram(binwidth = 1, position = "stack") +
    coord_cartesian(xlim = c(0, 12)) +
    theme(legend.position = "buttom") +
    labs(title = "All TX tested")
  
  
  
  
  
  
  
  return(relevent_dmsq_results)
}
