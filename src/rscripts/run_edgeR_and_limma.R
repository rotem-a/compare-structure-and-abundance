library(limma)
library(edgeR)
library(bambu)

#'  Obtain a SummarisedExperiment from bambu, add libsize and MDS dimensional reductions for each assay
#'
#' bambu returns its results as a SummarisedExperiment R object which includes transcript counts and full length read counts as separate assays.
#' Gene counts are obtained from it by using bambu::getGeneCounts() function, which effectively adds the transcript counts for each gene.
#'
#' This function uses edgeR to obtain library size and cpm counts matrix for each assay in the bambu output (SummarisedExperiment):
#' TX_counts, fullLength_counts and including gene counts.
#'
#' @param bambu_se bambu output: an object of class "SummarisedExperiment".
#' @param update_genes Boolean. Does not calculate for gene counts unless True.
#' @param counts_only Boolean. If TRUE, calculate only for transcript counts (exclude gene counts and full-length read counts).
#' @param calc_mds Boolean. If FALSE, only return library size without dimensional reduction.
#' @param top Integer. Number of features to include in the MDS reductions.
#'
#' @returns Returns the same bambu_output object with dimensional reduction and library sizes in column data
#' @export
#' @importFrom edgeR plotMDS
#' @examples
add_dimReductions <- function(bambu_se,
                              update_genes = TRUE,
                              counts_only = FALSE,
                              calc_mds = TRUE,
                              top = 500) {
  if (update_genes) {
    rel_assays <- list("counts", "fullLengthCounts", "gene_counts")
  } else {
    rel_assays <- list("counts", "fullLengthCounts")
  }
  
  if (counts_only) {
    rel_assays <- list("counts")
  }
  
  
  for (name in rel_assays) {
    print(name)
    
    if (name == "gene_counts") {
      genes_se <- transcriptToGeneExpression(bambu_se)
      current_counts <- SE2DGEList(genes_se)
    } else {
      current_assay <- assay(bambu_se, name)
      
      
      # (make sure this is doing what I think it is doing.....)
      current_counts <- DGEList(round(current_assay))
      
      cpm_str <- paste(name, "cpm", sep = "_")
      assay(bambu_se, cpm_str) <- edgeR::cpm(current_counts)
      
      lcpm_str <- paste(name, "lcpm", sep = "_")
      assay(bambu_se, lcpm_str) <- edgeR::cpm(current_counts, log = TRUE)
    }
    # otherwise, move back to round only counts explicitly:
    # current_counts$counts <- round(current_counts$counts )
    
    
    libsize_str <- paste(name, "libsize", sep = "_")
    colData(bambu_se)[, libsize_str] <- current_counts$samples$lib.size
    
    if (calc_mds) {
      mds_reduction <- edgeR_MDS(current_counts, top = top)
      
      mds_str <- paste(name, "mds", sep = "_")
      
      colData(bambu_se)[, paste(mds_str, "1", sep = "_")] <- mds_reduction$mds12$x
      colData(bambu_se)[, paste(mds_str, "2", sep = "_")] <- mds_reduction$mds12$y
      colData(bambu_se)[, paste(mds_str, "3", sep = "_")] <- mds_reduction$mds34$x
      colData(bambu_se)[, paste(mds_str, "4", sep = "_")] <- mds_reduction$mds34$y
    }
  }
  
  return(bambu_se)
}

edgeR_MDS <- function(counts_data, top = 500) {
  # Obtains a DGElist and returns MDS (1-4 first dim) dimensional reduction
  
  
  counts_data <- calcNormFactors(counts_data)
  
  MDS1_2 <- plotMDS(counts_data, top = top, plot = FALSE)
  MDS3_4 <- plotMDS(counts_data,
                    top = top,
                    dim = c(3, 4),
                    plot = FALSE)
  returned_list <- list(mds12 = MDS1_2, mds34 = MDS3_4)
  
  return(returned_list)
}


#' Streamline DE testing using limma-voom
#'
#' @param counts
#' @param design
#' @param contrasts
#' @param lfc
#' @param min_count
#' @param min_total_count
#' @param return_fit
#' @param run_single_test
#'
#' @returns
#' @export
#' @import limma
#' @importFrom edgeR filterByExpr
#' @examples
de_filter_fit_test <- function(counts,
                               design,
                               contrasts,
                               lfc = 0.5,
                               min_count = 10,
                               min_total_count = 20,
                               return_fit = FALSE,
                               run_single_test = "") {
  # Filter out lowly expressed features:
  keep.exprs <- filterByExpr(
    counts,
    design = design,
    min.count = min_count,
    min.total.count = min_total_count
  )
  
  counts <- counts[keep.exprs, ]
  
  # Fit a linear model with limma-voom
  counts_matrix_voom <- voomWithQualityWeights(counts, design, plot = FALSE)
  lm_fit <- lmFit(counts_matrix_voom, design)
  
  
  vfit <- contrasts.fit(lm_fit, contrasts = contrasts)
  efit <- eBayes(vfit)
  tfit <- treat(efit, lfc = lfc)
  
  
  
  dt_lfc <- decideTests(
    tfit,
    method = "separate",
    adjust.method = "BH",
    p.value = 0.05,
    
  )
  decide_results <- as.data.frame(dt_lfc)
  decide_results$feature_id <- rownames(decide_results)
  
  if (return_fit) {
    decide_results <- list(decide_results = decide_results, efit = efit)
    # print("=============de_filter_fit_test return ===============")
    # print(colnames(efit$coefficients))
  }
  
  
  
  return(decide_results)
}
