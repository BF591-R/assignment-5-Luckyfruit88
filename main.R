library('tidyverse')
library('SummarizedExperiment')
library('DESeq2')
library('biomaRt')
library('testthat')
library('fgsea')

#' Function to generate a SummarizedExperiment object with counts and coldata
#' to use in DESeq2
#'
#' @param csv_path (str): path to the file verse_counts.tsv
#' @param metafile (str): path to the metadata sample_metadata.csv
#' @param selected_times (list): list of sample timepoints to use
#' 
#'   
#' @return SummarizedExperiment object with subsetted counts matrix
#'   and sample data. Ensure that the timepoints column used as input 
#'   to the model design has 'vP0' set as the reference factor level. Your 
#'   colData dataframe should have columns named samplename and timepoint.
#' @export
#'
#' @examples se <- make_se('verse_counts.tsv', 'sample_metadata.csv', c('vP0', 'vAd'))
make_se <- function(counts_csv, metafile_csv, selected_times) {
  counts_df <- readr::read_tsv(counts_csv,
                               show_col_types = FALSE,
                               name_repair    = "minimal")
  meta_df <- readr::read_csv(metafile_csv, show_col_types = FALSE)
  
  gene_ids <- counts_df[[1]]
  count_mat <- counts_df[, -1] %>% as.data.frame()
  rownames(count_mat) <- gene_ids
  
  meta_sub <- meta_df %>%
    dplyr::filter(timepoint %in% selected_times) %>%
    dplyr::select(samplename, timepoint)
  
  sample_order <- intersect(meta_sub$samplename, colnames(count_mat))
  meta_sub <- meta_sub %>%
    dplyr::filter(samplename %in% sample_order) %>%
    dplyr::mutate(
      timepoint = factor(timepoint, levels = selected_times)
    )
  
  if (!"vP0" %in% levels(meta_sub$timepoint)) {
    stop("vP0 must be included in selected_times so it can be used as the reference level.")
  }
  
  meta_sub$timepoint <- stats::relevel(meta_sub$timepoint, ref = "vP0")
  count_mat <- count_mat[, sample_order, drop = FALSE]
  count_mat <- as.matrix(count_mat)
  storage.mode(count_mat) <- "integer"
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = count_mat),
    colData = S4Vectors::DataFrame(meta_sub)
  )
  
  return(se)
}


#' Function that runs DESeq2 and returns a named list containing the DESeq2
#' results as a dataframe and the dds object returned by DESeq2
#'
#' @param se (obj): SummarizedExperiment object containing counts matrix and
#' coldata
#' @param design: the design formula to be used in DESeq2
#'
#' @return list with DESeqDataSet object after running DESeq2 and results from
#'   DESeq2 as a dataframe
#' @export
#'
#' @examples results <- return_deseq_res(se, ~ timepoint)
return_deseq_res <- function(se, design) {
  dds <- DESeq2::DESeqDataSet(se, design = design)
  
  # Pre-filter: keep only rows with at least 10 reads total
  keep <- rowSums(DESeq2::counts(dds)) >= 10
  dds  <- dds[keep, ]
  
  dds <- DESeq2::DESeq(dds)
  
  res <- DESeq2::results(dds)
  
  res_df <- as.data.frame(res) %>%
    tibble::rownames_to_column("genes") %>%
    tibble::as_tibble()
  
  return(list(
    results = res_df,
    dds     = dds
  ))
}



#' Function that takes the DESeq2 results dataframe, converts it to a tibble and
#' adds a column to denote plotting status in volcano plot. Column should denote
#' whether gene is either 1. Significant at padj < .10 and has a positive log
#' fold change, 2. Significant at padj < .10 and has a negative log fold change,
#' 3. Not significant at padj < .10. Have the values for these labels be UP,
#' DOWN, NS, respectively. The column should be named `volc_plot_status`. Ensure
#' that the column name for your rownames is called "genes". 
#'
#' @param deseq2_res (df): results from DESeq2 
#' @param padj_threshold (float): threshold for considering significance (padj)
#'
#' @return Tibble with all columns from DESeq2 results and one additional column
#'   labeling genes by significant and up-regulated, significant and
#'   downregulated, and not significant at padj < .10.
#'   
#' @export
#'
#' @examples labeled_results <- label_res(res, .10)
label_res <- function(deseq2_res, padj_threshold) {
  if ("genes" %in% colnames(deseq2_res)) {
    tib <- tibble::as_tibble(deseq2_res)
  } else {
    tib <- tibble::as_tibble(deseq2_res, rownames = "genes")
  }
  
  tib$volc_plot_status <- ifelse(
    !is.na(tib$padj) & tib$padj < padj_threshold & tib$log2FoldChange > 0, "UP",
    ifelse(
      !is.na(tib$padj) & tib$padj < padj_threshold & tib$log2FoldChange < 0, "DOWN",
      "NS"
    )
  )
  
  return(tib)
}






#' Function to plot the unadjusted p-values as a histogram
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one additional
#' column denoting status in volcano plot
#'
#' @return ggplot: a histogram of the raw p-values from the DESeq2 results
#' @export
#'
#' @examples pval_plot <- plot_pvals(labeled_results)
plot_pvals <- function(labeled_results) {
  p <- ggplot2::ggplot(labeled_results,
                       ggplot2::aes(x = pvalue)) +
    ggplot2::geom_histogram(bins = 30,
                            fill  = "lightblue",
                            color = "white",
                            na.rm = TRUE) +
    ggplot2::labs(
      title = "Histogram of raw pvalues obtained from DE analysis (vP0 vs. vAd)",
      x     = "pvalue",
      y     = "count"
    ) +
    ggplot2::theme_bw()
  return(p)
}



#' Function to plot the log2foldchange from DESeq2 results in a histogram
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one additional
#' column denoting status in volcano plot
#' @param padj_threshold (float): threshold for considering significance (padj)
#'
#' @return ggplot: a histogram of log2FC values from genes significant at padj 
#' threshold of 0.1
#' @export
#'
#' @examples log2fc_plot <- plot_log2fc(labeled_results, .10)
plot_log2fc <- function(labeled_results, padj_threshold) {
  sig_genes <- labeled_results %>%
    dplyr::filter(!is.na(padj) & padj < padj_threshold)
  
  p <- ggplot2::ggplot(sig_genes,
                       ggplot2::aes(x = log2FoldChange)) +
    ggplot2::geom_histogram(bins = 30,
                            fill  = "lightblue",
                            color = "white") +
    ggplot2::labs(
      title = "Histogram of Log2FoldChanges for DE Genes (vP0 vs. vAd)",
      x     = "log2FoldChange",
      y     = "count"
    ) +
    ggplot2::theme_bw()
  return(p)
}


#' Function to make scatter plot of normalized counts for top ten genes ranked
#' by ascending padj
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one
#'   additional column denoting status in volcano plot
#' @param dds_obj (obj): The object returned by running DESeq (dds) containing
#' the updated DESeqDataSet object with test results
#' @param num_genes (int): Number of genes to plot
#'
#' @return ggplot: a scatter plot with the normalized counts for each sample for
#' each of the top ten genes ranked by ascending padj
#' @export
#'
#' @examples norm_counts_plot <- scatter_norm_counts(labeled_results, dds, 10)
scatter_norm_counts <- function(labeled_results, dds_obj, num_genes) {
  top_genes <- labeled_results %>%
    dplyr::filter(!is.na(padj)) %>%
    dplyr::arrange(padj) %>%
    dplyr::slice_head(n = num_genes) %>%
    dplyr::pull(genes)
  
  norm_counts <- DESeq2::counts(dds_obj, normalized = TRUE)
  norm_sub <- norm_counts[top_genes, , drop = FALSE] %>%
    as.data.frame() %>%
    tibble::rownames_to_column("genes") %>%
    tidyr::pivot_longer(-genes,
                        names_to  = "samplename",
                        values_to = "norm_counts")
  
  sample_info <- as.data.frame(SummarizedExperiment::colData(dds_obj)) %>%
    tibble::as_tibble() %>%
    dplyr::select(samplename, timepoint)
  
  norm_sub <- dplyr::left_join(norm_sub, sample_info, by = "samplename")
  
  p <- ggplot2::ggplot(norm_sub,
                       ggplot2::aes(x     = genes,
                                    y     = log10(norm_counts),
                                    color = samplename)) +
    ggplot2::geom_point(size = 2) +
    ggplot2::labs(
      title  = paste0("Plot of Log10(normalized counts) for top ten DE genes"),
      x      = NULL,
      y      = "log10(norm_counts)",
      color  = "samplenames"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, size = 7))
  return(p)
}




#' Function to generate volcano plot from DESeq2 results
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one
#'   additional column denoting status in volcano plot
#'
#' @return ggplot: a scatterplot (volcano plot) that displays log2foldchange vs
#'   -log10(padj) and labeled by status
#' @export
#'
#' @examples volcano_plot <- plot_volcano(labeled_results)
#' 
plot_volcano <- function(labeled_results) {
  p <- ggplot2::ggplot(labeled_results,
                       ggplot2::aes(x     = log2FoldChange,
                                    y     = -log10(padj),
                                    color = volc_plot_status)) +
    ggplot2::geom_point(alpha = 0.6, size = 1.5, na.rm = TRUE) +
    ggplot2::geom_hline(yintercept = -log10(0.1),
                        linetype   = "dashed",
                        color      = "black") +
    ggplot2::scale_color_manual(
      values = c("UP"   = "blue",
                 "DOWN" = "red",
                 "NS"   = "green3")
    ) +
    ggplot2::labs(
      title = "Volcano plot of DESeq2 differential expression results (vP0 vs. vAd)",
      x     = "log2FoldChange",
      y     = "-log10(padj)"
    ) +
    ggplot2::theme_bw()
  return(p)
}



#' Function to generate a named vector ranked by log2FC descending
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one
#'   additional column denoting status in volcano plot
#' @param id2gene_path (str): Path to the file containing the mapping of
#' ensembl IDs to MGI symbols
#'
#' @return Named vector with gene symbols as names, and log2FoldChange as values
#' ranked in descending order
#' @export
#'
#' @examples rnk_list <- make_ranked_log2fc(labeled_results, 'data/id2gene.txt')

make_ranked_log2fc <- function(labeled_results, id2gene_path) {
  id2gene <- readr::read_tsv(id2gene_path,
                             col_names      = c("ensembl_id", "gene_symbol"),
                             show_col_types = FALSE)
  
  ranked <- labeled_results %>%
    dplyr::inner_join(id2gene, by = c("genes" = "ensembl_id")) %>%
    dplyr::filter(!is.na(log2FoldChange)) %>%
    dplyr::arrange(dplyr::desc(log2FoldChange))
  
  rnk_list <- ranked$log2FoldChange
  names(rnk_list) <- ranked$gene_symbol
  
  return(rnk_list)
}




#' Function to run fgsea with arguments for min and max gene set size
#'
#' @param gmt_file_path (str): Path to the gene sets of interest in GMT format
#' @param rnk_list (named vector): Named vector generated previously with gene 
#' symbols and log2Fold Change values in descending order
#' @param min_size (int): Minimum number of genes in gene sets to be allowed
#' @param max_size (int): Maximum number of genes in gene sets to be allowed
#'
#' @return Tibble of results from running fgsea
#' @export
#'
#' @examples fgsea_results <- run_fgsea('data/m2.cp.v2023.1.Mm.symbols.gmt', rnk_list, 15, 500)
run_fgsea <- function(gmt_file_path, rnk_list, min_size, max_size) {
  pathways <- fgsea::gmtPathways(gmt_file_path)
  
  fgsea_res <- fgsea::fgsea(
    pathways = pathways,
    stats    = rnk_list,
    minSize  = min_size,
    maxSize  = max_size
  )
  
  fgsea_tib <- tibble::as_tibble(fgsea_res) %>%
    dplyr::arrange(padj)
  
  return(fgsea_tib)
}


#' Function to plot top ten positive NES and top ten negative NES pathways
#' in a barchart
#'
#' @param fgsea_results (tibble): the fgsea results in tibble format returned by
#'   the previous function
#' @param num_paths (int): the number of pathways for each direction (top or
#'   down) to include in the plot. Set this at 10.
#'
#' @return ggplot with a barchart showing the top twenty pathways ranked by positive
#' and negative NES
#' @export
#'
#' @examples fgsea_plot <- top_pathways(fgsea_results, 10)
top_pathways <- function(fgsea_results, num_paths) {
  top_pos <- fgsea_results %>%
    dplyr::filter(NES > 0) %>%
    dplyr::arrange(dplyr::desc(NES)) %>%
    dplyr::slice_head(n = num_paths)
  
  top_neg <- fgsea_results %>%
    dplyr::filter(NES < 0) %>%
    dplyr::arrange(NES) %>%
    dplyr::slice_head(n = num_paths)
  
  top_both <- dplyr::bind_rows(top_pos, top_neg) %>%
    dplyr::mutate(
      pathway_label = stringr::str_trunc(
        gsub("_", " ", pathway),
        width = 50, ellipsis = "..."),
      pathway_label = forcats::fct_reorder(pathway_label, NES),
      direction = ifelse(NES > 0, "positive", "negative")
    )
  
  p <- ggplot2::ggplot(top_both,
                       ggplot2::aes(x    = NES,
                                    y    = pathway_label,
                                    fill = direction)) +
    ggplot2::geom_col() +
    ggplot2::scale_fill_manual(
      values = c("positive" = "red", "negative" = "blue"),
      guide  = "none"
    ) +
    ggplot2::labs(
      title = "fgsea results for Hallmark MSigDB gene se",
      x     = "Normalized Enrichment Score (NES)",
      y     = NULL
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 7))
  
  return(p)
}




