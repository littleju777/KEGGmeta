#' Compare Metabolic Pathway Differences Between Two Conditions
#'
#' This function compares the metabolic pathway differences between two conditions
#' for a specified cell type within a Seurat object. The function performs
#' Wilcoxon tests for each pathway and saves the results to a CSV file.
#'
#' @param tam A Seurat object containing the single-cell RNA-seq data. 
#' Please make sure that the seurat@meta.data has "celltype" column stored cell type information, 
#' and "condition" column stored the condition information you want to compare with.
#' @param celltype A character string specifying the cell type to analyze.
#' @param cond1 A character string specifying the first condition. Defaults to "tumor".
#' @param cond2 A character string specifying the second condition. Defaults to "healthy".
#' @param output_dir A character string specifying the directory to save the results. 
#' Default to make a folder named "diffpath_cond1_vs_cond2" under current working directory.
#' @param make_volcano_plot A logical value indicating whether to generate a volcano plot. Defaults to TRUE.
#'
#' @return A data frame containing the results of the Wilcoxon tests, including p-values,
#' adjusted p-values, and differences in pathway scores.
#' @export
#'
#' @importFrom dplyr mutate case_when
#' @importFrom stringr str_c
#' @importFrom coin wilcox_test pvalue
#' @importFrom stats sd p.adjust
#' @importFrom ggplot2 ggplot geom_point scale_color_manual geom_vline geom_hline theme_classic
#' @importFrom ggrepel geom_text_repel
#'
#' @examples
#' for (j in unique(tam$celltype)) {DEpathway(tam, celltype = j, cond1 = "tumor", cond2 = "healthy")}

DEpathway <- function(tam, 
                      celltype = "T_cell", 
                      cond1 = "tumor", 
                      cond2 = "healthy", 
                      output_dir = getwd(),
                      make_volcano_plot = TRUE,
                      fig_wid = 7,
                      fig_wid = 5
) {
  # Subset the data for the specified cell type
  cell <- subset(tam, subset = celltype == celltype)
  tumor <- subset(cell, subset = condition == cond1)
  healthy <- subset(cell, subset = condition == cond2)
  
  # Combine meta and lipid scores for tumor and healthy conditions
  df_tumor <- rbind(as.data.frame(tumor@assays[["meta_UCell"]]@scale.data),
                    as.data.frame(tumor@assays[["lipid_UCell"]]@scale.data))
  
  df_healthy <- rbind(as.data.frame(healthy@assays[["meta_UCell"]]@scale.data),
                      as.data.frame(healthy@assays[["lipid_UCell"]]@scale.data))
  
  # Initialize vectors for storing results
  p_values <- numeric(nrow(df_tumor))
  changes <- numeric(nrow(df_tumor))
  
  # Perform Wilcoxon test for each pathway
  for (i in seq_len(nrow(df_tumor))) {
    sd_tumor <- sd(as.numeric(df_tumor[i, ]), na.rm = TRUE)
    sd_healthy <- sd(as.numeric(df_healthy[i, ]), na.rm = TRUE)
    
    if (!is.na(sd_tumor) && sd_tumor != 0 && !is.na(sd_healthy) && sd_healthy != 0) {
      combined_data <- data.frame(
        value = c(as.numeric(df_tumor[i, ]), as.numeric(df_healthy[i, ])),
        group = factor(c(rep("tumor", ncol(df_tumor)), rep("healthy", ncol(df_healthy))))
      )
      
      test_result <- wilcox_test(value ~ group, data = combined_data, distribution = "exact")
      p_values[i] <- pvalue(test_result)
      
      changes[i] <- (mean(as.numeric(df_tumor[i, ]), na.rm = TRUE) - mean(as.numeric(df_healthy[i, ]), na.rm = TRUE)) / mean(as.numeric(df_healthy[i, ]), na.rm = TRUE)
    } else {
      p_values[i] <- NA
      changes[i] <- NA
    }
  }
  
  adjusted_p_values <- p.adjust(p_values, method = "BH")
  
  # results to a data frame
  results <- data.frame(
    pathway = rownames(df_tumor),
    p_value = p_values,
    adj_p_value = adjusted_p_values,
    difference = changes
  )
  
  # Create the output directory if it doesn't exist
  folder_path <- file.path(output_dir, str_c("diffpath_", cond1, "_vs_", cond2))
  if (!dir.exists(folder_path)) {
    dir.create(folder_path, recursive = TRUE)
    message("Folder created: ", folder_path)
  }

  write.csv(results,str_c(folder_path,"/",celltype,".csv"))

  #make volcanoplots if make_volcano_plot = TRUE
  if (make_volcano_plot) {
    DE <- results %>% 
      mutate(group = case_when(
        adj_p_value < 0.05 & difference > 0.2 ~ "Upregulated",
        adj_p_value < 0.05 & difference < -0.2 ~ "Downregulated",
        TRUE ~ "No difference"
      ))
    
    topgenes <- DE %>% 
      filter(group %in% c("Upregulated", "Downregulated")) %>% 
      group_by(group) %>% 
      pull(pathway)
    
    p <- ggplot(DE, aes(x = difference, y = -log10(adj_p_value))) +
      geom_point(size = 1, aes(color = group)) +
      scale_color_manual(values = c("#4197D9", "grey", "orange")) +
      geom_vline(xintercept = c(-0.2, 0.2), linetype = 2) +
      geom_hline(yintercept = -log10(0.05), linetype = 2) +
      theme_classic() +
      geom_text_repel(data = subset(DE, pathway %in% topgenes), aes(label = pathway), 
                      min.segment.length = 0, seed = 42, box.padding = 0.5, max.overlaps = 25) +
      ggtitle(str_c("Difference of ", celltype, " between ", cond1, " and ", cond2))
    
    ggsave(p, filename = file.path(folder_path, str_c(celltype, ".pdf")), width = fig_wid, height = fig_ht)
  }
}
