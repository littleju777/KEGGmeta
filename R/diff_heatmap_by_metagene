#' Differential Heatmap by Gene
#'
#' Generates heatmaps of differential gene expression between two groups
#' for a given set of genes in a Seurat object.
#'
#' @param seurat_obj A Seurat object containing single-cell RNA-seq data.
#' @param geneset A named list of gene sets, where each element is a vector of gene names.
#' @param group_col The column in the metadata indicating the group (default: "condition").
#' @param group1 The name of the first group (default: "Fish").
#' @param group2 The name of the second group (default: "Soy").
#' @param celltype_col The column in the metadata indicating cell types (default: "cell_type").
#' @export
diff_heatmap_by_gene <- function(seurat_obj, geneset, group_col = "condition", group1 = "Fish", group2 = "Soy", celltype_col = "cell_type") {
  # Function body as given
  # Create the output directory if it does not exist
  output_dir <- paste0("/data/msun/projects/Stephen/PDAC_scRNA/04_metabolism/heatmap_by_genes/", group1, "_vs_", group2, "/")
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Loop over each pathway in the geneset
  for (pathway in names(geneset)) {
    
    DefaultAssay(seurat_obj) <- "RNA"
    Idents(seurat_obj) <- celltype_col
    genes_to_use <- geneset[[pathway]]
    mat <- seurat_obj[["RNA"]]@data[rownames(seurat_obj[["RNA"]]@data) %in% genes_to_use, ]
    genes_in_data <- rownames(mat)
    
    # Initialize a matrix to store log2FC results
    cell_types <- unique(seurat_obj@meta.data[[celltype_col]])
    mat_log2fc <- matrix(NA, nrow = length(genes_in_data), ncol = length(cell_types))
    rownames(mat_log2fc) <- genes_in_data
    colnames(mat_log2fc) <- cell_types
    
    # Initialize a matrix to store exact p-values
    p_values <- matrix(NA, nrow = length(genes_in_data), ncol = length(cell_types))
    rownames(p_values) <- genes_in_data
    colnames(p_values) <- cell_types
    
    # Loop over each cell type
    for (cell_type in cell_types) {
      # Subset data for the current cell type
      cell_subset <- WhichCells(seurat_obj, idents = cell_type)
      group1_cells <- cell_subset[seurat_obj@meta.data[cell_subset, group_col] == group1]
      group2_cells <- cell_subset[seurat_obj@meta.data[cell_subset, group_col] == group2]
      
      # Calculate mean expression for each group and cell type
      mean_group1 <- rowMeans(mat[, group1_cells, drop = FALSE], na.rm = TRUE)
      mean_group2 <- rowMeans(mat[, group2_cells, drop = FALSE], na.rm = TRUE)
      
      # Calculate log2 fold change
      log2fc <- log2((mean_group1 + 1) / (mean_group2 + 1))
      
      # Store log2FC in matrix
      mat_log2fc[, cell_type] <- log2fc
      
      # Perform Wilcoxon rank-sum test using the `coin` package for each gene and store exact p-values
      for (gene in genes_in_data) {
        if (gene %in% rownames(mat)) {
          group1_values <- mat[gene, group1_cells, drop = FALSE]
          group2_values <- mat[gene, group2_cells, drop = FALSE]
          
          # Remove NAs and check for zero variance
          group1_values <- group1_values[!is.na(group1_values)]
          group2_values <- group2_values[!is.na(group2_values)]
          
          if (length(group1_values) > 0 && length(group2_values) > 0 && 
              var(group1_values) > 0 && var(group2_values) > 0) {
            # Use the `coin` package to compute p-values (asymptotic distribution handles ties better)
            data <- data.frame(values = c(group1_values, group2_values),
                               group = factor(c(rep(group1, length(group1_values)), rep(group2, length(group2_values)))))
            
            # Assign p-value to correct location in the matrix using numeric indices
            row_idx <- which(rownames(p_values) == gene)
            col_idx <- which(colnames(p_values) == cell_type)
            p_values[row_idx, col_idx] <- pvalue(wilcox_test(values ~ group, data = data, distribution = "asymptotic"))
          }
        }
      }
    }
    
    # Define the color function for the heatmap
    col_fun <- circlize::colorRamp2(c(min(mat_log2fc, na.rm = TRUE), 0, max(mat_log2fc, na.rm = TRUE)), c("blue", "white", "red"))
    
    # Define the output file path
    filepath <- paste0(output_dir, pathway, "_heatmap.pdf")
    
    # Draw the heatmap
    pdf(file = filepath, width = 8, height = max(5, nrow(mat_log2fc) * 0.3))
    p <- Heatmap(
      mat_log2fc,
      name = "Log2FC",
      col = col_fun,
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      show_row_names = TRUE,
      show_column_names = TRUE,
      column_title = pathway,
      row_names_gp = grid::gpar(fontsize = 8),
      column_names_gp = grid::gpar(fontsize = 8),
      cell_fun = function(j, i, x, y, width, height, fill) {
        # Annotate significant genes with stars based on exact p-values
        if (!is.na(p_values[i, j])) {
          if (p_values[i, j] < 0.001) {
            grid.text("***", x, y, gp = grid::gpar(fontsize = 10))
          } else if (p_values[i, j] < 0.01) {
            grid.text("**", x, y, gp = grid::gpar(fontsize = 10))
          } else if (p_values[i, j] < 0.05) {
            grid.text("*", x, y, gp = grid::gpar(fontsize = 10))
          }
        }
      }
    )
    print(p)
    dev.off()
  }
}
