#' Read, Extract, and Merge Specified Columns from Multiple CSV Files
#'
#' This function reads CSV files from a specified folder, extracts a specific column from each file,
#' and merges these columns into a single data frame.
#'
#' @param compare A character string specifying the comparison identifier (e.g., "tumor_vs_healthy"). Defaults to "tumor_vs_healthy".
#' please make sure the comapre string is the same as the folder name!
#' @param column_name A character string specifying the column to extract from each CSV file.
#'
#' @return A data frame containing the merged data from the specified column across all CSV files.
#' @export
#'
#' @importFrom dplyr select all_of
#' @importFrom stringr str_c
#' @importFrom purrr map
#' @importFrom utils read.csv
#' @importFrom tools file_path_sans_ext
#'
#' @examples
#' merged_data <- read_extract_merge(compare = "Fish_vs_LFD", column_name = "difference")
read_extract_merge <- function(compare = "tumor_vs_healthy", column_name, df_folder = NULL) {
  
  # List all CSV files in the folder
  csv_files <- list.files(path = df_folder, pattern = "\\.csv$", full.names = TRUE)
  
  # Function to read a CSV file and extract the specified column
  read_and_extract_column <- function(file) {
    data <- read.csv(file)
    if (column_name %in% names(data)) {
      return(data %>% select(all_of(column_name)))
    } else {
      stop(paste("Column", column_name, "not found in", file))
    }
  }
  
  # Read, extract, and merge the specified column from all CSV files
  merged_data <- map(csv_files, read_and_extract_column) %>%
    bind_cols()
  
  # Extract pathways from the first file and set as rownames
  pathways <- read.csv(csv_files[[1]])$pathway
  rownames(merged_data) <- pathways
  
  # Rename the columns based on the filenames
  file_names <- basename(csv_files)
  file_names <- tools::file_path_sans_ext(file_names) # Get file names without extension
  names(merged_data) <- file_names
  
  return(merged_data)
}

#' Generate a Heatmap of Differential Expression Data
#'
#' This function generates a heatmap of differential expression data and highlights significant p-values.
#' The heatmap is saved as a PDF file.
#'
#' @param compare A character string specifying the comparison identifier (e.g., "Fish_vs_LFD"). Defaults to "Fish_vs_LFD".
#'
#' @return A heatmap is generated and saved as a PDF file. The function does not return any objects.
#' @export
#'
#' @importFrom dplyr mutate all_of
#' @importFrom ComplexHeatmap Heatmap draw
#' @importFrom circlize colorRamp2
#' @importFrom grid unit grid.text
#' @importFrom grDevices pdf dev.off
#' @importFrom stringr str_c
#'
#' @examples
#' diff_heatmap(compare = "Fish_vs_LFD")
diff_heatmap <- function(compare = "Fish_vs_LFD",df_folder = getwd()) {
  # Merge the difference and adjusted p-value data across cell types
  results_path = str_c(df_folder, "/diffpath_", compare)
  df.diff <- read_extract_merge(compare = compare, "difference", df_folder = results_path)
  df.adjp <- read_extract_merge(compare = compare, "adj_p_value", df_folder = results_path)
  
  # Check if the dimensions and row/column names match
  if (!all(dim(df.diff) == dim(df.adjp))) {
    stop("Dimensions of expression data and p-value data do not match.")
  }
  if (!all(rownames(df.diff) == rownames(df.adjp)) || !all(colnames(df.diff) == colnames(df.adjp))) {
    stop("Row and column names of expression data and p-value data do not match.")
  }
  
  # Replace NA and Inf values
  df.diff[is.na(df.diff)] <- 0
  df.adjp[is.na(df.adjp)] <- 1
  df.diff <- as.matrix(df.diff)
  df.adjp <- as.matrix(df.adjp)
  max_finite_value <- max(df.diff[is.finite(df.diff)], na.rm = TRUE)
  df.diff[is.infinite(df.diff)] <- max_finite_value
  
  # Order the pathways
  ordered_indices <- order(rownames(df.diff))
  ordered_diff <- df.diff[ordered_indices, ]
  ordered_adjp <- df.adjp[ordered_indices, ]
  
  # color for the heatmap
  col_fun <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
  
  # Define the file path for saving the heatmap
  filepath <- str_c(results_path, "/heatmap_", compare, ".pdf")
  
  # Generate and save the heatmap
  pdf(file = filepath, width = 8, height = 12)
  draw(Heatmap(ordered_diff,  
               cluster_rows = FALSE, 
               cluster_columns = FALSE, 
               show_row_dend = FALSE,
               na_col = "black", 
               col = col_fun,
               heatmap_legend_param = list(legend_direction = "horizontal"),
               row_names_gp = grid::gpar(fontsize = 10),
               cell_fun = function(j, i, x, y, w, h, fill) {
                 if (ordered_adjp[i, j] < 0.001) {
                   grid.text("***", x, y)
                 } else if (ordered_adjp[i, j] < 0.01) {
                   grid.text("**", x, y)
                 } else if (ordered_adjp[i, j] < 0.05) {
                   grid.text("*", x, y)
                 }
               }), 
        heatmap_legend_side = "bottom",
        padding = unit(c(2, 2, 2, 40), "mm"))
  dev.off()
}




#' Generate a Heatmap for a Subset of Pathways
#'
#' This function generates a heatmap for a subset of pathways based on differential expression
#' and adjusted p-values. The subset of pathways is provided by the user, and the heatmap is saved
#' as a PDF file.
#'
#' @param compare A character string specifying the comparison identifier (e.g., "Fish_vs_LFD"). Defaults to "Fish_vs_LFD".
#' @param subset_path A character vector containing the names of the pathways to include in the heatmap.
#' @param df_folder A character string specifying the directory containing the differential expression data.
#' Defaults to the current working directory.
#'
#' @return The heatmap is generated and saved as a PDF file. The function does not return any objects.
#' @export
#'
#' @importFrom dplyr select all_of
#' @importFrom ComplexHeatmap Heatmap draw
#' @importFrom circlize colorRamp2
#' @importFrom grid unit grid.text
#' @importFrom grDevices pdf dev.off
#' @importFrom stringr str_c
#'
#' @examples
#' subset_path <- c("Pathway1", "Pathway2", "Pathway3")
#' diff_heatmap_subset(compare = "Fish_vs_LFD", subset_path = subset_path)
diff_heatmap_subset <- function(compare = "Fish_vs_LFD", subset_path, df_folder = getwd(), fig_ht = 5, fig_wid = 7) {
  # Read and merge the difference and adjusted p-value data
  results_path = str_c(df_folder, "/diffpath_", compare)
  df.diff <- read_extract_merge(compare = compare, "difference", df_folder = results_path)
  df.adjp <- read_extract_merge(compare = compare, "adj_p_value", df_folder = results_path)
  
  # Check if the dimensions and row/column names match
  if (!all(dim(df.diff) == dim(df.adjp))) {
    stop("Dimensions of expression data and p-value data do not match.")
  }
  if (!all(rownames(df.diff) == rownames(df.adjp)) || !all(colnames(df.diff) == colnames(df.adjp))) {
    stop("Row and column names of expression data and p-value data do not match.")
  }
  
  # Subset the data to include only the specified pathways
  subset.diff <- df.diff[rownames(df.diff) %in% subset_df, ]
  subset.adjp <- df.adjp[rownames(df.adjp) %in% subset_df, ]
  
  # Replace NA and Inf values in the expression data and p-value data
  subset.diff[is.na(subset.diff)] <- 0
  subset.adjp[is.na(subset.adjp)] <- 1
  subset.diff <- as.matrix(subset.diff)
  subset.adjp <- as.matrix(subset.adjp)
  max_finite_value <- max(subset.diff[is.finite(subset.diff)], na.rm = TRUE)
  subset.diff[is.infinite(subset.diff)] <- max_finite_value
  
  # Order the pathways
  ordered_indices <- order(rownames(subset.diff))
  ordered_diff <- subset.diff[ordered_indices, ]
  ordered_adjp <- subset.adjp[ordered_indices, ]
  
  # Define the color function for the heatmap
  col_fun <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
  
  # Define the file path for saving the heatmap
  filepath <- str_c(results_path, "/subheatmap_", compare, ".pdf")
  
  # Generate and save the heatmap
  pdf(file = filepath, width = fig_wid, height = fig_ht)
  draw(Heatmap(ordered_diff,  
               cluster_rows = FALSE, 
               cluster_columns = FALSE, 
               show_row_dend = FALSE,
               na_col = "black", 
               col = col_fun,
               heatmap_legend_param = list(legend_direction = "horizontal"),
               row_names_gp = grid::gpar(fontsize = 10),
               cell_fun = function(j, i, x, y, w, h, fill) {
                 if (ordered_adjp[i, j] < 0.001) {
                   grid.text("***", x, y)
                 } else if (ordered_adjp[i, j] < 0.01) {
                   grid.text("**", x, y)
                 } else if (ordered_adjp[i, j] < 0.05) {
                   grid.text("*", x, y)
                 }
               }), 
       heatmap_legend_side = "bottom",
       padding = unit(c(2, 2, 2, 40), "mm"))
  dev.off()
}
