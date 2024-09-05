#' Compute Meta and Lipid Scores and Store Them in the Seurat Object
#'
#' This function computes meta and lipid scores for a given Seurat object using the irGSEA.score method.
#' The scores are based on the meta and fatty genesets included in the package.
#' The function then stores the computed scores in new assays within the Seurat object.
#'
#' @param tam A Seurat object containing the single-cell RNA-seq data.
#' @param species Optional. A character string specifying the name of the dataset containing the meta and fatty genesets. 
#' Defaults to "mice".
#' @param assay_name A character string specifying the assay to use. Defaults to "RNA".
#' @param slot_name A character string specifying the slot to use. Defaults to "data".
#'
#' @importFrom irGSEA irGSEA.score
#' @importFrom Seurat CreateAssayObject GetAssayData
#' @return A Seurat object with added assays for meta and lipid scores.
#' @export
#'
#' @examples
#' tam <- get_metascore(tam)
get_metascore <- function(tam, 
                          assay_name = "RNA", 
                          slot_name = "data",
                          species = "mice"
  ) {
  # Load the geneset data from the package
  data_name = str_c(species,"_meta_lipid_geneset")
  data(data_name, package = "KEGG")

  # Compute meta scores
  meta <- irGSEA.score(object = tam, assay = assay_name, slot = slot_name,
                       custom = TRUE, geneset = meta.geneset, method = method,
                       kcdf = 'Gaussian')

  # Compute lipid scores
  lipid <- irGSEA.score(object = tam, assay = assay_name, slot = slot_name,
                        custom = TRUE, geneset = fatty.geneset, method = method,
                        kcdf = 'Gaussian')

  # Store the meta_UCell data in the Seurat object
  tam[["meta_UCell"]] <- CreateAssayObject(counts = GetAssayData(object = meta@assays[["UCell"]], slot = 'scale.data'))
  tam[["meta_UCell"]]@scale.data <- meta@assays[["UCell"]]@scale.data

  # Store the lipid_UCell data in the Seurat object
  tam[["lipid_UCell"]] <- CreateAssayObject(counts = GetAssayData(object = lipid@assays[["UCell"]], slot = 'scale.data'))
  tam[["lipid_UCell"]]@scale.data <- lipid@assays[["UCell"]]@scale.data

  return(tam)
}
