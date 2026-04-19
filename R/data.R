#' Cell Name Mapping Table
#'
#' A data frame mapping internal CIBERSORT cell-type identifiers
#' to human-readable abbreviated cell type names.
#'
#' @format A data frame with 24 rows and 2 columns:
#' \describe{
#'   \item{Var2}{Internal cell-type identifier used in deconvolution
#'     reference matrices.}
#'   \item{Cell.Type}{Abbreviated human-readable cell type name
#'     (e.g., \code{"Macro M1"}, \code{"CD4 Tcm"}).}
#' }
#'
#' @source Custom curation for TCGA immune deconvolution analyses.
#'
#' @examples
#' data(cellname)
#' head(cellname)
"cellname"


#' Comprehensive Secreted Protein Gene List
#'
#' A curated list of genes encoding secreted proteins (secretome),
#' compiled from multiple databases including HPA and VerSeDa.
#'
#' @format A data frame with 4919 rows and 3 columns:
#' \describe{
#'   \item{symbol}{HGNC gene symbol.}
#'   \item{type}{Functional category of the secreted protein
#'     (e.g., \code{"predicted secreted protein"},
#'     \code{"growth_factors_&_receptors"}).}
#'   \item{database}{Source database
#'     (e.g., \code{"HPA"}, \code{"VerSeDa"}, \code{"Online resource"}).}
#' }
#'
#' @source Human Protein Atlas (HPA), VerSeDa, and curated online resources.
#'
#' @examples
#' data(secretome_genes)
#' head(secretome_genes)
"secretome_genes"


#' CIBERSORT Reference Matrix
#'
#' A reference expression matrix for immune cell deconvolution using
#' the CIBERSORT algorithm via \code{tidybulk::deconvolve_cellularity}.
#'
#' @format A matrix or data frame with genes in rows and reference
#'   cell types in columns.
#'
#' @source Custom reference matrix for TCGA analyses.
#'
#' @examples
#' data(my_ref)
#' dim(my_ref)
"my_ref"
