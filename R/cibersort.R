#' Run CIBERSORT Immune Cell Deconvolution
#'
#' Aggregates transcript counts and performs immune cell deconvolution
#' using \code{tidybulk::deconvolve_cellularity}.
#'
#' @param x A \code{tibble} of expression data.
#' @param c1 Character string. Name of the sample column.
#' @param c2 Character string. Name of the gene/transcript column.
#' @param c3 Character string. Name of the abundance column.
#' @param ref A reference matrix for deconvolution.
#'   Use \code{surVisualiser::my_ref} for the bundled reference.
#' @param meth Character string. Deconvolution method
#'   (passed to \code{tidybulk::deconvolve_cellularity}).
#' @param act Character string. Action parameter
#'   (passed to \code{tidybulk::deconvolve_cellularity}).
#'
#' @return A \code{tibble} with estimated cell-type proportions.
#'
#' @importFrom data.table setDT
#' @importFrom tibble as_tibble
#' @importFrom rlang sym !!
#'
#' @export
#'
#' @examples
#' \dontrun{
#' deconv_result <- CIBERSORT(
#'   expr_data,
#'   c1 = "sample",
#'   c2 = "symbol",
#'   c3 = "raw_count",
#'   ref = surVisualiser::my_ref,
#'   meth = "cibersort",
#'   act = "get"
#' )
#' }
CIBERSORT <- function(x, c1, c2, c3, ref, meth, act) {
    tibble::as_tibble(
        data.table::setDT(x)[,
            lapply(.SD, sum),
            by = c(c1, c2),
            .SDcols = c3
        ]
    ) |>
        tidybulk::deconvolve_cellularity(
            .sample = !!rlang::sym(c1),
            .transcript = !!rlang::sym(c2),
            .abundance = !!rlang::sym(c3),
            reference = ref,
            method = meth,
            action = act
        )
}
