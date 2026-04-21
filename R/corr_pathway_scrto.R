#' Pearson Correlation Between Hallmark Scores and Secretome Expression
#'
#' Computes the Pearson correlation between per-sample hallmark
#' enrichment scores (typically from
#' \code{\link{gsea_test_scrto_gene}}) and per-sample expression of a
#' user-specified gene set (typically secretome genes). Samples are
#' aligned by their \code{sample} identifier and only samples present
#' in both tables are retained.
#'
#' @details
#' Because \code{\link{gsea_test_scrto_gene}} removes the supplied
#' gene list from the expression matrix \emph{before} ranking, the
#' returned hallmark scores do not share genes with the gene set
#' correlated against here. Correlations produced by this function
#' are therefore not inflated by shared membership between the
#' hallmark gene sets and the queried gene list, and can be
#' interpreted as a circular-analysis-safe measure of association.
#'
#' @param hallmark_scores A \code{tibble} returned by
#'   \code{\link{gsea_test_scrto_gene}} (or any tibble with columns
#'   \code{sample}, \code{ID}, and the score column given by
#'   \code{score_col}).
#' @param expression A \code{tibble} of per-sample gene expression
#'   with at least columns \code{sample}, \code{symbol}, and the
#'   abundance column given by \code{expression_col}.
#' @param secretome_symbols Character vector of gene symbols to
#'   correlate against. Only rows of \code{expression} matching these
#'   symbols are used.
#' @param score_col Character string. Name of the score column in
#'   \code{hallmark_scores} to correlate. Defaults to \code{"NES"}
#'   (normalised enrichment score).
#' @param expression_col Character string. Name of the abundance
#'   column in \code{expression}. Defaults to
#'   \code{"raw_count_scaled_adjusted"}.
#'
#' @return A \code{tibble} with columns \code{pathway} (hallmark ID),
#'   \code{gene} (secretome gene symbol), \code{value} (Pearson
#'   correlation coefficient), and \code{p.value}.
#'
#' @seealso \code{\link{gsea_test_scrto_gene}},
#'   \code{\link{three_heatmaps_of_pth_scrto_TN}}.
#'
#' @importFrom dplyr filter select inner_join all_of
#' @importFrom tidyr pivot_wider pivot_longer
#' @importFrom tibble column_to_rownames rownames_to_column
#' @importFrom Hmisc rcorr
#' @importFrom rlang .data
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data(secretome_genes)
#' hallmarks <- gsea_test_scrto_gene(
#'     data = expr_tbl,
#'     rank_by = "raw_count_scaled_adjusted",
#'     excluded_genes = secretome_genes$symbol
#' )
#' corr <- corr_pathway_scrto(
#'     hallmark_scores = hallmarks,
#'     expression = expr_tbl,
#'     secretome_symbols = secretome_genes$symbol
#' )
#' }
corr_pathway_scrto <- function(hallmark_scores,
                               expression,
                               secretome_symbols,
                               score_col = "NES",
                               expression_col =
                                   "raw_count_scaled_adjusted") {

    ####*Input validation*####
    stopifnot(is.character(secretome_symbols))
    req_hallmark <- c("sample", "ID", score_col)
    req_expr <- c("sample", "symbol", expression_col)
    miss_h <- setdiff(req_hallmark, names(hallmark_scores))
    miss_e <- setdiff(req_expr, names(expression))
    if (length(miss_h) > 0L) {
        stop("Missing columns in 'hallmark_scores': ",
             paste(miss_h, collapse = ", "), call. = FALSE)
    }
    if (length(miss_e) > 0L) {
        stop("Missing columns in 'expression': ",
             paste(miss_e, collapse = ", "), call. = FALSE)
    }

    ####*Pivot hallmark scores to sample x pathway matrix*####
    hallmark_mat <- hallmark_scores |>
        dplyr::select("sample", "ID", dplyr::all_of(score_col)) |>
        tidyr::pivot_wider(
            names_from = "ID",
            values_from = dplyr::all_of(score_col)
        ) |>
        tibble::column_to_rownames("sample") |>
        as.matrix()

    ####*Pivot secretome expression to sample x gene matrix*####
    expr_mat <- expression |>
        dplyr::filter(.data$symbol %in% secretome_symbols) |>
        dplyr::select(
            "sample", "symbol", dplyr::all_of(expression_col)
        ) |>
        tidyr::pivot_wider(
            names_from = "symbol",
            values_from = dplyr::all_of(expression_col)
        ) |>
        tibble::column_to_rownames("sample") |>
        as.matrix()

    ####*Align samples between the two matrices*####
    common <- intersect(rownames(hallmark_mat), rownames(expr_mat))
    if (length(common) < 3L) {
        stop("Fewer than 3 samples shared between 'hallmark_scores' ",
             "and 'expression'; cannot compute correlation.",
             call. = FALSE)
    }
    hallmark_mat <- hallmark_mat[common, , drop = FALSE]
    expr_mat <- expr_mat[common, , drop = FALSE]

    ####*Pearson correlation via Hmisc::rcorr*####
    combined <- cbind(hallmark_mat, expr_mat)
    rc <- Hmisc::rcorr(combined, type = "pearson")
    path_names <- colnames(hallmark_mat)
    gene_names <- colnames(expr_mat)
    r_block <- rc$r[path_names, gene_names, drop = FALSE]
    p_block <- rc$P[path_names, gene_names, drop = FALSE]

    ####*Reshape r and p blocks to long tibble*####
    r_long <- as.data.frame(r_block) |>
        tibble::rownames_to_column("pathway") |>
        tidyr::pivot_longer(
            -"pathway",
            names_to = "gene",
            values_to = "value"
        )
    p_long <- as.data.frame(p_block) |>
        tibble::rownames_to_column("pathway") |>
        tidyr::pivot_longer(
            -"pathway",
            names_to = "gene",
            values_to = "p.value"
        )

    dplyr::inner_join(r_long, p_long, by = c("pathway", "gene"))
}
