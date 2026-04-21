#' Per-Sample Pre-Ranked GSEA with Secretome Gene Exclusion
#'
#' Performs per-sample pre-ranked gene set enrichment analysis (GSEA)
#' against MSigDB gene collections (hallmarks by default), optionally
#' excluding user-specified "confounding" genes (for example secretome
#' genes) from the input expression prior to ranking. The exclusion is
#' applied to the expression matrix, not to the gene sets themselves,
#' so that the excluded genes cannot contribute to the enrichment
#' statistic and downstream correlations between those genes and the
#' returned scores are not inflated by shared membership between the
#' hallmark sets and the excluded gene list.
#'
#' @details
#' This function wraps \code{clusterProfiler::GSEA} over each sample.
#' The per-sample workflow is:
#' \enumerate{
#'   \item Remove \code{excluded_genes} from \code{data} by gene
#'     symbol.
#'   \item Map symbols to Entrez IDs via \code{org.Hs.eg.db}.
#'   \item Rank genes by \code{rank_by} in descending order.
#'   \item For each requested MSigDB category, build a \code{TERM2GENE}
#'     map and call \code{clusterProfiler::GSEA}.
#'   \item Collect the resulting enrichment tables with sample and
#'     category identifiers.
#' }
#'
#' @note
#' This is per-sample pre-ranked GSEA and is mathematically distinct
#' from matrix-based per-sample scoring methods such as GSVA, ssGSEA,
#' singscore, AUCell, or simple gene-set mean. It returns an
#' enrichment statistic (normalised enrichment score, NES) per
#' pathway per sample, not a continuous score directly comparable to
#' the matrix-based methods. For workflows that correlate a
#' per-sample score matrix against external gene-set expression,
#' GSVA or ssGSEA is usually a more natural choice; \code{GSEA} is
#' retained here to preserve parity with historical analyses.
#'
#' @param data A \code{tibble} or \code{data.frame} with at least
#'   columns \code{sample}, \code{symbol}, and the ranking column.
#'   \code{sample} is the per-sample identifier, \code{symbol} is the
#'   HGNC gene symbol, and the ranking column is the numeric value
#'   used to order genes for pre-ranked GSEA (for example
#'   \code{"raw_count_scaled_adjusted"}).
#' @param rank_by Character string. Name of the numeric column in
#'   \code{data} used to rank genes for pre-ranked GSEA.
#' @param excluded_genes Character vector of gene symbols to remove
#'   from \code{data} before ranking. Pass \code{character(0)} (the
#'   default) to skip exclusion.
#' @param species Character string passed to \code{msigdbr::msigdbr}
#'   to select the species whose gene collections are queried.
#'   Defaults to \code{"Homo sapiens"}.
#' @param gene_collections Character vector of MSigDB category codes
#'   (case-insensitive). Defaults to \code{"h"} (hallmarks).
#' @param pvalueCutoff Numeric passed to \code{clusterProfiler::GSEA}.
#'   Defaults to \code{1}, which returns all pathways regardless of
#'   significance.
#'
#' @return A \code{tibble} with one row per pathway per sample per
#'   gene collection. Columns include \code{sample}, \code{gs_cat},
#'   \code{idx_for_plotting}, and all columns produced by
#'   \code{clusterProfiler::GSEA} (for example \code{ID},
#'   \code{Description}, \code{NES}, \code{pvalue}, \code{p.adjust}).
#'
#' @seealso \code{\link[clusterProfiler:GSEA]{clusterProfiler::GSEA}},
#'   \code{\link[msigdbr:msigdbr]{msigdbr::msigdbr}}.
#'
#' @importFrom dplyr filter mutate select arrange desc inner_join
#'   bind_rows all_of
#' @importFrom purrr map
#' @importFrom tibble as_tibble deframe rowid_to_column
#' @importFrom rlang .data
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data(secretome_genes)
#' result <- gsea_test_scrto_gene(
#'     data = expr_tbl,
#'     rank_by = "raw_count_scaled_adjusted",
#'     excluded_genes = secretome_genes$symbol,
#'     species = "Homo sapiens",
#'     gene_collections = "h"
#' )
#' }
gsea_test_scrto_gene <- function(data,
                                 rank_by,
                                 excluded_genes = character(0),
                                 species = "Homo sapiens",
                                 gene_collections = "h",
                                 pvalueCutoff = 1) {

    ####*Input validation*####
    stopifnot(is.character(excluded_genes))
    required_cols <- c("sample", "symbol", rank_by)
    missing_cols <- setdiff(required_cols, names(data))
    if (length(missing_cols) > 0L) {
        stop("Missing required columns in 'data': ",
             paste(missing_cols, collapse = ", "), call. = FALSE)
    }

    ####*Dependency checks*####
    if (!requireNamespace("msigdbr", quietly = TRUE)) {
        stop("Package 'msigdbr' is required for ",
             "gsea_test_scrto_gene(); install it to proceed.",
             call. = FALSE)
    }
    if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
        stop("Package 'clusterProfiler' is required for ",
             "gsea_test_scrto_gene(); install it to proceed.",
             call. = FALSE)
    }

    ####*Gene exclusion*####
    cleaned <- dplyr::filter(
        data, !.data$symbol %in% excluded_genes
    )

    ####*MSigDB gene collection lookup*####
    gene_collection <- msigdbr::msigdbr(species = species) |>
        dplyr::filter(
            tolower(.data$gs_cat) %in% tolower(gene_collections)
        )
    cat_levels <- unique(gene_collection$gs_cat)

    ####*Per-sample pre-ranked GSEA*####
    samples <- unique(cleaned$sample)
    per_sample <- purrr::map(samples, function(s) {
        ranks <- cleaned |>
            dplyr::filter(.data$sample == s) |>
            dplyr::inner_join(
                AnnotationDbi::toTable(
                    org.Hs.eg.db::org.Hs.egSYMBOL
                ),
                by = "symbol"
            ) |>
            dplyr::arrange(dplyr::desc(.data[[rank_by]])) |>
            dplyr::select("gene_id", dplyr::all_of(rank_by)) |>
            tibble::deframe()

        purrr::map(cat_levels, function(cat0) {
            term2gene <- gene_collection |>
                dplyr::filter(.data$gs_cat == cat0) |>
                dplyr::select("gs_name", "entrez_gene")

            fit <- clusterProfiler::GSEA(
                geneList = ranks,
                TERM2GENE = term2gene,
                pvalueCutoff = pvalueCutoff
            )

            tibble::as_tibble(fit) |>
                tibble::rowid_to_column(var = "idx_for_plotting") |>
                dplyr::mutate(sample = s, gs_cat = cat0)
        }) |>
            dplyr::bind_rows()
    })

    ####*Combine results*####
    dplyr::bind_rows(per_sample)
}
