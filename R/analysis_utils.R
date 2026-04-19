#' Extract Pairwise Survival Comparisons for a Gene Pair
#'
#' Filters pairwise comparison results for a specific gene pair
#' and returns a clean summary table.
#'
#' @param data A \code{data.frame} with columns \code{comp1},
#'   \code{comp2}, \code{strata}, \code{P}, and \code{adjP}.
#' @param genepair Character string. Gene pair identifier to filter on.
#'
#' @return A \code{data.frame} with columns \code{comp1}, \code{comp2},
#'   \code{P}, \code{adjP}, and \code{gene_combine}.
#'
#' @importFrom dplyr filter mutate select
#'
#' @export
#'
#' @examples
#' \dontrun{
#' pw_results <- pw_diff_selection(pairwise_data, "TP53_BRCA1")
#' }
pw_diff_selection <- function(data, genepair) {
    x <- data |>
        dplyr::filter(grepl(genepair, .data$comp2)) |>
        dplyr::filter(grepl(genepair, .data$comp1)) |>
        dplyr::filter(!grepl("HLLH", .data$strata)) |>
        dplyr::mutate(
            gene_combine = genepair,
            comp1 = substring(.data$comp1, nchar(.data$comp1) - 2),
            comp2 = substring(.data$comp2, nchar(.data$comp2) - 2)
        ) |>
        dplyr::select("comp1", "comp2", "P", "adjP", "gene_combine") |>
        stats::na.omit() |>
        unique()
    rownames(x) <- NULL
    x
}


#' Split a Continuous Variable at the Median Within Groups
#'
#' Dichotomises a continuous variable into High (H) and Low (L)
#' groups based on the median within each level of a grouping variable.
#'
#' @param x A \code{data.frame} or \code{tibble}.
#' @param cat Character string. Name of the grouping column.
#' @param item Character string. Name of the continuous column to split.
#'
#' @return The input data with an added \code{item} column coded as
#'   \code{"L"} (below median) or \code{"H"} (above median).
#'
#' @importFrom dplyr group_by mutate ungroup
#' @importFrom Hmisc cut2
#'
#' @export
#'
#' @examples
#' \dontrun{
#' dat <- split_median(my_data, cat = "cancer_type", item = "gene_expr")
#' }
split_median <- function(x, cat, item) {
    x |>
        dplyr::group_by(!!as.name(cat)) |>
        dplyr::mutate(
            item = factor(
                Hmisc::cut2(!!as.name(item), g = 2),
                labels = seq_len(
                    nlevels(Hmisc::cut2(!!as.name(item), g = 2))
                )
            )
        ) |>
        dplyr::ungroup() |>
        dplyr::mutate(item = gsub("1", "L", .data$item)) |>
        dplyr::mutate(item = gsub("2", "H", .data$item))
}
