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

    ####*Filter to the requested gene pair and clean rows*####
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

    ####*Reset row names and return*####
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

    ####*Per-group median split into two equal-size bins*####
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

        ####*Relabel bins as Low / High*####
        dplyr::mutate(item = gsub("1", "L", .data$item)) |>
        dplyr::mutate(item = gsub("2", "H", .data$item))
}


#' Stratified Pairwise Survival P-Values
#'
#' Runs a stratified pairwise log-rank test
#' (\code{\link[survival:pairwise_survdiff]{pairwise_survdiff}})
#' comparing survival across combinations of a categorical variable
#' (\code{cat}) and a group variable (\code{item}), then reshapes the
#' result into a tidy long table. Self-comparisons and symmetric
#' duplicates (e.g., H/L vs L/H) are removed.
#'
#' @param x A \code{data.frame} or \code{tibble} with at least
#'   columns \code{total_living_days}, \code{vital_status},
#'   \code{cat}, and \code{item}.
#' @param method0 Character string. P-value adjustment method passed
#'   to \code{\link[survival:pairwise_survdiff]{pairwise_survdiff}}
#'   (for example \code{"BH"}, \code{"bonferroni"}, or
#'   \code{"none"}).
#' @param p_name0 Character string. Name of the p-value column in
#'   the returned table (for example \code{"P"} or \code{"adjP"}).
#'
#' @return A \code{data.frame} with columns \code{cat} (stratum),
#'   \code{comb1} and \code{comb2} (the two groups being compared),
#'   and the p-value column named by \code{p_name0}.
#'
#' @importFrom survminer pairwise_survdiff
#' @importFrom survival Surv
#' @importFrom dplyr mutate filter select all_of
#' @importFrom tidyr pivot_longer separate unite
#' @importFrom rlang .data !!
#'
#' @export
#'
#' @examples
#' \dontrun{
#' pvals <- strata_p(surv_data,
#'     method0 = "BH",
#'     p_name0 = "adjP"
#' )
#' }
strata_p <- function(x, method0, p_name0) {

    ####*Run pairwise log-rank test across cat + item strata*####
    pw <- survminer::pairwise_survdiff(
        survival::Surv(total_living_days, vital_status) ~ cat + item,
        p.adjust.method = method0,
        data = x
    )

    ####*Extract and reshape the p-value matrix*####
    result <- as.data.frame(pw[[3]]) |>
        dplyr::mutate(cat = rownames(as.data.frame(pw[[3]]))) |>
        tidyr::pivot_longer(
            -"cat",
            names_to = "comb2",
            values_to = p_name0
        ) |>
        stats::na.omit()

    ####*Parse stratum and group labels from row/col names*####
    result <- result |>
        tidyr::separate(
            "cat", c("cat", "comb1"), sep = ","
        ) |>
        dplyr::mutate(
            cat   = gsub(".*=", "", .data$cat),
            comb1 = gsub(".*=", "", .data$comb1),
            comb2 = gsub(".*cat.", "", .data$comb2)
        ) |>
        dplyr::mutate(
            comb2 = gsub("\\.", "/", .data$comb2),
            cat1  = gsub(" ", "/", .data$cat)
        ) |>
        tidyr::separate(
            "comb2", c("cat2", "comb2"), sep = "//item/"
        )

    ####*Remove duplicate and symmetric comparisons*####
    result |>
        tidyr::unite(
            "test", c("comb1", "comb2"),
            sep = "/", remove = FALSE
        ) |>
        dplyr::filter(
            !.data$test %in% c(
                "H/L/L/H", "L/H/H/L",
                "L/L/H/H", "H/H/L/L"
            )
        ) |>
        dplyr::filter(.data$cat1 == .data$cat2) |>
        dplyr::select(
            "cat", "comb1", "comb2", dplyr::all_of(p_name0)
        )
}
