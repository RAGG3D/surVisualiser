#' Three-Panel Heatmap of Hallmark-Secretome Correlations
#'
#' Produces a side-by-side three-panel heatmap of pathway x secretome
#' gene correlations across tissue types: Tumor, Normal, and the
#' Tumor - Normal difference. Intended as the visual endpoint of the
#' workflow \code{\link{gsea_test_scrto_gene}} ->
#' \code{\link{corr_pathway_scrto}} -> this function, once per cancer
#' type.
#'
#' @details
#' Non-significant correlations (\code{p.value > p_cutoff}) in either
#' the tumor or the normal table are treated as zero before
#' differencing, matching the treatment in the original analysis
#' scripts. Hallmark-gene pairs present in only one table have the
#' missing side filled with zero.
#'
#' Rendering uses \code{\link[tidyHeatmap:heatmap]{tidyHeatmap::heatmap}},
#' a diverging \code{RdBu} palette via
#' \code{\link[circlize:colorRamp2]{circlize::colorRamp2}}, and
#' \code{\link[ComplexHeatmap:draw]{ComplexHeatmap::draw}} for the
#' side-by-side layout. These packages are declared under
#' \code{Suggests} and must be installed at call time.
#'
#' @param tumor_corr A \code{tibble} from
#'   \code{\link{corr_pathway_scrto}} computed on tumor samples.
#'   Must have columns \code{pathway}, \code{gene}, \code{value}, and
#'   \code{p.value}.
#' @param normal_corr A \code{tibble} from
#'   \code{\link{corr_pathway_scrto}} computed on normal samples.
#'   Same column requirements as \code{tumor_corr}.
#' @param cancer Character string used in panel titles and the
#'   default output file name (for example \code{"BRCA"}).
#' @param output_dir Character string. Directory in which the
#'   combined PDF is saved. Created if it does not already exist.
#' @param p_cutoff Numeric significance threshold applied to both
#'   \code{tumor_corr} and \code{normal_corr} before plotting.
#'   Defaults to \code{0.05}.
#' @param height,width Numeric PDF dimensions in inches. Default
#'   \code{8 x 15}, sized for three panels side by side.
#' @param file_name Character string. Optional custom PDF file name
#'   (without directory). Defaults to
#'   \code{paste0(cancer, "_heatmap_Tumor_Normal_Diff.pdf")}.
#'
#' @return Invisibly, a named list with elements \code{tumor},
#'   \code{normal}, and \code{diff} (the three \code{tidyHeatmap}
#'   objects) and \code{file} (the full path to the written PDF).
#'
#' @seealso \code{\link{corr_pathway_scrto}},
#'   \code{\link{gsea_test_scrto_gene}}.
#'
#' @importFrom dplyr filter select mutate left_join transmute
#' @importFrom tidyr pivot_wider pivot_longer replace_na
#' @importFrom rlang .data
#' @importFrom grDevices pdf dev.off
#'
#' @export
#'
#' @examples
#' \dontrun{
#' three_heatmaps_of_pth_scrto_TN(
#'     tumor_corr = corr_tumor,
#'     normal_corr = corr_normal,
#'     cancer = "BRCA",
#'     output_dir = "/path/to/figures"
#' )
#' }
three_heatmaps_of_pth_scrto_TN <- function(tumor_corr,
                                           normal_corr,
                                           cancer,
                                           output_dir,
                                           p_cutoff = 0.05,
                                           height = 8,
                                           width = 15,
                                           file_name = NULL) {

    ####*Dependency checks*####
    needed <- c("tidyHeatmap", "circlize", "RColorBrewer",
                "ComplexHeatmap")
    missing_pkgs <- needed[!vapply(
        needed, requireNamespace, logical(1), quietly = TRUE
    )]
    if (length(missing_pkgs) > 0L) {
        stop("The following packages are required for ",
             "three_heatmaps_of_pth_scrto_TN() but are not ",
             "installed: ", paste(missing_pkgs, collapse = ", "),
             call. = FALSE)
    }

    ####*Input validation*####
    req <- c("pathway", "gene", "value", "p.value")
    miss_t <- setdiff(req, names(tumor_corr))
    miss_n <- setdiff(req, names(normal_corr))
    if (length(miss_t) > 0L) {
        stop("Missing columns in 'tumor_corr': ",
             paste(miss_t, collapse = ", "), call. = FALSE)
    }
    if (length(miss_n) > 0L) {
        stop("Missing columns in 'normal_corr': ",
             paste(miss_n, collapse = ", "), call. = FALSE)
    }

    ####*Prepare output directory*####
    if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
    }
    if (is.null(file_name)) {
        file_name <- paste0(
            cancer, "_heatmap_Tumor_Normal_Diff.pdf"
        )
    }
    pdf_path <- file.path(output_dir, file_name)

    ####*Tumor panel data and heatmap*####
    x1 <- tumor_corr |>
        dplyr::filter(.data$p.value <= p_cutoff) |>
        dplyr::select("pathway", "gene", "value") |>
        tidyr::pivot_wider(
            names_from = "gene",
            values_from = "value",
            values_fill = 0
        ) |>
        tidyr::pivot_longer(
            -"pathway",
            names_to = "Secretome_Tumor",
            values_to = "value"
        )

    p1 <- x1 |>
        tidyHeatmap::heatmap(
            .column = .data$Secretome_Tumor,
            .row = .data$pathway,
            .value = .data$value,
            palette_value = circlize::colorRamp2(
                seq(min(x1$value), max(x1$value), length.out = 5L),
                rev(RColorBrewer::brewer.pal(5L, "RdBu"))
            ),
            clustering_distance_columns = "manhattan",
            clustering_method_columns = "ward.D",
            clustering_distance_rows = "manhattan",
            clustering_method_rows = "ward.D",
            column_title = paste0(cancer, " Tumour")
        )

    ####*Normal panel data and heatmap*####
    x2 <- normal_corr |>
        dplyr::filter(.data$p.value <= p_cutoff) |>
        dplyr::select("pathway", "gene", "value") |>
        tidyr::pivot_wider(
            names_from = "gene",
            values_from = "value",
            values_fill = 0
        ) |>
        tidyr::pivot_longer(
            -"pathway",
            names_to = "Secretome_Normal",
            values_to = "value"
        )

    p2 <- x2 |>
        tidyHeatmap::heatmap(
            .column = .data$Secretome_Normal,
            .row = .data$pathway,
            .value = .data$value,
            palette_value = circlize::colorRamp2(
                seq(min(x2$value), max(x2$value), length.out = 5L),
                rev(RColorBrewer::brewer.pal(5L, "RdBu"))
            ),
            clustering_distance_columns = "manhattan",
            clustering_method_columns = "ward.D",
            clustering_distance_rows = "manhattan",
            clustering_method_rows = "ward.D",
            column_title = paste0(cancer, " Normal")
        )

    ####*Tumor - Normal difference panel*####
    x3 <- x1 |>
        dplyr::transmute(
            pathway = .data$pathway,
            Secretome = .data$Secretome_Tumor,
            tumor = .data$value
        ) |>
        dplyr::left_join(
            x2 |>
                dplyr::transmute(
                    pathway = .data$pathway,
                    Secretome = .data$Secretome_Normal,
                    normal = .data$value
                ),
            by = c("pathway", "Secretome")
        ) |>
        dplyr::mutate(
            tumor  = tidyr::replace_na(.data$tumor, 0),
            normal = tidyr::replace_na(.data$normal, 0),
            value  = .data$tumor - .data$normal,
            Secretome_Complement = .data$Secretome
        ) |>
        dplyr::select("pathway", "Secretome_Complement", "value")

    p3 <- x3 |>
        tidyHeatmap::heatmap(
            .column = .data$Secretome_Complement,
            .row = .data$pathway,
            .value = .data$value,
            palette_value = circlize::colorRamp2(
                seq(min(x3$value), max(x3$value), length.out = 5L),
                rev(RColorBrewer::brewer.pal(5L, "RdBu"))
            ),
            clustering_distance_columns = "manhattan",
            clustering_method_columns = "ward.D",
            clustering_distance_rows = "manhattan",
            clustering_method_rows = "ward.D",
            column_title = paste0(cancer, " T-N")
        )

    ####*Render three-panel PDF*####
    grDevices::pdf(file = pdf_path, height = height, width = width)
    ComplexHeatmap::draw(p1 + p2 + p3)
    grDevices::dev.off()

    invisible(list(
        tumor  = p1,
        normal = p2,
        diff   = p3,
        file   = pdf_path
    ))
}
