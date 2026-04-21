#' Read and Aggregate TCGA Transcript-Level Expression Data
#'
#' Reads transcript-level count files from a TCGA cancer project directory,
#' maps Ensembl IDs to gene symbols, aggregates duplicate counts, and
#' optionally scales abundance. Supports filtering by sample type
#' (all, tumor, or normal).
#'
#' @param cancer Character string of the TCGA cancer type code
#'   (e.g., \code{"BRCA"}, \code{"LUAD"}).
#' @param data_dir Character string. Root directory containing the TCGA
#'   isoform data. Should contain a subfolder
#'   \code{TCGA-<cancer>/Transcript/} with count files and a
#'   corresponding sample sheet CSV.
#' @param sample_type Character string indicating which samples to include.
#'   One of \code{"all"} (default), \code{"tumor"}, or \code{"normal"}.
#'   When \code{"all"}, uses \code{gdc_sample_sheet.csv};
#'   when \code{"tumor"}, uses \code{tumor_gdc_sample_sheet.csv};
#'   when \code{"normal"}, uses \code{normal_gdc_sample_sheet.csv}.
#' @param scale Logical. Whether to apply \code{tidybulk::scale_abundance}
#'   to the aggregated counts. Default is \code{TRUE}.
#'
#' @return A \code{tibble} with columns \code{sample}, \code{symbol}, and
#'   \code{raw_count} (plus scaled columns if \code{scale = TRUE}).
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Reads all count files from the Transcript directory.
#'   \item Maps Ensembl IDs to gene symbols via
#'     \code{org.Hs.eg.db}.
#'   \item Joins with the appropriate GDC sample sheet to assign
#'     patient-level Case IDs.
#'   \item Sums raw counts for duplicate sample-gene pairs.
#'   \item Optionally scales abundance using \code{tidybulk::scale_abundance}.
#' }
#'
#' @importFrom data.table setDT fread
#' @importFrom purrr map_dfr
#' @importFrom dplyr mutate select inner_join filter
#' @importFrom readr read_csv
#' @importFrom tibble as_tibble
#' @importFrom rlang .data
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Read all samples for BRCA
#' brca_all <- TCGA_transcript("BRCA", data_dir = "/path/to/data")
#'
#' # Read only tumor samples without scaling
#' brca_tumor <- TCGA_transcript("BRCA",
#'   data_dir = "/path/to/data",
#'   sample_type = "tumor",
#'   scale = FALSE
#' )
#' }
TCGA_transcript <- function(cancer,
                            data_dir,
                            sample_type = c("all", "tumor", "normal"),
                            scale = TRUE) {

    ####*Argument resolution and path construction*####
    sample_type <- match.arg(sample_type)

    transcript_dir <- file.path(
        data_dir, paste0("TCGA-", cancer), "Transcript"
    )

    sheet_name <- switch(sample_type,
        all    = "gdc_sample_sheet.csv",
        tumor  = "tumor_gdc_sample_sheet.csv",
        normal = "normal_gdc_sample_sheet.csv"
    )
    sheet_path <- file.path(
        data_dir, paste0("TCGA-", cancer), sheet_name
    )

    ####*Read raw transcript-level counts*####
    count_files <- list.files(transcript_dir)

    counts <- purrr::map_dfr(as.list(count_files), function(i) {
        data.table::fread(file.path(transcript_dir, i)) |>
            dplyr::mutate(sample = i)
    })

    ####*Annotation mapping and sample sheet join*####
    counts <- counts |>
        tidybulk::rename(raw_count = "V2") |>
        dplyr::mutate(ensembl_id = gsub("\\..*", "", .data$V1)) |>
        dplyr::inner_join(
            AnnotationDbi::toTable(org.Hs.eg.db::org.Hs.egENSEMBL),
            by = "ensembl_id"
        ) |>
        dplyr::inner_join(
            AnnotationDbi::toTable(org.Hs.eg.db::org.Hs.egSYMBOL),
            by = "gene_id"
        ) |>
        dplyr::select("sample", "symbol", "raw_count") |>
        dplyr::mutate(sample = gsub("counts.*", "counts.gz", .data$sample)) |>
        dplyr::inner_join(
            readr::read_csv(sheet_path, show_col_types = FALSE) |>
                dplyr::mutate(sample = .data$`File Name`),
            by = "sample"
        ) |>
        dplyr::mutate(sample = .data$`Case ID`) |>
        dplyr::select("sample", "symbol", "raw_count")

    ####*Aggregate duplicated sample-gene rows*####
    result <- tibble::as_tibble(
        data.table::setDT(counts)[,
            list(raw_count = sum(raw_count)),
            by = c("sample", "symbol")
        ]
    )

    ####*Optional abundance scaling*####
    if (scale) {
        result <- tidybulk::scale_abundance(
            result,
            .sample = sample,
            .abundance = raw_count,
            .transcript = symbol
        )
    }

    result
}
