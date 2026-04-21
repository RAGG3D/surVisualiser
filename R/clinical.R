#' Combine Expression Data with TCGA Clinical Data
#'
#' Joins a sample-level expression tibble with TCGA clinical data,
#' computing survival time and status.
#'
#' @param x A \code{tibble} with a \code{sample} column containing
#'   TCGA patient barcodes.
#' @param cancer Character string of the TCGA cancer type code
#'   (e.g., \code{"BRCA"}).
#' @param clinical_dir Character string. Directory containing clinical
#'   CSV files named \code{clinical_<cancer>.csv}.
#'
#' @return A \code{tibble} with additional columns:
#'   \code{total_living_days} (survival time in days),
#'   \code{vital_status} (1 = Dead, 0 = Alive), and
#'   \code{age} (in years).
#'
#' @importFrom dplyr inner_join mutate
#' @importFrom utils read.csv
#'
#' @export
#'
#' @examples
#' \dontrun{
#' expr_clinical <- clinical_combine(expr_data, "BRCA",
#'   clinical_dir = "/path/to/clinical"
#' )
#' }
clinical_combine <- function(x, cancer, clinical_dir) {

    ####*Locate clinical CSV for the requested cancer type*####
    clinical_file <- file.path(
        clinical_dir, paste0("clinical_", tolower(cancer), ".csv")
    )

    ####*Join expression with clinical and derive survival fields*####
    x |>
        dplyr::inner_join(
            utils::read.csv(clinical_file),
            by = c("sample" = "bcr_patient_barcode")
        ) |>
        dplyr::mutate(
            total_living_days = as.numeric(as.character(.data$days_to_death)),
            age = -as.numeric(as.character(.data$days_to_birth)) / 365
        ) |>

        ####*Fill missing survival time with last contact days*####
        dplyr::mutate(na = is.na(.data$total_living_days)) |>
        dplyr::mutate(
            total_living_days = ifelse(
                .data$na == TRUE,
                as.numeric(as.character(.data$last_contact_days_to)),
                .data$total_living_days
            )
        ) |>

        ####*Encode vital status as 1 = Dead / 0 = Alive*####
        dplyr::mutate(
            vital_status = ifelse(.data$vital_status == "Dead", 1, 0)
        )
}
