# Synthetic fixtures shared across tests for the GSEA / correlation
# / heatmap functions. Keep small and deterministic so tests stay
# fast and independent of external databases.

.make_toy_expression <- function(n_samples = 20L,
                                 n_genes = 50L,
                                 seed = 42L) {
    set.seed(seed)
    samples <- paste0("S", seq_len(n_samples))
    genes   <- paste0("G", seq_len(n_genes))
    tidyr::expand_grid(sample = samples, symbol = genes) |>
        dplyr::mutate(
            raw_count_scaled_adjusted = stats::rnorm(
                dplyr::n(), mean = 5, sd = 2
            )
        )
}

.make_toy_hallmark_scores <- function(n_samples = 20L,
                                      n_hallmarks = 6L,
                                      seed = 7L) {
    set.seed(seed)
    samples <- paste0("S", seq_len(n_samples))
    ids     <- paste0("HALLMARK_TOY_", seq_len(n_hallmarks))
    tidyr::expand_grid(sample = samples, ID = ids) |>
        dplyr::mutate(
            NES = stats::rnorm(dplyr::n(), mean = 0, sd = 1)
        )
}

.make_toy_corr <- function(n_pathways = 4L,
                           n_genes = 8L,
                           seed = 11L) {
    set.seed(seed)
    pathways <- paste0("HALLMARK_TOY_", seq_len(n_pathways))
    genes    <- paste0("G", seq_len(n_genes))
    tidyr::expand_grid(pathway = pathways, gene = genes) |>
        dplyr::mutate(
            value   = stats::runif(dplyr::n(), -1, 1),
            p.value = stats::runif(dplyr::n(), 0, 0.1)
        )
}
