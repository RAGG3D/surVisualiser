test_that("corr_pathway_scrto returns long tibble with expected columns", {
    hs <- .make_toy_hallmark_scores()
    ex <- .make_toy_expression()
    secretome <- paste0("G", 1:5)

    out <- corr_pathway_scrto(
        hallmark_scores   = hs,
        expression        = ex,
        secretome_symbols = secretome
    )

    expect_s3_class(out, "tbl_df")
    expect_named(out, c("pathway", "gene", "value", "p.value"))
    expect_equal(
        nrow(out),
        length(unique(hs$ID)) * length(secretome)
    )
    expect_true(all(out$value >= -1 & out$value <= 1, na.rm = TRUE))
})

test_that("corr_pathway_scrto errors on column mismatches", {
    hs <- .make_toy_hallmark_scores()
    ex <- .make_toy_expression()
    bad_hs <- dplyr::rename(hs, Pathway = "ID")
    bad_ex <- dplyr::rename(ex, gene = "symbol")

    expect_error(
        corr_pathway_scrto(bad_hs, ex, "G1"),
        "Missing columns in 'hallmark_scores'"
    )
    expect_error(
        corr_pathway_scrto(hs, bad_ex, "G1"),
        "Missing columns in 'expression'"
    )
})

test_that("corr_pathway_scrto errors when fewer than 3 samples overlap", {
    hs <- .make_toy_hallmark_scores(n_samples = 20L)
    ex <- .make_toy_expression(n_samples = 2L)
    ex <- dplyr::mutate(ex, sample = paste0("X", .data$sample))
    expect_error(
        corr_pathway_scrto(hs, ex, secretome_symbols = "G1"),
        "Fewer than 3 samples"
    )
})

test_that("corr_pathway_scrto reproduces a known Pearson correlation", {
    set.seed(1)
    samples <- paste0("S", 1:30)
    x <- stats::rnorm(30)
    y <- x + stats::rnorm(30, sd = 0.1)

    hs <- tibble::tibble(sample = samples, ID = "HALLMARK_A", NES = x)
    ex <- tibble::tibble(
        sample = samples,
        symbol = "G1",
        raw_count_scaled_adjusted = y
    )
    out <- corr_pathway_scrto(hs, ex, secretome_symbols = "G1")

    expected <- stats::cor(x, y)
    expect_equal(out$value, expected, tolerance = 1e-6)
    expect_true(out$p.value < 1e-10)
})
