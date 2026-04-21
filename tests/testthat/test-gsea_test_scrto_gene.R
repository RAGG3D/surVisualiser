test_that("gsea_test_scrto_gene errors clearly on missing columns", {
    bad <- data.frame(sample = "S1", foo = 1)
    expect_error(
        gsea_test_scrto_gene(bad, rank_by = "raw_count_scaled_adjusted"),
        "Missing required columns"
    )
})

test_that("gsea_test_scrto_gene validates excluded_genes type", {
    ok_data <- .make_toy_expression(n_samples = 2L, n_genes = 3L)
    expect_error(
        gsea_test_scrto_gene(
            ok_data,
            rank_by = "raw_count_scaled_adjusted",
            excluded_genes = 1:3
        )
    )
})

test_that("gsea_test_scrto_gene signals missing Suggests packages", {
    skip_if(
        requireNamespace("msigdbr", quietly = TRUE) &&
            requireNamespace("clusterProfiler", quietly = TRUE),
        "msigdbr + clusterProfiler are installed; skipping guard test"
    )
    expect_error(
        gsea_test_scrto_gene(
            .make_toy_expression(),
            rank_by = "raw_count_scaled_adjusted"
        ),
        "required"
    )
})
