test_that("three_heatmaps_of_pth_scrto_TN errors on missing columns", {
    good <- .make_toy_corr()
    bad  <- dplyr::rename(good, p_val = "p.value")
    out_dir <- tempfile("heatmap_out_")

    expect_error(
        three_heatmaps_of_pth_scrto_TN(
            tumor_corr = bad,
            normal_corr = good,
            cancer = "TEST",
            output_dir = out_dir
        ),
        "Missing columns in 'tumor_corr'"
    )
    expect_error(
        three_heatmaps_of_pth_scrto_TN(
            tumor_corr = good,
            normal_corr = bad,
            cancer = "TEST",
            output_dir = out_dir
        ),
        "Missing columns in 'normal_corr'"
    )
})

test_that("three_heatmaps_of_pth_scrto_TN reports missing Suggests", {
    viz_installed <- all(vapply(
        c("tidyHeatmap", "circlize", "RColorBrewer", "ComplexHeatmap"),
        requireNamespace, logical(1), quietly = TRUE
    ))
    skip_if(viz_installed, "Viz packages present; skipping guard test")

    expect_error(
        three_heatmaps_of_pth_scrto_TN(
            tumor_corr  = .make_toy_corr(),
            normal_corr = .make_toy_corr(),
            cancer      = "TEST",
            output_dir  = tempfile("hm_")
        ),
        "not.*installed"
    )
})

test_that("three_heatmaps_of_pth_scrto_TN writes a PDF when viz is present", {
    viz_installed <- all(vapply(
        c("tidyHeatmap", "circlize", "RColorBrewer", "ComplexHeatmap"),
        requireNamespace, logical(1), quietly = TRUE
    ))
    skip_if_not(viz_installed, "Viz packages not installed")

    tumor  <- .make_toy_corr()
    normal <- .make_toy_corr(seed = 22L)
    out_dir <- tempfile("heatmap_out_")

    res <- three_heatmaps_of_pth_scrto_TN(
        tumor_corr  = tumor,
        normal_corr = normal,
        cancer      = "TEST",
        output_dir  = out_dir
    )

    expect_named(res, c("tumor", "normal", "diff", "file"))
    expect_true(file.exists(res$file))
    expect_gt(file.info(res$file)$size, 0)
})
