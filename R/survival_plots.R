#' Generate Faceted Kaplan-Meier Survival Plots
#'
#' Creates a Kaplan-Meier survival plot faceted by a categorical variable,
#' with groups defined by the \code{item} column.
#'
#' @param x A \code{data.frame} or \code{tibble} containing at least
#'   columns \code{total_living_days}, \code{vital_status},
#'   \code{item}, and \code{cat}.
#' @param p Logical. Whether to display p-values on the plot.
#' @param nrow Integer. Number of rows for facet layout.
#' @param palette Character vector of colors for the survival curves.
#' @param xlab Character string. X-axis label.
#' @param ylab Character string. Y-axis label.
#' @param title Character string. Plot title.
#' @param size0 Numeric. Base text size for the plot theme.
#'
#' @return A \code{ggsurvplot} object.
#'
#' @importFrom survival survfit Surv
#' @importFrom survminer ggsurvplot
#' @importFrom ggplot2 theme_bw theme element_text labs guides
#'
#' @export
#'
#' @examples
#' \dontrun{
#' ggsurvp_cat_item(
#'   surv_data,
#'   p = TRUE, nrow = 2,
#'   palette = c("#E41A1C", "#377EB8"),
#'   xlab = "Days", ylab = "Survival Probability",
#'   title = "Survival by Group", size0 = 12
#' )
#' }
ggsurvp_cat_item <- function(x, p, nrow, palette, xlab, ylab, title, size0) {
    fit <- survival::survfit(
        survival::Surv(total_living_days, vital_status) ~ item,
        data = x
    )

    survminer::ggsurvplot(
        fit,
        data = x,
        facet.by = "cat",
        scales = "free_x",
        conf.int = TRUE,
        risk.table = FALSE,
        conf.int.alpha = 0.15,
        pval = p,
        nrow = nrow,
        legend.title = " ",
        short.panel.labs = TRUE,
        palette = palette
    ) +
        ggplot2::theme_bw() +
        ggplot2::theme(
            text = ggplot2::element_text(size = size0, family = "sans"),
            legend.position = "bottom"
        ) +
        ggplot2::labs(x = xlab, y = ylab, title = title) +
        ggplot2::guides(linetype = "none")
}
