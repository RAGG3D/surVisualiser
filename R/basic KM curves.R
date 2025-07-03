# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

split_median <- function(x, cat, item){
  x %>%
    group_by(!!as.name(cat)) %>%
    mutate(item = factor(Hmisc::cut2(!!as.name(item), g = 2), labels = c(1:nlevels(Hmisc::cut2(!!as.name(item), g = 2))))) %>%
    ungroup() %>%
    mutate(item = gsub("1", "L", item)) %>%
    mutate(item = gsub("2", "H", item))
}

ggsurvp_cat_item <- function(x, p, nrow, palette, xlab, ylab, title, size0){
  ggsurvplot(
    survfit(
      Surv(total_living_days, vital_status) ~ item,
      data = x
    ),
    data = x,
    facet.by = c("cat"),
    scales = "free_x",
    conf.int = T,
    risk.table = F,
    conf.int.alpha = 0.15,
    pval = p,
    nrow = nrow,
    legend.title = " ",
    short.panel.labs = T,
    palette = palette
  ) + theme_bw() +
    theme(text=element_text(size=size0, family="sans"),
          legend.position = "bottom") +
    labs(x = xlab, y = ylab, title = title) +
    guides(linetype = FALSE)
}
