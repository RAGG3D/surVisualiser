
#Run CIBERSORT
CIBERSORT <- function(x, c1, c2, c3, ref, meth, act) {
  as_tibble(setDT(x)[, lapply(.SD, sum), by = c(c1, c2), .SDcols = c3]) %>%
    deconvolve_cellularity(
      .sample = !!sym(c1),
      .transcript = !!sym(c2),
      .abundance = !!sym(c3),
      reference = ref,
      method = meth,
      action = act
    )
}

########### Analysis

pw_diff_selection <- function(data, genepair){
  x <- data %>%
    filter(grepl(genepair, comp2)) %>%
    filter(grepl(genepair, comp1)) %>%
    filter(!grepl("HLLH", strata)) %>%
    mutate(gene_combine = genepair, 
           comp1 = substring(comp1, nchar(comp1)-2),
           comp2 = substring(comp2, nchar(comp2)-2)) %>%
    dplyr::select(comp1, comp2,P, adjP, gene_combine) %>%
    na.omit() %>%
    unique()
  rownames(x) <- NULL
  x
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


ggsurvp_cat_item_inf <- function(x, p, nrow, palette, xlab, ylab, title){
  ggsurvplot(
    survfit(
      Surv(total_living_days, vital_status) ~ item,
      data = x
    ),
    data = x,
    facet.by = c("cat", "infection"),
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
    theme(text=element_text(size=30, family="sans"),
          legend.position = "bottom") +
    labs(x = xlab, y = ylab, title = title) +
    guides(linetype = FALSE)
}

split_median <- function(x, cat, item){
  x %>%
    group_by(!!as.name(cat)) %>%
    mutate(item = factor(Hmisc::cut2(!!as.name(item), g = 2), labels = c(1:nlevels(Hmisc::cut2(!!as.name(item), g = 2))))) %>%
    ungroup() %>%
    mutate(item = gsub("1", "L", item)) %>%
    mutate(item = gsub("2", "H", item))
}
