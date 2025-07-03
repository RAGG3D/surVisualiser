# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

pw_diff_selection <- function(data, genepair){
  x <- data %>%
    dplyr::filter(grepl(genepair, comp2)) %>%
    dplyr::filter(grepl(genepair, comp1)) %>%
    dplyr::filter(!grepl("HLLH", strata)) %>%
    dplyr::mutate(gene_combine = genepair,
           comp1 = substring(comp1, nchar(comp1)-2),
           comp2 = substring(comp2, nchar(comp2)-2)) %>%
    dplyr::select(comp1, comp2,P, adjP, gene_combine) %>%
    na.omit() %>%
    unique()
  rownames(x) <- NULL
  x
}

strata_p <- function(x, method0, p_name0){

  as.data.frame(pairwise_survdiff(Surv(total_living_days, vital_status) ~ cat + item,
                                  p.adjust.method = method0,
                                  data = x)[3]) %>%
    dplyr::mutate(cat = rownames(.)) %>%
    pivot_longer(-cat, names_to = "comb2", values_to = p_name0) %>%
    na.omit() %>%
    separate("cat", c("cat", "comb1"), sep = ",") %>%
    dplyr::mutate(cat = gsub(".*=", "", cat),
           comb1 = gsub(".*=", "", comb1),
           comb2 = gsub(".*cat.", "", comb2)) %>%
    dplyr::mutate(comb2 = gsub("\\.", "/", comb2),
           cat1 = gsub(" ", "/", cat)) %>%
    separate("comb2", c("cat2", "comb2"), sep = "//item/") %>%   ### keep cat2 in complicated analysis for double check
    unite("test", c("comb1", "comb2"), sep = "/", remove = F) %>%
    dplyr::filter(!test %in% c("H/L/L/H", "L/H/H/L", "L/L/H/H", "H/H/L/L")) %>%
    dplyr::filter(cat1 == cat2) %>%
    dplyr::select(cat, comb1, comb2, !!as.name(p_name0))
}
