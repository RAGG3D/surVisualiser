ref = readRDS("C:/Rstudio/Chapter 3/data/my_ref.rds")[,-c(1,4,11,17,21,24,28,31,32,33)]
colnames(ref) <- (data.frame(Var2 = colnames(ref)) %>%
                    inner_join(read_csv("C:/Rstudio/Chapter 3/data/cellname.csv")))$`Cell Type`
scrto <- read_csv("C:/Rstudio/Chapter 3/data/comprehensive secreted protein gene list.csv")
h_df = msigdbr(species = "Homo sapiens", category = "H")

library(tidyverse)
library(tidybulk)
library(org.Hs.eg.db)
library(Hmisc)
library(msigdbr)
####*Functions*####
entrez_rank_to_gsea = function(my_entrez_rank, species, gene_collections  = NULL){
  my_gene_collection = msigdbr::msigdbr(species = species) %>%  filter( tolower(gs_cat) %in% tolower(gene_collections) )
  my_gene_collection |>
    nest(data = -gs_cat) |>
    mutate(fit =
             map(
               data,
               ~ 	clusterProfiler::GSEA(
                 my_entrez_rank,
                 TERM2GENE=.x %>% dplyr::select(gs_name, entrez_gene),
                 pvalueCutoff = 1
               )
               
             )) |>
    mutate(test =
             map(
               fit,
               ~ .x |>
                 # ggplot2::fortify(showCategory=Inf) %>%
                 as_tibble() |>
                 rowid_to_column(var = "idx_for_plotting")
               #%>%
               #	mutate(plot = future_imap(ID, ~ enrichplot::gseaplot2(fit, geneSetID = .y, title = .x)))
               
             )) |>
    dplyr::select(-data)
  
}

H_scrto <- h_df %>% filter(gene_symbol %in% scrto$symbol)
map(list(cancerlist),
    function(cancer0){
      a <- read_csv(
        paste0("C:/Rstudio/Chapter 3/data/Adjusted counts all gene/",cancer0,".csv")
      )
      sample0 = unique(a$sample_ct)
      
      #samplelist = sample0[(which(sample0 == "TCGA-67-6217-Tumor", arr.ind=TRUE)+1):length(sample0)]
      samplelist = sample0
      map(as.list(unique(scrto$symbol)), function(scrto1){
        dir.create(paste0("C:/Rstudio/Chapter 3/output/TCGA pathway scores/Hallmarks/Single scrto removal/",cancer0,"/"))
        dir.create(paste0("C:/Rstudio/Chapter 3/output/TCGA pathway scores/Hallmarks/Single scrto removal/",cancer0,"/", scrto1, "/"))
        map(as.list(samplelist), function(sample1){
          if(scrto1 %in% H_scrto$gene_symbol) {
            x <- (a) %>%
          filter(sample_ct == sample1) %>%
          filter(symbol != scrto1) %>%
          inner_join(toTable(org.Hs.egSYMBOL)) 
          } else {
            x <- (a) %>%
              filter(sample_ct == sample1) %>%
              inner_join(toTable(org.Hs.egSYMBOL)) 
          }
        

        .data = x
        .entrez = x$gene_id
        .arrange_desc = x$raw_count_scaled_adjusted
        species = "Homo sapiens"
        gene_sets = c("h")
        .sample = x$sample_ct
        
        a <- x %>%
          pivot_transcript(.transcript = gene_id) %>%
          arrange(desc(raw_count_scaled_adjusted)) %>%
          dplyr::select(gene_id, raw_count_scaled_adjusted) %>%
          deframe() %>%
          entrez_rank_to_gsea(species, gene_collections = gene_sets) 
        
        as.data.frame(a$test)%>%
          write.csv(paste0("C:/Rstudio/Chapter 3/output/TCGA pathway scores/Hallmarks/Single scrto removal/",cancer0,"/", scrto1, "/", sample1, ".csv"), row.names = F)
      })
      })
       
    })

three_heatmaps_of_pth_scrto_TN <- function(cancerlist){
  map(cancerlist, function(cancer0){
    dir.create(paste0("C:/Rstudio/Chapter 3/output/Figures/", cancer0, "/"))
    #### tumor hallmark score and scrtos
    x1 <- read_csv(paste0("C:/Rstudio/Chapter 3/output/TCGA pathway scores/Hallmarks/", cancer0, " whole gene pathway+scrto corr.csv")) %>%
      filter(!grepl("HALLMARK_", X2)) %>%
      filter(p.value <= 0.05) %>%
      dplyr::select(X1, X2, value) %>% 
      pivot_wider(names_from = "X2", values_from = value) %>% 
      replace(is.na(.), 0) %>%
      pivot_longer(-X1, names_to = "Secretome_Tumor", values_to = "value") 
    
    
    p1 <- x1 %>%
      tidybulk::rename("Pathway" = "X1") %>%
      tidyHeatmap::heatmap(
        .column = Secretome_Tumor,
        .row = Pathway,
        .value = value,   
        palette_value = circlize::colorRamp2(
          seq(min(x1$value), max(x1$value), length.out = 5), 
          rev(RColorBrewer::brewer.pal(5, "RdBu"))
        ),
        
        clustering_distance_columns = "manhattan",
        clustering_method_columns = "ward.D",
        clustering_distance_rows = "manhattan",
        clustering_method_rows = "ward.D",
        column_title = paste0(cancer0, " Tumour")
      ) 
    
    #### normal hallmark score and scrtos 
    x2 <- read_csv(paste0("C:/Rstudio/Chapter 3/output/TCGA pathway scores/Hallmarks/", cancer0, " whole gene normal pathway+scrto corr.csv")) %>%
      filter(!grepl("HALLMARK_", X2)) %>%
      filter(p.value <= 0.05) %>%
      dplyr::select(X1, X2, value) %>% 
      pivot_wider(names_from = "X2", values_from = value) %>% 
      replace(is.na(.), 0) %>%
      pivot_longer(-X1, names_to = "Secretome_Normal", values_to = "value") 
    
    
    p2 <- x2 %>%
      tidybulk::rename("Pathway" = "X1") %>%
      tidyHeatmap::heatmap(
        .column = Secretome_Normal,
        .row = Pathway,
        .value = value,   
        palette_value = circlize::colorRamp2(
          seq(min(x2$value), max(x2$value), length.out = 5), 
          rev(RColorBrewer::brewer.pal(5, "RdBu"))
        ),
        # Arguments passed to ComplexHeatmap 
        clustering_distance_columns = "manhattan",
        clustering_method_columns = "ward.D",
        clustering_distance_rows = "manhattan",
        clustering_method_rows = "ward.D",
        column_title = paste0(cancer0, " Normal")
      ) 
    
    #### tumor-normal hallmark score and scrtos
    x3 <- x1 %>% mutate(tumor = value, Secretome = Secretome_Tumor) %>%
      dplyr::select(X1, Secretome, tumor) %>%
      left_join(x2 %>% mutate(normal = value, Secretome = Secretome_Normal) %>%
                  dplyr::select(X1, Secretome, normal)) %>%
      mutate(value = tumor - normal, Secretome_Complement = Secretome)%>%
      dplyr::select(X1, Secretome_Complement, value) %>% 
      replace(is.na(.), 0)
    
    p3 <- x3 %>%
      tidybulk::rename("Pathway" = "X1") %>%
      tidyHeatmap::heatmap(
        .column = Secretome_Complement,
        .row = Pathway,
        .value = value,   
        palette_value = circlize::colorRamp2(
          seq(min(x3$value), max(x3$value),  length.out = 5), 
          rev(RColorBrewer::brewer.pal(5, "RdBu"))
        ),
        
        clustering_distance_columns = "manhattan",
        clustering_method_columns = "ward.D",
        clustering_distance_rows = "manhattan",
        clustering_method_rows = "ward.D",
        column_title = paste0(cancer0, " T-N")
      ) 
    
    pdf(file=paste0("C:/Rstudio/Chapter 3/output/Figures/", cancer0, "/heatmap Tumor Normal.pdf"), height = 8, width = 10)
    draw(p1+p2)
    dev.off()
    
    
    pdf(file=paste0("C:/Rstudio/Chapter 3/output/Figures/", cancer0, "/heatmap Tumor (Tumor-Normal).pdf"), height = 8, width = 10)
    draw(p1+p3)
    dev.off()
    
  })
}

