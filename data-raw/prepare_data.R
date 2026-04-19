## Script to prepare package data from raw files

# Cell name mapping
cellname <- read.csv("../surVisualizer/data/cellname.csv",
    stringsAsFactors = FALSE
)
usethis::use_data(cellname, overwrite = TRUE)

# Secretome gene list
secretome_genes <- read.csv(
    "../surVisualizer/data/comprehensive secreted protein gene list.csv",
    stringsAsFactors = FALSE
)
usethis::use_data(secretome_genes, overwrite = TRUE)

# CIBERSORT reference matrix
my_ref <- readRDS("../surVisualizer/data/my_ref.rds")
usethis::use_data(my_ref, overwrite = TRUE)
