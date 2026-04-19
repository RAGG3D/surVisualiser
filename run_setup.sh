#!/bin/bash
# ============================================================
# Run this script in your terminal to complete package setup:
#   cd /home/nazdaq_44sun/Rbioconductor/surVisualiser
#   bash run_setup.sh
# ============================================================

set -e

echo "=== Step 1: Install Bioconductor dependencies ==="
sudo Rscript -e 'options(download.file.method="wget"); BiocManager::install(c("AnnotationDbi","org.Hs.eg.db"), ask=FALSE, update=FALSE)'

echo "=== Step 2: Generate roxygen2 documentation ==="
cd /home/nazdaq_44sun/Rbioconductor/surVisualiser
Rscript -e 'roxygen2::roxygenise()'

echo "=== Step 3: Build package ==="
cd /home/nazdaq_44sun/Rbioconductor
R CMD build surVisualiser

echo "=== Step 4: Check package ==="
R CMD check surVisualiser_0.99.0.tar.gz

echo "===== ALL DONE ====="
