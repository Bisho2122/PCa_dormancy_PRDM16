# Create a vector of required packages
required_packages <- c(
  "DESeq2",
  "ggpubr",
  "ggrepel",
  "EnhancedVolcano",
  "pheatmap",
  "org.Hs.eg.db",
  "org.Mm.eg.db",
  "fplot",
  "clusterProfiler",
  "ggvenn",
  "RColorBrewer",
  "ReactomePA",
  "plyr",
  "GEOquery",
  "dplyr",
  "tidyr",
  "stringr"
)

# Install any packages that are not currently installed
installed_packages <- rownames(installed.packages())
packages_to_install <- setdiff(required_packages, installed_packages)

if (length(packages_to_install) > 0) {
  install.packages(packages_to_install, dependencies = TRUE)
}

# Install Bioconductor packages separately if needed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

bioc_packages <- c(
  "DESeq2",
  "EnhancedVolcano",
  "org.Hs.eg.db",
  "org.Mm.eg.db",
  "clusterProfiler",
  "ReactomePA",
  "GEOquery"
)

for (pkg in bioc_packages) {
  if (!pkg %in% installed_packages) {
    BiocManager::install(pkg)
  }
}
