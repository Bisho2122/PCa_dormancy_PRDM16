library(tidyverse)
library(decoupleR)
library(pheatmap)
library(DESeq2)
library(fplot)

source("All_functions.R")

# Create plots dir --------------------------------------------------------
if (!dir.exists(file.path(getwd(), "Plots"))){
  dir.create(file.path(getwd(), "Plots"))
}
main_path = file.path(getwd(), "Plots")
# Load data ---------------------------------------------------------------

all_DE_res = readRDS("../Data/all_DE_res.rds")


# Run TF inference --------------------------------------------------------
human_regulons = decoupleR::get_collectri(organism=9606L, genesymbols=TRUE, loops=TRUE)
mouse_regulons = decoupleR::get_collectri(organism=10090, genesymbols=TRUE, loops=TRUE)

cell_lines = c("RM1", "22RV1")

TF_inf_res = list()
for (c in cell_lines){
  if(c == "RM1"){
    bg = mouse_regulons
  }
  else{
    bg = human_regulons
  }
  TF_inf_res[[c]] = Run_TF_inf(DE_obj = all_DE_res[[c]],
                               regulon_bg = bg)
}


# Plotting ----------------------------------------------------------------
GOI_sel = c("IRF1", "IRF9", "FOXO3","TP53", "RB1", "STAT6","JUND",
            "E2F1","E2F2","E2F3","E2F4", "MYCN", "MYC","FOXM1")

p1 = plot_TF_inf_contrasts(TF_inf_res, top_n_TF = 30,
                           n_common_cell = 2,remove_LAPC4 = T,
                           sel_genes = GOI_sel)

export_graph_start(file.path(main_path, "S2E.pdf"),
                   pt = 10,
                   sideways = T)
print(p1)
export_graph_end()
