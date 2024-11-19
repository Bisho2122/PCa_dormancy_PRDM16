library(tidyr)
library(dplyr)
library(stringr)
library(DESeq2)
library(ggpubr)
library(ggrepel)
library(EnhancedVolcano)
library(pheatmap)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(fplot)
library(clusterProfiler)
library(ggvenn)
library(RColorBrewer)
library(ReactomePA)
library(plyr)

source("All_functions.R")
setFplot_page(page = "us")


# Create plots dir --------------------------------------------------------
if (!dir.exists(file.path(getwd(), "Plots"))){
  dir.create(file.path(getwd(), "Plots"))
}
main_path = file.path(getwd(), "Plots")

# Load data ---------------------------------------------------------------
metadata = readxl::read_xlsx("../Data/Sample_metadata.xlsx")
metadata$Condition = str_replace_all(metadata$Condition,
                                     c("dormant" = "Dormant",
                                       "reawaken" = "Reawakened",
                                       "ctrl" = "Active"))
all_DE_res = readRDS("../Data/all_DE_res.rds")

sel_enrichment = readxl::read_xlsx("../Data/Selected_enrichment.xlsx", sheet = 2)
# Main analysis -----------------------------------------------------------
cell_lines = c("RM1", "22RV1", "LAPC4")
all_DE_res = list()

for (c in cell_lines){
  count_file_path = switch(c, "RM1" = "../Data/Count_matrices/Lynch_3277_RM1_mouse_RNAseq_rsem_raw_count_20210920.txt",
                      "22RV1" = "../Data/Count_matrices/all_sample_RSEM_gene_count_geneNameType.20240105.txt",
                      "LAPC4" = "../Data/Count_matrices/all_sample_RSEM_gene_count_geneNameType.20240105.txt")
  counts = read.delim(count_file_path, check.names = F)
  obj = dds_from_counts(raw = counts, meta_df = metadata, combine_human = F,
                        cell_line = c, min_grp_size = 3,
                        cell_line_in_colnames = T,
                        min_count_per_sample = 10)

  all_DE_res[[c]] = deseq_workflow_wrapper(dds = obj)

}

saveRDS(all_DE_res, file = "../Data/all_DE_res.rds")


# Generating PCA plotting data ----------------------------------------------------
RM1_pca = plotPCA(all_DE_res$RM1$rlog, intgroup=c("Condition"))
RM1_pca_data = RM1_pca[["data"]]
RM1_pca_data$grp_name = paste0("Mouse\nRM1")
RM1_pca_data$x_title = RM1_pca$labels$x
RM1_pca_data$y_title = RM1_pca$labels$y

pca_22Rv1 = plotPCA(all_DE_res$`22RV1`$rlog, intgroup=c("Condition"))
pca_22Rv1_data = pca_22Rv1[["data"]]
pca_22Rv1_data$grp_name = paste0("Human\n22Rv1")
pca_22Rv1_data$x_title = pca_22Rv1$labels$x
pca_22Rv1_data$y_title = pca_22Rv1$labels$y

# PCA and volcano plots --------------------------------------------------------
RM1_pca = ggplot(data = RM1_pca_data, aes(x = PC1, y = PC2,
                                          color = Condition)) +
  geom_point(size = 5) +
  theme_pubr(legend = "right") +
  geom_text_repel(mapping = aes(label = name),show.legend = F,
                  max.overlaps = 20, size = 4) +
  ggtitle("Mouse",subtitle = "RM1") +
  theme(text = element_text(size = 20),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, size = 16),
        legend.title = element_blank()) +
  xlab(unique(RM1_pca_data$x_title)) +
  ylab(unique(RM1_pca_data$y_title)) +
  scale_color_manual(values = c("Dormant" = "#FB8514",
                                "Reawakened" = "#6666F6",
                                "Active" = "#0000FF"))

export_graph_start(file.path(main_path, "Fig_2A.pdf"), pt = 10,
                   sideways = T)
print(RM1_pca)
export_graph_end()

pca_22Rv1 = ggplot(data = pca_22Rv1_data, aes(x = PC1, y = PC2,
                                            color = Condition)) +
  geom_point(size = 5) +
  theme_pubr(legend = "none") +
  geom_text_repel(mapping = aes(label = name),show.legend = F,
                  max.overlaps = 20, size = 4) +
  ggtitle("Human",subtitle = "22Rv1") +
  theme(text = element_text(size = 20),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, size = 16),
        legend.title = element_blank()) +
  xlab(unique(pca_22Rv1_data$x_title)) +
  ylab(unique(pca_22Rv1_data$y_title))+
  scale_color_manual(values = c("Dormant" = "#FB8514",
                                "Reawakened" = "#6666F6",
                                "Active" = "#0000FF"))

export_graph_start(file.path(main_path, "Fig_S2A.pdf"), pt = 10,
                   sideways = T)
print(pca_22Rv1)
export_graph_end()

sel_genes = get_top_DE(dds = all_DE_res$RM1$dds,
                       cond_name = "Condition_Dormant_vs_Active",
                       top_x = 10, sort_by = "pval") %>%
  unlist() %>%
  as.character()
mouse_volcano = EnhancedVolcano(all_DE_res$RM1$DE_obj,
                lab = rownames(all_DE_res$RM1$DE_obj),
                selectLab = sel_genes,
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                title = "Mouse",
                subtitle = "RM1",
                drawConnectors = T,
                max.overlaps = Inf) +
  theme_pubr() +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

export_graph_start(file.path(main_path, "Fig_2B.pdf"), pt = 10,
                   sideways = T)
print(mouse_volcano)
export_graph_end()


sel_genes = get_top_DE(dds = all_DE_res$`22RV1`$dds,
                       cond_name = "Condition_Dormant_vs_Active",
                       top_x = 10, sort_by = "pval") %>%
  unlist() %>%
  as.character()
human_volcano = EnhancedVolcano(all_DE_res$`22RV1`$DE_obj,
                                lab = rownames(all_DE_res$`22RV1`$DE_obj),
                                selectLab = sel_genes,
                                x = 'log2FoldChange',
                                y = 'padj',
                                pCutoff = 0.05,
                                title = "Human",
                                subtitle = "22Rv1",
                                drawConnectors = T,
                                max.overlaps = Inf) +
  theme_pubr() +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

export_graph_start(file.path(main_path, "Fig_S2B.pdf"), pt = 10,
                   sideways = T)
print(p)
export_graph_end()

# Prepare data for enrichment ---------------------------------------------
mouse_gene_symbols = all_DE_res$RM1$DE_df$row %>% unique()
mouse_entrez_mapping = clusterProfiler::bitr(mouse_gene_symbols, fromType = "SYMBOL", toType = "ENTREZID",
                                             OrgDb = org.Mm.eg.db)
human_gene_symbols = all_DE_res$`22RV1`$DE_df$row %>% unique()
human_entrez_mapping = clusterProfiler::bitr(human_gene_symbols, fromType = "SYMBOL", toType = "ENTREZID",
                                             OrgDb = org.Hs.eg.db)

mouse_univ = get_topvariable_entrez(rlog_obj = all_DE_res$RM1$rlog,
                                    entrez_mapping = mouse_entrez_mapping,
                                    top = 8000)
human_univ = get_topvariable_entrez(rlog_obj = all_DE_res$`22RV1`$rlog,
                                          entrez_mapping = human_entrez_mapping,
                                          top = 8000)

human_DE = all_DE_res["22RV1"]
clust_comp_input_human = list()
for (c in names(human_DE)){
  DE_res = human_DE[[c]]$DE_df
  DE_res = DE_res %>%
    dplyr::left_join(human_entrez_mapping, by = c("row" = "SYMBOL")) %>%
    dplyr::filter(!is.na(ENTREZID), abs(log2FoldChange) > 1, padj < 0.05)
  down_reg = DE_res %>% dplyr::filter(log2FoldChange < 0) %>% dplyr::pull(ENTREZID)
  up_reg = DE_res %>% dplyr::filter(log2FoldChange > 0) %>% dplyr::pull(ENTREZID)
  clust_comp_input_human[[paste0(c, "_down")]] = intersect(down_reg, human_univ)
  clust_comp_input_human[[paste0(c, "_up")]] = intersect(up_reg, human_univ)
}

mouse_DE = all_DE_res["RM1"]
clust_comp_input_mouse = list()
for (c in names(mouse_DE)){
  DE_res = mouse_DE[[c]]$DE_df
  DE_res = DE_res %>%
    dplyr::left_join(mouse_entrez_mapping, by = c("row" = "SYMBOL")) %>%
    dplyr::filter(!is.na(ENTREZID), abs(log2FoldChange) > 1, padj < 0.05)
  down_reg = DE_res %>% dplyr::filter(log2FoldChange < 0) %>% dplyr::pull(ENTREZID)
  up_reg = DE_res %>% dplyr::filter(log2FoldChange > 0) %>% dplyr::pull(ENTREZID)
  clust_comp_input_mouse[[paste0(c, "_down")]] = intersect(down_reg, mouse_univ)
  clust_comp_input_mouse[[paste0(c, "_up")]] = intersect(up_reg, mouse_univ)
}


# Enrichment analysis -----------------------------------------------------
ck_reactome_22rv1 <- compareCluster(geneCluster = clust_comp_input_human[c("22RV1_down", "22RV1_up")],
                                    fun = enrichPathway,
                     organism = "human",qvalueCutoff = 1,
                     universe = human_univ)

ck_reactome_mouse <- compareCluster(geneCluster = clust_comp_input_mouse, fun = enrichPathway,
                                    organism = "mouse",
                                    universe = mouse_univ)
combined_reactome = ck_reactome_mouse
combined_reactome@compareClusterResult = rbind.data.frame(combined_reactome@compareClusterResult,
                                                          ck_reactome_22rv1@compareClusterResult)
combined_reactome@geneClusters = c(combined_reactome@geneClusters,
                                   ck_reactome_22rv1@geneClusters,)

combined_reactome@compareClusterResult = combined_reactome@compareClusterResult %>%
  tidyr::separate(Cluster, c("Group", "direction"),
                  sep = "_",remove = F)


# Select common enriched pathways -----------------------------------------
enrich_res_all = combined_reactome@compareClusterResult %>%
  tidyr::separate(Cluster, c("Group", "direction"),
                  sep = "_",remove = F)

common_enrich = enrich_res_all %>%
  dplyr::group_by(Description) %>%
  dplyr::summarise(n_clust = n_distinct(Cluster)) %>%
  dplyr::filter(n_clust > 1)

common_enrich_res = enrich_res_all[enrich_res_all$Description %in% common_enrich$Description,]
RM1_down = common_enrich_res$Description[common_enrich_res$Group == "RM1" & common_enrich_res$direction == "down"] %>%
  unique()
`22RV1_down` = common_enrich_res$Description[common_enrich_res$Group == "22RV1" & common_enrich_res$direction == "down"] %>%
  unique()

RM1_22RV1_down = common_enrich_res[common_enrich_res$Description %in% intersect(RM1_down, `22RV1_down`),]

#Redundant pathways were manually removed

curated_RM1_22Rv1_down = read.csv("../Data/Reactome_common_down_RM1_22rv1_curated.csv")
to_rm = setdiff(RM1_22RV1_down$Description, curated_RM1_22Rv1_down$x)

common_enrich_res = common_enrich_res[!common_enrich_res$Description %in% to_rm,]

final_enrich_res = combined_reactome
final_enrich_res@compareClusterResult = final_enrich_res@compareClusterResult[final_enrich_res@compareClusterResult$Description %in% common_enrich_res$Description,]

final_enrich_res@compareClusterResult$direction = str_replace_all(final_enrich_res@compareClusterResult$direction,
                                                       c("down" = "Downregulated", "up" = "Upregulated"))

saveRDS(final_enrich_res, "../Data/enrichment_results_clusterprofiler.rds")
# Plotting enrichment results ---------------------------------------------
final_enrich_res = readRDS("../Data/enrichment_results_clusterprofiler.rds")
final_enrich_res@compareClusterResult$Cluster = str_replace_all(final_enrich_res@compareClusterResult$Cluster,
                                                   c("22RV1_up" = "22Rv1_up",
                                                     "22RV1_down" = "22Rv1_down"))

final_enrich_res@compareClusterResult = final_enrich_res@compareClusterResult[final_enrich_res@compareClusterResult$Description %in% sel_enrichment$Description,]

S2C = dot_plot_custom(final_enrich_res, showCategory = 30,
                    facet_direction = T,x = "Cluster",
                    label_format = 50)
S2C

#PDF exported from R graph window US regular size 8x14 landscape

# Venn diagram ------------------------------------------------------------
venn_list = list()
for (c in names(all_DE_res)){
  DE_res = all_DE_res[[c]]$DE_df
  DE_res = DE_res %>%
    dplyr::filter(log2FoldChange > 1, padj < 0.05)
  venn_list[[c]] = unique(DE_res$row) %>% str_to_upper()
}

mycol = brewer.pal(4, "Pastel2")
p = ggvenn(venn_list,stroke_size = 0.5, set_name_size = 8,
           fill_color = mycol, stroke_color = NA, show_percentage = F,
           fill_alpha = 0.5, text_size = 10)
  # ggtitle(paste0("FDR ", i*100, " % ")) +
  # theme(title = element_text(size = 12, hjust = 0))

export_graph_start(file.path(main_path, "Fig_2C.pdf"), pt = 8,
                   sideways = T)
print(p)
export_graph_end()
# LFC Heatmap common up dormant -------------------------------------------
all_common = Reduce(intersect, venn_list)

dormant_hm_data = list()
for (c in names(all_DE_res)){
  DE_res = all_DE_res[[c]]$DE_df
  DE_res$row = str_to_upper(DE_res$row)
  DE_res = DE_res %>%
    dplyr::filter(row %in% all_common) %>%
    dplyr::select(row, log2FoldChange) %>%
    distinct()
  dormant_hm_data[[c]] = DE_res
}
dormant_hm_data = dormant_hm_data %>% dplyr::bind_rows(.id = "cell_line")
dormant_hm_data$cell_line = gsub("22RV1", "22Rv1", dormant_hm_data$cell_line)

dormant_hm_data = dormant_hm_data %>%
  dplyr::mutate(cell_line = paste0("D-", cell_line))
wide_mat = dormant_hm_data %>%
  spread(cell_line, log2FoldChange) %>%
  column_to_rownames("row") %>%
  as.matrix()

custom_order = c("PRDM16", "NR4A2",
                 setdiff(rownames(wide_mat), c("PRDM16", "NR4A2")))
wide_mat = wide_mat[custom_order,]

label_mat = dormant_hm_data %>%
  spread(cell_line, log2FoldChange) %>%
  column_to_rownames("row") %>%
  as.matrix()
label_mat = round(label_mat, digits = 1)
label_mat = label_mat[custom_order,]

p = pheatmap(wide_mat, cluster_rows = F, cluster_cols = T, show_rownames = T,
         color = colorRampPalette(brewer.pal(n = 6, name =
                                                   "Reds"))(10),
         display_numbers = label_mat,number_color = "black",angle_col = 0,
         fontsize = 12,border_color = "white")

export_graph_start(file.path(main_path, "Fig_2D.pdf"), pt = 8,
                   sideways = T)
print(p)
export_graph_end()


# Heatmap Known markers ----------------------------------------------
cell_lines = c("RM1", "22RV1")

GOI = c("CDKN1A", "CDKN1B","IRF7", "IRF9","SOX9",
        "TGFBR2", "TGFBR3","NDRG1","BMP2","BMP7",
        "WNT5A", "IL15","NR2F1", "MKI67")

plot_data = list()
for (c in cell_lines){
  dorm_cell_sign = subset_normalized_data(rlog_obj = all_DE_res[[c]]$rlog,
                                          metadata = metadata,
                                          genes_of_interest = GOI,
                                          scale_z = T)
  plot_data[[c]] = dorm_cell_sign %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    gather(-gene, key = "Sample", value = "zscore")
}
plot_data = plot_data %>% dplyr::bind_rows(.id = "cell_line")


wide_mat = plot_data %>%
  dplyr::select(-cell_line) %>%
  dplyr::distinct() %>%
  pivot_wider(names_from = Sample, values_from = zscore) %>%
  column_to_rownames("gene") %>%
  as.matrix()
wide_mat = wide_mat[GOI,]



annot_df = plot_data %>%
  dplyr::select(cell_line, Sample) %>%
  distinct() %>%
  dplyr::left_join(metadata) %>%
  dplyr::select(-Species) %>%
  column_to_rownames("Sample") %>%
  dplyr::arrange(Condition)

annot_df = annot_df[colnames(wide_mat),]
annot_df = annot_df %>% relocate(Condition, .before = cell_line)
annot_df$cell_line = str_replace_all(annot_df$cell_line, c("22RV1" = "22Rv1"))

annot_colors = list(
  cell_line = c("RM1" = "#8DD3C7", "22Rv1" = "#FB8072"),
  Condition = c("Dormant" = "#FB8514",
                "Active" = "#0000FF")
  )

p = ComplexHeatmap::pheatmap(wide_mat, cluster_rows = F, cluster_cols = F, show_colnames = F,
             color = colorRampPalette(rev(brewer.pal(n = 9, name =
                                                   "RdBu")))(10),
             annotation_col = annot_df, fontsize = 10,border_color = "white",
             na_col = "white", annotation_colors = annot_colors,annotation_legend = T,
             annotation_names_row = F,
             heatmap_legend_param = list(title = "Scaled Expression"))

export_graph_start(file.path(main_path, "Fig_S2D"), pt = 10,
                   sideways = T)
print(p)
export_graph_end()


