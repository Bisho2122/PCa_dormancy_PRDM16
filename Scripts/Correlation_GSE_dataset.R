library(tidyverse)
library(limma)
library(GEOquery)
library(DESeq2)
library(ggpubr)
library(ggvenn)
library(RColorBrewer)
library(fplot)

all_DE_res = readRDS("../Data/all_DE_res.rds")

source("All_functions.R")

# Calculate correlation with NED  ------------------------------------------------------
raw_non_ratio = read.csv("../Data/gene_sample_norm_expr_GSE48995.csv")
gset <- getGEO("GSE48995", GSEMatrix =TRUE, getGPL=T,AnnotGPL = T)

metadata = phenoData(gset$GSE48995_series_matrix.txt.gz)@data
metadata$patient_id = sub("_.*", "", metadata$source_name_ch2)

NED_samples = metadata$geo_accession[metadata$`patient status:ch2` == "NED"]

exprs_data =raw_non_ratio %>%
  column_to_rownames("X") %>%
  as.matrix()

GSE_corr = curated_pc_calc_corr(dataset = exprs_data[,NED_samples], gene = "PRDM16")
colnames(GSE_corr)[2] = "NED_ILSA"

Nasr_RM1 = as.data.frame(assay(all_DE_res$RM1$rlog))
Nasr_22Rv1 = as.data.frame(assay(all_DE_res$`22RV1`$rlog))

Nasr_22Rv1_corr = curated_pc_calc_corr(dataset = Nasr_22Rv1, "PRDM16")
Nasr_RM1_corr = curated_pc_calc_corr(dataset = Nasr_RM1, "Prdm16") %>%
  mutate(gene = str_to_upper(gene))

all_coefs = purrr:::reduce(list(Nasr_22Rv1_corr, Nasr_RM1_corr,
                                GSE_corr),
                           dplyr::left_join, by = "gene")

interesting_genes = all_coefs %>%
  na.omit() %>%
  gather(key = "Dataset", value = "coef", -gene) %>%
  dplyr::group_by(gene) %>%
  dplyr::summarise(mean = mean(coef), std = sd(coef)) %>%
  dplyr::mutate(std_z = as.numeric(scale(std))) %>%
  dplyr::filter(abs(mean) > 0.5, std_z < qnorm(0.05)) %>%
  pull(gene)

interesting_gene_coefs = all_coefs[all_coefs$gene %in% interesting_genes,]
interesting_gene_coefs$high_corr_common = "Yes"

all_coef_common = all_coefs %>%
  na.omit() %>%
  dplyr::left_join(interesting_gene_coefs)
all_coef_common$high_corr_common[is.na(all_coef_common$high_corr_common)] = "No"

write.csv(all_coef_common, "../Data/GSE_NED_NASR_coefs_raw_noratio.csv",
          row.names = F)


# Plot correlation ----------------------------------------------
all_coefs_no_NA = na.omit(all_coefs)
all_coef_gathered = gather(all_coefs_no_NA, key = "Dataset", value = "coef", -gene)
sel_genes = c("WNT5B", "MXI1", "ITGA5", "BMP2", "BMP4","IL15", "MCM10",
              "CDC42", "PCNA", "MCM6", "CDK4", "BUB1")
sel_ceof = all_coef_gathered[all_coef_gathered$gene %in% sel_genes,]
sel_ceof$gene = factor(sel_ceof$gene, levels = sel_genes)

p2 = sel_ceof %>%
  ggplot(aes(x = gene, y = coef, group = Dataset, color = Dataset)) +
  geom_line(aes(group = gene), color = "black", size = 1) +
  geom_point(aes(y = coef), size = 7, alpha = 0.7) +
  theme_pubr() +
  rotate_x_text(angle = 30) +
  geom_hline(yintercept = 0, linetype = "dotted", size = 1) +
  xlab("") +
  ylab("Spearman Coefficient") +
  scale_color_brewer(palette = "Set3") +
  scale_y_continuous(n.breaks = 13,limits = c(-1,1))

export_graph_start(file.path(main_path, "S5D.pdf"),
                   pt = 10,
                   sideways = T)
print(p2)
export_graph_end()


# Venn diagram ------------------------------------------------------------
venn_data = read.csv("../Data/GSE_NED_NASR_coefs_raw_noratio.csv")
venn_data = venn_data %>%
  gather(key = "Group", value = "Coef", -high_corr_common, -gene)

grps = unique(venn_data$Group)

pos_corr = list()
neg_corr = list()
for (grp in grps){
  filt = venn_data[venn_data$Group == grp,]
  pos_corr[[grp]] = filt$gene[filt$Coef >= 0.4] %>% unique()
  neg_corr[[grp]] = filt$gene[filt$Coef <= -0.4] %>% unique()
}

all_common_pos = Reduce(intersect, pos_corr)
all_common_neg = Reduce(intersect, neg_corr)

mycol = brewer.pal(4, "Pastel2")
pos_venn = ggvenn(pos_corr,stroke_size = 0.5, set_name_size = 8,
                  fill_color = mycol, stroke_color = NA, show_percentage = F,
                  fill_alpha = 0.5, text_size = 10)

export_graph_start(file.path(main_path, "S5A_1"), pt = 8,
                   sideways = T)
print(pos_venn)
export_graph_end()

neg_venn = ggvenn(neg_corr,stroke_size = 0.5, set_name_size = 8,
                  fill_color = mycol, stroke_color = NA, show_percentage = F,
                  fill_alpha = 0.5, text_size = 10)

export_graph_start(file.path(main_path, "S5A_2"), pt = 8,
                   sideways = T)
print(neg_venn)
export_graph_end()
