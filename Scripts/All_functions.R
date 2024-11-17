# Functions ---------------------------------------------------------------
get_top_DE = function(dds,cond_name, top_x = 100, alpha = 0.05,
                      sort_by = "LFC", LFC_thresh = 1){
  DE_genes = results(dds,
                     name = cond_name,
                     altHypothesis = "greaterAbs",
                     lfcThreshold = 0,
                     alpha = 0.1,
                     tidy = T) %>%
    dplyr::filter(padj < alpha) %>%
    dplyr::filter(abs(log2FoldChange) > LFC_thresh)

  n_down = length(which(DE_genes$log2FoldChange < 0))
  n_up = length(which(DE_genes$log2FoldChange > 0))

  if(sort_by == "LFC"){
    DE_genes$rank = abs(DE_genes$log2FoldChange)
  }
  else{
    DE_genes$rank = -log10(DE_genes$padj)
  }


  down_genes = DE_genes %>%
    dplyr::filter(log2FoldChange < 0) %>%
    dplyr::arrange(desc(abs(rank))) %>%
    slice_max(abs(rank),n=min(top_x, n_down)) %>%
    pull(row)

  up_genes = DE_genes %>%
    dplyr::filter(log2FoldChange > 0) %>%
    dplyr::arrange(desc(abs(rank))) %>%
    slice_max(abs(rank),n=min(top_x, n_up)) %>%
    pull(row)

  return(list("up" = up_genes, "down" = down_genes))
}
dds_from_counts = function(raw, meta_df, cell_line = c("RM1", "22RV1", "LAPC4"),
                           combine_human = F,
                           min_grp_size = 3,
                           min_count_per_sample = 10,
                           cell_line_in_colnames = F,
                           remove_reawaken = F){

  colnames(raw)[which(colnames(raw) == "geneName")] = "GeneName"
  sample_names = intersect(colnames(raw), meta_df$Sample)

  colnames(raw) = gsub("[.]", "-", colnames(raw))
  if(cell_line_in_colnames){
    if (!combine_human){
      raw = raw[,c(str_which(colnames(raw), match.arg(cell_line)),
                   which(colnames(raw) == "GeneName"))]
    }
  }
  else{
    raw = raw[,c(sample_names, "GeneName")]
  }

  if (remove_reawaken){
    meta_df = meta_df %>% dplyr::filter(Condition != "reawaken")
  }
  raw$count_sum = rowSums(raw[,sample_names])

  # Remove duplicated genes by keeping the one with the highest count sum
  message("Removing duplicate genes")
  raw = raw %>%
    group_by(GeneName) %>%
    slice_max(count_sum,n=1,with_ties = F) %>%
    ungroup()

  message("Removing genes with zero counts")
  counts = raw %>%
    dplyr::filter(count_sum > 0) %>%
    dplyr::select(GeneName, any_of(sample_names)) %>%
    mutate_if(is.numeric, round) %>%
    column_to_rownames("GeneName") %>%
    as.matrix()

  metadata = meta_df %>%
    dplyr::filter(Sample %in% colnames(counts)) %>%
    column_to_rownames("Sample")

  dds <- DESeqDataSetFromMatrix(countData = counts,
                                colData = metadata,
                                design = ~ Condition)


  smallestGroupSize <- min_grp_size
  keep <- rowSums(counts(dds) >= min_count_per_sample) >= smallestGroupSize
  dds <- dds[keep,]

  return(dds)

}

deseq_workflow_wrapper = function(dds,test = "Wald",fitType = "local",
                                  sfType = "ratio",betaPrior = F,useT = F,
                                  contrast_name = "Condition_Dormant_vs_Active",
                                  altHypothesis = "greaterAbs",
                                  lfcThreshold = 0,
                                  DE_alpha = 0.1){

  dds <- DESeq(object = dds,
               test = test,
               fitType = fitType,
               sfType = sfType,
               betaPrior = betaPrior,
               useT = useT)
  res_obj <- results(dds,
                     name = contrast_name,
                     altHypothesis = altHypothesis,
                     lfcThreshold =lfcThreshold,
                     alpha = DE_alpha,
                     tidy = F)
  res_df <- results(dds,
                    name = contrast_name,
                    altHypothesis = altHypothesis,
                    lfcThreshold =lfcThreshold,
                    alpha = DE_alpha,
                    tidy = T)
  rld <- rlog(dds, blind=FALSE)

  return(list("dds" = dds,
              "DE_obj" = res_obj,
              "DE_df" = res_df,
              "rlog" = rld))
}

fortify_custom = function (model, data, showCategory = 5, by = "geneRatio", split = NULL,
                           includeAll = TRUE) {
  clProf.df <- as.data.frame(model)
  clProf.df$Cluster = as.character(clProf.df$Cluster)
  .split <- split
  if (is.null(showCategory)) {
    result <- clProf.df
  }
  else if (is.numeric(showCategory)) {
    Cluster <- NULL
    topN <- function(res, showCategory) {
      plyr::ddply(.data = res, .variables = .(Cluster), .fun = function(df,
                                                                        N) {
        if (length(df$Count) > N) {
          if (any(colnames(df) == "pvalue")) {
            idx <- order(df$pvalue, decreasing = FALSE)[1:N]
          }
          else {
            idx <- order(df$Count, decreasing = T)[1:N]
          }
          return(df[idx, ])
        }
        else {
          return(df)
        }
      }, N = showCategory)
    }
    if (!is.null(.split) && .split %in% colnames(clProf.df)) {
      lres <- split(clProf.df, as.character(clProf.df[,
                                                      .split]))
      lres <- lapply(lres, topN, showCategory = showCategory)
      result <- do.call("rbind", lres)
    }
    else {
      result <- topN(clProf.df, showCategory)
    }
  }
  else {
    result <- subset(clProf.df, Description %in% showCategory)
  }
  ID <- NULL
  if (includeAll == TRUE) {
    result <- subset(clProf.df, ID %in% result$ID)
  }
  result$Description <- as.character(result$Description)
  GOlevel <- result[, c("ID", "Description")]
  GOlevel <- unique(GOlevel)
  result <- result[result$Count != 0, ]
  result$Description <- factor(result$Description, levels = unique(rev(GOlevel[,
                                                                               2])))
  if (by == "rowPercentage") {
    Description <- Count <- NULL
    result <- ddply(result, .(Description), transform, Percentage = Count/sum(Count),
                    Total = sum(Count))
    x <- mdply(result[, c("Description", "Total")], paste,
               sep = " (")
    y <- sapply(x[, 3], paste, ")", sep = "")
    result$Description <- y
    xx <- result[, c(2, 3)]
    xx <- unique(xx)
    rownames(xx) <- xx[, 1]
    Termlevel <- xx[as.character(GOlevel[, 1]), 2]
    result <- result[, colnames(result) != "Total"]
    result$Description <- factor(result$Description, levels = rev(Termlevel))
  }
  else if (by == "count") {
  }
  else if (by == "geneRatio") {
    gsize <- as.numeric(sub("/\\d+$", "", as.character(result$GeneRatio)))
    gcsize <- as.numeric(sub("^\\d+/", "", as.character(result$GeneRatio)))
    result$GeneRatio <- gsize/gcsize
    cluster <- paste(as.character(result$Cluster), "\n",
                     "(", gcsize, ")", sep = "")
    lv <- unique(cluster)[order(as.numeric(unique(result$Cluster)))]
    result$Cluster <- factor(cluster, levels = lv)
  }
  else {
  }
  return(result)
}

dot_plot_custom = function (object, x = "Cluster", colorBy = "p.adjust", showCategory = 5,
                            by = "geneRatio", size = "geneRatio", split = NULL, includeAll = TRUE,
                            font.size = 12, title = "", label_format = 30, group = FALSE,
                            shape = FALSE, facet_direction = F) {
  color <- NULL
  if (is.null(size))
    size <- by
  df <- fortify_custom(object, showCategory = showCategory, by = size,
                       includeAll = includeAll, split = split)
  if (by != "geneRatio")
    df$GeneRatio <- parse_ratio(df$GeneRatio)
  label_func <- enrichplot:::default_labeller(label_format)
  if (is.function(label_format)) {
    label_func <- label_format
  }
  if (size %in% c("rowPercentage", "count", "geneRatio")) {
    by2 <- switch(size, rowPercentage = "Percentage", count = "Count",
                  geneRatio = "GeneRatio")
  }
  else {
    by2 <- size
  }

  if(facet_direction){
    df$Cluster = gsub("_down", "", df$Cluster)
    df$Cluster = gsub("_up", "", df$Cluster)
  }

  df$Description <- gsub("\n", " ", df$Description)

  p <- ggplot(df, aes_string(x = x, y = "Description", size = by2))
  if (group) {
    p <- p + geom_line(aes_string(color = "Cluster", group = "Cluster"),
                       size = 0.3) + ggnewscale::new_scale_colour()
  }
  if (shape) {
    ggstar <- "ggstar"
    require(ggstar, character.only = TRUE)
    p <- p + ggstar::geom_star(aes_string(starshape = "Cluster",
                                          fill = colorBy)) + scale_fill_continuous(low = "red",
                                                                                   high = "blue", guide = guide_colorbar(reverse = TRUE))
  }
  else {
    p <- p + geom_point(aes_string(color = colorBy))
  }

  p = p + scale_color_continuous(low = "red",
                                 high = "blue",
                                 guide = guide_colorbar(reverse = TRUE)) +
    ylab(NULL) +
    ggtitle(title) +
    DOSE::theme_dose(font.size) +
    scale_size_continuous(range = c(3, 8)) +
    scale_y_discrete(labels = label_func) +
    guides(size = guide_legend(order = 1), color = guide_colorbar(order = 2))

  if (facet_direction) {
    p = p + facet_wrap(~direction, scales = "free_x")
  }
  return(p)
}

get_topvariable_entrez = function(rlog_obj, entrez_mapping, top = 5000){
  topVarGenes <- rowVars(assay(rlog_obj))
  names(topVarGenes) <- rownames(assay(rlog_obj))
  topVarGenes <- sort(topVarGenes, decreasing = TRUE)
  topVarGenes <- names(topVarGenes)[1:top]
  y = entrez_mapping$ENTREZID[entrez_mapping$SYMBOL %in% topVarGenes] %>%
    unique()
  return(y)
}

subset_normalized_data = function(rlog_obj,metadata,genes_of_interest,
                              scale_z = T, add_genes = NULL){
  rlog_data = assay(rlog_obj)
  rownames(rlog_data) = rownames(rlog_data) %>% str_to_upper()
  sel_samples =
    metadata %>%
    dplyr::filter(Sample %in% colnames(rlog_data)) %>%
    dplyr::filter(Condition != "Reawakened") %>%
    pull(Sample)

  dormant_samples = metadata %>%
    dplyr::filter(Sample %in% colnames(rlog_data)) %>%
    dplyr::filter(Condition == "Dormant") %>%
    pull(Sample)

  GOI = genes_of_interest %>% str_to_upper()
  if(!is.null(add_genes)){
    GOI = c(GOI, add_genes)
  }
  sel_genes = intersect(rownames(rlog_data),GOI)
  norm_data = rlog_data[sel_genes,sel_samples]
  if(scale_z){
    norm_data = t(scale(t(norm_data)))
  }
  return(norm_data)
}

curated_pc_calc_corr = function(dataset, gene, cor_method = "spearman") {
  name <- deparse(substitute(dataset))

  # data_subset <- as.data.frame(t(dataset[gene, ]))
  # data_subset$sample <- rownames(data_subset)
  # data_subset[, 1] <- as.numeric(data_subset[, 1])

  dataset <- as.data.frame(t(dataset))


  cor_values <- numeric(0)

  for (i in 1:ncol(dataset)) {
    indvidual_cor_values <- cor(dataset[, i], dataset[, gene], method = cor_method,
                                use = "pairwise.complete.obs")
    cor_values <- c(cor_values, indvidual_cor_values)
  }

  cor_matrix <- data.frame(matrix(NA, nrow = ncol(dataset), ncol = 2))
  colnames <- colnames(dataset)[1:ncol(dataset)]
  colnames <- as.vector(colnames)
  cor_matrix[, 1] <- colnames
  cor_matrix[, 2] <- cor_values
  # paste0("spearman_correlation_coeff_",deparse(substitute(dataset)))
  colnames(cor_matrix) <- c("gene", name)
  return(cor_matrix)
}

ensembl_API_id = function(ens_id, req_type = "archive"){
  server <- "https://rest.ensembl.org"
  ext <- paste0("/", req_type, "/id/", ens_id, "?")

  r <- httr::GET(paste(server, ext, sep = ""), content_type("application/json"))


  if(req_type == "archive"){
    ordered_colnames = c("id","version","is_current", "release", "assembly",
                         "type", "latest")
  }
  if(req_type == "lookup"){
    ordered_colnames = c("id","version","Parent","start", "end","length","biotype",
                         "is_canonical")
  }

  if(r$status_code == 200){
    res = fromJSON(toJSON(content(r)))

    if(str_detect(ens_id, "ENSG")){
      transc_id = sub("[.].*", "", res$canonical_transcript)
      return(ensembl_API_id(ens_id = transc_id,
                            req_type = "lookup"))
    }

    res_comps = unlist(res)
    res_df = matrix(res_comps, ncol=length(res_comps))
    colnames(res_df) = names(res_comps)
    return(res_df[,ordered_colnames])
  }
  else{
    return(NULL)
  }

  # use this if you get a simple nested list back, otherwise inspect its structure
  # head(data.frame(t(sapply(content(r),c))))

}

remove_redund_probes_GSE = function(gene_info, gene_sym, canonic_info = NULL){
  GOI = gene_info[gene_info$GENE_SYMBOL == gene_sym,]
  if(nrow(GOI) == 0){
    stop("Gene not found in info")
  }
  else if (nrow(GOI) == 1){
    return(GOI)
  }
  else{
    if(all(is.na(GOI$ENSEMBL_ID)) || sum(as.integer(GOI$is_current), na.rm = T) == 0){
      # chr_info = GOI %>%
      #   dplyr::select(ID, CHROMOSOMAL_LOCATION) %>%
      #   tidyr::separate(CHROMOSOMAL_LOCATION, into = c("chr", "coord"), sep = ":") %>%
      #   tidyr::separate(coord, into = c("start", "end"), sep = "-")
      # chr_info$length = abs(chr_info$end - chr_info$start)
      # sel_probe = chr_info$ID[which.max(chr_info$length)]
      # GOI = GOI[GOI$ID == sel_probe,]
      return(GOI)
    }
    else{
      if (sum(as.integer(GOI$is_current), na.rm = T) == 1){
        GOI = GOI[GOI$is_current == "1",]
      }
      else{
        if (sum(as.integer(GOI$is_canonical), na.rm = T) == 1){
          GOI = GOI[GOI$is_canonical == "1",]
        }
        else{
          GOI = GOI[GOI$is_current == "1",] %>% na.omit()
          if (nrow(GOI) > 1 & !is.null(canonic_info)){
            if (unique(GOI$Parent) %in% names(canonic_info)){
              GOI$Parent_canonic_len = canonic_info[[unique(GOI$Parent)]][,"length"] %>%
                as.integer()
              GOI$len_diff = abs(as.integer(GOI$length) - GOI$Parent_canonic_len)
              GOI = GOI[which.min(abs(GOI$len_diff)),] %>%
                dplyr::select(-Parent_canonic_len, -len_diff)
            }
          }
        }
      }
      GOI = GOI[!is.na(GOI$ID),]
      return(GOI)
    }
  }
}

