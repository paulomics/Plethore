#' @export
Generate_QC_SeuObj <- function(seurat_object,
                               clustering_tested,
                               sample_IDs,
                               level = 3) {
  if (!inherits(seurat_object, "Seurat")) {
    stop("'seurat_object' must be a Seurat object.")
  }

  if (!is.character(clustering_tested) ||
      length(clustering_tested) != 1 ||
      is.na(clustering_tested)) {
    stop("'clustering_tested' must be a character string of length 1.")
  }

  if (!clustering_tested %in% colnames(seurat_object@meta.data)) {
    stop(sprintf(
      "'%s' was not found in seurat_object@meta.data.",
      clustering_tested
    ))
  }

  if (all(is.na(seurat_object@meta.data[[clustering_tested]]))) {
    stop(sprintf("Column '%s' contains only NA values.", clustering_tested))
  }

  if (!is.character(sample_IDs) ||
      length(sample_IDs) != 1 ||
      is.na(sample_IDs)) {
    stop("'sample_IDs' must be a character string of length 1.")
  }

  if (!sample_IDs %in% colnames(seurat_object@meta.data)) {
    stop(sprintf("'%s' was not found in seurat_object@meta.data.", sample_IDs))
  }

  if (all(is.na(seurat_object@meta.data[[sample_IDs]]))) {
    stop(sprintf("Column '%s' contains only NA values.", sample_IDs))
  }

  if (!is.numeric(level) ||
      length(level) != 1 ||
      is.na(level) ||
      level <= 0 ||
      level %% 1 != 0) {
    stop("'level' must be a strictly positive integer.")
  }

  level_1 <- strrep("#", level)
  level_2 <- strrep("#", level + 1)
  cat("\n\n::: {.panel-tabset}\n\n")
  cat("\n\n", level_1, " By Clusters\n\n", sep = "")
  cat("\n\n::: {.panel-tabset}\n\n")
  cat("\n\n", level_2, " Feature plot\n\n", sep = "")
  p <- Seurat::DimPlot(seurat_object, reduction = "umap", label = TRUE)
  print(p)
  cat("\n\n", level_2, " Cell count by clusters\n\n", sep = "")
  table_count_cell <- table(seurat_object@meta.data[[clustering_tested]])
  df <- as.data.frame(table_count_cell)
  colnames(df) <- c("Cluster", "Count")

  p <- ggplot2::ggplot(df, aes(x = factor(Cluster), y = Count)) +
    geom_col(fill = "#FCA192") +
    geom_text(aes(label = Count), vjust = -0.3, size = 4) +
    labs(x = "Cluster", y = "Number of cells", title = "Cell counts per cluster")  +
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(
        size = 18,
        face = "bold",
        hjust = 0.5
      ),
      axis.title.y = element_text(size = 15, face = "bold"),
      axis.title.x = element_text(size = 15, face = "bold"),
      axis.text.y = element_text(
        size = 12,
        angle = 0,
        hjust = 1,
        face = "bold",
        colour = "black"
      ),
      axis.text.x = element_text(
        size = 12,
        angle = 45,
        hjust = 1,
        face = "bold",
        colour = "black"
      ),
      axis.line = element_line(linewidth = 0.8),
      axis.ticks = element_line(linewidth = 0.7),
      legend.position = "none"
    )  +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x = "", y = "Number of cells")
  print(p)

  cat("\n\n::: \n\n")


  cat("\n\n", level_1, " By sample\n\n", sep = "")
  cat("\n\n::: {.panel-tabset}\n\n")
  cat("\n\n", level_2, " Feature plot\n\n", sep = "")
  p <- Seurat::DimPlot(
    seurat_object,
    reduction = "umap",
    group.by = sample_IDs,
    label.box = TRUE
  )
  print(p)
  cat("\n\n", level_2, " Cell count by samples\n\n", sep = "")
  table_count_cell <- table(seurat_object@meta.data[[sample_IDs]], seurat_object@meta.data[[clustering_tested]])
  p <- pheatmap::pheatmap(
    table_count_cell,
    border_color = "white",
    scale = "none",
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    legend = FALSE,
    show_rownames = TRUE,
    show_colnames = TRUE,
    display_numbers = TRUE,
    fontsize = 12,
    fontsize_number = 0.8 * 12,
    number_format = "%.0f",
    angle_col = "0"
  )
  print(p)
  cat("\n\n::: \n\n")

  stat.list <- c("nFeature_RNA",
                 "nCount_RNA",
                 "percent.mt",
                 "S.Score",
                 "G2M.Score")
  title.list <- c(
    "By Gene count",
    "By UMI count",
    "By % MT genes expression",
    "By Cell cylce S",
    "By Cell cylce G2M"
  )

  for (i in seq_along(stat.list)) {
    cat("\n\n", level_1, " ", title.list[i], "\n\n", sep = "")
    p <- Seurat::FeaturePlot(seurat_object,
                             features = stat.list[i],
                             reduction = "umap",
    )
    v <- Seurat::VlnPlot(
      seurat_object,
      features = stat.list[i],
      group.by = clustering_tested,
      slot = "data",
      assay = "RNA"
    )
    cat("\n\n::: {.panel-tabset}\n\n")
    cat("\n\n", level_2, " Violin plot\n\n", sep = "")
    print(v)
    cat("\n\n", level_2, " Feature plot\n\n", sep = "")
    print(p)
    cat("\n\n::: \n\n")
  }

  cat("\n\n::: \n\n")
}
