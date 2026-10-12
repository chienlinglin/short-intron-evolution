library(dplyr)
library(igraph)

# ============================================================
# Figure 6C, D: Network visualization of intron-retention KD correlations
#
# Input:
#   class_colors_2020.txt
#   IR_LongIntrons_edgelist.tab
#   IR_shortIntrons_edgelist.tab
#
# Output:
#   Longer introns IR envents.pdf
#   Short introns IR envents.pdf
# ============================================================


# ------------------------------------------------------------
# 1. Load functional classification
# ------------------------------------------------------------

class <- read.delim(
  "class_colors_2020.txt",
  header = TRUE
)

# Consolidate selected classification colors
rows_to_change <- grepl("^#c5dd[5-6]", class$color)
class$color[rows_to_change] <- "#c5dd51"

legend_data <- class %>%
  dplyr::select(`Class...family.small`, color) %>%
  distinct() %>%
  na.omit()

# ------------------------------------------------------------
# 2. Reusable network plotting function
# ------------------------------------------------------------

plot_ir_network <- function(
    edgelist_file,
    output_file,
    plot_title,
    flip_x = FALSE,
    legend_position = "topright",
    abs_cor_cutoff = 0.5,
    fdr_cutoff = 0.2
) {
  
  # Load and filter correlation edges
  Q1 <- read.table(edgelist_file, header = TRUE)
  
  Q1top <- Q1 %>%
    filter(
      absCor >= abs_cor_cutoff,
      fdr < fdr_cutoff
    )
  
  # Preserve the original source/target column selection
  edge_df <- Q1top[, c(
    7, 8,
    which(colnames(Q1top) == "absCor"),
    which(colnames(Q1top) == "cor")
  )]
  
  colnames(edge_df)[1:2] <- c("from", "to")
  
  # Build undirected graph
  g <- graph_from_data_frame(
    d = edge_df,
    directed = FALSE
  )
  
  # ----------------------------------------------------------
  # Community detection
  # ----------------------------------------------------------
  
  cfg <- cluster_fast_greedy(
    g,
    modularity = TRUE,
    weights = E(g)$weight
  )
  
  # ----------------------------------------------------------
  # Network layout
  # ----------------------------------------------------------
  
  memb <- membership(cfg)
  el <- as_edgelist(g, names = FALSE)
  same <- memb[el[, 1]] == memb[el[, 2]]
  
  # Amplify within-cluster weights and attenuate between-cluster weights
  w_lay <- ifelse(same, 5, 0.5)
  
  set.seed(42)
  lay <- layout_with_fr(
    g,
    weights = w_lay,
    niter = 5000
  )
  
  rownames(lay) <- V(g)$name
  
  if (flip_x) {
    lay[, 1] <- -lay[, 1]
  }
  
  # ----------------------------------------------------------
  # Node colors
  # ----------------------------------------------------------
  
  ii <- match(V(g)$name, class$Gene.Symbol)
  
  V(g)$label_color <- class$color[ii]
  V(g)$label_color[is.na(V(g)$label_color)] <- "grey"
  
  # ----------------------------------------------------------
  # Edge widths
  # ----------------------------------------------------------
  
  min_width <- 1
  max_width <- 10
  
  abs_cor_values <- E(g)$absCor
  cor_values <- E(g)$cor
  
  cor_min <- min(abs_cor_values)
  cor_max <- max(abs_cor_values)
  
  E(g)$width <- min_width +
    (abs_cor_values - cor_min) /
    (cor_max - cor_min) *
    (max_width - min_width)
  
  # ----------------------------------------------------------
  # Edge colors
  # Positive: pink to red
  # Negative: light blue to dark blue
  # ----------------------------------------------------------
  
  pos_palette <- colorRampPalette(
    c("#FFB6C1", "#B22222")
  )(100)
  
  neg_palette <- colorRampPalette(
    c("#ADD8E6", "#00008B")
  )(100)
  
  color_indices <- as.integer(
    (abs_cor_values - cor_min) /
      (cor_max - cor_min) * 99
  ) + 1
  
  edge_colors <- character(ecount(g))
  
  pos_indices <- which(cor_values > 0)
  neg_indices <- which(cor_values < 0)
  
  edge_colors[pos_indices] <-
    pos_palette[color_indices[pos_indices]]
  
  edge_colors[neg_indices] <-
    neg_palette[color_indices[neg_indices]]
  
  E(g)$color <- adjustcolor(
    edge_colors,
    alpha.f = 0.2
  )
  
  # ----------------------------------------------------------
  # Community background colors
  # ----------------------------------------------------------
  
  colorColsLinks <- rep(
    "lightsalmon",
    length(table(cfg$membership))
  )
  
  colorColsLinks <- adjustcolor(
    colorColsLinks,
    alpha.f = 0.15
  )
  
  # ----------------------------------------------------------
  # Plot and save
  # ----------------------------------------------------------
  
  pdf(output_file, width = 12, height = 12)
  
  plot(
    cfg,
    g,
    layout = lay,
    
    vertex.size = 3,
    col = V(g)$label_color,
    vertex.label.cex = 1,
    vertex.label.color = "black",
    vertex.label.dist = 0.5,
    
    edge.width = E(g)$width,
    edge.color = E(g)$color,
    
    main = paste0(
      plot_title,
      "\nabsCor >= ", abs_cor_cutoff,
      " and fdr < ", fdr_cutoff
    ),
    
    mark.col = colorColsLinks[-1],
    mark.border = NA
  )
  
  legend(
    legend_position,
    legend = legend_data$`Class...family.small`,
    pch = 19,
    col = legend_data$color,
    pt.cex = 1,
    bty = "n",
    ncol = 2,
    title = "Gene Functional Class",
    cex = 1,
    x.intersp = 0.5,
    y.intersp = 1.0,
    text.width = 0.3
  )
  
  dev.off()
  
  invisible(list(
    graph = g,
    community = cfg,
    layout = lay
  ))
}



# ------------------------------------------------------------
# 3. Plot Figure 6C: Long-intron network
# ------------------------------------------------------------

long_network <- plot_ir_network(
  edgelist_file = "IR_LongIntrons_edgelist.tab",
  output_file = "Longer introns IR envents.pdf",
  plot_title = "Longer introns IR envents",
  flip_x = FALSE,
  legend_position = "topright"
)

# ------------------------------------------------------------
# 4. Plot Figure 6D: Short-intron network
# ------------------------------------------------------------

short_network <- plot_ir_network(
  edgelist_file = "IR_shortIntrons_edgelist.tab",
  output_file = "Short introns IR envents.pdf",
  plot_title = "Short introns IR envents",
  flip_x = TRUE,
  legend_position = "topleft"
)
