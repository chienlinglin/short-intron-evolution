library(dplyr)
library(igraph)

# ============================================================
# Figure 6B: Protein KD correlation network
#
# Input:
#   fig6B_proteon_class.csv
#   both_80_allEvent_edgelist.tab
#
# Output:
#   160 protein all events.pdf
# ============================================================

# ------------------------------------------------------------
# 1. Load protein classification and correlation data
# ------------------------------------------------------------

class <- read.csv("fig6B_proteon_class.csv")

Q1 <- read.table(
  "both_80_allEvent_edgelist.tab",
  header = TRUE
)

# Retain strong and statistically significant correlations
Q1top <- Q1 %>%
  filter(absCor >= 0.5, fdr < 0.2)

# Prepare edge attributes
edge_df <- Q1top[, c(
  7, 8,
  which(colnames(Q1top) == "absCor"),
  which(colnames(Q1top) == "cor")
)]

colnames(edge_df)[1:2] <- c("from", "to")

# Build undirected network
g <- graph_from_data_frame(
  d = edge_df,
  directed = FALSE
)

# ------------------------------------------------------------
# 2. Community detection and network layout
# ------------------------------------------------------------

cfg <- cluster_fast_greedy(
  g,
  modularity = TRUE,
  weights = E(g)$weight
)

# Amplify within-cluster edge weights and attenuate
# between-cluster edge weights for layout visualization
memb <- membership(cfg)
el <- as_edgelist(g, names = FALSE)
same <- memb[el[, 1]] == memb[el[, 2]]

w_lay <- ifelse(same, 5, 0.5)

set.seed(42)
lay <- layout_with_fr(
  g,
  weights = w_lay,
  niter = 5000
)

rownames(lay) <- V(g)$name

# ------------------------------------------------------------
# 3. Node colors and legend
# ------------------------------------------------------------

ii <- match(V(g)$name, class$protein)

V(g)$label_color <- class$color[ii]
V(g)$label_color[is.na(V(g)$label_color)] <- "grey"

legend_data <- class %>%
  select(type, color) %>%
  distinct() %>%
  na.omit()

# ------------------------------------------------------------
# 4. Edge widths and colors
# ------------------------------------------------------------

abs_cor_values <- E(g)$absCor
cor_values <- E(g)$cor

cor_min <- min(abs_cor_values)
cor_max <- max(abs_cor_values)

# Map correlation strength to edge width (1–5)
E(g)$width <- 1 +
  (abs_cor_values - cor_min) /
  (cor_max - cor_min) * 4

# Map correlation strength to color intensity
color_indices <- as.integer(
  (abs_cor_values - cor_min) /
    (cor_max - cor_min) * 99
) + 1

pos_palette <- colorRampPalette(
  c("#FFB6C1", "#B22222")
)(100)

neg_palette <- colorRampPalette(
  c("#ADD8E6", "#00008B")
)(100)

edge_colors <- character(ecount(g))

pos_indices <- which(cor_values > 0)
neg_indices <- which(cor_values < 0)

edge_colors[pos_indices] <-
  pos_palette[color_indices[pos_indices]]

edge_colors[neg_indices] <-
  neg_palette[color_indices[neg_indices]]

E(g)$color <- adjustcolor(
  edge_colors,
  alpha.f = 0.6
)

# ------------------------------------------------------------
# 5. Community background colors
# ------------------------------------------------------------

community_colors <- rep(
  "lightsalmon",
  length(table(cfg$membership))
)

community_colors <- adjustcolor(
  community_colors,
  alpha.f = 0.15
)

# ------------------------------------------------------------
# 6. Plot Figure 6B
# ------------------------------------------------------------

dset <- "160 protein all events"

pdf(file = file.path(paste0(dset, ".pdf")),
    width = 8,height = 7)

plot(
  cfg,
  g,
  layout = lay,
  
  vertex.size = 5,
  col = V(g)$label_color,
  vertex.label.cex = 0.7,
  vertex.label.color = "black",
  vertex.label.dist = 1,
  
  edge.width = E(g)$width,
  edge.color = E(g)$color,
  
  main = paste0(dset, "\nabsCor >= 0.5 and fdr < 0.2"),
  
  mark.col = community_colors[-1],
  mark.border = NA
)

legend(
  "bottomleft",
  legend = legend_data$type,
  pch = 19,
  col = legend_data$color,
  pt.cex = 1,
  bty = "n",
  ncol = 1,
  title = "Gene Functional Class",
  cex = 0.8
)

dev.off()
