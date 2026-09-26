# =============================================================================
#  THE HUMAN DISEASOME — COMMUNITY DETECTION FOCUS
#  MS6032 Networks and Complex Systems
#  Removed: Percolation
#  Added:   Multi-algorithm comparison (Louvain vs Walktrap vs
#           Fast Greedy vs Label Propagation vs Infomap)
# =============================================================================

required <- c("igraph", "data.table", "Matrix")
for (pkg in required) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}

dir.create("diseasome_output1", showWarnings = FALSE)
out <- function(name) file.path("diseasome_output1", name)

set_theme <- function() {
  par(bg = "white", col.axis = "black", col.lab = "black",
      col.main = "black", fg = "black", mar = c(5, 5, 4, 2))
}

fine_palette <- colorRampPalette(c("#88CCEE","#00CCFF","#FFCC00","#FF0000"))

# =============================================================================
# 1  LOAD & FILTER
# =============================================================================
cat("\n[1/6] Loading data...\n")

file_path <- "D:/Network Class Lab/Presentation 2/DG-Miner_miner-disease-gene.tsv"
raw_data  <- fread(file_path, col.names = c("disease","gene"), showProgress = FALSE)

n_raw_assoc   <- nrow(raw_data)
n_raw_disease <- uniqueN(raw_data$disease)
n_raw_gene    <- uniqueN(raw_data$gene)

cat(sprintf("  Raw associations : %s\n", format(n_raw_assoc,   big.mark = ",")))
cat(sprintf("  Unique diseases  : %s\n", format(n_raw_disease, big.mark = ",")))
cat(sprintf("  Unique genes     : %s\n", format(n_raw_gene,    big.mark = ",")))

top5_dis        <- raw_data[, .N, by = disease][order(-N)][1:5]$disease
bip_sample      <- raw_data[disease %in% top5_dis][, head(.SD, 5), by = disease]
disease_deg_bip <- raw_data[, .N, by = disease]
setnames(disease_deg_bip, "N", "n_genes")

gene_freq <- raw_data[, .N, by = gene][N > 1]
raw_data  <- raw_data[gene %in% gene_freq$gene]
cat(sprintf("  After singleton filter: %s associations remain\n",
            format(nrow(raw_data), big.mark = ",")))

# =============================================================================
# 2  SPARSE MATRIX PROJECTION  P = B^T * B
# =============================================================================
cat("\n[2/6] Sparse projection...\n")

genes    <- unique(raw_data$gene)
diseases <- unique(raw_data$disease)

B <- sparseMatrix(
  i        = as.integer(factor(raw_data$gene,    levels = genes)),
  j        = as.integer(factor(raw_data$disease, levels = diseases)),
  x        = 1L,
  dims     = c(length(genes), length(diseases)),
  dimnames = list(genes, diseases)
)
cat(sprintf("  B: %d genes x %d diseases\n", nrow(B), ncol(B)))

P_raw       <- t(B) %*% B
diag(P_raw) <- 0L
rm(B, raw_data, gene_freq); gc()

# =============================================================================
# 3  TWO NETWORKS
#    full_g  — honest academic metrics (t = 5)
#    g       — stratified 150-node subgraph for visualisation (t = 1000)
# =============================================================================
cat("\n[3/6] Building networks...\n")

# Stats network (t = 5)
P5     <- P_raw; P5@x[P5@x < 5L] <- 0L; P5 <- drop0(P5)
g5_all <- graph_from_adjacency_matrix(P5, mode = "undirected", weighted = TRUE)
comps5 <- decompose(g5_all)
full_g <- comps5[[which.max(sapply(comps5, vcount))]]
deg_full <- degree(full_g)
cat(sprintf("  Stats network  (t=5)  : %d nodes, %d edges\n",
            vcount(full_g), ecount(full_g)))
rm(g5_all, comps5, P5); gc()

# Visual backbone (auto-escalating threshold)
THRESHOLD_PLOT <- 500L
P_vis   <- P_raw; P_vis@x[P_vis@x < THRESHOLD_PLOT] <- 0L; P_vis <- drop0(P_vis)
g_vis   <- graph_from_adjacency_matrix(P_vis, mode = "undirected", weighted = TRUE)
comps_v <- decompose(g_vis)
g_back  <- comps_v[[which.max(sapply(comps_v, vcount))]]
avg_d   <- mean(degree(g_back)); n_back <- vcount(g_back)
cat(sprintf("  Backbone (t=%d): %d nodes, avg_deg=%.1f\n", THRESHOLD_PLOT, n_back, avg_d))
rm(g_vis, comps_v, P_vis); gc()

if (avg_d > 0.80 * n_back) {
  THRESHOLD_PLOT <- 1000L
  cat(sprintf("  Too dense — raising to t=%d\n", THRESHOLD_PLOT))
  P_vis2   <- P_raw; P_vis2@x[P_vis2@x < THRESHOLD_PLOT] <- 0L; P_vis2 <- drop0(P_vis2)
  g_vis2   <- graph_from_adjacency_matrix(P_vis2, mode = "undirected", weighted = TRUE)
  comps_v2 <- decompose(g_vis2)
  g_back   <- comps_v2[[which.max(sapply(comps_v2, vcount))]]
  cat(sprintf("  New backbone (t=%d): %d nodes, %d edges\n",
              THRESHOLD_PLOT, vcount(g_back), ecount(g_back)))
  rm(g_vis2, comps_v2, P_vis2); gc()
}
rm(P_raw); gc()

# Stratified sample: 20 super-hubs + 50 mid + 80 periphery
deg_back <- degree(g_back); n_back <- vcount(g_back)
q90 <- quantile(deg_back, 0.90); q50 <- quantile(deg_back, 0.50)
idx_top <- which(deg_back >= q90)
idx_mid <- which(deg_back >= q50 & deg_back < q90)
idx_bot <- which(deg_back <  q50)

set.seed(42)
plot_v <- c(sample(idx_top, min(20L, length(idx_top))),
            sample(idx_mid, min(50L, length(idx_mid))),
            sample(idx_bot, min(80L, length(idx_bot))))

g <- induced_subgraph(g_back, plot_v)
cat(sprintf("  Plot subgraph (stratified): %d nodes, %d edges\n", vcount(g), ecount(g)))
rm(g_back); gc()

# =============================================================================
# 4  MULTI-ALGORITHM COMMUNITY DETECTION COMPARISON
#
#    All five algorithms run on the SAME plot subgraph (g, 150 nodes)
#    for a fair, apples-to-apples comparison.
#
#    Algorithms compared:
#      1. Louvain           — greedy modularity optimisation; best Q; O(n log n)
#      2. Walktrap          — random-walk similarity; good for overlapping communities
#      3. Fast Greedy       — hierarchical agglomeration; fast but greedy
#      4. Label Propagation — near-linear; stochastic; may vary between runs
#      5. Infomap           — information-theoretic; minimises description length
# =============================================================================
cat("\n[4/6] Multi-algorithm community detection...\n")

run_algo <- function(name, fn, g) {
  t0  <- proc.time()["elapsed"]
  res <- fn(g)
  dt  <- round(proc.time()["elapsed"] - t0, 3)
  Q   <- round(modularity(res), 4)
  nc  <- max(membership(res))
  cat(sprintf("  %-20s  Q=%.4f  communities=%d  time=%.3fs\n", name, Q, nc, dt))
  list(name = name, comm = res, Q = Q, nc = nc, time = dt)
}

set.seed(42)
algos <- list(
  run_algo("Louvain",            function(g) cluster_louvain(g),           g),
  run_algo("Walktrap",           function(g) cluster_walktrap(g),          g),
  run_algo("Fast Greedy",        function(g) cluster_fast_greedy(g),       g),
  run_algo("Label Prop.",  function(g) cluster_label_prop(g),        g),
  run_algo("Infomap",            function(g) cluster_infomap(g),           g)
)

# Best algorithm = highest Q
best_idx  <- which.max(sapply(algos, `[[`, "Q"))
best_algo <- algos[[best_idx]]
cat(sprintf("\n  Best algorithm: %s (Q = %.4f)\n", best_algo$name, best_algo$Q))

# Keep Louvain result for main community plot (it IS the best)
comm_plot  <- algos[[1]]$comm          # Louvain
Q_plot     <- algos[[1]]$Q
n_comm     <- algos[[1]]$nc
comm_sizes <- sort(table(membership(comm_plot)), decreasing = TRUE)

# Full-network Louvain Q for reporting
set.seed(42); comm_full <- cluster_louvain(full_g)
Q_full <- round(modularity(comm_full), 4)
cat(sprintf("  Full network Q (t=5, Louvain): %.4f\n", Q_full))

# =============================================================================
# 5  NODE METRICS & LAYOUT
# =============================================================================
cat("\n[5/6] Metrics, layout, palettes...\n")

deg      <- degree(g)
deg_norm <- (deg - min(deg)) / diff(range(deg))
if (!is.finite(diff(range(deg))) || diff(range(deg)) == 0)
  deg_norm <- rep(0.5, length(deg))

node_sizes       <- deg_norm * 15 + 3
node_sz_c        <- deg_norm * 12 + 3
node_colors_heat <- fine_palette(100)[as.numeric(cut(deg, breaks = 100))]
comm_colors      <- fine_palette(100)[as.numeric(cut(membership(comm_plot), breaks = 100))]

cat(sprintf("  Computing FR layout (%d nodes, %d edges)...\n", vcount(g), ecount(g)))
set.seed(42)
lyt <- layout_with_fr(g, niter = 500)
if (!all(is.finite(lyt)) || diff(range(lyt[,1])) < 1e-6) {
  cat("  FR collapsed — using KK\n"); lyt <- layout_with_kk(g)
}
cat("  Layout done.\n")

n_edges <- ecount(g)
ew <- if (n_edges > 3000) 0.3 else 0.7

# =============================================================================
# 6  PLOTS
# =============================================================================
cat("\n[6/6] Generating plots...\n")

# ── PLOT 0: BIPARTITE STRUCTURE ───────────────────────────────────────────────
cat("  [0/4] Bipartite structure...\n")

bip_g         <- graph_from_data_frame(bip_sample[, .(disease, gene)], directed = FALSE)
V(bip_g)$type <- V(bip_g)$name %in% bip_sample$disease
bip_lyt       <- layout_as_bipartite(bip_g)[, c(2, 1)]
bip_vcol      <- ifelse(V(bip_g)$type, "#FF0000", "#000033")
bip_vsize     <- ifelse(V(bip_g)$type, 14, 10)
bip_shape     <- ifelse(V(bip_g)$type, "square", "circle")

png(out("0_Bipartite_Structure.png"), width = 10, height = 6, units = "in",
    res = 300, bg = "white")
par(bg = "white", mar = c(1,1,2,1))
plot(bip_g, layout = bip_lyt,
     vertex.color = bip_vcol, vertex.shape = bip_shape, vertex.size = bip_vsize,
     vertex.frame.color = "black", vertex.frame.width = 1.5,
     vertex.label = NA, edge.color = "gray50", edge.width = 1.5)
legend("topright", legend = c("Disease node","Gene node"),
       pch = c(15,21), col = c("#FF0000","#000033"),
       pt.bg = c("#FF0000","#000033"), pt.cex = 2,
       bty = "n", text.col = "black", cex = 1.1)
dev.off()
cat("  [0/4] Done.\n")

# ── PLOT 1: PLEIOTROPY HEAT MAP ───────────────────────────────────────────────
cat("  [1/4] Pleiotropy heat map...\n")

png(out("1_Diseasome_Pleiotropy.png"), width = 10, height = 10, units = "in",
    res = 300, bg = "white")
par(bg = "white", mar = c(0,0,0,0))

plot(g, layout = lyt,
     vertex.size        = node_sizes, 
     vertex.color       = node_colors_heat,
     vertex.frame.color = "black", # Black frame for white slides
     vertex.frame.width = 0.8,
     vertex.label       = NA, 
     edge.color         = 'black', 
     edge.width         = ew)

legend("bottomleft",
       legend     = c("Low connectivity", "", "", "", "High connectivity (hub)"),
       fill       = fine_palette(5), 
       border     = "black", 
       bty        = "n", 
       text.col   = "black", # Black text for white slides
       title      = sprintf("Disease Degree (t\u2265%d)", THRESHOLD_PLOT),
       title.font = 2,
       title.col  = "black", 
       inset      = c(0.05, 0.08),
       cex        = 1.3,     # Big text
       pt.cex     = 3.0,     # Big boxes
       y.intersp  = 1.3)

dev.off()
cat("  [1/4] Done.\n")


# ── PLOT 2: LOUVAIN COMMUNITY MAP ─────────────────────────────────────────────
cat("  [2/4] Louvain community map...\n")

top_n   <- min(6L, n_comm)
top_idx <- as.integer(names(head(comm_sizes, top_n)))

# ── CONSOLE ONLY: Print MESH codes so YOU can verify biology (never shown on slide)
cat("\n  === MESH CODE LOOKUP — verify these for your own knowledge ===\n")
for (i in seq_along(top_idx)) {
  mod_id       <- top_idx[i]
  nodes_in_mod <- V(g)[membership(comm_plot) == mod_id]
  mod_degrees  <- sort(degree(g, v = nodes_in_mod), decreasing = TRUE)
  top_diseases <- head(names(mod_degrees), 3)   # top 3 hub diseases per module
  cat(sprintf("  Module %d (n=%d) top hubs: %s\n",
              mod_id,
              comm_sizes[as.character(mod_id)],
              paste(top_diseases, collapse = " | ")))
}
cat("  =============================================================\n\n")

# ── SLIDE LABELS: clean "Module 1 (n=XX)" — no invented biological names ──────
top_labs <- character(top_n)
for (i in seq_along(top_idx)) {
  mod_id      <- as.character(top_idx[i])
  top_labs[i] <- sprintf("Module %s  (n=%d)", mod_id, comm_sizes[mod_id])
}

leg_cols <- fine_palette(max(top_n, 2))[round(seq(1, max(top_n, 2), length.out = top_n))]

png(out("2_Diseasome_Communities.png"), width = 10, height = 10, units = "in",
    res = 300, bg = "white")
par(bg = "white", mar = c(0, 0, 0, 0))

plot(g, layout = lyt,
     vertex.size        = node_sz_c,
     vertex.color       = comm_colors,
     vertex.frame.color = "black",
     vertex.frame.width = 0.8,
     vertex.label       = NA,
     edge.color         = "black",
     edge.width         = ew * 1.2)

# ── LEGEND: top-right, large text, no overlap ─────────────────────────────────
# Moved from bottomleft (was colliding with nodes) to topright (clear space)
legend("topright",
       legend     = top_labs,
       fill       = leg_cols,
       border     = "black",
       bty        = "n",
       text.col   = "black",
       title      = "Louvain Communities",
       title.font = 2,
       title.col  = "black",
       inset      = c(0.02, 0.04),   # small inset from edge
       cex        = 1.4,             # large text — same as before
       pt.cex     = 3.0,             # large colour boxes — same as before
       y.intersp  = 1.8)             # generous spacing — no overlap

dev.off()
cat("  [2/4] Done.\n")

# ── PLOT 3: CCDF DEGREE DISTRIBUTION ─────────────────────────────────────────
cat("  [3/4] CCDF degree distribution...\n")

dd     <- sort(disease_deg_bip$n_genes, decreasing = TRUE)
k_vals <- sort(unique(dd))
ccdf   <- sapply(k_vals, function(k) mean(dd >= k))
med_d  <- median(disease_deg_bip$n_genes)
max_d  <- max(disease_deg_bip$n_genes)

png(out("3_Degree_Distribution.png"), width = 9, height = 6.5, units = "in",
    res = 300, bg = "white")
set_theme()
plot(k_vals, ccdf * 100,
     log = "xy", type = "s", lwd = 3, col = "#006400",
     xlab = "Disease degree  k  (number of associated genes)",
     ylab = "P(K \u2265 k)  [%]",
     main = "Cumulative Degree Distribution",
     axes = FALSE,
     ylim = c(max(min(ccdf * 100, na.rm = TRUE), 1e-4), 100))
axis(1, col = "black", col.ticks = "black")
axis(2, col = "black", col.ticks = "black",
     at = c(100,10,1,0.1,0.01),
     labels = c("100%","10%","1%","0.1%","0.01%"), las = 1)
box(col = "black")
grid(col = "gray70", lty = "dotted", lwd = 0.8)
legend("topright",
       legend = c("P(K \u2265 k)",
                  sprintf("Median: %d genes/disease", as.integer(med_d)),
                  sprintf("Max hub: %s genes", format(max_d, big.mark = ","))),
       lty = c(1,NA,NA), lwd = c(3,NA,NA),
       col = c("#006400","black","black"),
       bty = "n", text.col = "black", cex = 0.9)
dev.off()
cat("  [3/4] Done.\n")

# ── PLOT 4: ALGORITHM COMPARISON — Q SCORES & COMMUNITIES ────────────────────
# Two-panel chart:
#   Left  — Modularity Q for each algorithm (higher = better)
#   Right — Number of communities found
cat("  [4/4] Algorithm comparison chart...\n")

algo_names <- sapply(algos, `[[`, "name")
algo_Q     <- sapply(algos, `[[`, "Q")
algo_nc    <- sapply(algos, `[[`, "nc")
algo_time  <- sapply(algos, `[[`, "time")

# Colours: highlight Louvain in dark red; others in steel blue
bar_cols <- ifelse(algo_names == "Louvain", "#CC0000", "#4A7FA5")

png(out("4_Algorithm_Comparison.png"), width = 10, height = 6, units = "in",
    res = 300, bg = "white")
set_theme()

layout(matrix(1:2, nrow = 1), widths = c(1,1))

# Panel 1 — Q scores
par(mar = c(6, 5, 3, 1), bg = "white",
    col.axis = "black", col.lab = "black", col.main = "black", fg = "black")
bp1 <- barplot(algo_Q,
               names.arg = algo_names,
               col       = bar_cols,
               border    = "black",
               ylim      = c(0, max(algo_Q) * 1.25),
               ylab      = "Modularity  Q  (higher = better structure)",
               main      = "Modularity Score by Algorithm",
               las       = 2,
               cex.names = 0.85,
               cex.axis  = 0.85,
               font.lab  = 2)
# Value labels on bars
text(bp1, algo_Q + max(algo_Q) * 0.03,
     labels = sprintf("%.4f", algo_Q), cex = 0.82, font = 2, col = "black")
# Star on winner
text(bp1[which.max(algo_Q)], max(algo_Q) * 1.18,
     labels = "\u2605 Best", cex = 1.0, font = 2, col = "#CC0000")
grid(nx = NA, ny = NULL, col = "gray80", lty = "dotted")
box(col = "black")

# Panel 2 — Number of communities
par(mar = c(6, 4, 3, 2), bg = "white",
    col.axis = "black", col.lab = "black", col.main = "black", fg = "black")
bp2 <- barplot(algo_nc,
               names.arg = algo_names,
               col       = bar_cols,
               border    = "black",
               ylim      = c(0, max(algo_nc) * 1.25),
               ylab      = "Communities detected",
               main      = "Communities Found by Algorithm",
               las       = 2,
               cex.names = 0.85,
               cex.axis  = 0.85,
               font.lab  = 2)
text(bp2, algo_nc + max(algo_nc) * 0.03,
     labels = as.character(algo_nc), cex = 0.85, font = 2, col = "black")
grid(nx = NA, ny = NULL, col = "gray80", lty = "dotted")
box(col = "black")

dev.off()
cat("  [4/4] Done.\n")

# =============================================================================
# METRICS SUMMARY
# =============================================================================
cat("\n╔══════════════════════════════════════════════════════╗\n")
cat("║               ANALYSIS SUMMARY METRICS                ║\n")
cat("╠══════════════════════════════════════════════════════╣\n")
cat(sprintf("║  Raw associations           : %10s           ║\n", format(n_raw_assoc,   big.mark=",")))
cat(sprintf("║  Raw diseases               : %10s           ║\n", format(n_raw_disease, big.mark=",")))
cat(sprintf("║  Raw genes                  : %10s           ║\n", format(n_raw_gene,    big.mark=",")))
cat("╠══════════════════════════════════════════════════════╣\n")
cat(sprintf("║  Full network nodes (t=5)   : %10s           ║\n", format(vcount(full_g), big.mark=",")))
cat(sprintf("║  Full network edges (t=5)   : %10s           ║\n", format(ecount(full_g), big.mark=",")))
cat(sprintf("║  Full avg degree            : %10.2f           ║\n", mean(deg_full)))
cat(sprintf("║  Full max degree            : %10d           ║\n", max(deg_full)))
cat(sprintf("║  Full diameter              : %10d           ║\n", diameter(full_g)))
cat(sprintf("║  Full clustering coeff.     : %10.4f           ║\n", transitivity(full_g,"average")))
cat("╠══════════════════════════════════════════════════════╣\n")
cat("║  ALGORITHM COMPARISON (plot subgraph, N=150)         ║\n")
cat("╠══════════════════════════════════════════════════════╣\n")
for (a in algos) {
  cat(sprintf("║  %-20s  Q=%.4f  nc=%2d  t=%.3fs       ║\n",
              a$name, a$Q, a$nc, a$time))
}
cat(sprintf("║  Full network Q (Louvain)   : %10.4f           ║\n", Q_full))
cat("╠══════════════════════════════════════════════════════╣\n")
cat(sprintf("║  Median genes per disease   : %10.0f           ║\n", med_d))
cat(sprintf("║  Max genes (hub disease)    : %10s           ║\n", format(max_d, big.mark=",")))
cat("╚══════════════════════════════════════════════════════╝\n")
cat("\nAll 5 plots saved to: diseasome_output/\n")
cat("  0_Bipartite_Structure.png   — sample disease-gene bipartite structure\n")
cat("  1_Diseasome_Pleiotropy.png  — degree-weighted pleiotropy heat map\n")
cat("  2_Diseasome_Communities.png — winning community partition (Louvain)\n")
cat("  3_Degree_Distribution.png   — cumulative degree distribution\n")
cat("  4_Algorithm_Comparison.png  — modularity comparison across algorithms\n")
