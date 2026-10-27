# Author: Yang Liu
# Date: 2026-10-27 | R 4.5.2
# Packages: clusterProfiler, ReactomePA, org.Hs.eg.db, GOSemSim, igraph, ggplot2.
# Input: full Olink screening results and optional SYMBOL/ENTREZID mapping.

# 1. GO, KEGG and Reactome enrichment
# e <- enrichment_run(screen)
# Database-wide annotated backgrounds; raw P<0.05, with no multiple-testing correction.
enrichment_run <- function(screen, symbol_to_entrez = NULL) {
  candidates <- enrichment_candidates(screen)
  mapping <- enrichment_map_assays(candidates, symbol_to_entrez)
  if (!length(mapping$entrez)) stop("No candidate Entrez IDs after mapping.")
  # Omitting universe retains each database-wide background.
  go <- lapply(c("BP", "CC", "MF"), function(ontology) enrichment_tidy(
    clusterProfiler::enrichGO(gene = mapping$entrez, OrgDb = org.Hs.eg.db::org.Hs.eg.db,
      keyType = "ENTREZID", ont = ontology, pAdjustMethod = "none",
      pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE), "GO", ontology))
  kegg <- enrichment_tidy(clusterProfiler::enrichKEGG(gene = mapping$entrez,
    organism = "hsa", keyType = "kegg", pAdjustMethod = "none",
    pvalueCutoff = 1, qvalueCutoff = 1), "KEGG")
  reactome <- enrichment_tidy(ReactomePA::enrichPathway(gene = mapping$entrez,
    organism = "human", pAdjustMethod = "none", pvalueCutoff = 1,
    qvalueCutoff = 1, readable = TRUE), "Reactome")
  tables <- c(go, list(kegg, reactome)); tables <- tables[vapply(tables, nrow, integer(1)) > 0L]
  full <- if (length(tables)) do.call(rbind, tables) else data.frame()
  list(full = full, candidate_mapping = mapping, candidates = candidates)
}

# 2. Hierarchical pathway reduction
# modules <- enrichment_reduce(e)
# GO Wang>=0.50; KEGG/Reactome Jaccard>=0.30/0.50; merged Jaccard>=0.30.
enrichment_reduce <- function(result) {
  d <- result$full
  if (!nrow(d)) return(list(database = list(), merged = NULL, display_nodes = d, display_edges = data.frame()))
  d <- d[!is.na(d$pvalue) & d$pvalue < 0.05, , drop = FALSE]
  d$fold_enrichment <- enrichment_ratio(d$GeneRatio)/enrichment_ratio(d$BgRatio)
  mapping <- result$candidate_mapping$mapping
  lookup <- split(mapping$ENTREZID, mapping$SYMBOL)
  d$genes <- lapply(strsplit(d$geneID, "/", fixed = TRUE), function(ids) {
    numeric_ids <- ids[grepl("^[0-9]+$", ids)]
    symbols <- ids[!grepl("^[0-9]+$", ids)]
    sort(unique(stats::na.omit(c(numeric_ids, unlist(lookup[symbols], use.names = FALSE)))))
  })
  d <- d[is.finite(d$fold_enrichment) & lengths(d$genes) > 0L, , drop = FALSE]
  d$key <- paste(d$Database, d$Ontology, d$ID, sep = "::")
  d <- d[order(d$pvalue, -d$Count, d$Description), , drop = FALSE]
  go <- d[d$Database == "GO", , drop = FALSE]
  go_edges <- lapply(c("BP", "CC", "MF"), function(ontology) {
    x <- go[go$Ontology == ontology, , drop = FALSE]
    if (nrow(x) < 2L) return(NULL)
    sem <- GOSemSim::godata(annoDb = "org.Hs.eg.db", ont = ontology, computeIC = FALSE)
    sim <- GOSemSim::mgoSim(x$ID, x$ID, semData = sem, measure = "Wang", combine = NULL)
    ij <- which(upper.tri(sim) & is.finite(sim) & sim >= 0.50, arr.ind = TRUE)
    data.frame(from = x$key[ij[, 1]], to = x$key[ij[, 2]], similarity = sim[ij])
  })
  go_edges <- Filter(Negate(is.null), go_edges)
  go_edges <- if (length(go_edges)) do.call(rbind, go_edges) else
    data.frame(from = character(), to = character(), similarity = numeric())
  kegg <- d[d$Database == "KEGG", , drop = FALSE]
  reactome <- d[d$Database == "Reactome", , drop = FALSE]
  database <- list(
    GO = enrichment_reduce_network(go, go_edges, "edge_betweenness", "GO"),
    KEGG = enrichment_reduce_network(kegg, enrichment_jaccard_edges(kegg, 0.30), "fast_greedy", "KEGG"),
    Reactome = enrichment_reduce_network(reactome, enrichment_jaccard_edges(reactome, 0.50), "fast_greedy", "Reactome"))
  reps <- do.call(rbind, lapply(database, `[[`, "representatives"))
  merged <- enrichment_reduce_network(reps, enrichment_jaccard_edges(reps, 0.30), "fast_greedy", "FM")
  display_nodes <- utils::head(merged$representatives, 10L)
  display_edges <- enrichment_jaccard_edges(display_nodes, 0.10)
  display_nodes <- display_nodes[display_nodes$key %in% unique(c(display_edges$from, display_edges$to)), , drop = FALSE]
  display <- enrichment_reduce_network(display_nodes, display_edges, "edge_betweenness", "Display")
  list(database = database, merged = merged, display_nodes = display$annotated, display_edges = display_edges)
}

# 3. Nominal enrichment plot
# enrichment_plot(e$full)
enrichment_plot <- function(full, n = 15L) {
  if (!nrow(full)) stop("No enrichment rows to plot.")
  d <- full[!is.na(full$pvalue) & full$pvalue < 0.05, , drop = FALSE]
  d <- utils::head(d[order(d$pvalue, -d$Count, d$Description), , drop = FALSE], n)
  d$term <- factor(paste(d$Database, d$Ontology, d$Description, sep = ": "),
                   levels = rev(unique(paste(d$Database, d$Ontology, d$Description, sep = ": "))))
  d$neg_log10_p <- -log10(pmax(d$pvalue, .Machine$double.xmin))
  ggplot2::ggplot(d, ggplot2::aes(neg_log10_p, term, size = Count, colour = Database)) +
    ggplot2::geom_point() + ggplot2::labs(x = "-log10(nominal P)", y = NULL) + ggplot2::theme_classic()
}

# Helpers
enrichment_candidates <- function(screen) {
  a <- screen$path_a; b <- screen$path_b
  stopifnot(!anyDuplicated(a$feature), !anyDuplicated(b$feature))
  candidates <- sort(intersect(
    a$feature[!is.na(a$p_bonferroni) & a$p_bonferroni < 0.05],
    b$feature[!is.na(b$cox_p_bonferroni) & b$cox_p_bonferroni < 0.05]))
  if (!setequal(candidates, screen$candidates)) stop("Candidate list differs from the double-Bonferroni intersection.")
  if (!length(candidates)) stop("No intersection candidates.")
  candidates
}

enrichment_map_assays <- function(assays, symbol_to_entrez = NULL) {
  # Enrichment expands composite assays; STRING uses its own mapping.
  special <- list(NTproBNP = "NPPB", CERT = "CERT1", WARS = "WARS1", GBA = "GBA1",
    DEFA1_DEFA1B = c("DEFA1", "DEFA1B"), MICB_MICA = c("MICB", "MICA"),
    SPACA5_SPACA5B = c("SPACA5", "SPACA5B"), ERVV_1 = "ERVV-1",
    HLA_A = "HLA-A", HLA_DRA = "HLA-DRA", HLA_E = "HLA-E", MYLPF = "MYL11")
  assays <- unique(as.character(assays))
  if (!length(assays) || anyNA(assays) || any(!nzchar(assays))) stop("Assay labels must be nonempty.")
  expanded <- do.call(rbind, lapply(assays, function(a) {
    symbols <- if (a %in% names(special)) special[[a]] else a
    data.frame(assay = a, SYMBOL = symbols, stringsAsFactors = FALSE)
  }))
  if (is.null(symbol_to_entrez)) {
    symbol_to_entrez <- suppressMessages(clusterProfiler::bitr(unique(expanded$SYMBOL),
      fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db::org.Hs.eg.db))
  }
  if (!all(c("SYMBOL", "ENTREZID") %in% names(symbol_to_entrez))) stop("Mapping requires SYMBOL and ENTREZID.")
  mapping <- unique(symbol_to_entrez[c("SYMBOL", "ENTREZID")])
  mapping$ENTREZID <- as.character(mapping$ENTREZID)
  mapping <- mapping[!is.na(mapping$ENTREZID) & nzchar(mapping$ENTREZID), , drop = FALSE]
  full <- merge(expanded, mapping, by = "SYMBOL", all.x = TRUE, sort = FALSE)
  list(mapping = full, entrez = sort(unique(stats::na.omit(full$ENTREZID))),
       unmapped = full[is.na(full$ENTREZID), c("assay", "SYMBOL"), drop = FALSE])
}

enrichment_tidy <- function(object, database, ontology = NA_character_) {
  d <- as.data.frame(object)
  if (!nrow(d)) return(data.frame())
  d <- d[setdiff(names(d), c("p.adjust", "qvalue"))]
  d$Database <- database; d$Ontology <- ontology
  d$Significant_nominal <- !is.na(d$pvalue) & d$pvalue < 0.05
  d
}

enrichment_ratio <- function(x) vapply(strsplit(as.character(x), "/", fixed = TRUE), function(y) {
  if (length(y) != 2L) return(NA_real_)
  z <- suppressWarnings(as.numeric(y))
  if (all(is.finite(z)) && z[2] > 0) z[1]/z[2] else NA_real_
}, numeric(1))

enrichment_jaccard_edges <- function(nodes, cutoff) {
  empty <- data.frame(from = character(), to = character(), similarity = numeric())
  if (nrow(nodes) < 2L) return(empty)
  pairs <- utils::combn(seq_len(nrow(nodes)), 2)
  rows <- lapply(seq_len(ncol(pairs)), function(i) {
    a <- pairs[1, i]; b <- pairs[2, i]
    u <- union(nodes$genes[[a]], nodes$genes[[b]])
    score <- if (length(u)) length(intersect(nodes$genes[[a]], nodes$genes[[b]]))/length(u) else 0
    if (is.finite(score) && score >= cutoff)
      data.frame(from = nodes$key[a], to = nodes$key[b], similarity = score) else NULL
  })
  rows <- Filter(Negate(is.null), rows)
  if (length(rows)) do.call(rbind, rows) else empty
}

enrichment_reduce_network <- function(nodes, edges, algorithm, prefix) {
  if (!"constituent_ids" %in% names(nodes)) nodes$constituent_ids <- nodes$ID
  nodes$module <- if (nrow(nodes)) paste0(prefix, "-isolated-", seq_len(nrow(nodes))) else character()
  nodes$representative <- rep(FALSE, nrow(nodes))
  if (!nrow(nodes)) return(list(annotated = nodes, representatives = nodes, edges = edges))
  if (anyDuplicated(nodes$key)) stop("Duplicate pathway keys.")
  connected <- unique(c(edges$from, edges$to))
  if (length(connected)) {
    g <- igraph::graph_from_data_frame(edges, directed = FALSE,
      vertices = data.frame(name = connected))
    membership <- if (igraph::vcount(g) <= 2L) setNames(rep(1L, length(connected)), connected) else
      if (algorithm == "edge_betweenness") igraph::membership(igraph::cluster_edge_betweenness(g, directed = FALSE)) else
        igraph::membership(igraph::cluster_fast_greedy(g))
    nodes$module[match(names(membership), nodes$key)] <- paste0(prefix, "-M", unname(membership))
  }
  representative <- vapply(split(seq_len(nrow(nodes)), nodes$module), function(ii) {
    ii[order(nodes$pvalue[ii], -nodes$Count[ii], nodes$Description[ii])][1]
  }, integer(1))
  nodes$representative <- seq_len(nrow(nodes)) %in% representative
  reps <- nodes[representative, , drop = FALSE]
  for (i in seq_len(nrow(reps))) {
    ids <- nodes$constituent_ids[nodes$module == reps$module[i]]
    reps$constituent_ids[i] <- paste(unique(unlist(strsplit(ids, ";", fixed = TRUE))), collapse = ";")
  }
  reps <- reps[order(reps$pvalue, -reps$Count, reps$Description), , drop = FALSE]
  list(annotated = nodes, representatives = reps, edges = edges)
}
