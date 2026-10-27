# Author: Yang Liu
# Date: 2026-10-27 | R 4.5.2
# Packages: igraph, clusterProfiler, AnnotationDbi, org.Hs.eg.db, GOSemSim.
# Input: candidate assays and STRING v12.0 human info, aliases and links tables.

# 1. Primary and sensitivity networks
# p <- ppi_run(candidates, protein_info, aliases, links)
ppi_run <- function(candidates, protein_info, aliases, links) {
  mapping <- ppi_map_assays(candidates, protein_info, aliases)
  primary <- ppi_build_network(mapping, links, 400L)
  sensitivity <- ppi_build_network(mapping, links, 700L)
  list(mapping = mapping, primary = primary, sensitivity = sensitivity,
       threshold_statistics = rbind(primary$statistics, sensitivity$statistics), hubs = ppi_hubs(primary))
}

# 2. Network construction
# primary <- ppi_build_network(p$mapping, links, threshold = 400)
# Candidate nodes, weighted fast-greedy modules and full-network centralities.
ppi_build_network <- function(mapping, links, threshold) {
  links <- ppi_columns(links)
  mapped <- mapping[!is.na(mapping$STRING_ID), , drop = FALSE]
  if (!nrow(mapped)) stop("No STRING-mapped candidate proteins.")
  if (anyDuplicated(mapped$STRING_ID) || anyDuplicated(mapped$Protein))
    stop("Multiple assays map to one network node; resolve explicitly before graph construction.")
  links <- links[grepl("^9606\\.", links$protein1) & grepl("^9606\\.", links$protein2), , drop = FALSE]
  links$combined_score <- suppressWarnings(as.numeric(links$combined_score))
  if (any(!is.finite(links$combined_score)) || any(links$combined_score < 0 | links$combined_score > 1000))
    stop("STRING scores must be finite values in [0,1000].")
  x <- links[links$protein1 %in% mapped$STRING_ID & links$protein2 %in% mapped$STRING_ID &
    links$protein1 != links$protein2 & links$combined_score >= threshold, , drop = FALSE]
  if (nrow(x)) {
    edges <- stats::aggregate(x$combined_score,
      by = list(from = pmin(x$protein1, x$protein2), to = pmax(x$protein1, x$protein2)), FUN = max)
    names(edges)[3] <- "combined_score"
  } else edges <- data.frame(from = character(), to = character(), combined_score = numeric())
  vertices <- data.frame(name = mapped$STRING_ID, Protein = mapped$Protein, Assay = mapped$Assay)
  g <- igraph::graph_from_data_frame(edges, directed = FALSE, vertices = vertices)
  degree <- igraph::degree(g)
  connected <- igraph::induced_subgraph(g, which(degree > 0))
  memberships <- setNames(rep(0L, igraph::vcount(g)), igraph::V(g)$name)
  modularity <- NA_real_; clustering <- NULL
  if (igraph::vcount(connected) > 1L && igraph::ecount(connected) > 0L) {
    clustering <- igraph::cluster_fast_greedy(connected, weights = igraph::E(connected)$combined_score/1000)
    z <- igraph::membership(clustering); memberships[names(z)] <- as.integer(z)
    modularity <- igraph::modularity(clustering)
  }
  eigenvector <- if (igraph::vcount(g) > 1L && igraph::ecount(g) > 0L)
    as.numeric(igraph::eigen_centrality(g, weights = igraph::E(g)$combined_score/1000)$vector) else
      rep(0, igraph::vcount(g))
  nodes <- data.frame(Protein = igraph::V(g)$Protein, Assay = igraph::V(g)$Assay,
    STRING_ID = igraph::V(g)$name, Cluster = unname(memberships[igraph::V(g)$name]),
    Degree = unname(degree), Betweenness = as.numeric(igraph::betweenness(g, weights = NA, normalized = TRUE)),
    Eigenvector_centrality = eigenvector, stringsAsFactors = FALSE)
  nodes <- nodes[order(-nodes$Degree, -nodes$Betweenness, nodes$Protein), , drop = FALSE]
  nodes$Hub_rank <- seq_len(nrow(nodes))
  edge_table <- igraph::as_data_frame(g, what = "edges")
  edge_table$Protein_1 <- nodes$Protein[match(edge_table$from, nodes$STRING_ID)]
  edge_table$Protein_2 <- nodes$Protein[match(edge_table$to, nodes$STRING_ID)]
  edge_table$normalized_score <- edge_table$combined_score/1000
  stats <- data.frame(threshold = threshold, input_assays = nrow(mapping), mapped_assays = nrow(mapped),
    nodes = nrow(nodes), connected_nodes = sum(nodes$Degree > 0), isolated_nodes = sum(nodes$Degree == 0),
    edges = nrow(edges), mean_degree = mean(nodes$Degree), density = igraph::edge_density(g),
    clusters = length(unique(nodes$Cluster[nodes$Cluster > 0])), modularity = modularity)
  list(graph = g, connected_graph = connected, nodes = nodes, edges = edge_table,
       clustering = clustering, statistics = stats)
}

# 3. Core-subgraph MCC hubs
# hubs <- ppi_hubs(p$primary)
# Core: degree, betweenness and closeness all above their full-network medians.
ppi_hubs <- function(network) {
  g <- network$graph
  degree <- igraph::degree(g)
  between <- igraph::betweenness(g, directed = FALSE, weights = NA, normalized = TRUE)
  close <- igraph::closeness(g, weights = NA, normalized = TRUE)
  close[!is.finite(close)] <- 0
  core <- degree > stats::median(degree) & between > stats::median(between) & close > stats::median(close)
  core_graph <- igraph::induced_subgraph(g, names(degree)[core])
  mcc <- setNames(rep(0, igraph::vcount(core_graph)), igraph::V(core_graph)$name)
  # MCC sums (clique size - 1)! over maximal cliques in the core graph.
  for (clique in igraph::max_cliques(core_graph, min = 2)) {
    ids <- igraph::as_ids(clique)
    mcc[ids] <- mcc[ids] + factorial(length(ids) - 1)
  }
  if (any(!is.finite(mcc))) stop("MCC score overflow.")
  table <- data.frame(STRING_ID = names(degree), Protein = igraph::V(g)$Protein,
    Degree = unname(degree), Betweenness = unname(between), Closeness = unname(close),
    Degree_median = stats::median(degree), Betweenness_median = stats::median(between),
    Closeness_median = stats::median(close), Core = unname(core), MCC = unname(mcc[names(degree)]),
    MCC_rank = NA_integer_, Hub_display = FALSE)
  if (length(mcc)) {
    table$MCC_rank[table$Core] <- rank(-table$MCC[table$Core], ties.method = "min")
    # Retain top-10 positive ranks, including cutoff ties.
    table$Hub_display <- table$Core & !is.na(table$MCC_rank) & table$MCC_rank <= 10 &
      !is.na(table$MCC) & table$MCC > 0
  }
  list(table = table, core_graph = core_graph,
       hub_graph = igraph::induced_subgraph(g, table$STRING_ID[table$Hub_display]))
}

# 4. Module GO biological processes
# go <- ppi_module_go(p$primary)
# Formal and broad annotation families are separate; top labels may be nonsignificant.
ppi_module_go <- function(primary, symbol_to_entrez = NULL) {
  nodes <- primary$nodes; db <- org.Hs.eg.db::org.Hs.eg.db
  if (is.null(symbol_to_entrez)) symbol_to_entrez <- ppi_gene_mapping(nodes$Protein)
  if (!all(c("SYMBOL", "ENTREZID") %in% names(symbol_to_entrez)) || anyDuplicated(symbol_to_entrez$SYMBOL))
    stop("Provide one pre-resolved ENTREZID per SYMBOL for PPI GO mapping.")
  universe <- sort(unique(as.character(AnnotationDbi::keys(db, keytype = "ENTREZID"))))
  formal <- annotation <- list()
  for (cluster in sort(unique(nodes$Cluster[nodes$Cluster > 0]))) {
    symbols <- nodes$Protein[nodes$Cluster == cluster]
    genes <- unique(stats::na.omit(as.character(symbol_to_entrez$ENTREZID[match(symbols, symbol_to_entrez$SYMBOL)])))
    # Formal GO: >=3 genes, full human Entrez universe, default term-size bounds.
    if (length(genes) >= 3L) {
      fit <- clusterProfiler::enrichGO(gene = genes, universe = universe, OrgDb = db,
        keyType = "ENTREZID", ont = "BP", pvalueCutoff = 1, qvalueCutoff = 1,
        pAdjustMethod = "BH", readable = TRUE)
      d <- as.data.frame(fit)
      if (nrow(d)) {
        d$Cluster <- cluster; d$Cluster_size <- length(symbols)
        d$Fold_enrichment <- ppi_fold_enrichment(d)
        formal[[length(formal) + 1L]] <- d
      }
    }
    # Broad annotation: all terms of size 1..5000, its own BH family.
    if (length(genes)) {
      fit <- clusterProfiler::enrichGO(gene = genes, OrgDb = db, keyType = "ENTREZID", ont = "BP",
        pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1,
        minGSSize = 1, maxGSSize = 5000, readable = TRUE)
      d <- as.data.frame(fit)
      if (nrow(d)) {
        d <- d[order(d$pvalue, d$p.adjust, -d$Count, d$Description), , drop = FALSE]
        d$Cluster <- cluster; d$Cluster_size <- length(symbols)
        d$Cluster_mapped_gene_count <- length(genes); d$Rank_within_cluster <- seq_len(nrow(d))
        d$Significant_BH_0_05 <- d$p.adjust < 0.05
        d$Significant_formal_threshold <- d$pvalue < 0.05 & d$p.adjust < 0.05 & d$qvalue < 0.20
        d$Fold_enrichment <- ppi_fold_enrichment(d)
        annotation[[length(annotation) + 1L]] <- d
      }
    }
  }
  pool <- function(rows) {
    if (!length(rows)) return(data.frame())
    d <- do.call(rbind, rows)
    d$Global_BH_FDR <- stats::p.adjust(d$pvalue, "BH", n = nrow(d))
    d
  }
  formal <- pool(formal); annotation <- pool(annotation)
  top <- if (nrow(annotation)) annotation[!duplicated(annotation$Cluster), , drop = FALSE] else annotation
  list(formal = formal, annotation = annotation, top_annotation = top,
       mapping = symbol_to_entrez, universe = universe)
}

# 5. Within-module GO similarity
# annotation <- ppi_go_similarity(go$annotation)
# Wang>=0.25 annotates similar terms without removing them.
ppi_go_similarity <- function(annotation) {
  if (!nrow(annotation)) return(annotation)
  sem <- GOSemSim::godata(annoDb = "org.Hs.eg.db", ont = "BP", computeIC = FALSE)
  annotation$Redundant_with_term <- NA_character_
  for (cluster in unique(annotation$Cluster)) {
    ii <- which(annotation$Cluster == cluster)
    if (length(ii) < 2L) next
    sim <- GOSemSim::mgoSim(annotation$ID[ii], annotation$ID[ii], semData = sem, measure = "Wang", combine = NULL)
    for (j in seq_along(ii)) {
      linked <- which(is.finite(sim[j, ]) & sim[j, ] >= 0.25 & seq_along(ii) != j)
      if (length(linked)) annotation$Redundant_with_term[ii[j]] <- paste(annotation$Description[ii[linked]], collapse = "|")
    }
  }
  annotation
}

# 6. Network plot
# ppi_plot(p$primary)
ppi_plot <- function(network, seed = 2026L) {
  g <- network$graph
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit(if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv) else
    if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv))
  set.seed(seed)
  pos <- igraph::layout_with_fr(g, weights = NA)
  plot(g, layout = pos, vertex.label = igraph::V(g)$Protein,
       vertex.size = 4 + sqrt(igraph::degree(g)), edge.width = 0.5 + 2*igraph::E(g)$combined_score/1000)
  invisible(pos)
}

# Helpers
ppi_columns <- function(data) {
  names(data) <- sub("^X\\.", "", names(data))
  names(data) <- sub("^#", "", names(data))
  data
}

ppi_map_assays <- function(assays, protein_info, aliases) {
  protein_info <- ppi_columns(protein_info); aliases <- ppi_columns(aliases)
  assays <- unique(toupper(as.character(assays)))
  if (!length(assays) || anyNA(assays) || any(!nzchar(assays))) stop("Supply nonempty candidate assays.")
  special <- c(NTPROBNP = "NPPB", DEFA1_DEFA1B = "DEFA1",
               DEFB4A_DEFB4B = "DEFB4A", LGALS7_LGALS7B = "LGALS7")
  input <- data.frame(Assay = assays, Protein = ifelse(assays %in% names(special),
    unname(special[assays]), assays), stringsAsFactors = FALSE)
  info <- protein_info[grepl("^9606\\.", protein_info$string_protein_id), , drop = FALSE]
  info <- info[!duplicated(info$string_protein_id), , drop = FALSE]
  info$key <- toupper(trimws(info$preferred_name))
  aliases <- aliases[grepl("^9606\\.", aliases$string_protein_id), , drop = FALSE]
  aliases$key <- toupper(trimws(aliases$alias))
  direct <- info[info$key %in% input$Protein, , drop = FALSE]
  direct <- data.frame(Protein = direct$key, STRING_ID = direct$string_protein_id,
    preferred_name = direct$preferred_name, mapping_priority = rep(1L, nrow(direct)),
    mapping_method = rep("STRING preferred name", nrow(direct)))
  a <- aliases[aliases$key %in% input$Protein, , drop = FALSE]
  preferred <- info$preferred_name[match(a$string_protein_id, info$string_protein_id)]
  hgnc <- !is.na(a$source) & grepl("HGNC", a$source, ignore.case = TRUE)
  same_name <- !is.na(preferred) & toupper(trimws(preferred)) == a$key
  fallback <- data.frame(Protein = a$key, STRING_ID = a$string_protein_id,
    preferred_name = preferred, mapping_priority = ifelse(hgnc, 2L, ifelse(same_name, 3L, 4L)),
    mapping_method = ifelse(hgnc, "HGNC-supported STRING alias", paste0("STRING alias: ", a$source)))
  hits <- rbind(direct, fallback)
  hits <- hits[order(hits$Protein, hits$mapping_priority, hits$STRING_ID), , drop = FALSE]
  hits <- hits[!duplicated(hits$Protein), , drop = FALSE]
  mapping <- merge(input, hits, by = "Protein", all.x = TRUE, sort = FALSE)
  mapping$Mapping_status <- ifelse(is.na(mapping$STRING_ID), "Unmapped", "Mapped")
  mapping
}

ppi_gene_mapping <- function(symbols) {
  db <- org.Hs.eg.db::org.Hs.eg.db
  official <- AnnotationDbi::keys(db, keytype = "SYMBOL")
  canonical <- official[match(toupper(symbols), toupper(official))]
  valid <- unique(stats::na.omit(canonical))
  if (!length(valid)) stop("No candidate symbols map to the human annotation database.")
  lookup <- suppressMessages(AnnotationDbi::select(db, keys = valid, keytype = "SYMBOL", columns = "ENTREZID"))
  lookup <- lookup[!is.na(lookup$ENTREZID), , drop = FALSE]
  lookup <- lookup[!duplicated(lookup$SYMBOL), , drop = FALSE]
  data.frame(SYMBOL = symbols, ENTREZID = lookup$ENTREZID[match(canonical, lookup$SYMBOL)])
}

ppi_fold_enrichment <- function(data) {
  ratio <- function(x) vapply(strsplit(as.character(x), "/", fixed = TRUE),
    function(y) as.numeric(y[1])/as.numeric(y[2]), numeric(1))
  ratio(data$GeneRatio)/ratio(data$BgRatio)
}
