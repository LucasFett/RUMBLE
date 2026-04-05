#'
#' Internal functions for TreeSHAP extraction, S-R-I decomposition and
#' ecological visualization assets.
#'
#' These helpers extend the classical RUMBLE interpretability layer with a
#' tree-specific stage based on exact TreeSHAP interaction values. The goal is
#' to decompose the global contribution of each taxon into synergy,
#' redundancy and independence components and to prepare ecological summaries
#' and topology graphs.
#'
#' @name tree-sri-internal
#' @keywords internal
#' @return No return value. This object is internal and serves only for organization.
NULL


#' @noRd
.projectVector <- function(x, onto, eps = 1e-12) {
  denom <- sum(onto^2, na.rm = TRUE)
  if (is.na(denom) || denom < eps) {
    return(rep(0, length(x)))
  }
  coef <- sum(x * onto, na.rm = TRUE) / denom
  coef * onto
}


#' @noRd
.orthogonalizeInteraction <- function(pij, pii, pjj, eps = 1e-12) {
  corrected <- pij
  corrected <- corrected - .projectVector(corrected, pii, eps = eps)
  corrected <- corrected - .projectVector(corrected, pjj, eps = eps)
  corrected
}


#' @noRd
.preprocessForTreeShap <- function(fitted_workflow, new_data, target_var) {
  mold <- workflows::extract_mold(fitted_workflow)
  predictors_only <- new_data[, !colnames(new_data) %in% target_var, drop = FALSE]

  processed <- hardhat::forge(
    predictors_only,
    blueprint = mold$blueprint
  )$predictors

  list(
    reference = as.data.frame(mold$predictors),
    processed = as.data.frame(processed)
  )
}


#' @noRd
.rangerUnifyRobust <- function(rf_model, data, positive_class = NULL) {
  if (!inherits(rf_model, "ranger")) {
    stop("Object rf_model was not of class 'ranger'.")
  }

  ranger_common <- getFromNamespace("ranger_unify.common", "treeshap")
  n_trees <- rf_model$num.trees

  tree_list <- lapply(seq_len(n_trees), function(tree_idx) {
    tree_data <- as.data.frame(
      ranger::treeInfo(rf_model, tree = tree_idx),
      stringsAsFactors = FALSE
    )

    pred_cols <- grep("^pred\\.", colnames(tree_data), value = TRUE)

    if (!("prediction" %in% colnames(tree_data))) {
      if (length(pred_cols) == 0L) {
        stop(
          "Could not identify a compatible prediction column in ranger::treeInfo()."
        )
      }

      preferred_col <- NULL
      if (!is.null(positive_class)) {
        candidate <- paste0("pred.", positive_class)
        if (candidate %in% pred_cols) {
          preferred_col <- candidate
        }
      }
      if (is.null(preferred_col)) {
        preferred_col <- if (length(pred_cols) >= 2L) pred_cols[2L] else pred_cols[1L]
      }

      tree_data$prediction <- tree_data[[preferred_col]]
    }

    keep_cols <- c(
      "nodeID", "leftChild", "rightChild",
      "splitvarName", "splitval", "prediction"
    )
    tree_data[, keep_cols, drop = FALSE]
  })

  feature_names <- rf_model$forest$independent.variable.names
  ranger_common(x = tree_list, n = n_trees, data = data, feature_names = feature_names)
}


#' @noRd
.unifyTreeModel <- function(fitted_workflow, reference_data, positive_class = NULL) {
  engine_fit <- workflows::extract_fit_engine(fitted_workflow)

  if (inherits(engine_fit, "ranger")) {
    unified_model <- tryCatch(
      treeshap::ranger.unify(engine_fit, data = reference_data),
      error = function(e) {
        .rangerUnifyRobust(
          rf_model = engine_fit,
          data = reference_data,
          positive_class = positive_class
        )
      }
    )

    return(list(
      unified_model = unified_model,
      engine = "RF"
    ))
  }

  if (inherits(engine_fit, "xgb.Booster")) {
    return(list(
      unified_model = treeshap::xgboost.unify(engine_fit, data = as.matrix(reference_data)),
      engine = "XGB"
    ))
  }

  stop(
    "TreeSHAP is currently implemented only for ranger and xgboost engines."
  )
}


#' @noRd
.computeSriMatrices <- function(shap_df, interactions, eps = 1e-12) {
  feature_names <- colnames(shap_df)
  n_features <- length(feature_names)

  synergy <- matrix(0, nrow = n_features, ncol = n_features,
                    dimnames = list(feature_names, feature_names))
  redundancy <- matrix(0, nrow = n_features, ncol = n_features,
                       dimnames = list(feature_names, feature_names))
  independence <- matrix(0, nrow = n_features, ncol = n_features,
                         dimnames = list(feature_names, feature_names))

  shap_mat <- as.matrix(shap_df)

  for (i in seq_len(n_features)) {
    pi <- shap_mat[, i]
    denom_i <- sum(pi^2, na.rm = TRUE)

    if (is.na(denom_i) || denom_i < eps) {
      independence[i, i] <- 1
      next
    }

    for (j in seq_len(n_features)) {
      if (i == j) {
        independence[i, j] <- 1
        next
      }

      pj <- shap_mat[, j]
      pij <- interactions[i, j, ]
      pii <- interactions[i, i, ]
      pjj <- interactions[j, j, ]

      pij_corr <- .orthogonalizeInteraction(pij, pii, pjj, eps = eps)

      s_ij <- .projectVector(pi, pij_corr, eps = eps)
      s_ji <- .projectVector(pj, pij_corr, eps = eps)

      a_ij <- pi - s_ij
      a_ji <- pj - s_ji

      r_ij <- .projectVector(a_ij, a_ji, eps = eps)
      i_ij <- a_ij - r_ij

      s_val <- sum(s_ij^2, na.rm = TRUE) / denom_i
      r_val <- sum(r_ij^2, na.rm = TRUE) / denom_i
      i_val <- sum(i_ij^2, na.rm = TRUE) / denom_i

      s_val <- max(0, min(1, s_val))
      r_val <- max(0, min(1, r_val))
      i_val <- max(0, min(1, i_val))

      total <- s_val + r_val + i_val
      if (is.finite(total) && total > eps) {
        s_val <- s_val / total
        r_val <- r_val / total
        i_val <- i_val / total
      }

      synergy[i, j] <- s_val
      redundancy[i, j] <- r_val
      independence[i, j] <- i_val
    }
  }

  list(
    synergy = synergy,
    redundancy = redundancy,
    independence = independence
  )
}


#' @noRd
.aggregateMatrixList <- function(matrices, component) {
  valid <- matrices[!vapply(matrices, is.null, logical(1L))]

  if (length(valid) == 0L) {
    return(NULL)
  }

  feature_names <- rownames(valid[[1L]][[component]])
  acc <- matrix(0, nrow = length(feature_names), ncol = length(feature_names),
                dimnames = list(feature_names, feature_names))

  for (obj in valid) {
    acc <- acc + obj[[component]][feature_names, feature_names, drop = FALSE]
  }

  acc / length(valid)
}


#' @noRd
.buildEcologicalProfiles <- function(sri_agg, global_importance, top_n = 20L) {
  if (is.null(sri_agg$synergy) || nrow(sri_agg$synergy) == 0L) {
    return(data.frame())
  }

  feature_rank <- head(global_importance$variable, top_n)
  feature_rank <- feature_rank[feature_rank %in% rownames(sri_agg$synergy)]

  if (length(feature_rank) == 0L) {
    return(data.frame())
  }

  out <- lapply(feature_rank, function(feature_name) {
    partner_names <- setdiff(feature_rank, feature_name)

    if (length(partner_names) == 0L) {
      synergy_frac <- 0
      redundancy_frac <- 0
      independence_frac <- 1
    } else {
      synergy_frac <- mean(sri_agg$synergy[feature_name, partner_names], na.rm = TRUE)
      redundancy_frac <- mean(sri_agg$redundancy[feature_name, partner_names], na.rm = TRUE)
      independence_frac <- mean(sri_agg$independence[feature_name, partner_names], na.rm = TRUE)
      total <- synergy_frac + redundancy_frac + independence_frac
      if (is.finite(total) && total > 0) {
        synergy_frac <- synergy_frac / total
        redundancy_frac <- redundancy_frac / total
        independence_frac <- independence_frac / total
      }
    }

    imp_row <- global_importance[global_importance$variable == feature_name, , drop = FALSE]
    total_importance <- imp_row$mean_shap[1]
    direction <- imp_row$direction[1]

    data.frame(
      variable = feature_name,
      total_importance = total_importance,
      synergy_frac = synergy_frac,
      redundancy_frac = redundancy_frac,
      independence_frac = independence_frac,
      synergy_value = total_importance * synergy_frac,
      redundancy_value = total_importance * redundancy_frac,
      independence_value = total_importance * independence_frac,
      direction = direction,
      ecological_role = ifelse(direction >= 0, "Pathogen", "Protector"),
      stringsAsFactors = FALSE
    )
  })

  dplyr::bind_rows(out) %>%
    dplyr::arrange(dplyr::desc(total_importance))
}


#' @noRd
.buildEdgeTable <- function(component_matrix, top_features) {
  if (is.null(component_matrix) || length(top_features) < 2L) {
    return(data.frame())
  }

  edge_rows <- list()
  idx <- 1L

  for (i in seq_len(length(top_features) - 1L)) {
    for (j in (i + 1L):length(top_features)) {
      a <- top_features[i]
      b <- top_features[j]
      weight <- mean(
        c(component_matrix[a, b], component_matrix[b, a]),
        na.rm = TRUE
      )

      edge_rows[[idx]] <- data.frame(
        source = a,
        target = b,
        weight = weight,
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }
  }

  dplyr::bind_rows(edge_rows) %>%
    dplyr::arrange(dplyr::desc(weight))
}


#' @noRd
.renderSRINetworkPython <- function(nodes, edges, output_path, title) {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    warning("Package 'reticulate' is not available; skipping Python network rendering.")
    return(NULL)
  }

  py <- reticulate::import_main(convert = FALSE)

  py_code <- paste(
    "import csv",
    "import math",
    "import os",
    "import matplotlib",
    "matplotlib.use('Agg')",
    "",
    "def rumble_render_sri_network(nodes_csv, edges_csv, output_path, title):",
    "    try:",
    "        import networkx as nx",
    "        import matplotlib.pyplot as plt",
    "    except Exception as exc:",
    "        raise RuntimeError('Python modules networkx and matplotlib are required.') from exc",
    "",
    "    G = nx.Graph()",
    "    with open(nodes_csv, newline='', encoding='utf-8') as f:",
    "        reader = csv.DictReader(f)",
    "        for row in reader:",
    "            node = row['variable']",
    "            size = float(row.get('node_size', '1') or '1')",
    "            role = row.get('ecological_role', 'Neutral')",
    "            G.add_node(node, size=size, role=role)",
    "",
    "    with open(edges_csv, newline='', encoding='utf-8') as f:",
    "        reader = csv.DictReader(f)",
    "        for row in reader:",
    "            weight = float(row.get('weight', '0') or '0')",
    "            if weight > 0:",
    "                G.add_edge(row['source'], row['target'], weight=weight)",
    "",
    "    if len(G.nodes) == 0:",
    "        raise RuntimeError('No nodes available for network rendering.')",
    "",
    "    pos = nx.spring_layout(G, seed=42, weight='weight')",
    "    color_map = {'Pathogen': '#d62728', 'Protector': '#1f77b4', 'Neutral': '#7f7f7f'}",
    "    node_colors = [color_map.get(G.nodes[n].get('role', 'Neutral'), '#7f7f7f') for n in G.nodes]",
    "    node_sizes = [300 + 1800 * math.sqrt(max(G.nodes[n].get('size', 0), 0)) for n in G.nodes]",
    "    edge_widths = [1 + 8 * G[u][v].get('weight', 0) for u, v in G.edges]",
    "    edge_alphas = [min(0.9, 0.2 + G[u][v].get('weight', 0)) for u, v in G.edges]",
    "",
    "    plt.figure(figsize=(10, 8))",
    "    for (u, v), width, alpha in zip(G.edges, edge_widths, edge_alphas):",
    "        nx.draw_networkx_edges(G, pos, edgelist=[(u, v)], width=width, alpha=alpha, edge_color='#7f7f7f')",
    "",
    "    nx.draw_networkx_nodes(G, pos, node_size=node_sizes, node_color=node_colors, linewidths=1.2, edgecolors='black')",
    "    nx.draw_networkx_labels(G, pos, font_size=9)",
    "    plt.title(title)",
    "    plt.axis('off')",
    "    plt.tight_layout()",
    "    os.makedirs(os.path.dirname(output_path), exist_ok=True)",
    "    plt.savefig(output_path, dpi=300, bbox_inches='tight')",
    "    plt.close()",
    sep = "\n"
  )

  reticulate::py_run_string(py_code)

  nodes_csv <- tempfile(fileext = ".csv")
  edges_csv <- tempfile(fileext = ".csv")
  utils::write.csv(nodes, nodes_csv, row.names = FALSE)
  utils::write.csv(edges, edges_csv, row.names = FALSE)

  py$rumble_render_sri_network(nodes_csv, edges_csv, output_path, title)
  output_path
}


#' @noRd
computeTreeShapEcology <- function(final_models,
                                   train_data,
                                   prediction_data,
                                   target_var,
                                   global_importance,
                                   top_n = 20L,
                                   output_dir = NULL,
                                   project_name = "RUMBLE_analysis",
                                   verbose = TRUE) {
  msg <- function(...) {
    if (verbose) message(...)
  }

  if (!requireNamespace("treeshap", quietly = TRUE)) {
    stop("The 'treeshap' package is required for the TreeSHAP ecology stage.")
  }

  msg("Computing TreeSHAP interaction ecology...")

  per_model <- purrr::imap(final_models, function(model_obj, model_name) {
    fitted_wf <- model_obj$model_fit
    processed_train <- .preprocessForTreeShap(fitted_wf, train_data, target_var)
    processed_pred <- .preprocessForTreeShap(fitted_wf, prediction_data, target_var)
    positive_class <- levels(train_data[[target_var]])[2L]
    unified <- .unifyTreeModel(
      fitted_wf,
      processed_train$reference,
      positive_class = positive_class
    )


    ts_obj <- treeshap::treeshap(
      unified_model = unified$unified_model,
      x = processed_pred$processed,
      interactions = TRUE,
      verbose = FALSE
    )

    sri_mats <- .computeSriMatrices(ts_obj$shaps, ts_obj$interactions)

    list(
      model = model_name,
      engine = unified$engine,
      shaps = ts_obj$shaps,
      interactions = ts_obj$interactions,
      sri = sri_mats
    )
  })

  sri_agg <- list(
    synergy = .aggregateMatrixList(purrr::map(per_model, "sri"), "synergy"),
    redundancy = .aggregateMatrixList(purrr::map(per_model, "sri"), "redundancy"),
    independence = .aggregateMatrixList(purrr::map(per_model, "sri"), "independence")
  )

  profile <- .buildEcologicalProfiles(sri_agg, global_importance, top_n = top_n)
  top_features <- profile$variable

  synergy_edges <- .buildEdgeTable(sri_agg$synergy, top_features)
  redundancy_edges <- .buildEdgeTable(sri_agg$redundancy, top_features)

  node_table <- profile %>%
    dplyr::transmute(
      variable = variable,
      node_size = independence_value,
      ecological_role = ecological_role,
      total_importance = total_importance,
      independence_value = independence_value,
      synergy_value = synergy_value,
      redundancy_value = redundancy_value
    )

  network_files <- list(synergy = NULL, redundancy = NULL)
  render_dir <- if (is.null(output_dir)) tempdir() else output_dir

  if (nrow(node_table) > 0L) {
    synergy_path <- file.path(render_dir, paste0(project_name, "_synergy_network.png"))
    redundancy_path <- file.path(render_dir, paste0(project_name, "_redundancy_network.png"))

    try({
      network_files$synergy <- .renderSRINetworkPython(
        nodes = node_table,
        edges = synergy_edges,
        output_path = synergy_path,
        title = "Synergy Topology"
      )
    }, silent = TRUE)

    try({
      network_files$redundancy <- .renderSRINetworkPython(
        nodes = node_table,
        edges = redundancy_edges,
        output_path = redundancy_path,
        title = "Redundancy Topology"
      )
    }, silent = TRUE)
  }

  list(
    per_model = per_model,
    aggregate = sri_agg,
    profile = profile,
    nodes = node_table,
    synergy_edges = synergy_edges,
    redundancy_edges = redundancy_edges,
    network_files = network_files
  )
}
