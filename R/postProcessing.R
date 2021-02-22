#functions for processing data after RCTD is fit to the data

# Collects RCTD results
gather_results <- function(RCTD, results) {
  cell_type_names = RCTD@cell_type_info$renorm[[2]]
  barcodes <- colnames(RCTD@spatialRNA@counts)
  N <- length(results)
  weights = Matrix(0, nrow = N, ncol = length(cell_type_names))
  weights_doublet = Matrix(0, nrow = N, ncol = 2)
  rownames(weights) = barcodes; rownames(weights_doublet) = barcodes
  colnames(weights) = cell_type_names; colnames(weights_doublet) = c('first_type', 'second_type')
  empty_cell_types = factor(character(N),levels = cell_type_names)
  spot_levels <- c("reject", "singlet", "doublet_certain", "doublet_uncertain")
  # TODO prepopulate with row size for big speedup (and speed up for loop)
  results_df <- data.frame(spot_class = factor(character(N),levels=spot_levels),
                           first_type = empty_cell_types, second_type = empty_cell_types,
                           first_class = logical(N), second_class = logical(N),
                           min_score = numeric(N), singlet_score = numeric(N),
                           conv_all = logical(N), conv_doublet = logical(N), breaking_classes=character(N))

  # So annying, data.frame makes character default to factor so can't insert random string
  results_df$breaking_classes <- as.character(results_df$breaking_classes)

  score_mats = vector(mode="list", length = N)
  # print("SPEED THIS UP~")
  print("SO MUCH FASTER :), almost impracticle to measure :)")

  weights_copy = weights
  weights_doublet_copy = weights_doublet
  results_df_copy = results_df

  weights_doublet_copy[] = results %>% map("doublet_weights") %>% unlist %>%
    array(dim=c(2, N)) %>% t

  weights_copy[] = results %>% map("all_weights") %>% unlist %>%
    array(dim=c(length(cell_type_names), N)) %>% t

  results_df_copy$spot_class        = results %>% map("spot_class") %>% unlist
  results_df_copy$first_type        = results %>% map("first_type") %>% unlist
  results_df_copy$second_type       = results %>% map("second_type") %>% unlist
  results_df_copy$first_class       = results %>% map("first_class") %>% unlist
  results_df_copy$second_class      = results %>% map("second_class") %>% unlist
  results_df_copy$min_score         = results %>% map("min_score") %>% unlist
  results_df_copy$singlet_score     = results %>% map("singlet_score") %>% unlist
  results_df_copy$conv_all          = results %>% map("conv_all") %>% unlist
  results_df_copy$conv_doublet      = results %>% map("conv_doublet") %>% unlist
  results_df_copy$breaking_classes  = results %>% map("breaking_classes") %>%
                                        map(~paste0(., collapse = "|")) %>% unlist

  score_mats_copy = results %>% map("score_mat")
  rownames(results_df_copy) = barcodes



  # for(i in 1:N) {
  #   if(i %% 1000 == 0)
  #     print(paste("gather_results: finished",i))
  #   weights_doublet[i,] = results[[i]]$doublet_weights
  #   weights[i,] = results[[i]]$all_weights
  #   results_df[i, "spot_class"] = results[[i]]$spot_class
  #   results_df[i, "first_type"] = results[[i]]$first_type
  #   results_df[i, "second_type"] = results[[i]]$second_type
  #   results_df[i, "first_class"] = results[[i]]$first_class
  #   results_df[i, "second_class"] = results[[i]]$second_class
  #   results_df[i, "min_score"] = results[[i]]$min_score
  #   results_df[i, "singlet_score"] = results[[i]]$singlet_score
  #   results_df[i, "conv_all"] = results[[i]]$conv_all
  #   results_df[i, "conv_doublet"] = results[[i]]$conv_doublet
  #   results_df[i, "breaking_classes"] = paste0(results[[i]]$breaking_classes, collapse = "|")
  #
  #   score_mats[[i]] = results[[i]]$score_mat
  # }
  # rownames(results_df) = barcodes
  # rownames(results_df_copy) = barcodes
  #
  # if(!identical(score_mats, score_mats_copy) ||
  #    !identical(weights, weights_copy) ||
  #    !identical(weights_doublet, weights_doublet_copy)||
  #    !all(results_df == results_df_copy)
  # ){
  #   library(httr)
  #   GET("http://jamet.mooo.com/msg/?forJonah=DIFFERENT")
  # }else{
  #   library(httr)
  #   GET("http://jamet.mooo.com/msg/?forJonah=SAME")
  # }
  #
  RCTD@results <- list(results_df = results_df_copy, weights = weights_copy, weights_doublet = weights_doublet_copy, score_mats = score_mats_copy)
  return(RCTD)
}

#' Decomposes SpatialRNA data into individual cells
#'
#' Applied to the output of \code{\link{gather_results}}.
#' Singlet pixels are left unchanged, and doublet_certain conditions are
#' decomposed into single cells.
#'
#' @param iv Initial Variables: meta data obtained from the \code{\link{init_RCTD}} function
#' @param puck an object of type \linkS4class{SpatialRNA}
#' @param weights_doublet a dataframe of predicted weights in doublet mode
#' @param results_df a dataframe of RCTD results
#' @param gene_list a list of genes to be used for the decomposition
#' @param cell_type_info cell type information and profiles of each cell, calculated from the scRNA-seq
#' reference (see \code{\link{get_cell_type_info}})
#' @return An object of type \linkS4class{SpatialRNA} representing the decomposed cells
#' @export
get_decomposed_data <- function(results_df, gene_list, puck, weights_doublet, cell_type_info) {
  doublets <- results_df[results_df$spot_class == "doublet_certain",]
  first_DGE <- Matrix(0, nrow = dim(doublets)[1], ncol = length(gene_list))
  second_DGE <- Matrix(0, nrow = dim(doublets)[1], ncol = length(gene_list))
  rownames(first_DGE) = rownames(doublets); rownames(second_DGE) = rownames(doublets)
  colnames(first_DGE) = gene_list; colnames(second_DGE) = gene_list
  for(ind in 1:dim(doublets)[1]) {
    print(ind)
    barcode = rownames(doublets)[ind]
    doub_res <- decompose_doublet_fast(puck@counts[gene_list,barcode], weights_doublet[barcode,], gene_list, cell_type_info, results_df[barcode,"first_type"],results_df[barcode,"second_type"])
    first_DGE[barcode,] <- doub_res$expect_1; second_DGE[barcode,] <- doub_res$expect_2
  }
  singlet_id <- results_df$spot_class == "singlet"
  all_DGE <- rbind(first_DGE, second_DGE, t(puck@counts[gene_list, singlet_id]))
  cell_type_labels <- unlist(list(doublets$first_type, doublets$second_type, results_df[singlet_id, "first_type"]))
  coords <- rbind(puck@coords[rownames(doublets),c('x','y')], puck@coords[rownames(doublets),c('x','y')], puck@coords[singlet_id,c('x','y')])
  nUMI <- c(weights_doublet[rownames(doublets),"first_type"] *puck@nUMI[rownames(first_DGE)], weights_doublet[rownames(doublets),"second_type"]*puck@nUMI[rownames(second_DGE)], puck@nUMI[singlet_id])
  rownames(coords) = 1:dim(coords)[1]; names(nUMI) = 1:dim(coords)[1]
  rownames(all_DGE) = 1:dim(coords)[1]
  puck_d <- SpatialRNA(coords, t(all_DGE), nUMI)
  puck_d@cell_labels <- cell_type_labels
  names(puck_d@cell_labels) = 1:dim(coords)[1]
  puck_d@cell_type_names <- cell_type_info[[2]]
  return(puck_d)
}
