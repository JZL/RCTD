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

