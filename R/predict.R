#' Cell type prediction with Perigee
#' @param pg_model A Perigee model created with pg_train
#' @param counts Matrix of single-cell RNA-seq count data
#' @param return_dcs Should the output contain the full matrix of pairwise
#' differences of cosine similarities? Default is FALSE
#' @param verbosity 0 is quiet; larger integers may lead to more print
#' output; default is 1
#'
#' @return A data frame of predictions with columns
#'
#' Output columns:
#' \describe{
#' \item{label}{The predicted cell type label}
#' \item{n_comp}{The number of pairwise comparisons that support the predicted label}
#' \item{score}{The prediction score (based on model specificity)}
#' }
#'
#' @importFrom stats approxfun
#'
#' @export
pg_predict <- function(pg_model, counts, return_dcs = FALSE, verbosity = 1) {
  if (!inherits(x = counts, what = 'dgCMatrix')) {
    stop('Problem with input: counts must be a dgCMatrix')
  }
  start_time = Sys.time()

  missing_features <- setdiff(rownames(pg_model$mean_diff), rownames(counts))
  if (length(missing_features) > 0) {
    w <- sprintf(paste0('%d features out of %d used by the classifier are missing in ',
                        'the input data. This could be simply because these features are ',
                        'not detected in your data, but it could also be due to a different ',
                        'annotation than used for training.'),
                 length(missing_features), nrow(pg_model$mean_diff))
    warning(w)
    missing_counts <- Matrix(data = 0, nrow = length(missing_features), ncol = ncol(counts),
                             dimnames = list(missing_features, colnames(counts)))
    counts <- rbind(counts, missing_counts)
  }

  counts2 <- counts[rownames(pg_model$mean_diff), ]
  counts2@x <- log(counts2@x + 1)
  tmp <- crossprod(counts2, pg_model$mean_diff)
  c_mag <- sqrt(crossprod(counts2^2, pg_model$mean_diff != 0))
  diff_cos_sim <- tmp / c_mag

  # extract the upper and lower bound for the threshold
  th1 <- pg_model$thresholds[2, ]
  th2 <- pg_model$thresholds[1, ]
  # create functions to turn delta sim into prediction score
  pred_score1 <- lapply(pg_model$pred_score1, function(df) {
    approxfun(x = df$dcs, y = df$value, rule = 2)
  })
  pred_score2 <- lapply(pg_model$pred_score2, function(df) {
    approxfun(x = df$dcs, y = df$value, rule = 2)
  })

  pred <- apply(diff_cos_sim, 1, function(y) {
    class1 <- y > th1
    class2 <- y < th2
    pred_int <- c(pg_model$comparisons[class1, 1],
                  pg_model$comparisons[class2, 2])
    tab <- tabulate(pred_int)
    m <- which.max(tab)
    sup <- tab[m]

    # handle ties
    candidates <- which(tab == sup)
    pred_scores <- vapply(X = candidates, FUN = function(mc) {
      idx1 <- which(pg_model$comparisons[, 1] == mc)  # where we will use pred_score1
      idx2 <- which(pg_model$comparisons[, 2] == mc)  # where we will use pred_score2
      score1 <- vapply(X = idx1, FUN = function(i) pred_score1[[i]](y[i]), FUN.VALUE = 1)
      score2 <- vapply(X = idx2, FUN = function(i) pred_score2[[i]](y[i]), FUN.VALUE = 1)
      prod(score1, score2)
    }, FUN.VALUE = 1)
    j <- which.max(pred_scores)
    m <- candidates[j]
    if (length(m) == 0) {
      return(c(NA_real_, NA_real_ , NA_real_))
    }
    return(c(m, tab[m], pred_scores[j]))
  })
  pred_df <- data.frame(t(pred))
  colnames(pred_df) <- c('label', 'n_comp', 'score')
  pred_df$label <- factor(pg_model$cell_labels[pred_df$label], levels = pg_model$cell_labels)

  if (!return_dcs) {
    ret <- pred_df
  } else {
    ret <- list(prediction = pred_df,
                dcs = diff_cos_sim)
  }
  if (verbosity > 0) {
    message('Time passed: ', format(round(difftime(Sys.time(), start_time), 2)))
  }
  return(ret)
}
