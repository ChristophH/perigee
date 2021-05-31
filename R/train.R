
#' Training a Perigee model
#'
#' @param counts Matrix of single-cell RNA-seq count data; columns are cells
#' @param cell_labels Cell type labels
#' @param de_res Optional; results of differential expression tests
#' (output of sctransform::diff_mean_test with compare argument set to 'all_vs_all');
#' if NULL, DE testing will be done during training; default is NULL
#' @param max_cells Maximum number of cells to use per type; default uses all cells
#' @param exclude Cell type labels to exclude from training and the final model
#' @param verbosity Level of printed output; 0 is quiet; 1 provides progress information;
#' 2 provide more details; default is 1
#'
#' @return A Perigee model (named list)
#'
#' @importFrom Matrix Matrix crossprod colSums
#' @importFrom stats quantile approx
#' @importFrom utils object.size
#'
#' @export
pg_train <- function(counts, cell_labels, de_res = NULL, max_cells = NULL,
                     exclude = NULL, verbosity = 1) {
  if (!inherits(x = counts, what = 'dgCMatrix')) {
    stop('Problem with input: counts must be a dgCMatrix')
  }
  if (length(cell_labels) != ncol(counts)) {
    stop('Problem with input: length of cell_labels must be equal to the number of columns in counts')
  }
  if (!is.null(max_cells)) {
    if(!is.numeric(max_cells) || length(max_cells) > 1) {
      stop('Problem with input: max_cells must be NULL or a numeric of length 1 that can be converted to integer')
    }
  }
  times <- list(start_time = Sys.time())

  cell_labels <- droplevels(factor(cell_labels))
  sorted_levels <- names(sort(table(cell_labels), decreasing = TRUE))
  sorted_levels <- setdiff(sorted_levels, exclude)
  l <- length(sorted_levels)
  L <- (l^2-l)/2
  # for each combination of cell labels
  # create binary classifier
  classifiers <- list()
  cntr <- 0
  for (i in 1:(l-1)) {
    label1 <- sorted_levels[i]
    sel1 <- which(cell_labels == label1)
    # down-sample if max_cells is set
    if (!is.null(max_cells) && length(sel1) > max_cells) {
      sel1 <- sort(sample(x = sel1, size = max_cells))
    }
    counts1 <- counts[, sel1]

    for (j in (i+1):l) {
      label2 <- sorted_levels[j]
      comp <- paste(label1, label2, sep = ' vs ')
      cntr <- cntr + 1
      if (verbosity > 0) {
        message(cntr, ' of ', L, ': ', comp)
      }

      sel2 <- which(cell_labels == label2)
      # down-sample if max_cells is set
      if (!is.null(max_cells) && length(sel2) > max_cells) {
        sel2 <- sort(sample(x = sel2, size = max_cells))
      }
      counts2 <- counts[, sel2]

      disc_feats <- get_discriminating_features(counts1, counts2, label1, label2, de_res, verbosity)
      disc_feats_comb <- unlist(disc_feats)

      counts1_goi <- counts1[disc_feats_comb, ]
      mean1 <- sctransform:::row_gmean(counts1_goi)
      counts2_goi <- counts2[disc_feats_comb, ]
      mean2 <- sctransform:::row_gmean(counts2_goi)

      # log-transform and length-normalize the mean vectors
      mean1 <- log1p(mean1)
      mean1 <- mean1 / sqrt(sum(mean1^2))
      mean2 <- log1p(mean2)
      mean2 <- mean2 / sqrt(sum(mean2^2))
      mean_diff <- mean1 - mean2

      # get the cos sim diff
      counts_log <- cbind(counts1_goi, counts2_goi)
      counts_log@x <- log(counts_log@x + 1)
      tmp <- crossprod(counts_log, Matrix(mean_diff))
      c_mag <- sqrt(colSums(counts_log^2))
      diff_cos_sim <- tmp / c_mag

      # get BAC-based threshold - bootstrap the process and keep 95% confidence interval
      thv <- get_th(type = rep(x = c(label1, label2), times = c(ncol(counts1_goi), ncol(counts2_goi))),
                    value = diff_cos_sim,
                    label1 = label1,
                    B = 111)
      th_range <- quantile(thv[1, ], probs = c(0.025, 0.975))
      if (verbosity > 1) {
        message('  Threshold is between ',
                round(th_range[1], 3), ' and ', round(th_range[2], 3), ' (95% CI)')
        message('  BAC is between ',
                round(quantile(thv[2, ], 0.025), 3), ' and ', round(quantile(thv[2, ], 0.975), 3), ' (95% CI)')
      }

      # we need sensitivity (TPR) and specificity (TNR) wrt label1
      # this is what the prediction score will be based on
      o <- order(-diff_cos_sim)
      labels_sorted = rep(x = c(label1, label2), times = c(ncol(counts1), ncol(counts2)))[o]
      spec <- rep(NA_real_, length(diff_cos_sim))
      spec[o] <- (length(sel2) - cumsum(labels_sorted == label2)) / length(sel2)
      sens <- rep(NA_real_, length(diff_cos_sim))
      sens[o] <- cumsum(labels_sorted == label1) / length(sel1)

      length_out <- 101
      # specificity will be score of group1
      spec_approx <- data.frame(dcs = seq(from = min(diff_cos_sim), to = max(diff_cos_sim), length.out = length_out))
      spec_approx$value <- approx(x = diff_cos_sim, y = spec, xout = spec_approx$dcs)$y
      # sensitivity will be score of group2
      sens_approx <- data.frame(dcs = seq(from = min(diff_cos_sim), to = max(diff_cos_sim), length.out = length_out))
      sens_approx$value <- approx(x = diff_cos_sim, y = sens, xout = sens_approx$dcs)$y

      # potential optimization: cut off repetitive values from head and tail

      classifiers[[comp]] <- list(label1 = label1, label2 = label2, mean_diff = mean_diff,
                                  th_range = th_range, spec_approx = spec_approx,
                                  sens_approx = sens_approx)


    }
  }
  # all we really need is the sparse matrix of differences of mean of the discriminating genes
  # the names of the classes
  # and the thresholds
  disc_feats <- unlist(lapply(classifiers, function(x) names(x$mean_diff)), use.names = FALSE)
  unique_disc_feats <- unique(disc_feats)
  mean_diff <- matrix(0, length(unique_disc_feats), length(classifiers))
  rownames(mean_diff) <- unique_disc_feats
  for (i in 1:length(classifiers)) {
    mean_diff[names(classifiers[[i]]$mean_diff), i] <- classifiers[[i]]$mean_diff
  }
  mean_diff <- Matrix(mean_diff)
  final_cell_labels <- setdiff(levels(cell_labels), exclude)
  comparisons <- cbind(match(x = sapply(classifiers, function(x) x$label1), table = final_cell_labels),
                       match(x = sapply(classifiers, function(x) x$label2), table = final_cell_labels))
  thresholds <- sapply(classifiers, function(x) x$th_range)
  dimnames(thresholds) <- NULL
  pred_score1 <- lapply(classifiers, function(x) x$spec_approx)
  pred_score2 <- lapply(classifiers, function(x) x$sens_approx)
  res <- list(mean_diff = mean_diff, comparisons = comparisons, thresholds = thresholds,
              cell_labels = final_cell_labels, pred_score1 = pred_score1,
              pred_score2 = pred_score2)
  if (verbosity > 0) {
    message('Model size in memory: ', format(object.size(x = res), standard = 'SI', unit = 'auto'))
    message('Time passed: ', format(round(difftime(Sys.time(), times$start_time), 2)))
  }
  return(res)
}


#' Helper function to extract the discriminating features
#'
#' @param counts1 Count matrix for group 1 (columns are cells)
#' @param counts2 Count matrix for group 2 (columns are cells)
#' @param label1 The label of group 1
#' @param label2 The label of group 2
#' @param de_res Optional; results of differential expression tests
#' (output of sctransform::diff_mean_test with compare argument set to 'all_vs_all');
#' if not provided, DE testing will be done during training; default is NULL
#' @param verbosity Level of printed output; 0 is quiet; default is 0
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr filter mutate arrange slice_head pull
#' @importFrom rlang .data
#' @importFrom sctransform diff_mean_test
get_discriminating_features <- function(counts1, counts2, label1, label2, de_res = NULL,
                                        verbosity = 0) {
  if (!is.null(de_res)) {
    tmp <- filter(de_res, .data$group1 == label1, .data$group2 == label2)
    if (nrow(tmp) < 1) {
      # maybe the groups are switched
      tmp <- filter(de_res, .data$group1 == label2, .data$group2 == label1)
      if (nrow(tmp) < 1) {
        stop('Could not find the requested group comparison in provided DE results: ', label1, ' ', label2)
      } else {
        de_res <- tmp
        de_res$group1 <- de_res$group2
        de_res$mean1 <- de_res$mean2
        de_res$cells1 <- de_res$cells2
        de_res$group2 <- tmp$group1
        de_res$mean2 <- tmp$mean1
        de_res$cells2 <- tmp$cells1
        de_res <- mutate(de_res, mean_diff = -.data$mean_diff, log2FC = -.data$log2FC, zscore = -.data$zscore)
      }
    } else {
      de_res <- tmp
    }
  } else {
    de_res <- diff_mean_test(y = cbind(counts1, counts2),
                             group_labels = rep(x = c(label1, label2),
                                                times = c(ncol(counts1), ncol(counts2))),
                             R = 499, log2FC_th = log2(1.5), mean_th = 0.05,
                             compare = c(label1, label2), verbosity = verbosity)
  }
  set1 <- filter(de_res, .data$mean1 > 0.05, .data$log2FC > 0, .data$emp_pval < 0.01) %>%
    mutate(rm = rank(.data$mean1), rfc = rank(abs(.data$log2FC))) %>%
    arrange(-(.data$rm+.data$rfc)) %>% slice_head(n = 50) %>% pull(.data$gene)
  set2 <- filter(de_res, .data$mean2 > 0.05, .data$log2FC < 0, .data$emp_pval < 0.01) %>%
    mutate(rm = rank(.data$mean2), rfc = rank(abs(.data$log2FC))) %>%
    arrange(-(.data$rm+.data$rfc)) %>% slice_head(n = 50) %>% pull(.data$gene)
  return(list(set1, set2))
}

# helper function to get score threshold based on BAC
get_th <- function(type, value, label1, B = 1) {
  n <- length(value)
  if (B == 1) {
    p <- sum(type == label1)
    o <- order(-value)
    sens <- rep(NA_real_, n)
    sens[o] <- cumsum(type[o] == label1) / p
    spec <- rep(NA_real_, n)
    spec[o] <- (n - p - cumsum(type[o] != label1)) / (n - p)
    bac <- (sens + spec) / 2
    m <- which.max(bac)
    return(c(value[m], bac[m]))
  }
  # use sampling with replacement to get threshold distribution
  th <- sapply(1:B, function(i) {
    idx <- sample(x = 1:n, size = n, replace = TRUE)
    return(get_th(type[idx], value[idx], label1, B = 1))
  })
  return(th)
}
