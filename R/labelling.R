library('progress')

otsu_threshold_multiple <- function(hist_counts, hist_breaks, n_thresholds) {
  total_pixels <- sum(hist_counts)
  thresholds <- numeric(n_thresholds)
  class_intervals <- list()

  # Create a loop to find multiple thresholds
  for (i in 1:n_thresholds) {
    max_variance <- 0
    threshold <- 0
    best_split <- NULL

    # Try all possible splits
    for (t in 2:(length(hist_counts) - 1)) {
      # Split histogram into background and foreground
      weight_b <- sum(hist_counts[1:t]) / total_pixels
      weight_f <- 1 - weight_b

      # Skip invalid splits
      if (weight_b == 0 || weight_f == 0) next  # If either class has zero weight, skip this split

      # Calculate class means
      mean_b <- sum(hist_counts[1:t] * hist_breaks[1:t]) / sum(hist_counts[1:t])
      mean_f <- sum(hist_counts[(t+1):length(hist_counts)] * hist_breaks[(t+1):length(hist_breaks)]) /
        sum(hist_counts[(t+1):length(hist_counts)])

      # Calculate between-class variance
      between_class_variance <- weight_b * weight_f * (mean_b - mean_f)^2

      # Check if between-class variance is a valid number
      if (!is.na(between_class_variance) && between_class_variance > max_variance) {
        max_variance <- between_class_variance
        threshold <- hist_breaks[t]
        best_split <- t
      }
    }

    # Save the threshold and update the histogram for the next iteration
    thresholds[i] <- threshold

    # Split the histogram data for the next iteration (class-wise segmentation)
    hist_counts <- hist_counts[(best_split+1):length(hist_counts)]
    hist_breaks <- hist_breaks[(best_split+1):length(hist_breaks)]
  }

  return(thresholds)
}

multi_threshold_with_nan_by_column <- function(matrix, num_thresholds = 2, verbose = FALSE, algorithm = "otsu") {
  segmented_matrix <- matrix(0, nrow = nrow(matrix), ncol = ncol(matrix))
  thresholds <- list()
  n <- ncol(matrix)

  if (verbose) {
    pb <- txtProgressBar(min = 1, max = n, style = 3)
  }
  for (col_idx in seq_len(n)) {
    if (verbose) {
      setTxtProgressBar(pb, col_idx)
    }
    col <- matrix[, col_idx]
    valid_values <- col[!is.na(col)]

    # Skip empty columns
    if (length(valid_values) == 0) {
      segmented_matrix[, col_idx] <- col
      next
    }

    # Thresholding logic
    if (algorithm == "otsu") {
      if (length(unique(valid_values)) < num_thresholds) {
        # Too few unique values, create equally spaced thresholds
        thresh <- seq(min(valid_values), max(valid_values), length.out = num_thresholds + 1)[-c(1, num_thresholds + 1)]
      } else {
        # Otsu's thresholding (manual for histograms)
        hist_info <- hist(valid_values, breaks = seq(min(valid_values), max(valid_values), length.out = 60), plot = FALSE)
        hist_counts <- hist_info$counts
        hist_breaks <- hist_info$mids

        # Apply custom Otsu's method
        thresh <- otsu_threshold_multiple(hist_counts, hist_breaks, n_thresholds = num_thresholds-1)
      }
    } else if (algorithm == "linspace") {
      # Use linspace method for thresholds
      thresh <- seq(min(valid_values), max(valid_values), length.out = num_thresholds + 1)[-c(1, num_thresholds + 1)]
    } else {
      stop("Thresholding method not supported.")
    }

    # Save thresholds
    thresholds[[col_idx]] <- thresh

    # Adjust the breaks and labels to match
    breaks_with_inf <- c(-Inf, thresh, Inf)  # Add -Inf and Inf to breaks
    labels <- seq_along(thresh)  # Number of labels should match number of thresholds plus one
    col_digitize <- cut(col, breaks = breaks_with_inf, right = FALSE)
    segmented_matrix[, col_idx] <- as.numeric(col_digitize)
  }
  if (verbose)
    close(pb)

  return(list(segmented_matrix = segmented_matrix, thresholds = thresholds))
}

# Manual implementation of Otsu's method for histogram
otsu_histogram <- function(hist_counts, hist_breaks) {
  total_pixels <- sum(hist_counts)
  sum_total <- sum(hist_counts * hist_breaks)

  sum_b <- 0
  w_b <- 0
  max_variance <- 0
  threshold <- 0

  for (i in 1:(length(hist_counts) - 1)) {
    w_f <- total_pixels - w_b
    if (w_b == 0 || w_f == 0) next

    sum_f <- sum(hist_counts[(i+1):length(hist_counts)] * hist_breaks[(i+1):length(hist_breaks)])
    m_b <- sum_b / w_b
    m_f <- sum_f / w_f

    between_class_variance <- w_b * w_f * (m_b - m_f)^2

    if (between_class_variance > max_variance) {
      max_variance <- between_class_variance
      threshold <- hist_breaks[i]
    }

    w_b <- w_b + hist_counts[i]
    sum_b <- sum_b + hist_counts[i] * hist_breaks[i]
  }

  return(threshold)
}

rows_with_all_nan <- function(df) {
  idx <- which(apply(df, 1, function(row) all(is.na(row))))
  return(idx)
}

modemax_nan <- function(a, reducefoo = max) {
  # Ensure `a` is a matrix and check if it's empty
  if (!is.matrix(a) && !is.data.frame(a)) {
    stop("Input `a` must be a matrix or data frame.")
  }
  # If `a` is a data frame, convert to a matrix
  if (is.data.frame(a)) {
    a <- as.matrix(a)
  }
  res <- rep(NA, nrow(a))  # Initialize a vector with NAs, size is the number of rows in the matrix
  # Iterate over rows of the matrix
  for (i in 1:nrow(a)) {
    row <- a[i, , drop = FALSE]  # Extract the ith row

    # Remove missing values (NA) and find the mode
    non_missing_values <- row[!is.na(row)]

    if (length(non_missing_values) > 0) {
      # Calculate the frequency of each unique value in the row
      freq_table <- table(non_missing_values)

      # Get the mode(s), which are the most frequent value(s)
      max_freq <- max(freq_table)
      mode_list <- as.numeric(names(freq_table[freq_table == max_freq]))  # Get values with the maximum frequency

      # If there are multiple modes, apply the reduce function (default is max)
      res[i] <- reducefoo(mode_list)
    } else {
      # If no non-missing values, return NA
      res[i] <- NA
    }
  }

  return(res)  # Return the result as a vector
}

modemax_nan_my <- function(a, reducefoo = max) {
  # Ensure `a` is a matrix and check if it's empty
  if (!is.matrix(a) && !is.data.frame(a)) {
    stop("Input `a` must be a matrix or data frame.")
  }
  # If `a` is a data frame, convert to a matrix
  if (is.data.frame(a)) {
    a <- as.matrix(a)
  }
  res <- rep(NA, nrow(a))  # Initialize a vector with NAs, size is the number of rows in the matrix
  # Iterate over rows of the matrix
  for (i in 1:nrow(a)) {
    row <- a[i, , drop = FALSE]  # Extract the ith row

    # Remove missing values (NA) and find the mode
    non_missing_values <- row[!is.na(row)]

    if (length(non_missing_values) > 0) {
      # Calculate the frequency of each unique value in the row
      freq_table <- table(non_missing_values)

      # Get the mode(s), which are the most frequent value(s)
      max_freq <- max(freq_table)
      mode_list <- as.numeric(names(freq_table[freq_table == max_freq]))  # Get values with the maximum frequency
      # If there are multiple modes, apply the reduce function (default is max)
      res[i] <- reducefoo(mode_list)
    } else {
      # If no non-missing values, return NA
      res[i] <- NA
    }
    print(res[i])
  }

  return(res)  # Return the result as a vector
}

labelling_core <- function(df, columns = character(0), n_classes = 2, show_progress = FALSE,
                           verbose = FALSE, labelnames = list(0 <- "E", 1 <- "NE"), mode = "flat-multi",
                           algorithm = "otsu", rowname = "gene", colname = "label",
                           reducefoo = max) {

  # Get numerical columns (threshold will work on them)
  df_str <- df[, sapply(df, is.character)]
  dfnum <- df[, setdiff(names(df), names(df_str))]

  # Create T matrix from selected columns
  T <- if (length(columns) == 0) as.matrix(dfnum) else as.matrix(dfnum[, columns])

  # Mode: two-by-two
  if (mode == "two-by-two") {
    # Placeholder for multi_threshold_with_nan_by_column (you need to define it)
    result <- multi_threshold_with_nan_by_column(T, 2, algorithm = algorithm, verbose = show_progress)
    Q <- result$segmented_matrix
    Thr <- result$thresholds
    Labels2 <- modemax_nan(Q, reducefoo = reducefoo)

    dfout <- cbind(dfnum[, columns], data.frame(colname = Labels2))

    # Genes where label is 1 (NE genes)
    dff <- subset(dfout, colname == 2)
    dff_str <- dff[, sapply(dff, is.character)]
    dffnum <- dff[setdiff(names(dff), names(dff_str))]

    TNE <- if (length(columns) == 0) as.matrix(dffnum) else as.matrix(dffnum[, columns])

    # Run again on NE genes
    resultNE <- multi_threshold_with_nan_by_column(TNE, 2, algorithm = algorithm, verbose = show_progress)
    QNE <- resultNE$segmented_matrix
    ThrNE <- resultNE$thresholds
    Labels <- modemax_nan(QNE, reducefoo = reducefoo)

    dfout2 <- data.frame(colname = Labels)

    # Re-label the data
    dfout2$colname <- ifelse(dfout2$colname == 2, 3, ifelse(dfout2$colname == 1, 2, dfout2$colname))

    # Merge back to dfout
    NE_genes <- which(dfout$colname == 2)
    dfout$colname[NE_genes] <- dfout2$colname
    return(dfout)

  } else if (mode == "flat-multi") {
    # Placeholder for multi_threshold_with_nan_by_column (you need to define it)
    result <- multi_threshold_with_nan_by_column(T, n_classes, algorithm = algorithm, verbose = show_progress)
    Q <- result$segmented_matrix
    Thr <- result$thresholds
    Labels <- modemax_nan(Q, reducefoo = reducefoo)

    dfout <- cbind(df[, names(df_str)], data.frame(colname = Labels))
    return(dfout)

  } else {
    stop(paste("Labelling mode", mode, "not supported"))
  }
}

labelling <- function(
    df,
    columns,
    n_classes = 2,
    show_progress = FALSE,
    verbose = FALSE,
    labelnames = list("2" = "NE", "1" = "E"),
    mode = "flat-multi",
    rowname = "gene",
    colname = "label",
    algorithm = "otsu",
    reducefoo = max,
    warning = FALSE
) {
  # Adjust n_classes for "two-by-two" mode
  if (mode == "two-by-two") {
    n_classes <- 3
  }

  # Ensure labelnames length matches n_classes
  if (length(labelnames) != n_classes) {
    stop("Label dictionary not the same size as the number of classes!")
  }

  # If columns are a list of vectors
  if (all(sapply(columns, is.character)) && is.list(columns)) {
    if (verbose) {
      message(sprintf("Performing mode of mode on %d-class labelling (%s) on %d contexts", n_classes, mode,length(columns)))
    }

    # Initialize matrix to hold results
    L_tot <- matrix(NA, nrow = nrow(df), ncol = 0)

    # Loop through column groups
    for (lines in columns) {
      # Get labels for the current group of columns
      labels <- labelling_core(
        df,
        columns = lines,
        verbose = verbose,
        show_progress = show_progress,
        mode = mode,
        labelnames = labelnames,
        rowname = rowname,
        colname = colname,
        n_classes = n_classes,
        algorithm = algorithm,
        reducefoo = reducefoo
      )
      # Append current labels to result matrix
      L_tot <- cbind(L_tot, data.frame(colname = labels))
    }
    # Get the mode of the modes
    mode_of_mode <- modemax_nan(L_tot, reducefoo = reducefoo)
    dfout <- data.frame(colname = mode_of_mode)
    rownames(dfout) <- rownames(df)
  } else if (any(sapply(columns, is.character)) && is.list(columns)) {
    stop("Wrong columns partition format.")
  } else {
    if (verbose) {
      message(sprintf("Performing flat mode on %d-class labelling (%s) on %d cell lines.",
                      n_classes, mode, length(columns)))
    }

    # Identify rows with all NaNs (if needed)
    nanrows <- rows_with_all_nan(df[columns])

    if (length(nanrows) > 0 && warning) {
      warning("There are rows with all NaNs. Please remove them using the function 'rows_with_all_nan()' and reapply the labelling. Otherwise, you will have NaN labels in your output.")
    }

    # Perform labelling
    dfout <- labelling_core(
      df,
      columns = columns,
      verbose = verbose,
      show_progress = show_progress,
      mode = mode,
      labelnames = labelnames,
      rowname = rowname,
      colname = colname,
      n_classes = n_classes,
      algorithm = algorithm,
      reducefoo = reducefoo
    )
  }

  # Replace numerical labels with their names
  dfout$colname <- unname(labelnames[as.character(dfout$colname)])
  colnames(dfout)[colnames(dfout) == "colname"] <- colname
  #dfout[[colname]] <- sapply(dfout[[colname]], function(x) labelnames[[as.character(x)]])

  if (verbose) {
    unique_values <- unique(dfout[[colname]])
    counts <- sapply(unique_values, function(x) sum(dfout[[colname]] == x))
    count_df <- data.frame(Value = unlist(unique_values), Count = counts)
    print(count_df)
  }
  return(dfout)
}

