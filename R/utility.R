# HELPr utility functions
#

library(dplyr)
#' Filter a CRISPR DataFrame by Model
#' This function filters a CRISPR dataset (`df`) based on model-related criteria provided in a mapping data frame (`df_map`).
#'
#' @param df A data frame containing CRISPR data to be filtered.
#' @param df_map A mapping data frame containing model-related information used for filtering.
#' @param minlines An integer specifying the minimum number of lines required for a model to be included in the filtered dataset. Default is `1`.
#' @param line_colname A string specifying the column name in `df_map` that contains the model identifiers. Default is `"ModelID"`.
#' @param line_group A string specifying the grouping column name in `df_map`, typically representing lineage or category. Default is `"OncotreeLineage"`.
#'
#' @return A filtered data frame containing rows from `df` filtered based on the selected cell lines and conditions.
#'
#' @examples
#' # Example usage:
#' df <- data.frame(ModelID = c("A", "B", "C"), Value = c(1.2, 3.4, 5.6))
#' df_map <- data.frame(ModelID = c("A", "B", "D"), OncotreeLineage = c("Lineage1", "Lineage2", "Lineage3"))
#' filtered_df <- filter_crispr_by_model(df, df_map, minlines = 1)
#'
#' @export
filter_crispr_by_model <- function(df, df_map,
                                   minlines = 1,
                                   line_colname = "ModelID",
                                   line_group = "OncotreeLineage") {

  # Get cell lines from the mapping DataFrame
  map_cell_lines <- df_map %>%
    filter(!is.na(.data[[line_group]])) %>%
    pull(.data[[line_colname]])

  # Intersect cell lines in the CRISPR DataFrame with cell lines in the mapping DataFrame
  dep_cell_lines <- intersect(colnames(df), map_cell_lines)
  rest_columns <- setdiff(colnames(df), map_cell_lines)

  # Filter mapping DataFrame based on common cell lines
  df_map_filtered <- df_map %>%
    filter(.data[[line_colname]] %in% dep_cell_lines)

  # Select tissue models with lines greater than or equal to minlines
  sel_cell_lines <- filter_cellmap(df_map_filtered, minlines = minlines, line_group = line_group)

  # Return filtered CRISPR DataFrame based on selected cell lines
  filtered_df <- df %>%
    select(all_of(c(rest_columns, intersect(colnames(df), sel_cell_lines[[line_colname]]))))

  return(filtered_df)
}
filter_cellmap <- function(df_map, minlines = 1, line_group = "OncotreeLineage") {
  # Filters a cell map DataFrame based on the minimum number of lines per group.
  #
  # Arguments:
  # df_map: The input DataFrame containing cell map information.
  # minlines: The minimum number of lines required to retain a group. Default is 1.
  # line_group: Column name for the grouping information in the cell map DataFrame. Default: "OncotreeLineage".
  #
  # Returns:
  # A filtered DataFrame containing only the groups that meet the minimum lines criteria.

  # Calculate group counts and filter groups meeting the minlines criteria
  tissue_list <- df_map %>%
    group_by(.data[[line_group]]) %>%
    summarize(count = n(), .groups = "drop") %>%
    filter(count >= minlines) %>%
    pull(.data[[line_group]])

  # Filter the original DataFrame to include only rows in the selected groups
  filtered_df <- df_map %>%
    filter(.data[[line_group]] %in% tissue_list)

  return(filtered_df)
}

#' Selects cell lines based on tissue and mapping information.
#'
#' This function filters a data frame of cell line data (`df`) based on a list of tissues and a mapping data frame (`df_map`) containing metadata about the cell lines.
#'
#' @param df A data frame containing cell line data to be filtered.
#' @param df_map A mapping data frame that includes metadata for cell lines, such as tissue or lineage information.
#' @param tissue_list A character vector specifying the list of tissues to use for filtering.
#' @param line_group A string specifying the column in `df_map` representing the group or lineage for each cell line. Default is `"OncotreeLineage"`.
#' @param line_col A string specifying the column in `df_map` that contains the unique identifiers for cell lines. Default is `"ModelID"`.
#' @param nested Logical, if `TRUE`, the function assumes the cell lines are nested and applies additional grouping logic. Default is `FALSE`.
#' @param verbose Logical, if `TRUE`, the function prints additional details about the filtering process. Default is `FALSE`.
#'
#' @return A list of selected cell lines, either flattened or nested based on the `nested` parameter..
#'
#' @examples
#' # Example usage:
#' df <- data.frame(ModelID = c("A", "B", "C"), Value = c(1.2, 3.4, 5.6))
#' df_map <- data.frame(ModelID = c("A", "B", "C"), OncotreeLineage = c("Tissue1", "Tissue2", "Tissue3"))
#' tissue_list <- c("Tissue1", "Tissue3")
#' filtered_df <- select_cell_lines(df, df_map, tissue_list, verbose = TRUE)
#'
#' @export
select_cell_lines <- function(df, df_map, tissue_list,
                              line_group = "OncotreeLineage",
                              line_col = "ModelID",
                              nested = FALSE,
                              verbose = FALSE) {

  # Ensure tissue_list is valid
  if (is.character(tissue_list) && tissue_list != "all") {
    stop('tissue_list argument must be "all" or a list of strings!')
  }

  # Initialize the list of cell lines
  lines <- if (nested) list() else character()

  # Handle "all" case for tissue_list
  if (identical(tissue_list, "all")) {
    map_cell_lines <- df_map[[line_col]]
    dep_cell_lines <- intersect(colnames(df), map_cell_lines)
    tissue_list <- df_map %>%
      filter(.data[[line_col]] %in% dep_cell_lines) %>%
      pull(.data[[line_group]]) %>%
      unique()
  }

  # Iterate over each tissue in tissue_list
  for (tissue in tissue_list) {
    map_cell_lines <- df_map %>%
      filter(coalesce(.data[[line_group]] == tissue, FALSE)) %>%
      pull(.data[[line_col]])

    dep_cell_lines <- intersect(colnames(df), map_cell_lines)

    if (length(dep_cell_lines) == 0) {
      stop(sprintf('Empty list of lines ... the tissue "%s" may not be in the model.', tissue))
    }

    if (nested) {
      lines[[tissue]] <- dep_cell_lines
    } else {
      lines <- unique(c(lines, dep_cell_lines))
    }

    if (verbose) {
      message(sprintf(
        'There are %d "%s" cell-lines in the CRISPR file in common with the %d cell-lines in DepMap',
        length(dep_cell_lines), tissue, length(map_cell_lines)
      ))
    }
  }

  # Print summary if verbose
  if (verbose) {
    total_selected <- if (nested) sum(sapply(lines, length)) else length(lines)
    message(sprintf(
      'A total of %d cell lines have been selected for tissues: %s',
      total_selected, paste(tissue_list, collapse = ", ")
    ))
  }

  return(lines)
}

#' Filter rows in a DataFrame based on the percentage of NaN values.
#'
#' This function removes rows from a data frame where the percentage of `NaN` (or `NA`) values exceeds a specified threshold.
#'
#' @param df A data frame from which rows will be removed.
#' @param perc A numeric value specifying the maximum allowed percentage of `NaN` (or `NA`) values in a row. Rows exceeding this threshold will be removed. Default is `100.0`.
#' @param verbose Logical, if `TRUE`, the function will print the number of rows removed. Default is `FALSE`.
#'
#' @return A data frame with rows exceeding the specified percentage of `NaN` values removed.
#'
#' @examples
#' # Example usage:
#' df <- data.frame(a = c(1, NA, 3), b = c(NA, NA, 6), c = c(7, 8, 9))
#' df_clean <- delrows_with_nan_percentage(df, perc = 50.0, verbose = TRUE)
#'
#' @export
delrows_with_nan_percentage <- function(df, perc = 100.0, verbose = FALSE) {

  # Calculate the minimum count of non-missing values required
  mincount <- round(((100 - perc) / 100) * ncol(df)) + 1

  # Filter rows based on the condition
  df_filtered <- df[rowSums(is.na(df)) <= (ncol(df) - mincount), ]

  # Verbose logging
  if (verbose) {
    message(sprintf(
      "Removed %d rows from %d with at least %.2f%% NaN values",
      nrow(df) - nrow(df_filtered), nrow(df), perc
    ))
  }
  return(df_filtered)
}

