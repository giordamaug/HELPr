# HELPr utility functions
#

library(dplyr)
filter_crispr_by_model <- function(df, df_map,
                                   minlines = 1,
                                   line_colname = "ModelID",
                                   line_group = "OncotreeLineage") {
  # Filter a CRISPR DataFrame based on a mapping DataFrame and specified
  # conditions.
  #
  # Parameters:
  # df: The CRISPR DataFrame to be filtered.
  # df_map: The mapping DataFrame containing information about cell lines and models.
  # minlines: The minimum number of lines required for a tissue in the model. Default is 1.
  # line_colname: The column name in both DataFrames representing the cell line ID. Default is "ModelID".
  # line_group: The column name in the mapping DataFrame representing the tissue/lineage group. Default is "OncotreeLineage".
  #
  # Returns:
  # A new DataFrame with CRISPR data filtered based on the selected cell lines and conditions.


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

select_cell_lines <- function(df, df_map, tissue_list,
                              line_group = "OncotreeLineage",
                              line_col = "ModelID",
                              nested = FALSE,
                              verbose = FALSE) {
  # Selects cell lines based on tissue and mapping information.
  #
  # Arguments:
  # df: DataFrame containing cell line information.
  # df_map: DataFrame containing mapping information.
  # tissue_list: List of tissues for which cell lines need to be selected, or "all".
  # line_group: Column in `df_map` to use for line selection (default is "OncotreeLineage").
  # line_col: Column in `df_map` to use for tissue selection (default is "ModelID").
  # nested: Whether to return cell lines as nested lists (lists for each tissue).
  # verbose: Verbosity level for printing information.
  #
  # Returns:
  # A list of selected cell lines, either flattened or nested based on the `nested` parameter.

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

delrows_with_nan_percentage <- function(df, perc = 100.0, verbose = FALSE) {
  # Filter rows in a DataFrame based on the percentage of NaN values.
  #
  # Arguments:
  # df: The input DataFrame.
  # perc: The percentage of NaN values allowed in each row. Default is 100.0.
  # verbose: Whether to print information about the filtering process.
  #
  # Returns:
  # A new DataFrame with rows filtered based on the specified percentage of NaN values.

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

