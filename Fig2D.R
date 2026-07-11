############################################################
# Figure 2D
# Odds ratio heatmap of B-cell distribution across
# treatment response groups
############################################################

library(plyr)
library(data.table)
library(grid)

#-----------------------------------------------------------
# Calculate odds ratios using Fisher's exact test
#-----------------------------------------------------------
calculate_odds_ratio <- function(count_matrix, min_row_sum = 0) {

  count_matrix <- count_matrix[rowSums(count_matrix) >= min_row_sum, , drop = FALSE]

  column_totals <- colSums(count_matrix)
  row_totals <- rowSums(count_matrix)

  count_table <- as.data.table(
    as.data.frame(count_matrix),
    keep.rownames = "celltype"
  )

  count_long <- melt(
    count_table,
    id.vars = "celltype",
    variable.name = "response_group",
    value.name = "count"
  )

  result_list <- lapply(seq_len(nrow(count_long)), function(i) {

    current_celltype <- count_long$celltype[i]
    current_group <- count_long$response_group[i]
    observed_count <- count_long$count[i]

    remaining_group_count <- column_totals[current_group] - observed_count

    contingency_table <- matrix(
      c(
        observed_count,
        row_totals[current_celltype] - observed_count,
        remaining_group_count,
        sum(column_totals) -
          row_totals[current_celltype] -
          remaining_group_count
      ),
      ncol = 2
    )

    fisher_result <- fisher.test(contingency_table)

    data.table(
      celltype = current_celltype,
      response_group = current_group,
      p_value = fisher_result$p.value,
      odds_ratio = unname(fisher_result$estimate)
    )
  })

  odds_ratio_results <- rbindlist(result_list)

  odds_ratio_results <- merge(
    count_long,
    odds_ratio_results,
    by = c("celltype", "response_group")
  )

  odds_ratio_results <- as.data.table(odds_ratio_results)

  odds_ratio_results[
    ,
    adjusted_p_value := p.adjust(p_value, method = "BH")
  ]

  return(odds_ratio_results)
}

#-----------------------------------------------------------
# Plot odds ratio heatmap
#-----------------------------------------------------------
plot_odds_ratio_heatmap <- function(
    metadata,
    celltype,
    group,
    output_prefix,
    figure_title,
    pdf_width = 5,
    pdf_height = 5,
    return_results = FALSE
){

  dir.create(dirname(output_prefix), recursive = TRUE, showWarnings = FALSE)

  metadata <- data.table(metadata)
  metadata$celltype <- as.character(celltype)

  if (is.factor(group)) {
    metadata$response_group <- group
  } else {
    metadata$response_group <- factor(group)
  }

  group_levels <- levels(metadata$response_group)

  count_matrix <- unclass(
    metadata[, table(celltype, response_group)]
  )[, group_levels]

  proportion_matrix <- sweep(
    count_matrix,
    1,
    rowSums(count_matrix),
    "/"
  )

  print(floor(proportion_matrix * 10))

  odds_ratio_results <- calculate_odds_ratio(count_matrix)

  pvalue_table <- dcast(
    odds_ratio_results,
    celltype ~ response_group,
    value.var = "p_value"
  )

  odds_ratio_table <- dcast(
    odds_ratio_results,
    celltype ~ response_group,
    value.var = "odds_ratio"
  )

  odds_ratio_matrix <- as.matrix(odds_ratio_table[, -1])

  rownames(odds_ratio_matrix) <- odds_ratio_table[[1]]

  sscVis::plotMatrix.simple(
    odds_ratio_matrix,
    out.prefix = sprintf("%s_OR", output_prefix),
    show.number = TRUE,
    exp.name = expression(italic(OR)),
    z.hi = 2,
    palatte = viridis::viridis(7),
    pdf.width = pdf_width,
    pdf.height = pdf_height,
    mytitle = figure_title,
    my.cell_fun = function(j, i, x, y, width, height, fill) {
      grid.text(round(odds_ratio_matrix[i, j], 2), x, y)
    }
  )

  if (return_results) {
    return(
      list(
        odds_ratio_results = odds_ratio_results,
        pvalue_table = pvalue_table,
        odds_ratio_table = odds_ratio_table,
        odds_ratio_matrix = odds_ratio_matrix
      )
    )
  }

  return(odds_ratio_matrix)
}

#-----------------------------------------------------------
# Load metadata
#-----------------------------------------------------------

metadata <- seu_bcell@meta.data

celltype_order <- c(
  "B1_TCL1A_naiveB",
  "B2_SPP1_B",
  "B3_EGR1_AcB",
  "B4_NR4A3_AcB",
  "B5_non_SwBm",
  "B6_SwBm",
  "B7_ITGAX_AtM",
  "B8_Plasma_cell"
)

metadata$celltype_fine <- factor(
  metadata$celltype_fine,
  levels = celltype_order
)

#-----------------------------------------------------------
# Generate Figure 2D
#-----------------------------------------------------------

odds_ratio_matrix <- plot_odds_ratio_heatmap(
  metadata = metadata,
  celltype = metadata$celltype_fine,
  group = metadata$response_group,
  output_prefix = "Figure2D",
  figure_title = "Odds Ratio",
  pdf_width = 5,
  pdf_height = 5
)
