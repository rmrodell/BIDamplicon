# ------------------------------------------------------------------------------
#
# Script: plot.R
#
# Description:
#   Plots in cellulo Nano-BID-Amp data.
#   The script is designed for interactive use.
#
# Dependencies:
#   - R version >= 4.3.0
#   - R packages: dplyr, readr, car, tidyr, purrr, ggplot2,
#                 furrr, forcats
#
# Usage:
#   1. Set all parameters in the "USER CONFIGURATION" section.
#   2. Run the entire script in RStudio / other IDE.
#
# ------------------------------------------------------------------------------

# ---- USER CONFIGURATION ----

# --- Global Settings ---
input_file      <- "/scratch/users/rodell/20260305_endoBID/counts_delpos/BIDdetect_data.txt"
output_dir      <- "/scratch/users/rodell/20260305_endoBID/counts_delpos"

# --- Analysis-Specific Parameters & Palettes ---

# Conditions to iterate over. "Both" pools data from all specified cell types.
analysis_conditions <- c("Both")

pus7_config <- list(
  factor_col = "PUS7level",
  comparisons = list(
    WT_v_KO = c("WT", "PUS7KO"),
    # WT_v_KD = c("P102", "P101"),
    # WT_v_KD = c("P36", "P97"),
    # WT_v_KD = c("pLKO", "shPUS7"),
    # WT_v_OE = c("P4", "P3")
  ),
  colors = list(input = "#eee8aa", KO = "#f3ccdf", KD = "#ed93c0", WT = "#e25098", OE = "#a21b5d")
)

celltype_config <- list(
  factor_col = "celltype",
  levels = c("HepG2", "293T"),
  wt_vectors = c("WT", "P102", "P4", "P36", "pLKO"),
  colors = list(input = "#eee8aa", "293T" = "#32cd32", "HepG2" = "#009acd")
)

# ---- Environment & Library Setup ----
cat("--- Setting up Environment ---\n")

# 1. Check R Version
if (getRversion() < "4.3.0") {
  stop("ERROR: R version 4.3.0 or higher is required.")
}

# 2. Install and Load Libraries
required_packages <- c("dplyr", "readr", "tidyr", "purrr", "ggplot2", "forcats")
missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
if (length(missing_packages) > 0) {
  cat("Installing missing packages:", paste(missing_packages, collapse = ", "), "\n")
  install.packages(missing_packages, repos = "http://cran.us.r-project.org")
}
suppressPackageStartupMessages({lapply(required_packages, library, character.only = TRUE)})
cat("✓ Libraries loaded successfully.\n\n")


# ---- Directory and File Setup ----
plot_dir <- file.path(output_dir, "plots")
dir.create(plot_dir, showWarnings=FALSE, recursive=TRUE)

# --- Plotting Functions ---

# Function: Generate individual bar plots for all sites
generate_allsites_plot <- function(plot_data, num_facets, plot_title, filename, fill_palette, color_by = "replicate", color_palette_jitter = NULL) {
  # Generates and saves the complex, multi-panel bar plot for a given dataset.
  
  cat("   Generating allsites plot:", filename, "\n")
  
  # 1. Create summary for bars and error bars
  summary_for_bars <- plot_data %>%
    group_by(facet_label, group) %>%
    summarise(
      mean_delrate = mean(delrate, na.rm = TRUE),
      sd_delrate = sd(delrate, na.rm = TRUE),
      .groups = "drop"
    )
  
  # 2. Build the base plot
  p <- ggplot() +
    geom_bar(data = summary_for_bars, aes(x = group, y = mean_delrate, fill = group), stat = "identity") +
    geom_errorbar(data = summary_for_bars, aes(x = group, ymin = pmax(0, mean_delrate - sd_delrate), 
                                               ymax = mean_delrate + sd_delrate), width = 0.3)
  
  # 3. Add jitter points with conditional coloring
  if (color_by == "celltype") {
    p <- p + 
      geom_jitter(data = plot_data, aes(x = group, y = delrate, color = celltype), width = 0.2, size = 2, alpha = 0.8) +
      scale_color_manual(name = "Cell Type", values = color_palette_jitter)
  } else { # Default to coloring by replicate
    p <- p + 
      geom_jitter(data = plot_data, aes(x = group, y = delrate, color = factor(rep)), width = 0.2, size = 2, alpha = 0.8) +
      scale_color_brewer(name = "Replicate", palette = "Set2")
  }
  
  # 4. Add remaining layers and save
  p <- p +
    facet_wrap(~ fct_inorder(facet_label), ncol = 8) +
    labs(title = plot_title, x = "", y = "Deletion Fraction") +
    theme_boxplot(base_size = 10) +
    theme(
      strip.background = element_rect(fill = "gray90"),
      strip.text = element_text(size = rel(0.8))
    ) +
    scale_fill_manual(name = "Group", values = fill_palette) +
    coord_cartesian(ylim = c(0, NA))
  
  plot_height <- max(3, ceiling(n_distinct(plot_data$facet_label) / 8) * 3)
  plot_width <- (num_facets * 1.5)
  ggsave(filename, p, width = plot_width, height = plot_height, units = "in", limitsize = FALSE)
  cat("      ... plot saved.\n")
}


# ---- Plot Themes ----
theme_base_custom <- function(base_size=16) { 
  theme_minimal(base_size=base_size, base_family="sans") %+replace% 
    theme(plot.title=element_text(hjust=0.5, size=rel(1.2)), 
          plot.subtitle=element_text(hjust=0.5, size=rel(0.8)), 
          axis.title=element_text(size=rel(1.0)), 
          axis.text=element_text(size=rel(1.0)), 
          plot.margin=margin(5,5,5,5)) }

theme_boxplot <- function(base_size=16) { 
  theme_base_custom(base_size) + 
    theme(legend.position="none", 
          panel.grid.major.x=element_blank(), 
          panel.grid.minor=element_blank()) }

# ---- MAIN WORKFLOW ----
cat("--- Starting Plotting ---\n")

# ---- Global Setup ----
if (!file.exists(input_file)) stop(paste("Error: Input file not found at", input_file))
main_data <- readr::read_tsv(input_file, show_col_types = FALSE) %>%
  mutate(treat = case_when(
    treat %in% c("in", "input") ~ "input",
    treat %in% c("BID", "BS") ~ "BS",
    TRUE ~ treat
  ))

# ---- 2. PUS7 Dependency Analysis ----
  
# Define vectors for grouping once
KO_vectors <- pus7_config$comparisons$WT_v_KO[[2]]
KD_vectors <- pus7_config$comparisons$WT_v_KD[[2]]
WT_vectors <- unique(c(pus7_config$comparisons$WT_v_KO[[1]], pus7_config$comparisons$WT_v_OE[[1]], pus7_config$comparisons$WT_v_KD[[1]]))
OE_vectors <- pus7_config$comparisons$WT_v_OE[[2]]

# --- Plot All Sites ---

cat("\n\n--- Generating All Sites Plots for PUS7 Expression ---\n")

for (condition in analysis_conditions) {
  
  # 1. Filter data and prepare for plotting
  if (condition == "Both") {
    plot_data <- main_data
  } else {
    plot_data <- main_data %>% filter(celltype == condition)
  }
  
  plot_data_prepared <- plot_data %>%
    mutate(group = case_when(
      treat == "input" ~ "input",
      .data[[pus7_config$factor_col]] %in% KO_vectors ~ "KO",
      .data[[pus7_config$factor_col]] %in% KD_vectors ~ "KD",
      .data[[pus7_config$factor_col]] %in% WT_vectors ~ "WT",
      .data[[pus7_config$factor_col]] %in% OE_vectors ~ "OE",
      TRUE ~ NA_character_
    )) %>%
    filter(!is.na(group)) %>%
    mutate(group = factor(group, levels = c("input", "KO", "KD", "WT", "OE")),
            facet_label = paste0(chr, ":", pos))
  
  num_sites <- n_distinct(plot_data_prepared$facet_label)

  # 2. Generate and save the plot using the helper function
  generate_allsites_plot(
    plot_data = plot_data_prepared,
    num_facets = num_sites,
    plot_title = paste("All Sites -", condition, "Data"),
    filename = file.path(plot_dir, paste0(condition, "_allsites.pdf")),
    fill_palette = pus7_config$colors,
    color_by = if (condition == "Both") "celltype" else "replicate",
    color_palette_jitter = if (condition == "Both") celltype_config$colors else NULL
  )

  if (num_sites < 20) {
    generate_allsites_plot(
      plot_data = plot_data_prepared,
      num_facets = num_sites,
      plot_title = paste("All Sites -", condition, "Data"),
      filename = file.path(plot_dir, paste0(condition, "_allsites.png")),
      fill_palette = pus7_config$colors,
      color_by = if (condition == "Both") "celltype" else "replicate",
      color_palette_jitter = if (condition == "Both") celltype_config$colors else NULL
    )
  }
}
