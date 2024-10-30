# scripts/analysis.R

# ============================================
# MSMuSig2 - Comprehensive Driver Mutation Analysis
# ============================================

# Set CRAN mirror to avoid interactive prompts during package installation
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Load Necessary Packages
required_packages <- c("magrittr", "minpack.lm", "ggplot2", "gridExtra", "MASS", 
                       "dplyr", "data.table", "fitdistrplus", "survival", 
                       "car", "ggrepel", "ggpubr", "optparse")
installed_packages <- rownames(installed.packages())
for (pkg in required_packages) {
  if (!pkg %in% installed_packages) {
    install.packages(pkg, dependencies = TRUE)
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

# Define Command-Line Options
option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="Path to input CSV file", metavar="character"),
  make_option(c("-o", "--output"), type="character", default="results", 
              help="Path to output directory [default= %default]", metavar="character"),
  make_option(c("-m", "--models"), type="character", default="lognormal,exponential", 
              help="Comma-separated list of models to apply (options: exponential, gaussian, weibull, lognormal) [default= %default]", metavar="character"),
  make_option(c("-t", "--threshold"), type="integer", default=45, 
              help="Minimum occurrence count to select motifs [default= %default]", metavar="integer")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Validate Input Arguments
if (is.null(opt$input)){
  print_help(opt_parser)
  stop("Input file must be supplied (--input).", call.=FALSE)
}

# Parse Selected Models
selected_models <- unlist(strsplit(opt$models, split = ","))
selected_models <- trimws(selected_models)  # Remove any extra whitespace

# Validate Selected Models
valid_models <- c("exponential", "gaussian", "weibull", "lognormal")
if (!all(selected_models %in% valid_models)) {
  stop(paste("Invalid model(s) selected. Choose from:", paste(valid_models, collapse = ", ")))
}

# Retrieve Threshold Value
threshold_val <- opt$threshold

# Create Output Directory (Single Folder)
output_dir <- opt$output
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Define Paths for Outputs
qqplot_dir <- output_dir
qqplot_paths <- paste0(qqplot_dir, "/QQplot_", selected_models, ".svg")

# Read Mutation Data
mutation_data <- read.csv(opt$input)

# Ensure Required Columns Exist
required_columns <- c("PATTERN", "REFERENCE_REPEATS", "COUNTS", "hgnc_symbol")
missing_cols <- setdiff(required_columns, names(mutation_data))
if (length(missing_cols) > 0) {
  stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
}

# Subset Unique Motifs
MOTIFS <- mutation_data[, c("PATTERN", "REFERENCE_REPEATS")] %>% .[!duplicated(.), ]

# Fix: Order MOTIFS and Calculate Occurrence
MOTIFS <- MOTIFS[order(MOTIFS$PATTERN, MOTIFS$REFERENCE_REPEATS), ]
MOTIFS$occurence <- as.vector(table(paste0(mutation_data$PATTERN, "-", mutation_data$REFERENCE_REPEATS)))

# Select Motifs Based on Occurrence Threshold
selected_motifs_df <- MOTIFS %>% filter(occurence > threshold_val)
selected_motifs <- paste0(selected_motifs_df$PATTERN, "-", selected_motifs_df$REFERENCE_REPEATS)

# Inform User About Selected Motifs
cat(paste0("Selected Motifs with occurrence > ", threshold_val, ": ", 
           length(selected_motifs), "\n"))

# Initialize Lists to Store Results
rez.allMotifs <- list()
all_p_values_per_model <- list()
all_hgnc_symbols_per_model <- list()
all_original_p_values_per_model <- list()

# Initialize Lists for Each Model
for (model in selected_models) {
  all_p_values_per_model[[model]] <- c()
  all_original_p_values_per_model[[model]] <- c()
  all_hgnc_symbols_per_model[[model]] <- c()
}

# Define the plot_pvalues Function
plot_pvalues <- function(pvalues, significant, gene_names, model_name, output_path) {
  # Ensure pvalues and gene_names have the same length
  if (length(pvalues) != length(gene_names)) {
    warning(paste("Mismatch in lengths: pvalues =", length(pvalues), 
                  "gene_names =", length(gene_names), ". Truncating to the minimum length."))
    min_len <- min(length(pvalues), length(gene_names))
    pvalues <- pvalues[1:min_len]
    gene_names <- gene_names[1:min_len]
    significant <- significant[1:min_len]
  }
  
  # Create a data frame and sort by p-values
  df <- data.frame(
    pvalue = pvalues,
    gene = gene_names,
    significant = significant
  ) %>%
    dplyr::arrange(pvalue) %>%
    dplyr::mutate(
      expected = (rank(pvalue) - 0.5) / length(pvalue)
    )
  
  # Determine the maximum value for both axes
  max_val <- max(-log10(df$expected), -log10(df$pvalue))
  
  # Generate the QQ plot
  qq_plot <- ggplot(df, aes(x = -log10(expected), y = -log10(pvalue))) +
    geom_point(aes(color = significant), size = 1) +
    scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    coord_fixed(ratio = 1, xlim = c(0, max_val), ylim = c(0, max_val)) +
    geom_text_repel(
      data = subset(df, significant),
      aes(label = gene),
      vjust = -1, 
      size = 3
    ) +
    theme_minimal() +
    labs(
      title = paste0("QQ Plot of Adjusted P-Values - ", toupper(substr(model_name, 1,1)), substr(model_name, 2, nchar(model_name))),
      x = "Expected -log10(p)",
      y = "Observed -log10(p)",
      color = "Significant"
    )
  
  # Save the QQ Plot
  ggsave(filename = output_path, plot = qq_plot, width = 8, height = 6)
}

# Function to Fit Models and Return Fit Objects
fit_models <- function(x_data, models) {
  fits <- list()
  for (model in models) {
    if (model == "exponential") {
      fits$exponential <- tryCatch({
        fitdistr(x_data, "exponential")
      }, error = function(e) { NULL })
    }
    if (model == "gaussian") {
      fits$gaussian <- tryCatch({
        fitdistr(x_data, "normal")
      }, error = function(e) { NULL })
    }
    if (model == "weibull") {
      fits$weibull <- tryCatch({
        fitdistr(x_data, "weibull")
      }, error = function(e) { NULL })
    }
    if (model == "lognormal") {
      fits$lognormal <- tryCatch({
        fitdistr(x_data, "lognormal")
      }, error = function(e) { NULL })
    }
  }
  return(fits)
}

# Loop Through Each Selected Motif and Perform Analysis
for (MOTIF in selected_motifs) {
  df <- subset(mutation_data, paste0(PATTERN, "-", REFERENCE_REPEATS) == MOTIF)
  x_data <- df$COUNTS
  x_data <- x_data[x_data > 0]
  
  # Fit Selected Models to Raw Data
  model_fits <- fit_models(x_data, selected_models)
  
  # Check if at least one model was fitted successfully
  if (all(sapply(model_fits, is.null))) {
    warning(paste("No valid models fitted for motif:", MOTIF))
    next
  }
  
  # Fit Models to Histogram Counts (RAW only, as per requirement)
  # Create Histogram Data
  h <- hist(x_data, breaks = 50, plot = FALSE)
  x_mids <- h$mids
  y_counts <- h$counts
  hist_data <- data.frame(x = x_mids, y = y_counts)
  
  # Define Starting Values for Models
  start_values <- list(
    exponential = list(rate = 1/mean(x_data)),
    gaussian = list(mean = mean(x_data), sd = sd(x_data)),
    weibull = list(shape = 1, scale = mean(x_data)),
    lognormal = list(meanlog = mean(log(x_data)), sdlog = sd(log(x_data)))
  )
  
  # Fit Models to Histogram Counts
  hist_fits <- list()
  for (model in selected_models) {
    if (model == "exponential") {
      hist_fits$exponential <- tryCatch({
        nlsLM(y ~ a * exp(-rate * x), data = hist_data, 
              start = list(a = max(y_counts), rate = start_values$exponential$rate))
      }, error = function(e) { NULL })
    }
    if (model == "gaussian") {
      hist_fits$gaussian <- tryCatch({
        nlsLM(y ~ a * exp(-((x - mean)^2) / (2 * sd^2)), data = hist_data, 
              start = list(a = max(y_counts), mean = start_values$gaussian$mean, sd = start_values$gaussian$sd))
      }, error = function(e) { NULL })
    }
    if (model == "weibull") {
      hist_fits$weibull <- tryCatch({
        nlsLM(y ~ a * (shape / scale) * (x / scale)^(shape - 1) * exp(- (x / scale)^shape), 
              data = hist_data, 
              start = list(a = max(y_counts), shape = start_values$weibull$shape, scale = start_values$weibull$scale))
      }, error = function(e) { NULL })
    }
    if (model == "lognormal") {
      hist_fits$lognormal <- tryCatch({
        nlsLM(y ~ a * (1 / (x * sigma * sqrt(2 * pi))) * 
                exp(- (log(x) - mu)^2 / (2 * sigma^2)), 
              data = hist_data, 
              start = list(a = max(y_counts), mu = start_values$lognormal$meanlog, sigma = start_values$lognormal$sdlog))
      }, error = function(e) { NULL })
    }
  }
  
  # Visualization of Fits Over Histogram with Consistent Axes
  p1 <- local({
    temp_x_data <- x_data
    temp_plot_data <- data.frame(
      x = seq(min(temp_x_data), max(temp_x_data), length.out = 1000)
    )
    
    for (model in selected_models) {
      if (!is.null(model_fits[[model]])) {
        if (model == "exponential") {
          temp_plot_data$Exponential <- dexp(temp_plot_data$x, rate = model_fits$exponential$estimate["rate"])
        }
        if (model == "gaussian") {
          temp_plot_data$Gaussian <- dnorm(temp_plot_data$x, 
                                          mean = model_fits$gaussian$estimate["mean"], 
                                          sd = model_fits$gaussian$estimate["sd"])
        }
        if (model == "weibull") {
          temp_plot_data$Weibull <- dweibull(temp_plot_data$x, 
                                            shape = model_fits$weibull$estimate["shape"], 
                                            scale = model_fits$weibull$estimate["scale"])
        }
        if (model == "lognormal") {
          temp_plot_data$LogNormal <- dlnorm(temp_plot_data$x, 
                                            meanlog = model_fits$lognormal$estimate["meanlog"], 
                                            sdlog = model_fits$lognormal$estimate["sdlog"])
        }
      }
    }
    
    # Calculate scaling factor to match histogram density
    hist_max <- max(ggplot_build(ggplot() + 
                                   geom_histogram(aes(x = temp_x_data, y = after_stat(density)), 
                                                 bins = ceiling(length(temp_x_data) / 3)) + 
                                   theme_minimal())$data[[1]]$density)
    density_values <- sapply(selected_models, function(model) {
      if (model == "exponential") return(max(temp_plot_data$Exponential, na.rm = TRUE))
      if (model == "gaussian") return(max(temp_plot_data$Gaussian, na.rm = TRUE))
      if (model == "weibull") return(max(temp_plot_data$Weibull, na.rm = TRUE))
      if (model == "lognormal") return(max(temp_plot_data$LogNormal, na.rm = TRUE))
    })
    density_max <- max(density_values, na.rm = TRUE)
    scaling_factor <- density_max / hist_max
    
    # Generate Plot
    plot <- ggplot() +
      geom_histogram(aes(x = temp_x_data, y = after_stat(density) * scaling_factor), 
                     bins = ceiling(length(temp_x_data) / 3), 
                     fill = "grey", color = "black", alpha = 0.6) +
      labs(title = MOTIF,
           x = paste0("COUNTS (n=", length(temp_x_data), ")"), y = "Density")
    
    # Add model lines
    color_mapping <- c()
    for (model in selected_models) {
      if (model == "exponential" && !is.null(temp_plot_data$Exponential)) {
        plot <- plot + 
          geom_line(data = temp_plot_data, 
                    aes(x = x, y = Exponential, color = "Exponential"), 
                    linewidth = 1)
        color_mapping <- c(color_mapping, "Exponential" = "#004D40")
      }
      if (model == "gaussian" && !is.null(temp_plot_data$Gaussian)) {
        plot <- plot + 
          geom_line(data = temp_plot_data, 
                    aes(x = x, y = Gaussian, color = "Gaussian"), 
                    linewidth = 1)
        color_mapping <- c(color_mapping, "Gaussian" = "#FFC107")
      }
      if (model == "weibull" && !is.null(temp_plot_data$Weibull)) {
        plot <- plot + 
          geom_line(data = temp_plot_data, 
                    aes(x = x, y = Weibull, color = "Weibull"), 
                    linewidth = 1)
        color_mapping <- c(color_mapping, "Weibull" = "#1E88E5")
      }
      if (model == "lognormal" && !is.null(temp_plot_data$LogNormal)) {
        plot <- plot + 
          geom_line(data = temp_plot_data, 
                    aes(x = x, y = LogNormal, color = "LogNormal"), 
                    linewidth = 1)
        color_mapping <- c(color_mapping, "LogNormal" = "#D81B60")
      }
    }
    
    plot <- plot +
      scale_color_manual(name = "Models", 
                         values = color_mapping) +
      theme_minimal() +
      theme(legend.position = "bottom")
    
    return(plot)
  })
  
  # Store Comparison Metrics for RAW Models
  comparison_raw <- data.frame(
    Model = character(),
    AIC = numeric(),
    BIC = numeric(),
    LogLikelihood = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (model in selected_models) {
    fit <- model_fits[[model]]
    if (!is.null(fit)) {
      logLik <- fit$loglik
      k <- length(fit$estimate)  # Number of parameters
      n <- length(x_data)        # Number of observations
      
      AIC_value <- -2 * logLik + 2 * k
      BIC_value <- -2 * logLik + k * log(n)
      
      comparison_raw <- rbind(comparison_raw, data.frame(
        Model = model,
        AIC = AIC_value,
        BIC = BIC_value,
        LogLikelihood = logLik
      ))
    }
  }
  
  # Best Model Based on AIC
  best_model_raw <- comparison_raw$Model[which.min(comparison_raw$AIC)]
  
  # Store Results
  rez.allMotifs[[MOTIF]] <- list(
    best_model_raw = best_model_raw,
    comparison_raw = comparison_raw,
    plot = p1
  )
  
  # Collect P-Values and HGNCSymbols for Each Model
  for (model in selected_models) {
    if (!is.null(model_fits[[model]])) {
      if (model %in% c("weibull", "lognormal", "exponential", "gaussian")) {
        # Calculate p-values based on the fitted model
        if (model == "exponential") {
          p_vals <- pexp(x_data, rate = model_fits$exponential$estimate["rate"], lower.tail = FALSE)
        }
        if (model == "gaussian") {
          p_vals <- pnorm(x_data, mean = model_fits$gaussian$estimate["mean"], 
                         sd = model_fits$gaussian$estimate["sd"], lower.tail = FALSE)
        }
        if (model == "weibull") {
          p_vals <- pweibull(x_data, shape = model_fits$weibull$estimate["shape"], 
                             scale = model_fits$weibull$estimate["scale"], lower.tail = FALSE)
        }
        if (model == "lognormal") {
          p_vals <- plnorm(x_data, meanlog = model_fits$lognormal$estimate["meanlog"], 
                          sdlog = model_fits$lognormal$estimate["sdlog"], lower.tail = FALSE)
        }
        # Remove NAs and Infinities
        valid_indices <- !is.na(p_vals) & is.finite(p_vals)
        p_vals <- p_vals[valid_indices]
        gene_names <- df$hgnc_symbol[x_data > 0][valid_indices]
        
        # Ensure alignment before storing
        if (length(p_vals) != length(gene_names)) {
          warning(paste("Mismatch in lengths for model", model, 
                        ": p_vals =", length(p_vals), 
                        "gene_names =", length(gene_names), ". Truncating to the minimum length."))
          min_len <- min(length(p_vals), length(gene_names))
          p_vals <- p_vals[1:min_len]
          gene_names <- gene_names[1:min_len]
        }
        
        # Store p-values and hgnc_symbols per model
        all_p_values_per_model[[model]] <- c(all_p_values_per_model[[model]], p_vals)
        all_original_p_values_per_model[[model]] <- c(all_original_p_values_per_model[[model]], p_vals)
        all_hgnc_symbols_per_model[[model]] <- c(all_hgnc_symbols_per_model[[model]], gene_names)
      }
    }
  }
}

# Adjust P-Values Globally for Each Model and Generate QQ Plots
for (model in selected_models) {
  p_vals <- all_p_values_per_model[[model]]
  hgnc_symbols <- all_hgnc_symbols_per_model[[model]]
  original_p_vals <- all_original_p_values_per_model[[model]]
  
  if (length(p_vals) > 0) {
    adjusted_p_values <- p.adjust(p_vals, method = "BH")
    significant <- adjusted_p_values < 0.1
    
    # Generate QQ Plot using the provided plot_pvalues function
    plot_pvalues(
      pvalues = p_vals,
      significant = significant,
      gene_names = hgnc_symbols,
      model_name = model,
      output_path = paste0(qqplot_dir, "/QQplot_", model, ".svg")
    )
    
    cat(paste0("QQ plot for model ", model, " saved to QQplot_", model, ".svg\n"))
    
    # Ensure p_vals and hgnc_symbols have the same length
    if (length(p_vals) != length(hgnc_symbols)) {
      warning(paste("Mismatch in lengths when creating all_results_df for model", model, 
                    ": p_vals =", length(p_vals), 
                    "hgnc_symbols =", length(hgnc_symbols), ". Truncating to the minimum length."))
      min_len <- min(length(p_vals), length(hgnc_symbols))
      p_vals <- p_vals[1:min_len]
      hgnc_symbols <- hgnc_symbols[1:min_len]
      adjusted_p_values <- adjusted_p_values[1:min_len]
    }
    
    # Create a Data Frame for All Results
    all_results_df <- data.frame(
      hgnc_symbol = hgnc_symbols,
      p_value = p_vals,
      adjusted_p_value = adjusted_p_values
    )
    
    # Save All Results to <model>.csv
    write.csv(all_results_df, paste0(output_dir, "/", model, ".csv"), row.names = FALSE)
    cat(paste0("All results for model ", model, " saved to ", model, ".csv\n"))
  } else {
    cat(paste0("No p-values available for model: ", model, ". QQ plot not generated.\n"))
  }
}

# Save Combined Model Fit Plots (if any motifs selected)
if (length(selected_motifs) > 0) {
  plots_to_combine <- lapply(selected_motifs, function(motif) {
    if (!is.null(rez.allMotifs[[motif]])) {
      return(rez.allMotifs[[motif]]$plot)
    } else {
      warning(paste("Motif", motif, "not found in the results."))
      return(NULL)
    }
  })
  # Remove NULLs
  plots_to_combine <- plots_to_combine[!sapply(plots_to_combine, is.null)]
  
  if (length(plots_to_combine) > 0) {
    combined_plot <- ggpubr::ggarrange(plotlist = plots_to_combine, 
                                       common.legend = TRUE, 
                                       legend = "bottom")
    ggsave(filename = paste0(output_dir, "/Indels_Models_fit.svg"), 
           plot = combined_plot, width = 15, height = 10)
    cat("Combined model fit plot saved to Indels_Models_fit.svg\n")
  } else {
    cat("No plots to combine for selected motifs.\n")
  }
} else {
  cat(paste0("No motifs selected with occurrence > ", threshold_val, ". Skipping combined plotting.\n"))
}

# Merge Comparison Tables for RAW Models
if (length(rez.allMotifs) > 0) {
  comparison_list <- lapply(names(rez.allMotifs), function(motif) { 
    df <- rez.allMotifs[[motif]]$comparison_raw
    df$Motif <- motif
    return(df)
  })
  
  Merged_t <- do.call(rbind, comparison_list)
  Merged_t <- Merged_t[order(Merged_t$Motif, Merged_t$AIC), ]
  
  write.csv(Merged_t, paste0(output_dir, "/Indels_Models_fit.csv"), row.names = FALSE)
  cat("Comparison table saved to Indels_Models_fit.csv\n")
} else {
  cat("No comparison tables to merge.\n")
}

# Output Summary of Significant Mutations for Each Model
for (model in selected_models) {
  p_vals <- all_original_p_values_per_model[[model]]
  adjusted_p_vals <- p.adjust(p_vals, method = "BH")
  hgnc_symbols <- all_hgnc_symbols_per_model[[model]]
  
  if (length(p_vals) > 0) {
    # Ensure p_vals and hgnc_symbols have the same length
    if (length(p_vals) != length(hgnc_symbols)) {
      warning(paste("Mismatch in lengths for model", model, 
                    ": p_vals =", length(p_vals), 
                    "hgnc_symbols =", length(hgnc_symbols), ". Truncating to the minimum length."))
      min_len <- min(length(p_vals), length(hgnc_symbols))
      p_vals <- p_vals[1:min_len]
      hgnc_symbols <- hgnc_symbols[1:min_len]
      adjusted_p_vals <- adjusted_p_vals[1:min_len]
    }
    
    significant_df <- data.frame(
      hgnc_symbol = hgnc_symbols,
      p_value = p_vals,
      adjusted_p_value = adjusted_p_vals
    )
    significant_df <- significant_df[significant_df$adjusted_p_value < 0.1, ]
    
    if (nrow(significant_df) > 0) {
      write.csv(significant_df, paste0(output_dir, "/Significant_", model, ".csv"), row.names = FALSE)
      cat(paste0("Significant mutations for model ", model, " saved to Significant_", model, ".csv\n"))
    } else {
      cat(paste0("No significant mutations detected for model ", model, ".\n"))
    }
  } else {
    cat(paste0("No p-values available for model: ", model, ". Significant mutations file not created.\n"))
  }
}

# Print Final Summary
cat("====================================\n")
cat("Analysis Completed Successfully.\n")
cat("====================================\n")
cat("Results saved in:", output_dir, "\n")
