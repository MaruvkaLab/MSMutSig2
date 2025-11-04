# run_analysis/run_analysis.R

# Example Script to Run MSMuSig2 Analysis

# Define Input and Output Paths
input_file <- "../data/Loci.table.ALL_Only.2025.csv"
output_dir <- "../results/All_2025/"
selected_models <- "exponential,gaussian,weibull,lognormal" # Running all models (default)

# Create Output Directories if They Don't Exist
if (!dir.exists(output_dir)) {  dir.create(output_dir, recursive = TRUE) }

# Run the Analysis Script with Specified Parameters
system(paste("Rscript", 
             "../scripts/analysis.R", 
             "--input", input_file, 
             "--output", output_dir, 
             "--models", selected_models,
             "--threshold", 40))
