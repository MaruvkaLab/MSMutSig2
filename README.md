# MSMuSig2

**MSMuSig2** is a tool designed to identify driver genes in Microsatellite-Unstable (MSI) cancers by analyzing recurrent MS-indel mutation counts. It employs various statistical models to determine significant driver mutations.

## Table of Contents

- Features
- Installation
- Usage
- Understanding the Output
- Repository Structure
- Methods
- Contributing
- License
- Authors
- Contact

## Features

- Fits multiple statistical models (**Exponential**, **Gaussian**, **Weibull**, **Log-Normal**) to mutation count data.
- Allows users to select which models to apply during analysis.
- Compares models using **AIC** and **BIC** to determine the best fit.
- Identifies significant driver mutations based on adjusted p-values across all models.
- Outputs results with detailed information including gene symbols, motifs, and genomic locations (**CHR**, **START**).
- Aggregates results to provide gene-level significance by combining p-values from multiple motifs per gene.
- Generates comprehensive visualizations, including model fits and QQ plots per model.

## Installation

1. **Clone the Repository**

   ```
   git clone https://github.com/yourusername/MSMuSig2.git
   cd MSMuSig2
   ```

2. **Install Required R Packages**

   Ensure you have the necessary R packages installed. You can install them using the following command in R:

   ```
   install.packages(c("magrittr", "minpack.lm", "ggplot2", "gridExtra", "MASS",
                      "dplyr", "data.table", "fitdistrplus", "survival", "car",
                      "ggrepel", "ggpubr", "optparse", "metap"))
   ```

   Alternatively, you can install all required packages at once using the `requirements.txt`:

   ```
   install.packages(readLines("requirements.txt"), dependencies = TRUE)
   ```

## Usage

### Running the Analysis

You can run the analysis using the provided `run_analysis.R` script or directly via the command line.

**Input**:

The required input is a CSV file containing the following columns:

- **PATTERN**
- **REFERENCE_REPEATS**
- **COUNTS**
- **hgnc_symbol**
- **CHR**
- **START**

#### 1. Using the `run_analysis.R` Script

The `run_analysis/run_analysis.R` script demonstrates how to execute the analysis with predefined parameters.

```
Rscript run_analysis/run_analysis.R
```

**Parameters to Modify**:

- **Input File**: Modify the `input_file` variable to point to your dataset.
- **Output Directory**: Change the `output_dir` variable if you prefer a different output location.
- **Selected Models**: Adjust the `selected_models` variable to choose different models (e.g., `"lognormal,exponential"`).

#### 2. Directly via Command Line

You can also run the `analysis.R` script directly with custom parameters.

```
Rscript scripts/analysis.R --input data/your_data.csv --output results --models exponential,gaussian,weibull,lognormal --threshold 50
```

**Command-Line Arguments**:

- `-i` or `--input`: **(Required)** Path to the input CSV file containing mutation data.
- `-o` or `--output`: Path to the output directory. Default is `results`.
- `-m` or `--models`: Comma-separated list of models to apply. Options include `exponential`, `gaussian`, `weibull`, `lognormal`. Default is `lognormal,exponential`.
- `-t` or `--threshold`: **(Optional)** Minimum motif occurrence threshold. Motifs with occurrence counts above this threshold will be included in the analysis. Default is `45`. Adjust according to your dataset size; smaller datasets may require a lower threshold.

**Examples**:

1. **Run All Models**:

   ```
   Rscript scripts/analysis.R --input data/your_data.csv --output results --models exponential,gaussian,weibull,lognormal
   ```

2. **Run Only Log-Normal and Weibull Models**:

   ```
   Rscript scripts/analysis.R --input data/your_data.csv --output results --models lognormal,weibull
   ```

3. **Run Only Gaussian Model**:

   ```
   Rscript scripts/analysis.R --input data/your_data.csv --output results --models gaussian
   ```

## Understanding the Output

All results are saved in the specified output directory (default is `results/`).

### Result Files per Model

For each selected model, the following files are generated:

- `<model>.csv`: Contains all results for that model, including:
  - **hgnc_symbol**: Gene symbol.
  - **CHR**: Chromosome.
  - **START**: Genomic start position.
  - **motif**: The microsatellite motif (combination of PATTERN and REFERENCE_REPEATS).
  - **p_value**: Raw p-value for each mutation.
  - **adjusted_p_value**: P-value adjusted for multiple testing using the Benjamini-Hochberg method.

- `Significant_<model>.csv`: Lists significant mutations for that model (adjusted p-value < 0.1).

- `<model>_aggregated.csv`: Contains aggregated results per gene for that model, including:
  - **hgnc_symbol**: Gene symbol.
  - **CHR**: Chromosome.
  - **START**: Genomic start position.
  - **num_p_values**: Number of p-values combined for the gene.
  - **combined_p_value**: Combined p-value using Fisher's method.
  - **adjusted_combined_p_value**: Combined p-value adjusted for multiple testing.

- `Significant_<model>_aggregated.csv`: Lists significant genes for that model after aggregation (adjusted combined p-value < 0.1).

### Plots

- **QQ Plots per Model**:

  - `QQplot_<model>.svg`: QQ plot of adjusted p-values for each model, highlighting significant mutations.

- **Combined Model Fit Plot**:

  - `Indels_Models_fit.svg`: Combined plot of model fits over histograms for selected motifs.

### Comparison Table

- `Indels_Models_fit.csv`: Merged comparison table of models per motif, including AIC, BIC, and LogLikelihood values.

### Directories and Files

```
results/
├── exponential.csv
├── exponential_aggregated.csv
├── Significant_exponential.csv
├── Significant_exponential_aggregated.csv
├── gaussian.csv
├── gaussian_aggregated.csv
├── Significant_gaussian.csv
├── Significant_gaussian_aggregated.csv
├── weibull.csv
├── weibull_aggregated.csv
├── Significant_weibull.csv
├── Significant_weibull_aggregated.csv
├── lognormal.csv
├── lognormal_aggregated.csv
├── Significant_lognormal.csv
├── Significant_lognormal_aggregated.csv
├── QQplot_exponential.svg
├── QQplot_gaussian.svg
├── QQplot_weibull.svg
├── QQplot_lognormal.svg
├── Indels_Models_fit.svg
└── Indels_Models_fit.csv
```

### Notes:

- **Significant Mutations**: Mutations with adjusted p-values less than 0.1 are considered significant.

- **Aggregated Results**: P-values are combined per gene using Fisher's method. If a gene is associated with multiple motifs, their p-values are aggregated to provide a gene-level significance.

- **Genomic Location**: The outputs now include genomic coordinates (**CHR** and **START**) for each mutation or gene.

## Repository Structure

```
MSMuSig2/
├── data/
│   └── example_data.csv
├── scripts/
│   └── analysis.R
├── results/
│   ├── exponential.csv
│   ├── exponential_aggregated.csv
│   ├── Significant_exponential.csv
│   ├── Significant_exponential_aggregated.csv
│   ├── ... (similar files for other models)
│   ├── QQplot_exponential.svg
│   ├── ... (QQ plots for other models)
│   ├── Indels_Models_fit.svg
│   └── Indels_Models_fit.csv
├── run_analysis/
│   └── run_analysis.R
├── README.md
├── LICENSE
├── .gitignore
└── requirements.txt
```

## Methods

### Driver MS-loci Analysis

To identify potential driver mutations, we conduct a statistical analysis on recurrent MS-indel data comprising mutation counts. To accurately model the distribution of these mutation counts, we evaluated several statistical models: **Exponential**, **Gaussian (Normal)**, **Weibull**, and **Log-Normal** distributions. We fitted each selected model to the data and compared their performance using the Akaike Information Criterion (**AIC**) and Bayesian Information Criterion (**BIC**). The model with the lowest AIC and BIC values is deemed the best fit.

Based on this evaluation, we selected the most appropriate model(s) for our data. We grouped the mutation counts according to specific microsatellite motifs and repeat lengths, fitting the chosen model(s) to each group using maximum likelihood estimation. This approach allowed us to estimate the distribution parameters for each microsatellite group.

To detect significant mutations that deviated from the expected distribution within each group, we calculated the cumulative probability (**p-value**) of observing a mutation count equal to or greater than each observed value under the fitted model(s). We then adjusted these p-values for multiple comparisons using the Benjamini-Hochberg false discovery rate method. Mutations with adjusted p-values below a threshold of **0.1** are considered statistically significant and identified as potential driver events.

We further aggregated the p-values per gene using Fisher's method to combine p-values from multiple motifs associated with the same gene. This provides a gene-level significance measure, highlighting genes that are significantly mutated across different motifs.

All significant p-values across these methods are adjusted globally to control the false discovery rate. QQ plots are generated for each model to visualize the distribution of these adjusted p-values.

## Contributing

Contributions are welcome! Please follow these steps to contribute:

1. **Fork the Repository**

2. **Create a New Branch**

   ```
   git checkout -b feature/YourFeatureName
   ```

3. **Make Your Changes**

4. **Commit Your Changes**

   ```
   git commit -m "Add your message here"
   ```

5. **Push to Your Fork**

   ```
   git push origin feature/YourFeatureName
   ```

6. **Open a Pull Request**

## License

This project is licensed under the [MIT License](LICENSE).

## Authors

Hagay Ladany, Dr. Yosef Maruvka.

## Contact

For any questions or suggestions, please open an issue or contact [MaruvkaLab@gmail.com](mailto:maruvkalab@gmail.com).
