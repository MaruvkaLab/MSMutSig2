# MSMuSig2

**MSMuSig2** is a tool designed to identify driver genes in Microsatellite-Unstable (MSI) cancers by analyzing recurrent MS-indel mutation counts.
It employs various statistical models to determine significant driver mutations.

## Table of Contents

- [Features](#features)
- [Installation](#installation)
- [Usage](#usage)
- [Repository Structure](#repository-structure)
- [Methods](#methods)
- [Contributing](#contributing)
- [License](#license)

## Features

- Fits multiple statistical models (**Exponential**, **Gaussian**, **Weibull**, **Log-Normal**) to mutation count data.
- Allows users to select which models to apply during analysis.
- Compares models using **AIC** and **BIC** to determine the best fit.
- Identifies significant driver mutations based on adjusted p-values across all models.
- Generates comprehensive visualizations, including model fits and a single QQ plot.

## Installation

1. **Clone the Repository**

   ```bash
   git clone https://github.com/yourusername/MSMuSig2.git
   cd MSMuSig2
   ```

2. **Install Required R Packages**

   Ensure you have the necessary R packages installed. You can install them using the following command in R:

   ```r
   install.packages(c("magrittr", "minpack.lm", "ggplot2", "gridExtra", "MASS", 
                      "dplyr", "data.table", "fitdistrplus", "survival", "car", 
                      "ggrepel", "ggpubr", "optparse"))
   ```

   Alternatively, you can install all required packages at once using the `requirements.txt`:

   ```r
   install.packages(readLines("requirements.txt"), dependencies = TRUE)
   ```

## Usage

### Running the Analysis

You can run the analysis using the provided run_analysis script or directly via the command line.

Input:
The required input is a CSV file containing the columns: PATTERN, REFERENCE_REPEATS, COUNTS, and hgnc_symbol.

#### 1. **Using the run_analysis Script**

The `run_analysis/run_analysis.R` script demonstrates how to execute the analysis with predefined parameters.

```bash
Rscript run_analysis/run_analysis.R
```

**Parameters to Modify**:

- **Input File**: Modify the `input_file` variable to point to your dataset.
- **Output Directory**: Change the `output_dir` variable if you prefer a different output location.
- **Selected Models**: Adjust the `selected_models` variable to choose different models (e.g., `"lognormal,exponential"`).

#### 2. **Directly via Command Line**

You can also run the `analysis.R` script directly with custom parameters.

```bash
Rscript scripts/analysis.R --input data/your_data.csv --output results --models exponential,gaussian,weibull,lognormal --threshold 50
```

**Command-Line Arguments**:

- `-i` or `--input`: **(Required)** Path to the input CSV file containing mutation data.
- `-o` or `--output`: Path to the output directory. Default is `results`.
- `-m` or `--models`: Comma-separated list of models to apply. Options include `exponential`, `gaussian`, `weibull`, `lognormal`. Default is `lognormal,exponential`.
- `-t` or `--threshold`: **(Optional)** Minimum `REFERENCE_REPEATS` length to select motifs. Default is `45`.

**Examples**:

1. **Run All Models**:

   ```bash
   Rscript scripts/analysis.R --input data/your_data.csv --output results --models exponential,gaussian,weibull,lognormal
   ```

2. **Run Only Log-Normal and Weibull Models**:

   ```bash
   Rscript scripts/analysis.R --input data/your_data.csv --output results --models lognormal,weibull
   ```

3. **Run Only Gaussian Model**:

   ```bash
   Rscript scripts/analysis.R --input data/your_data.csv --output results --models gaussian
   ```

### Understanding the Output

- **Results Directory (`results/`)**:

  - **driver_analysis/**:
    - `comparison_raw-Model_VS_RAW.csv`: Comparison metrics (**AIC**, **BIC**, **LogLikelihood**) for each selected RAW model per motif.
    - `Indels_Models_fit.csv`: Merged comparison table sorted by Motif and AIC.

  - **normality_check/**:
    - `Significant.csv`: Consolidated list of significant mutations across all selected models with adjusted p-values < 0.1.

- **Plots Directory (`plots/`)**:

  - `Indels_Models_fit.png` and `.svg`: Combined plots of model fits for selected motifs.
  - `QQplot.svg`: Single comprehensive QQ plot displaying adjusted p-values across all models and motifs.

## Repository Structure

```
MSMuSig2/
├── data/
│   └── example_data.csv
├── scripts/
│   └── analysis.R
├── results/
│   ├── driver_analysis/
│   │   ├── comparison_raw-Model_VS_RAW.csv
│   │   └── Indels_Models_fit.csv
│   └── normality_check/
│       └── Significant.csv
├── plots/
│   ├── Indels_Models_fit.png
│   ├── Indels_Models_fit.svg
│   └── QQplot.svg
├── run_analysis/
│   └── run_analysis.R
├── README.md
├── LICENSE
├── .gitignore
└── requirements.txt
```

## Methods

### Driver MS-loci

To identify potential driver mutations, we conduct a statistical analysis on recurrent MS-indel data comprising mutation counts.
To accurately model the distribution of these mutation counts, we evaluated several statistical models: **Exponential**, **Gaussian (Normal)**, **Weibull**, and **Log-Normal** distributions.
We fitted each selected model to the data and compared their performance using the Akaike Information Criterion (**AIC**) and Bayesian Information Criterion (**BIC**).
The model with the lowest AIC and BIC values is deemed the best fit.

Based on this evaluation, we selected the most appropriate model(s) for our data.
We grouped the mutation counts according to specific microsatellite motifs and repeat lengths, fitting the chosen model(s) to each group using maximum likelihood estimation.
This approach allowed us to estimate the distribution parameters—the mean (μ) and standard deviation (σ) of the logarithm of the mutation counts—for each microsatellite group.

To detect significant mutations that deviated from the expected distribution within each group, we calculated the cumulative probability (**p-value**) of observing a mutation count equal to or greater than each observed value under the fitted model(s).
We then adjusted these p-values for multiple comparisons using the Benjamini-Hochberg false discovery rate method.
Mutations with adjusted p-values below a threshold of **0.1** are considered statistically significant and identified as potential driver events.

All significant p-values across these methods are adjusted globally to control the false discovery rate, and a single QQ plot is generated to visualize the distribution of these adjusted p-values.

## Contributing

Contributions are welcome! Please follow these steps to contribute:

1. **Fork the Repository**

2. **Create a New Branch**

   ```bash
   git checkout -b feature/YourFeatureName
   ```

3. **Make Your Changes**

4. **Commit Your Changes**

   ```bash
   git commit -m "Add your message here"
   ```

5. **Push to Your Fork**

   ```bash
   git push origin feature/YourFeatureName
   ```

6. **Open a Pull Request**

## License

This project is licensed under the [MIT License] (LICENSE) haga.

## Authors
Hagay Ladany, Dr. Yosef Maruvka.

## Contact

For any questions or suggestions, please open an issue or contact [MaruvkaLab@gmail.com](mailto:maruvkalab@gmail.com).
