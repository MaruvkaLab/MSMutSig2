# MSMuSig2

**MSMuSig2** is a tool to identify driver genes in Microsatellite-Unstable (MSI) cancers. It analyzes recurrent MS-indel mutation counts using statistical models (Exponential, Gaussian, Weibull, Log-Normal) to find significant mutations.

---

## 💾 Installation

1.  **Clone the repo:**
    ```
    git clone [https://github.com/yourusername/MSMuSig2.git](https://github.com/yourusername/MSMuSig2.git)
    cd MSMuSig2
    ```

2.  **Install R packages:**
    Run this command in your R session:
    ```R
    install.packages(c("magrittr", "minpack.lm", "ggplot2", "gridExtra", "MASS",
                       "dplyr", "data.table", "fitdistrplus", "survival", "car",
                       "ggrepel", "ggpubr", "optparse", "metap"))
    ```

---

## 🚀 Usage

The main script is `scripts/analysis.R`.

### 1. Input File

Your input CSV file **must** contain the following columns:

* `CHR`
* `START`
* `PATTERN`
* `REFERENCE_REPEATS`
* `hgnc_symbol`
* `COUNTS`

### 2. Run from Command Line
Rscript scripts/analysis.R --input data/your_data.csv --output results/

### 3. Command-Line Arguments

* `-i` or `--input`: **(Required)** Path to your input CSV file.
* `-o` or `--output`: **(Optional)** Path to the output directory. (Default: `results`)
* `-m` or `--models`: **(Optional)** Comma-separated list of models (e.g., `lognormal,weibull`). (Default: `lognormal,exponential`)
* `-t` or `--threshold`: **(Optional)** Minimum motif occurrence count to analyze. (Default: `45`)

---

## 📊 Understanding the Output

All results are saved to the output directory (default: `results/`).

**Key Outputs:**

* **Raw Results:** `<model>.csv` (e.g., `lognormal.csv`)
    * Contains all p-values for every mutation.
* **Aggregated Results:** `<model>_aggregated.csv`
    * Contains gene-level combined p-values (using Fisher's method).
* **Significant Hits:** `Significant_<model>.csv` and `Significant_<model>_aggregated.csv`
    * Filtered lists of significant mutations and genes (adjusted p-value < 0.1).
* **Plots:** `QQplot_<model>.svg` and `Indels_Models_fit.svg`
    * Visualizations of p-value distributions and model fits.
* **Model Comparison:** `Indels_Models_fit.csv`
    * AIC/BIC scores for each model fit.

---

## 📄 License

This project is licensed under the [MIT License](LICENSE).

## 🧑‍💻 Authors

Hagay Ladany, Dr. Yosef Maruvka.

## 📧 Contact

For questions, please open an issue or contact [MaruvkaLab@gmail.com](mailto:maruvkalab@gmail.com).
