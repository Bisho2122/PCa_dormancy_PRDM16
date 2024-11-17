# PCa Dormancy PRDM16

Authors: Bishoy Wadie [bishoy.wadie\@embl.de](mailto:bishoy.wadie@embl.de) and Mostafa Nasr [mostafanasr\@usf.edu](mailto:mostafanasr@usf.edu)

This repository contains the necessary data and scripts to reproduce the plots and analyses from the paper "PRDM16 induces prostate cancer cell dormancy and prevents bone metastatic outgrowth" for reproducibility purposes.

## Contents

1.  Source data for main and supplementary figures are provided in `Data` folder. Count matrices are in `Data/Count_matrices`

2.  Scripts :

    1.  Plotting and utility functions can be found in `Scripts/All_functions.R`

    2.  Figure reproducibility scripts include

        1.  `Scripts/Correlation_GSE_data` for correlation against NED samples in Figure S5

        2.  `Scripts/Original_RNASeq_analysis` for DESeq based RNA-Seq analysis and other figure panels in Figures 2 and S2.

## Steps to reproduce the figures

Follow the steps below to clone the repository, install the required dependencies, and run the scripts to reproduce the figures.

#### Step 1: Clone the repository

Start by cloning the repository to your local machine. You can do this using the following command:

``` bash
git clone https://github.com/Bisho2122/PCa_dormancy_PRDM16.git
```

Navigate to the repository directory

``` bash
cd PCa_dormancy_PRDM16
```

#### Step 2: Install R Dependencies

Inside the repository, you'll find a script named `Install_R_Deps.R`. This script installs all the necessary R packages required to run the analysis scripts. Run it from your R environment:

1.  Open R or RStudio

2.  Set your working directory to the cloned repository

3.  Run the following command to install the dependencies

``` r
source("Install_R_Deps.R")
```

This will automatically install the necessary R packages, such as `ggplot2`, `dplyr`, `tidyr`, and any other dependencies needed to run the analysis scripts.

#### Step 3: Run the scripts to generate figures

After installing the dependencies, you can run the scripts that generate the figures for the paper.

For example, to run correlation analysis and generate correlation plots in figures S5

1.  Open the script `Correlation_GSE_dataset.R` located in the `/Scripts` folder.

2.  Execute the script in your R environment. If you're using RStudio, you can open the script and click **Run**.

3.  The generated plots will be found in `/Scripts/Plots`

## License

This repository is licensed under the MIT License.
