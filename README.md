# Marius Wank Masterthesis Repository

This repository contains the Python and R scripts that were used during the master thesis of Marius Wank. Additionally, the Images.zip contains all images created during analysis of the data while the coreABM-main-hyphal-branch-marius.zip contains the code base.

## Repository Structure

### 1. `DataFiles/`
- **Description**: This folder contains an `.xlsx` file with processed data obtained from Hoang et al. through their Imaris analysis.
- **Contents**: 
  - `ProcessedData.xlsx`: The main data file used for the analysis of the raw Imaris data and was used to determine certain parameter values like branching angles.

### 2. `DataRankScripts/`
- **Description**: Python scripts used to rank the results of the various simulation runs. The ranking is based on the relative error of output measurements when compared to the in-vitro data from Hoang et al..
- **Contents**: 
  - `Optimized_rank_data.py`: An updated and optimized version of the rank_data.py script. Utilizes parallelism. 
  - `rank_data.py`: Algorithm that was used to rank simulation results. Outdated, use Optimized_rank_data.py instead.

### 3. `DataVisualizationScripts/`
- **Description**: Python scripts for visualizing the simulation results. This folder also includes the implementation of the Welzl algorithm for visualization and an example video showcasing the Welzl process.
- **Contents**: 
  - `Optimised_VisualizeBoxplots.py`: Updated and optimized Script to visualize the top-performing parameter combinations. Utilizes parallelism. 
  - `Visualize_Boxplots_Growth_Branching.py`: Script that was used to visualize the top-performing parameter combinations. Outdated, use Optimised_VisualizeBoxplots.py instead.
  - `Visualize_Data_over_Time.py`: Script that visualizes the value of certain output measurements over time.
  - `Welzl.py`: Implementation of the Welzl algorithm for determining the smallest circle enclosing the hyphal growth.
  - `welzl_algorithm.avi`: Example video demonstrating the Welzl process.
  - `WelzlVis.py`: Python script used for generating the example video. 
  
### 4. `ImageJMacro/`
- **Description**: Custom ImageJ macro script used to calculate the average fractal dimension of individual conidia.
- **Contents**: 
  - `CustomMacro.ijm`: The ImageJ macro script.

### 5. `InteractionAnalysis/`
- **Description**: Python scripts used to visualize and analyze the interactions between certain parameters.
- **Contents**: 
  - `interaction_checker.py`: Script for generating interaction plots between 3 key parameters.
  - `interaction_checker_for_2.py`: Script for generating interaction plots between 2 key parameters.

### 6. `KolmogorovSmirnov/`
- **Description**: Contains an R script used for statistical evaluation using the Kolmogorov-Smirnov test.
- **Contents**: 
  - `kolmogorovSmirnov_Test.R`: R script to statistically compare an in-vitro data distribution to another distributions.

### 7. `coreABM-main-hyphal-branch-marius.zip`
- **Description**: ZIP file that contains the entire code base that was used in this thesis.

### 8. `Images.zip`
- **Description**: ZIP file that contains the entirety of images that were created during analysis. They were added here to save space in the appendix of the thesis.

### 9. `Pseudocode`
- **Description**: A file containing the pseudocode used in the thesis, outlining the basic agent-based approach used in the thesis.

## Requirements

- **Python 3.9**
- **R**
- **ImageJ/Fiji**
- Python packages: `matplotlib`, `numpy`, `scipy`, `math`, `pandas`

## Contributions

Feel free to open issues or contribute with pull requests if you find bugs or have suggestions for improvements.
