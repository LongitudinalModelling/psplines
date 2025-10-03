# Capturing infant and child growth dynamics with P-splines mixed effects models

This page contains the R analysis code and the example (synthetic) dataset from the above paper (authors: María Alejandra Hernandez, Zheyuan Li, Tim J Cole, Yi Ying Ong, Kate Tilling, Ahmed Elhakeem). The files are described below:

- **1_bspline_example.R**: Generates a figure demonstrating cubic B-spline basis functions.

- **2_setup.R**: Loads and cleans GUSTO growth data, then creates and saves a plot of the observed data. It also saves the cleaned data for later use.

- **3_lme_models_features.R**: Fits P-spline LME models to the growth data, calculates predictions from these models, and then extracts  growth features (BMI and age at peak BMI and rebound BMI, infant peak growth velocity and growth velocity at ages 1, 6, 12, and 24 months, and height/weight differences compared to WHO standards up to age 5 years). Uses helper functions from the next script (3_1).

- **3_1_features_helper_functions.R**: Defines a set of functions for estimating growth features from the P-spline LME predictions.

- **4_dist_vel_plots**: Generates figures based on the models and data prepared in the previous scripts. Plots the predicted growth curves, predicted velocity curves, and the predicted vs. observed values for 3 boys and 3 girls.

- **5_curves_features_plots.R**: Generates figures based on the models and data prepared in the previous scripts. Creates boxplots and correlation plots of the main growth features (peak/rebound BMI and infant peak growth velocity), and boxplots and correlation plots of infant growth velocity at ages 1, 6, 12, and 24 months. The file also Fits LME models to examine height and weight differences from WHO standards up to age 5 years and creates a plot of the results.

- **6_prenatal_assoc.R**: Performs linear regression analyses to examine the association between prenatal factors and the main derived growth features. Generates a figure to visualise the results. 

- **7_double_penalty.R**: Uses the helper functions from the next script (7_1) to fit double penalty P-spline LME models and compares them to standard psme models.

- **7_1_dp_helper_functions.R**: Defines a set of functions for fitting double penalty P-spline LME models.

- **synth_dat_F.csv**: Synthetic copy of the GUSTO girls dataset.

- **synth_dat_M.csv**: Synthetic copy of the GUSTO boys dataset.
