# Capturing infant and child growth dynamics with P-splines mixed effects models

This page contains the R analysis code and the example analysis dataset from the above paper (authors: María Alejandra Hernandez, Zheyuan Li, Tim J Cole, Yi Ying Ong, Kate Tilling, Ahmed Elhakeem). The files are described below:

**1_bspline_example.R**: Generates a figure demonstrating cubic B-spline basis functions.

**2_setup.R**: Loads and cleans GUSTO growth data, then creates and saves a plot of the observed data. It also saves the cleaned data for later use.

**3_lme_models_features.R**: Fits P-spline LME models to the growth data, calculates predictions from these models, and then extracts  growth features (BMI and age at peak BMI and rebound BMI, infant peak growth velocity, and height/weight differences compared to WHO standards). Uses helper functions from the next script (3_1).

**3_1_features_helper_functions.R**: Defines a set of functions for estimating growth features from the P-spline LME predictions.

**4_curves_features_plots.R**: Generates figures based on the models and data prepared in the previous scripts. It plots predicted growth curves, predicted velocity curves, compares predicted vs. observed values, and creates boxplots and correlation plots of the growth features.

**5_prenatal_assoc.R**: Performs a linear regression analysis to examine the association between prenatal factors and the growth features calculated in the previous scripts. It then generates a figure to visualize these associations. The file also Fits LME models to examine height and weight differences from WHO standards up to age 5 years by maternal ethnicity group, and creates a plot of the results.

**6_double_penalty.R**: Uses the helper functions from the next script (6_1) to fit double penalty P-spline LME models and compares them to standard psme models.

**6_1_dp_helper_functions.R**: Defines a set of functions for fitting double penalty P-spline LME models.

**synth_dat_F.csv**: Synthetic copy of the GUSTO girls dataset.

**synth_dat_M.csv**: Synthetic copy of the GUSTO boys dataset.
