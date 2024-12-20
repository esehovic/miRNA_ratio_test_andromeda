# miRNA_ratio_test_andromeda
The code used for analysing the plasma miRNA ratios for early breast cancer detection explained in the article "Validation of previously identified plasma microRNA ratios for breast cancer detection in a nested case-control study within a screening setting" (https://doi.org/10.1002/ctm2.70068). 

The dataset used in the code can be downloaded from the Zenodo data repository (https://doi.org/10.5281/zenodo.11234225).

In this R code we applied the complete model (miRNA ratios and non-molecular variables) and the miRNA ratio-only model obtained in the discovery set (published here: https://doi.org/10.1038/s41598-023-38886-0) in a validation set. Both sets of samples were obtained from a large prospectively sampled cohort of women undergoing mammography screening (ANDROMEDA). The model in the discovery paper was a penalised lasso logistic regression and the intercept and coefficient were applied using the sigmoid function for logistic regression. We assessed the discriminatory ability of the models in the validation set by investigating the ROC AUC, sensitivity and specificity at Youden's cutoff. Additionally, we assessed the calibration of the predicted probabilities by generating a calibration curve (https://cran.r-project.org/web/packages/CalibrationCurves/vignettes/CalibrationCurves.html) and calculating the Brier score. 

Importantly, in this study our primary goal was not to create a diagnostic model (as the breast cancer event size was low), but to identify most promising biomarkers associated with breast cancer in a screening setting.
