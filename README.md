# HeckmanSelect

The package carries out variable selection in binary Heckman selection model and use bootstrap validation technique to estimate optimisim in the metrics for predictive accuracy. AUROC (Area under the receiver operating curve), 
AUCROC (Area under the Receiver Operating Characteristic curve), AUPRC (Area under the precision-recall curve), BS (Brier Score), ECE (Expected Calibration Error) and MCE (Maximum Calibration Error) are implemented.
Lasso and Adaptive Lasso are implemented for variable selection. Normal error and AMH (Ali-Mikhail-Haq) copula errors are implemented. We also implemented the 
bootstrap approach for models developed via variable selection using P-values in the functions "HeckPval" and "bootValidate_Pval".

We also implemented Probit Lasso regression in our "ProbitLasso"" function. It is similar to the GLMNET package as they both implemented the coordinate descent
algorithm. The main addition is that the model can be validated using bootstrap validation method via the function "boot_ProbitLasso".
