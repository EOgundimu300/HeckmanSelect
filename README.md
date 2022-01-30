# HeckmanSelect

The package carries out variable selection in binary Heckman selection model and use bootstrap validation technique to estimate optimisim in the metrics for predictive accuracy. AUROC (Area under the receiver operating curve), 
AUCROC (Area under the Receiver Operating Characteristic curve), AUPRC (Area under the precision-recall curve), BS (Brier Score), ECE (Expected Calibration Error) and MCE (Maximum Calibration Error) are implemented.
Lasso and Adaptive Lasso are implemented for variable selection. Normal error and AMH (Ali-Mikhail-Haq) copula errors are implemented. We also implemented the 
bootstrap approach for models developed via variable selection using P-values in the functions "HeckPval" and "bootValidate_Pval".

We also implemented Probit Lasso regression in our "ProbitLasso"" function. It is similar to the GLMNET package as they both implemented the coordinate descent
algorithm. The main addition is that the model can be validated using bootstrap validation method via the function "boot_ProbitLasso".

Functions implemented (use help to read more about the functions)
### (1) HeckSelect 
 Function for binary Heckman model with variable selection. Adaptive Lasso and Lasso are implemented. Normal error and AMH copula (with probit marginals) based approach is implemented.
 
### (2) bootValidate
 Bootstrap internal validation technique to correct for overoptimism in predictions - mboot is the number of bootstrap samples.
### (3) HeckPval
  This function is based on the use of P-value to select variables in Binary Heckman selection model. Default P-value = 0.05. 
### (4) bootValidate_Pval
  Bootstrap internal validation technique to correct for overoptimism in predictions - the alpha value is inherited from the object HeckPval.
  Note that this is different from the "bootValidate" as this is based on dropping variables whose values are greater than the alpha value from the model. If no variable
  selection is required, please set alpha =1 in HeckPval object.
