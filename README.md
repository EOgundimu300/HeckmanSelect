# HeckmanSelect

The package carries out variable selection in binary Heckman selection model and use bootstrap validation technique to estimate optimisim in the metrics for predictive accuracy. AUROC (Area under the receiver operating curve), 
AUCROC (Area under the Receiver Operating Characteristic curve), AUPRC (Area under the precision-recall curve), BS (Brier Score), ECE (Expected Calibration Error) and MCE (Maximum Calibration Error) are implemented.
Lasso and Adaptive Lasso are implemented for variable selection. Normal error and AMH (Ali-Mikhail-Haq) copula errors are implemented. We implemented the 
bootstrap approach for models developed via variable selection using P-values in the functions "HeckPval" and "bootValidate_Pval".

We also implemented Probit Lasso regression in our "ProbitLasso"" function. It is similar to the GLMNET package as they both implemented the coordinate descent
algorithm. The main addition is that the model can be validated using bootstrap validation method via the function "boot_ProbitLasso".

Functions implemented (use help (e.g. ?HeckSelect) to read more about the functions)
#### (1) HeckSelect 
  Function for binary outcome with sample selection and variable selection. Adaptive Lasso and Lasso are implemented. Normal error and AMH copula (with probit marginals)     based approach is implemented at the moment.
 
#### (2) bootValidate
  Bootstrap internal validation technique to correct for overoptimism in predictions - "mboot" is the number of bootstrap samples. The function takes the object created by
  HeckSelect and use non-parametric bootstrap method to compute optimism corrected predictive accuracy measures.
 
#### (3) HeckPval
  This function is based on the use of P-value to select variables in Binary Heckman selection model. Default P-value = 0.05. 
  
#### (4) bootValidate_Pval
  Bootstrap internal validation technique to correct for overoptimism in predictions - the alpha value is inherited from the object HeckPval.
  Note that this is different from the "bootValidate" as this is based on dropping variables whose values are greater than the alpha value from the model. If no variable
  selection is required, please set alpha =1 in HeckPval object.
  

The package also contain functions for regularized probit regression. The results are similar to GLMNET package as they both implemented the coordinate descent algorithm.

#### (5) ProbitLasso
  This is probit regression. The missing data is delected to fit the model to complete data
    
#### (6) boot_ProbitLasso
  Bootstrap internal validation technique to correct for overoptimism in predictions - mboot is the number of bootstrap samples. The function takes the object created by 
  ProbitLasso.
    
## How to use the Package
 library(HeckmanSelect)
### Example simulated data
data()
### Data sets in package HeckmanSelect:
#### AmEx: American Express Credit Card data (see: Greene WH (1998) Sample selection in credit-scoring models )
#### binHeckman is a Simulated data



###### selection <- uu~ X1+ X2 + X3+ X4+ X5+ X6+ X7+ X8 + X9 + X10 +X11+X12
###### outcome <- yobs ~ X1+ X2 + X3+ X4+ X5+ X6+ X7+ X8 + X9 + X10 +X11

###### pp <- HeckSelect(selection, outcome, data=binHeckman, allowParallel = TRUE, penalty="ALASSO", Model="AMH",crit="bic")
###### names(pp)# to see objects created
###### options(scipen=999)
###### coef.HeckSelect(pp)# coefficients
###### aa <- bootValidate(pp, data=binHeckman, mboot=100, seed=1)# bootstrap validation
###### aa$resu # optimism corrected metrics



### We can also use the package for Probit regression with Lasso (complete case analysis). Here, we use the American Express Credit Card data for illustration.

##### datt <- AmEx
##### dat <- subset(datt,  !(CARDHLDR ==0))# we selected complete data set
##### default_eq <- DEFAULT~AGE+ACADMOS+ADEPCNT+AEMPMOS+MAJORDRG+ MINORDRG+OWNRENT+APADMOS+AMAMIND+INCOME+SELFEMPL+ TRADACCT+ INCPER+ EXP_INC+CPTOPNB+ CPTOPNG+ CPT30C+CPTF30+CPTAVRV+CBURDEN

##### aa <- ProbitLasso(formula=default_eq, data=dat, allowParallel = TRUE, penalty="ALASSO", crit="bic")
##### Nadlasso <- boot_ProbitLasso(aa, data=dat, mboot=100, seed=1)
