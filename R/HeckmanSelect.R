

###################################################################
# Function for binary Heckman model with variable selection
# Adaptive Lasso and Lasso are implemented. Normal error and AMH
# copula (with probit marginals) based approach is implemented.
####################################################################

#' Title Regularization method for binary Heckman selection model
#'
#' @param selection selection equation
#' @param outcome outcome equation
#' @param data data matrix containing both the outcome and selection variables
#' @param lambda shrinkage parameter, both scalar and vector are acceptable
#' When lambda=NULL, the internal vector of Lambda is used
#' @param allowParallel If true, the "doParallel" package is invoked
#' @param penalty can be ALASSO (for adaptive lasso) or LASSO (for Lasso) penalty
#' @param Model can either be Normal error of AMH (Ali-Mikhail-Haq) copula function
#' @param crit can be BIC, AIC or GCV, default is BIC
#' @param ...
#'
#' @return class HeckSelect containing penalized coefficients and
#' the parameters supplied for the creation of the object. Function call is
#' also returned
#' @export
#'
#' @examples
#' HeckSelect(selection, outcome, data=data, allowParallel = TRUE, penalty="ALASSO", Model="AMH",crit="bic")
#'
HeckSelect <- function(selection, outcome, data = sys.frame(sys.parent()), lambda = NULL, allowParallel = FALSE,
             penalty=c("LASSO","ALASSO"), Model = c("Normal","AMH"), crit=c("bic","aic","gcv"),...)
{
    if (match("pbivnorm",.packages(),0)==0) require(pbivnorm)
    if (match("copula",.packages(),0)==0) require(copula)
    if (match("foreach",.packages(),0)==0) require(foreach)
    if (allowParallel) library(doParallel)

    if (!missing(data)) {
        if (!inherits(data, "environment") & !inherits(data,
            "data.frame") & !inherits(data, "list")) {
            stop("'data' must be either environment, data.frame, or list (currently a ",
                class(data), ")")
        }
    }

         if (!Model %in% c("Normal", "AMH"))
            stop("Model must be Normal, AMH")
        funcCall <- match.call(expand.dots = FALSE)

    mf <- model.frame(selection, data=data)
    YS <- model.response(mf, "numeric")
    XS <- model.matrix(selection, data = data)
    NXS <- ncol(XS)

    mf2 <- model.frame(outcome, data)
    YO <- model.response(mf2, "numeric")
    XO <- model.matrix(outcome, data = data)
    NXO <- ncol(XO)
    ibetaS <- 1:NXS
    ibetaO <- seq(tail(ibetaS, 1) + 1, length = NXO)
    irho <- tail(ibetaO, 1) + 1

    size_s <- NXS
    size_o <- NXO
     Model <- match.arg(Model)
    pseudo <- pseudo_data(selection, outcome, data, Model)
    xstar <- pseudo$xstar
    ystar <- pseudo$ystar
    crit <- crits <- match.arg(crit)
    penalty  <-  match.arg(penalty)
    n=length(YS)

         if (allowParallel) {
            `%op%` <- `%dopar%`
            cl <- makeCluster(detectCores())
            clusterExport(cl, c("coord_select","softrule",
               "loglik","log_amhcop","makePD"))
            registerDoParallel(cl)
        } else {
            `%op%` <- `%do%`
        }

     if(is.null(lambda)) {

     lambda <- seq(0,10, length=50)
    }


      fitter <- function(size_s, size_o, xstar, ystar, lambda, penalty){
     fit <-  coord_select(size_s, size_o, xstar, ystar, lambda, penalty)

   if(Model == "Normal"){
     beta <- fit$beta
     df  <- fit$df
     fn <- loglik(as.matrix(beta),YS,XS,YO,XO)
      bic <- 2*fn + log(n)*df
      aic <- 2*fn +2*df
      gcv <- fn/(n*(1-df/n)^2)
     out <- list(beta = beta, df = df, bic= bic, aic = aic, gcv = gcv)
     }

   if(Model == "AMH"){
     beta <- fit$beta
     df  <- fit$df
     fn <- log_amhcop(as.matrix(beta),YS,XS,YO,XO)
      bic <- 2*fn + log(n)*df
      aic <- 2*fn +2*df
      gcv <- fn/(n*(1-df/n)^2)
     out <- list(beta = beta, df = df, bic= bic, aic = aic, gcv = gcv)
     }
      return(out)
    }

           comb <- function(...) {
                   mapply('cbind', ..., SIMPLIFY=FALSE)
         }

   H  <- foreach(i = lambda, .combine='comb', .multicombine=TRUE,
        .verbose = FALSE, .errorhandling = "stop") %op% fitter(size_s, size_o, xstar, ystar, i, penalty)

       if (allowParallel) stopCluster(cl)

       beta <- H$beta
       df <- as.vector(H$df)
       bic <- as.vector(H$bic)
       aic <-  as.vector(H$aic)
       gcv <-  as.vector(H$gcv)


   crit <- switch(crit, bic=bic, aic=aic, gcv=gcv)
    selected <- best.model <- which(crit == min(crit,na.rm=TRUE))
   ic <- c(bic=bic[selected],aic=aic[selected], gcv=gcv[selected])
    final_coeff <- beta[,selected]
   rename <- c(paste("S:",colnames(XS), sep=""),paste("O:",colnames(XO), sep=""),"theta")
   names(final_coeff) <- rename

   lambda_f <- min(lambda[selected])

     beta_t <- final_coeff
     i11 <- !(YS == 0)
     betaS <-  beta_t[ibetaS]
     betaO <-  beta_t[ibetaO]
     rho <- beta_t[irho]

      n.rho <- tanh(rho)
   n.rho2 <- ifelse(n.rho > 0.99, 0.99, ifelse(n.rho < -0.99,-0.99,n.rho))
   XS11.b <- drop(XS[i11, , drop = FALSE] %*% betaS)
   XO11.b <- drop(XO[i11, , drop = FALSE] %*% betaO)

    if(Model == "Normal"){
     rho_rep <- rep(n.rho2, length(XS11.b))
     p.pred <- pbivnorm(XO11.b, XS11.b, rho_rep)/pnorm(XS11.b)
    }

   if(Model == "AMH"){
     XXS11.b <- pnorm(XS11.b)
     XXO11.b <- pnorm(XO11.b)
     amh.cop <- amhCopula(n.rho2)
     p.pred <-  pCopula(cbind(XXS11.b, XXO11.b), amh.cop)/(1-pnorm(-XS11.b))
    }


     yyf <- ifelse(i11==1, YO,NA)
     yselect <- as.vector(na.omit(yyf))
     result <- structure(list(call = funcCall,coef=beta_t, betaS=betaS, betaO=betaO, lpS=XS11.b,lpO=XO11.b, pred=p.pred, yout=yselect,
                 rho=n.rho2,lambda=lambda_f, selection=selection, outcome = outcome,
               crit=crits, penalty = penalty, allowParallel = allowParallel, Model = Model),
               class = "HeckSelect")

   class(result) <- "HeckSelect"
    return(result)
}

#pp <- HeckSelect(selection, outcome, data=data, allowParallel = TRUE, penalty="ALASSO", Model="AMH",crit="bic")



#########################################################################
# Coefficients of object of class HeckSelect
#########################################################################

#' Title Coefficients from objects returned by HeckSelect function
#'
#' @param object a fitted object of class inheriting from "HeckSelect"
#' @param ...
#'
#' @return vector of coefficients for the selection and outcome models
#' @export
#'
#' @examples
#' coef(HeckSelect(selection, outcome, data=data, allowParallel = TRUE, penalty="ALASSO", Model="AMH",crit="bic"))

coef.HeckSelect<- function(object,...){
 x <- object$coef
return(x)
}


#####################################################################
# Predict function for object of class HeckSelect
##########################################################################
#' Title Predict Method for HeckSelect Fits
#'
#' @param object a fitted object of class inheriting from "HeckSelect"
#' @param newdata optionally, a data frame in which to look for variables with which to predict. If omitted, the fitted linear predictors are used.
#' @param na.action function determining what should be done with missing values in newdata. The default is to predict NA
#' @param ...
#'
#' @return vector of predicted values for P(outcome=1|selection=1)
#' @export
#'
#' @examples
#' pp <- HeckSelect(selection, outcome, data=data, allowParallel = TRUE, penalty="ALASSO", Model="AMH",crit="bic")
#' predict(pp)

predict.HeckSelect <- function(object, newdata=NULL,na.action = na.pass,...){
  if (match("pbivnorm",.packages(),0)==0) require(pbivnorm)
  if (match("copula",.packages(),0)==0) require(copula)
  if (missing(newdata)){
    pred <- object$pred
  } else {
    selection <- object$selection
    outcome <- object$outcome
    betaS <- object$betaS
    betaO <- object$betaO
    rho <- object$rho
    Model <- object$Model
    mf <- model.frame(selection, newdata)
    YS <- model.response(mf, "numeric")
    XS <- model.matrix(selection, newdata)

    mf2 <- model.frame(outcome, newdata)
    YO <- model.response(mf2, "numeric")
    XO <- model.matrix(outcome, newdata)

    i11 <- !(YS == 0)
    yyf <- ifelse(i11==1, YO,NA)
    yy <- as.vector(na.omit(yyf))

    XS11.b <- drop(XS[i11, , drop = FALSE] %*% betaS)
    XO11.b <- drop(XO[i11, , drop = FALSE] %*% betaO)

    if(Model == "Normal"){
      rho_rep <- rep(rho, length(XS11.b))
      pred <- pbivnorm(XO11.b,XS11.b, rho_rep)/pnorm(XS11.b)
    }
    if(Model == "AMH"){
      XXS11.b <- pnorm(XS11.b)
      XXO11.b <- pnorm(XO11.b)
      amh.cop <- amhCopula(rho)
      pred <-  pCopula(cbind(XXS11.b, XXO11.b), amh.cop)/(1-pnorm(-XS11.b))
    }
  }
  pred
}


#options(scipen=999)
#pcomb <- predict(pp)

#######################################################################
# Soft threshold operator
########################################################################
softrule=function(beta, lambda)
{
  (lambda>=abs(beta))*0+
    ((lambda<abs(beta)) & (beta>0))*(beta-lambda)+
    ((lambda<abs(beta)) & (beta<0))*(beta+lambda)
}





###########################################################################
# Coordinate descent algorithm
##########################################################################

coord_select <- function(size_s, size_o, x, y, lambda, penalty=c("LASSO","ALASSO"),...)
{
  s <- size_s; o <- size_o
  x <- x
  y <- y
  n <- dim(x)[1]
  p <- dim(x)[2]
  #lll <- s+o+1
  if(penalty=="LASSO"){
    weight <- rep(1, each= p)
  }
  if(penalty=="ALASSO"){
    weight <- as.vector(1/abs(glm(y~x-1)$coef))
  }

  lambda1 <- lambda
  lambda2 <- lambda1*weight

  maxstep <- min(n, p)*500
  beta <- matrix(NA, maxstep+1, p)
  beta[1, ] <- glm(y~x-1)$coef
  delta <- 1
  K <- 1
  while((delta>1e-10) & (K<maxstep))
  {
    K <- K+1
    beta.temp <- beta[K-1, ]  #working downwards through rows, starting with lm coeffs

    for (j in 1:p)
    {
      xminusj <- x[, -j]         #eliminate jth column from data matrix
      bminusj <- beta.temp[-j]
      yminusj <- xminusj%*%bminusj #fitted values (for each nrow)
      rminusj <- y-yminusj         #target - fitted values

      a11 <- softrule(sum(x[, j]*rminusj), lambda=0)/sum(x[, j]^2)
      a12 <- softrule(sum(x[, j]*rminusj), lambda=lambda2[j])/sum(x[, j]^2)
      bj <- ifelse((j==1) | (j==(s+1)) | (j==s+o+1), a11, a12)
      beta.temp[j] <- bj
    }
    beta[K, ] <- beta.temp                        #save new coeffs
    delta <- max(abs(beta[K, ]-beta[K-1, ]))      #compare to previous to check convergence
  }

  beta <- beta[K, ]

  beta <- as.matrix(beta)
  df <- sum(beta !=0)-3
  if(is.null(colnames(x))){colnames(x)=c(paste("x", 1:p, sep=""))}
  rownames(beta) <- colnames(x)

  object <- list(beta=beta, lambda=lambda,df = df, K=K, delta=delta)
  return(object)
}



####################################################################
#Likelihood function for binary sample selection (Normal error)
#########################################################################

loglik <- function(beta,YS,XS,YO,XO)
{
  if (match("pbivnorm",.packages(),0)==0) require(pbivnorm)
  nObs <- length(YS)
  #NO <- length(YS[YS > 0])
  i00 <- YS == 0
  i10 <- (YS == 1) & (YO == 0)
  i11 <- (YS == 1) & (YO == 1)
  NXS <- ncol(XS)
  NXO <- ncol(XO)
  ibetaS <- 1:NXS
  ibetaO <- seq(tail(ibetaS, 1) + 1, length = NXO)
  irho <- tail(ibetaO, 1) + 1
  betaS <- beta[ibetaS]
  betaO <- beta[ibetaO]
  rho <- tanh(beta[irho])
  XS00.b <- drop(XS[i00, , drop = FALSE] %*% betaS)
  XS10.b <- drop(XS[i10, , drop = FALSE] %*% betaS)
  XS11.b <- drop(XS[i11, , drop = FALSE] %*% betaS)
  XO10.b <- drop(XO[i10, , drop = FALSE] %*% betaO)
  XO11.b <- drop(XO[i11, , drop = FALSE] %*% betaO)
  loglik <- numeric(nObs)
  loglik[i00] <-  pnorm(-XS00.b, log.p = TRUE)
  rho_rep <- rep(rho, length(XS10.b))
  f2 <- pbivnorm(XS10.b, -XO10.b, -rho_rep)
  loglik[i10] <- log(f2)
  rho_rep2 <- rep(rho, length(XS11.b))
  f3 <- pbivnorm(XS11.b, XO11.b,  rho_rep2)
  loglik[i11] <-  log(f3)
  return(-sum(loglik))
}




####################################################################
#Likelihood function for AMH copula with Normal marginals
#########################################################################

log_amhcop <- function(beta,YS,XS,YO,XO)
{
  if (match("copula",.packages(),0)==0) require(copula)
  nObs <- length(YS)
  #NO <- length(YS[YS > 0])
  i00 <- YS == 0
  i10 <- (YS == 1) & (YO == 0)
  i11 <- (YS == 1) & (YO == 1)
  NXS <- ncol(XS)
  NXO <- ncol(XO)
  ibetaS <- 1:NXS
  ibetaO <- seq(tail(ibetaS, 1) + 1, length = NXO)
  irho <- tail(ibetaO, 1) + 1
  betaS <- beta[ibetaS]
  betaO <- beta[ibetaO]
  rho <- tanh(beta[irho])

  XS00.b <- drop(XS[i00, , drop = FALSE] %*% betaS)
  XS10.b <- drop(XS[i10, , drop = FALSE] %*% betaS)
  XS11.b <- drop(XS[i11, , drop = FALSE] %*% betaS)
  XO10.b <- drop(XO[i10, , drop = FALSE] %*% betaO)
  XO11.b <- drop(XO[i11, , drop = FALSE] %*% betaO)
  loglik <- numeric(nObs)
  loglik[i00] <-  pnorm(-XS00.b, log.p = TRUE)

  XXS10.b <- pnorm(XS10.b)
  XXO10.b <- pnorm(XO10.b)
  amh.cop <- amhCopula(rho)
  f2 <- (1-pnorm(-XS10.b)) - pCopula(cbind(XXS10.b, XXO10.b), amh.cop)

  loglik[i10] <- log(f2)
  XXS11.b <- pnorm(XS11.b)
  XXO11.b <- pnorm(XO11.b)
  f3 <- pCopula(cbind(XXS11.b, XXO11.b), amh.cop)
  loglik[i11] <-  log(f3)
  return(-sum(loglik))
}



#####################################################################
# Function to make the hessian/covariance matrix positive definite
# Alternatively use nearPD function from the package Matrix
# Function is taken from Chatterjee et al. (2018) to stabilize the computation
##########################################################################
makePD = function(mat){
  N = nrow(mat)
  HC = mat
  D = eigen(mat)
  E = D$values
  U = D$vectors
  v = as.numeric(E < 0)
  idx = which(v==1)
  m = sum(v) # number of negative values
  if(m > 0){
    S = sum(v*E)*2
    W = (S*S*100)+1
    P = min(abs(E[idx])) # smallest positive value
    for(i in idx){
      C = E[i]
      E[i] = P * (S-C)*(S-C)/W
    }
  }
  return(E)
}


#######################################################################
#Pseudo Data Generation for both the Normal and AMH model error
#####################################################################
#' Title Pseudo data generation for the Least Square Approximation Method
#'
#' @param selection selection equation
#' @param outcome outcome equation
#' @param data data data matrix containing both the outcome and selection variables
#' @param Model can either be Normal error of AMH (Ali-Mikhail-Haq) copula function
#'
#' @return xstar: a square matrix containing all the predictors along with
#' the correlation and the intercepts
#'
#' ystar: corresponding pseudo response
#' @export
#'
#' @examples
#' pseudo_data(selection, outcome,data=data, Model="AMH")

pseudo_data <- function(selection, outcome,data, Model = c("Normal","AMH")){
  require(GJRM)
  options(warn=-1)

  if (!Model %in% c("Normal", "AMH"))
    stop("Model must be Normal, AMH")
  funcCall <- match.call(expand.dots = FALSE)

  Model <- match.arg(Model)
  if(Model == "Normal"){
    m2 <- gjrm(list(selection, outcome), data = data, BivD = "N",
               margins = c("probit", "probit"), Model = "BSS")
  }
  if(Model == "AMH"){
    m2 <- gjrm(list(selection, outcome), data = data, BivD = "AMH",
               margins = c("probit", "probit"), Model = "BSS")
  }

  coefs <- coef(m2, part = "full")
  covb <- m2$Vb

  e <- eigen(covb)
  aa <- e$vectors %*% diag(1/sqrt(abs(e$values))) %*% t(e$vectors)
  bb <- e$vectors %*% diag(1/sqrt(makePD(covb))) %*% t(e$vectors)
  sigma2 <-  ifelse((det(covb) > 0),list(aa),list(bb))
  xstar <- sigma2[[1]]
  colnames(xstar) <- colnames(covb)
  ystar <- xstar%*%as.vector(coefs)
  XY <- list(xstar,ystar)
  names(XY) <- c('xstar','ystar')
  return(XY)
}

#hh <- pseudo_data(selection, outcome, dat)




#########################################################################
# Expected Calibration Error (ECE) and Maximum Calibration Error
#######################################################################
ece_mcefun <- function(y, prob, g){
  mtx = cbind(y, y_not = 1- y, prob, prob_not = 1-prob)
  mtx = as.data.frame(mtx)
  mtx = mtx[order(mtx$prob),]
  n <- length(prob)/g
  nr <- nrow(mtx)
  split_mtx = split(mtx, rep(1:ceiling(nr/n), each=n, length.out=nr))

  H_stat = c()
  for (i in 1:length(split_mtx)){
    obs = mean(split_mtx[[i]]$y == 1)
    exp = mean(split_mtx[[i]]$prob)

    H_stat = c(H_stat, abs(obs - exp))

  }
  return(list(ece = sum(H_stat)/length(split_mtx), mce = max(H_stat)))
}

#aab <- ece_mcefun(yobs, prob, g=10)




#######################################################################
## Bootstrap internal validation technique to correct for overoptimism
# in predictions - mboot is the number of bootstrap samples.
########################################################################
#' Title Bootstrap validation for object of model "HeckSelect"
#'
#' @param object a fitted object of class inheriting from "HeckSelect"
#' @param data data matrix containing both the outcome and selection variables
#' @param mboot number of bootstrap samples. 50 bootstrap samples are
#' used if not specified
#' @param seed integer value set for reproducibility
#' @param ...
#'
#' @return coefficients: same as in coef
#'
#' nonconvergence: the number of samples that did not converge
#' in the bootstrap sample
#'
#' lambda: optimal shrinkage parameter
#'
#' resu: data matrix containing apparent, bootstrap, test and optimism
#' corrected performance measures
#'
#' Vectors of performance measures for each bootstrap sample are also returned
#'
#' @export
#'
#' @examples
#' pp <- HeckSelect(selection, outcome, data=data, allowParallel = TRUE, penalty="ALASSO", Model="AMH",crit="bic")
#' bootValidate(pp, data=data, mboot=2, seed=1)

bootValidate <- function(object, data, mboot=50, seed,...){
  if (match("PRROC",.packages(),0)==0) require(PRROC)
  set.seed(seed)
  selection <- object$selection; outcome <- object$outcome; penalty <- object$penalty
  crit <- object$crit; allowParallel <- object$allowParallel
  Model <- object$Model

  auc <- pauc <- brier   <- ece <- mce <- c()
  auctest <- pauctest <- briertest <- ecetest <- mcetest <- c()

  # Looping to generate m bootstrap data
  data1 <- data
  bootdata <- list()
  for(i in 1:mboot)
  {

    iboot <- sample(1:nrow(data1), replace=TRUE)
    bootdata[[i]] <- data1[iboot,]

    options(warn=2)
    res <- try(HeckSelect(selection, outcome, data=bootdata[[i]],allowParallel, penalty=penalty, Model=Model, crit=crit),silent=TRUE)
    if (any(class(res)=="try-error"))
    {
      auc[i] <- pauc[i] <-  brier[i] <- ece[i] <- mce[i] <-  NA
      auctest[i] <- pauctest[i] <- briertest[i] <- ecetest[i] <- mcetest[i] <- NA
    } else{

      p <- res$pred
      yboot <- res$yout
      fg <- p[yboot == 1]
      bg <- p[yboot == 0]
      roc <- roc.curve(scores.class0 = fg, scores.class1 = bg)
      pr <- pr.curve(scores.class0 = fg, scores.class1 = bg)
      auc[i] <- roc$auc
      pauc[i] <- pr$auc.integral
      brier[i] <- mean((p-yboot)^2)
      b_boot <- ece_mcefun(yboot, p, g=10)
      ece[i] <- b_boot$ece
      mce[i] <- b_boot$mce

      # predicting test data back in the original data

      yy <- object$yout
      prtest <- predict(res, data)
      fgt <- prtest[yy == 1]
      bgt <- prtest[yy == 0]

      roct <- roc.curve(scores.class0 = fgt, scores.class1 = bgt)
      prt <- pr.curve(scores.class0 = fgt, scores.class1 = bgt)
      auctest[i] <- roct$auc
      pauctest[i] <- prt$auc.integral
      briertest[i] <- mean((prtest-yy)^2)
      b_test <- ece_mcefun(yy, prtest, g=10)
      ecetest[i] <- b_test$ece
      mcetest[i] <- b_test$mce
    }
  }

  auc_r <- auc; auctest_r <- auctest; pauc_r <- pauc; pauctest_r <- pauctest
  brier_r <- brier; briertest_r <- briertest; ece_r <- ece; ecetest_r <- ecetest
  mce_r <- mce; mcetest_r <- mcetest

  ck <- sum(is.na(auc_r))
  auc <- mean(auc, na.rm=T); pauc <- mean(pauc, na.rm=T); brier <- mean(brier, na.rm=T)
  ece <- mean(ece, na.rm=T); mce <- mean(mce, na.rm=T)

  auct <- mean(auctest, na.rm=T); pauct <- mean(pauctest, na.rm=T)
  briert <- mean(briertest, na.rm=T); ecet <- mean(ecetest, na.rm=T)
  mcet <- mean(mcetest, na.rm=T)

  # Index measures from the orginal data set
  options(warn=2)
  orig <- try(object, silent=TRUE)

  if (any(class(orig)=="try-error"))
  {
    Oauc <- Opauc <-  Obrier <- Oece <- Omce <- NA
  } else{
    prorig <- orig$pred
    yy <- orig$yout
    fgo <- prorig[yy == 1]
    bgo <- prorig[yy == 0]
    oroc <- roc.curve(scores.class0 = fgo, scores.class1 = bgo)
    opr <- pr.curve(scores.class0 = fgo, scores.class1 = bgo)
    Oauc <- oroc$auc
    Opauc <- opr$auc.integral
    Obrier <- mean((prorig-yy)^2)
    b_orig <- ece_mcefun(yy, prorig, g=10)
    Oece <- b_orig$ece
    Omce <- b_orig$mce
  }
  selected <- orig$coef
  lambda <- orig$lambda
  index.orig <- c(Oauc,Opauc,Obrier, Oece, Omce)
  training <- c(auc,pauc,brier, ece, mce)
  test <- c(auct,pauct,briert, ecet, mcet)
  data.a <- data.frame(index.orig,training,test)
  data.a$optimism <- training-test
  data.a$index.corrected <- index.orig-data.a$optimism
  data.a$n <- rep(mboot-ck, length(training))
  data.a <- round(data.a,digits=4)
  xx <- c("auc","pauc","brier", "ece", "mce")
  rownames(data.a) <- xx
  boot <- list(coef = selected, lambda=lambda, nonconvergence=ck,resu=data.a,auc_r =auc_r,auctest_r =auctest_r,
               pauc_r =pauc_r,pauctest_r=pauctest_r,brier_r=brier_r,briertest_r=briertest_r,
               ece_r = ece_r, ecetest_r = ecetest_r, mce_r=mce_r, mcetest_r=mcetest_r)

  return(boot)
}


#pp <- HeckSelect(selection, outcome, data=AmEx, allowParallel = TRUE, penalty="ALASSO", Model="AMH",crit="bic")
#comb <- bootValidate(pp, data=AmEx, mboot=2, seed=1)





#########################################################################
# This function is based on the use of P-value to select variables
# in Binary Heckman selection model. Default P-value = 0.05
####################################################################

#' Title Variable selection using P-value in binary sample selection model
#'
#' @param selection selection equation
#' @param outcome outcome equation
#' @param data data matrix containing both the outcome and selection variables
#' @param alpha significance level: default is 0.05
#' @param ...
#'
#' @return class HeckPval containing variables selected via P-value and
#' the parameters supplied for the creation of the object. Function call is
#' also returned
#' @export
#'
#' @examples
#' res <- HeckPval(selection, outcome, data=data,alpha=0.05))
#'
HeckPval <- function(selection, outcome, data = sys.frame(sys.parent()),alpha = 0.05,...){
  if (match("pbivnorm",.packages(),0)==0) require(pbivnorm)
  if (match("GJRM",.packages(),0)==0) require(GJRM)

  if (!missing(data)) {
    if (!inherits(data, "environment") & !inherits(data,
                                                   "data.frame") & !inherits(data, "list")) {
      stop("'data' must be either environment, data.frame, or list (currently a ",
           class(data), ")")
    }
  }


  funcCall <- match.call(expand.dots = FALSE)

  mf <- model.frame(selection, data = data)
  YS <- model.response(mf, "numeric")
  XS <- model.matrix(selection, data = data)

  mf2 <- model.frame(outcome, data = data)
  YO <- model.response(mf2, "numeric")
  XO <- model.matrix(outcome, data = data)

  options(warn=-1)


  m2 <- gjrm(list(selection, outcome), data = data, BivD = "N",
             margins = c("probit", "probit"), Model = "BSS")

  bpar <- coef(m2, part = "full")
  dff <- length(bpar)
  df <- nrow(data)-dff
  covb <- m2$Vb
  se <- sqrt(diag(covb))
  tval <- bpar/ se
  p.value <- round(2*pt(-abs(tval), df=df), digits=5)
  kk <- ifelse(p.value > alpha, 0, bpar)

  NXS <- ncol(XS)
  NXO <- ncol(XO)
  ibetaS <- 1:NXS
  ibetaO <- seq(tail(ibetaS, 1) + 1, length = NXO)
  irho <- tail(ibetaO, 1) + 1
  NXSS <- NXS+1
  kk[1] <- bpar[1]; kk[NXSS] <- bpar[NXSS]; kk[irho] <- bpar[irho]
  betaS <-  kk[ibetaS]
  betaO <-  kk[ibetaO]
  rho <- kk[irho]

  n.rho <- tanh(rho)

  n.rho2 <- ifelse(n.rho > 0.99, 0.99, ifelse(n.rho < -0.99,-0.99,n.rho))

  i11 <- !(YS == 0)

  XS11.b <- drop(XS[i11, , drop = FALSE] %*% betaS)
  XO11.b <- drop(XO[i11, , drop = FALSE] %*% betaO)
  rho_rep <- rep(n.rho2, length(XS11.b))
  p.pred <- pbivnorm(XO11.b,XS11.b, rho_rep)/pnorm(XS11.b)
  yyf <- ifelse(i11==1, YO,NA)
  yselect <- as.vector(na.omit(yyf))

  result <- structure(list(call = funcCall, coef=kk, betaS=betaS, betaO=betaO, lpS=XS11.b,lpO=XO11.b,
                           pred=p.pred, yout=yselect,rho=n.rho2,alpha = alpha,selection=selection, outcome = outcome),class = "HeckPval")
  return(result)
}


#res <- HeckPval(selection, outcome, data=data,alpha=0.05))


#########################################################################
# Coefficient
#########################################################################

coef.HeckPval <- function(object,...){
  x <- object$coef
  return(x)
}



#####################################################################
# Predict function for object of class HeckPval
##########################################################################
#' Title Predict Method for HeckPval Fits
#'
#' @param object a fitted object of class inheriting from "HeckPval"
#' @param newdata optionally, a data frame in which to look for variables with which to predict. If omitted, the fitted linear predictors are used.
#' @param na.action function determining what should be done with missing values in newdata. The default is to predict NA
#' @param ...
#'
#' @return vector of predicted values for P(outcome=1|selection=1)
#' @export
#'
#' @examples
#' res <- HeckPval(selection, outcome, data=data, alpha=alpha)
#' predict(res)
#'
predict.HeckPval <- function(object, newdata=NULL, na.action = na.pass,...){
  if (match("pbivnorm",.packages(),0)==0) require(pbivnorm)
  if (missing(newdata)){
    pred <- object$pred
  } else {
    selection <- object$selection
    outcome <- object$outcome
    betaS <- object$betaS
    betaO <- object$betaO
    rho <- object$rho
    mf <- model.frame(selection, newdata)
    YS <- model.response(mf, "numeric")
    XS <- model.matrix(selection, newdata)

    mf2 <- model.frame(outcome, newdata)
    YO <- model.response(mf2, "numeric")
    XO <- model.matrix(outcome, newdata)

    i11 <- !(YS == 0)
    yyf <- ifelse(i11==1, YO,NA)
    yy <- as.vector(na.omit(yyf))

    XS11.b <- drop(XS[i11, , drop = FALSE] %*% betaS)
    XO11.b <- drop(XO[i11, , drop = FALSE] %*% betaO)
    rho_rep <- rep(rho, length(XS11.b))
    pred <- pbivnorm(XO11.b,XS11.b, rho_rep)/pnorm(XS11.b)
  }
  pred
}

#options(scipen=999)
#aa <- predict(res)



#######################################################################
## Bootstrap internal validation technique to correct for overoptimism
# in predictions - the alpha value is inherited from the object HeckPval.
# Note that this is different from the "bootValidate" as this is based
# on dropping variables whose values are greater than 5% from the model
########################################################################

#' Title Bootstrap validation for object of model "HeckPval"(only implemented for Normal error distribution)
#'
#' @param object a fitted object of class inheriting from "HeckPval"
#' @param data data matrix containing both the outcome and selection variables
#' @param mboot  number of bootstrap samples. 50 bootstrap samples are
#' used if not specified
#' @param seed integer value set for reproducibility
#' @param ...
#'
#' @return coefficients: same as in coef
#'
#' nonconvergence: the number of samples that did not converge
#' in the bootstrap sample
#'
#' alpha: significance level: validation is based on model fitting without variable selection when alpha=1
#'
#' resu: data matrix containing apparent, bootstrap, test and optimism
#' corrected performance measures
#'
#' Vectors of performance measures for each bootstrap sample are also returned
#'
#' @export
#'
#' @examples
#' res <- HeckPval(selection, outcome, data=data,alpha=0.05)
#' bootValidate_Pval(res, data=data, mboot=2, seed=1)
#'

bootValidate_Pval <- function(object, data, mboot=50, seed,...){
  if (match("PRROC",.packages(),0)==0) require(PRROC)
  set.seed(seed)
  selection <- object$selection; outcome <- object$outcome; alpha <- object$alpha

  auc <- pauc <- brier   <- ece <- mce <- c()
  auctest <- pauctest <- briertest <- ecetest <- mcetest <- c()

  # Looping to generate mboot bootstrap data
  data1 <- data
  bootdata <- list()
  for(i in 1:mboot)
  {

    iboot <- sample(1:nrow(data1), replace=TRUE)
    bootdata[[i]] <- data1[iboot,]

    options(warn=2)
    res <- try(HeckPval(selection, outcome, data=bootdata[[i]], alpha = alpha),silent=TRUE)
    if (any(class(res)=="try-error"))
    {
      auc[i] <- pauc[i] <-  brier[i] <- ece[i] <- mce[i] <-  NA
      auctest[i] <- pauctest[i] <- briertest[i] <- ecetest[i] <- mcetest[i] <- NA
    } else{

      p <- res$pred
      yboot <- res$yout
      fg <- p[yboot == 1]
      bg <- p[yboot == 0]
      roc <- roc.curve(scores.class0 = fg, scores.class1 = bg)
      pr <- pr.curve(scores.class0 = fg, scores.class1 = bg)
      auc[i] <- roc$auc
      pauc[i] <- pr$auc.integral
      brier[i] <- mean((p-yboot)^2)
      b_boot <- ece_mcefun(yboot, p, g=10)
      ece[i] <- b_boot$ece
      mce[i] <- b_boot$mce

      # predicting test data back in the original data

      yy <- object$yout
      prtest <- predict(res, data)
      fgt <- prtest[yy == 1]
      bgt <- prtest[yy == 0]

      roct <- roc.curve(scores.class0 = fgt, scores.class1 = bgt)
      prt <- pr.curve(scores.class0 = fgt, scores.class1 = bgt)
      auctest[i] <- roct$auc
      pauctest[i] <- prt$auc.integral
      briertest[i] <- mean((prtest-yy)^2)
      b_test <- ece_mcefun(yy, prtest, g=10)
      ecetest[i] <- b_test$ece
      mcetest[i] <- b_test$mce
    }
  }

  auc_r <- auc; auctest_r <- auctest; pauc_r <- pauc; pauctest_r <- pauctest
  brier_r <- brier; briertest_r <- briertest; ece_r <- ece; ecetest_r <- ecetest
  mce_r <- mce; mcetest_r <- mcetest

  ck <- sum(is.na(auc_r))
  auc <- mean(auc, na.rm=T); pauc <- mean(pauc, na.rm=T); brier <- mean(brier, na.rm=T)
  ece <- mean(ece, na.rm=T); mce <- mean(mce, na.rm=T)

  auct <- mean(auctest, na.rm=T); pauct <- mean(pauctest, na.rm=T)
  briert <- mean(briertest, na.rm=T); ecet <- mean(ecetest, na.rm=T)
  mcet <- mean(mcetest, na.rm=T)

  # Index measures from the orginal data set
  options(warn=2)
  orig <- try(object, silent=TRUE)

  if (any(class(orig)=="try-error"))
  {
    Oauc <- Opauc <-  Obrier <- Oece <- Omce <- NA
  } else{
    prorig <- orig$pred
    yy <- orig$yout
    fgo <- prorig[yy == 1]
    bgo <- prorig[yy == 0]
    oroc <- roc.curve(scores.class0 = fgo, scores.class1 = bgo)
    opr <- pr.curve(scores.class0 = fgo, scores.class1 = bgo)
    Oauc <- oroc$auc
    Opauc <- opr$auc.integral
    Obrier <- mean((prorig-yy)^2)
    b_orig <- ece_mcefun(yy, prorig, g=10)
    Oece <- b_orig$ece
    Omce <- b_orig$mce
  }
  selected <- orig$coef
  lambda <- orig$lambda
  index.orig <- c(Oauc,Opauc,Obrier, Oece, Omce)
  training <- c(auc,pauc,brier, ece, mce)
  test <- c(auct,pauct,briert, ecet, mcet)
  data.a <- data.frame(index.orig,training,test)
  data.a$optimism <- training-test
  data.a$index.corrected <- index.orig-data.a$optimism
  data.a$n <- rep(mboot-ck, length(training))
  data.a <- round(data.a,digits=4)
  xx <- c("auc","pauc","brier", "ece", "mce")
  rownames(data.a) <- xx
  boot <- list(coef = selected, alpha=alpha, nonconvergence=ck,resu=data.a,auc_r =auc_r,auctest_r =auctest_r,
               pauc_r =pauc_r,pauctest_r=pauctest_r,brier_r=brier_r,briertest_r=briertest_r,
               ece_r = ece_r, ecetest_r = ecetest_r, mce_r=mce_r, mcetest_r=mcetest_r)

  return(boot)
}

#res <- HeckPval(selection, outcome, data=dat,alpha=0.05)
#ss <- bootValidate_Pval(res, data=dat, mboot=2, seed=1)



#######################################################################
# Complete cases uses Probit regression. Varibale selection is similar
# to the GLMNET package
######################################################################


###########################################################################
# Coordinate descent algorithm  lasso/adlasso complete cases
##########################################################################



coord_select_complete <- function(x, y, lambda, penalty=c("LASSO","ALASSO"),...)
{
  x <- x; y <- y
  n <- dim(x)[1]; p <- dim(x)[2]

  if(penalty=="LASSO"){
    weight <- rep(1, each= p)
  }
  if(penalty=="ALASSO"){
    weight <- as.vector(1/abs(glm(y~x-1)$coef))
  }

  lambda1 <- lambda
  lambda2 <- lambda1*weight

  maxstep <- min(n, p)*500
  beta <- matrix(NA, maxstep+1, p)
  beta[1, ] <- glm(y~x-1)$coef
  delta <- 1
  K <- 1
  while((delta>1e-10) & (K<maxstep))
  {
    K <- K+1
    beta.temp <- beta[K-1, ]  #working downwards through rows, starting with lm coeffs

    for (j in 1:p)
    {
      xminusj <- x[, -j]         #eliminate jth column from data matrix
      bminusj <- beta.temp[-j]
      yminusj <- xminusj%*%bminusj #fitted values (for each nrow)
      rminusj <- y-yminusj         #target - fitted values

      a11 <- softrule(sum(x[, j]*rminusj), lambda=0)/sum(x[, j]^2)
      a12 <- softrule(sum(x[, j]*rminusj), lambda=lambda2[j])/sum(x[, j]^2)
      bj <- ifelse((j==1),a11,a12)

      beta.temp[j] <- bj
    }
    beta[K, ] <- beta.temp
    delta <- max(abs(beta[K, ]-beta[K-1, ]))
  }
  beta <- beta[K, ]
  beta <- as.matrix(beta)
  df <- sum(beta !=0)-1
  if(is.null(colnames(x))){colnames(x)=c(paste("x", 1:p, sep=""))}
  rownames(beta) <- colnames(x)

  object <- list(beta=beta, lambda=lambda,df = df, K=K, delta=delta)
  return(object)
}




#######################################################################
#Pseudo Data Generation complete
#####################################################################
pseudo_data_complete <- function(formula, data){
  options(warn=-1)
  mle <- glm(formula,family = binomial(link = "probit"),data = data)
  beta0 <- coef(mle)
  names0 <- names(beta0)
  p <- length(beta0)
  covb <- vcov(mle)
  colnames(covb) <- rownames(covb) <- names0

  e <- eigen(covb)
  aa <- e$vectors %*% diag(1/sqrt(abs(e$values))) %*% t(e$vectors)
  bb <- e$vectors %*% diag(1/sqrt(makePD(covb))) %*% t(e$vectors)
  sigma2 <-  ifelse((det(covb) > 0),list(aa),list(bb))
  xstar <- sigma2[[1]]
  colnames(xstar) <- colnames(covb)
  ystar <- xstar%*%as.vector(beta0)
  XY <- list(xstar,ystar)
  names(XY) <- c('xstar','ystar')
  return(XY)
}

#hh <- pseudo_data_complete(formula, dat)



####################################################################
#Likelihood function probit regression
#######################################################################

logLik <- function(beta, XS, YS){
  -sum(YS*log(pnorm(XS%*%beta)) + (1-YS)*log(pnorm(-XS%*%beta)))
}




###########################################################################
# The model for predictions. Note that the probability of interest
# is P(Y=1|S=1), which can be derived using the ratio of the joint distribution
# and marginal distribution. This is Regularized probit regression
# for complete data
######################################################################

#' Title Lasso and Adaptive Lasso Probit Regression
#'
#' @param formula  model formula for probit regression. This is the outcome model when the missing data is deleted
#' @param data data matrix containing the outcome and covariates for probit regression
#' @param lambda shrinkage parameter, both scalar and vector are acceptable
#' When lambda=NULL, the internal vector of Lambda is used
#' @param allowParallel If true, the "doParallel" package is invoked
#' @param penalty can be ALASSO (for adaptive lasso) or LASSO (for Lasso) penalty
#' @param crit can be BIC, AIC or GCV, default is BIC
#' @param ...
#'
#' @return class ProbitLasso containing penalized probit coefficients
#' similar to the GLMNET package results. Linear predictor and fitted values are returned

#' @export
#'
#' @examples
#' ProbitLasso(formula, data=data, allowParallel = TRUE, penalty="ALASSO", crit="bic")
#'


ProbitLasso <- function(formula, data = sys.frame(sys.parent()), lambda = NULL, allowParallel = FALSE,
                    penalty=c("LASSO","ALASSO"), crit=c("bic","aic","gcv"),...)
{
  if (match("foreach",.packages(),0)==0) require(foreach)
  if (allowParallel) library(doParallel)

  if (!missing(data)) {
    if (!inherits(data, "environment") & !inherits(data,
                                                   "data.frame") & !inherits(data, "list")) {
      stop("'data' must be either environment, data.frame, or list (currently a ",
           class(data), ")")
    }
  }

  mf <- model.frame(formula, data=data)
  YS <- model.response(mf, "numeric")
  XS <- model.matrix(formula, data = data)

  pseudo <- pseudo_data_complete(formula, data)
  xstar <- pseudo$xstar
  ystar <- pseudo$ystar
  crit <- crits <- match.arg(crit)
  penalty  <-  match.arg(penalty)
  n <- length(YS)

  if (allowParallel) {
    `%op%` <- `%dopar%`
    cl <- makeCluster(detectCores())
    clusterExport(cl, c("coord_select_complete","softrule",
                        "logLik","makePD"))
    registerDoParallel(cl)
  } else {
    `%op%` <- `%do%`
  }

  if(is.null(lambda)) {

    lambda <- seq(0,10, length=50)
  }

  fitter <- function(xstar, ystar, lambda, penalty){
    fit <-  coord_select_complete(xstar, ystar, lambda, penalty)
    beta <- fit$beta
    df  <- fit$df
    fn <- logLik(as.matrix(beta), XS = XS, YS = YS)
    bic <- 2*fn + log(n)*df
    aic <- 2*fn +2*df
    gcv <- fn/(n*(1-df/n)^2)

    out <- list(beta = beta, df = df, bic= bic, aic = aic, gcv = gcv)
    return(out)
  }

  comb <- function(...) {
    mapply('cbind', ..., SIMPLIFY=FALSE)
  }

  H  <- foreach(i = lambda, .combine='comb', .multicombine=TRUE,
                .verbose = FALSE, .errorhandling = "stop") %op% fitter(xstar, ystar, i, penalty)

  if (allowParallel) stopCluster(cl)

  beta <- H$beta
  df <- as.vector(H$df)
  bic <- as.vector(H$bic)
  aic <-  as.vector(H$aic)
  gcv <-  as.vector(H$gcv)


  crit <- switch(crit, bic=bic, aic=aic, gcv=gcv)
  selected <- best.model <- which(crit == min(crit,na.rm=TRUE))
  ic <- c(bic=bic[selected],aic=aic[selected], gcv=gcv[selected])
  coef <- beta[,selected]

  lambda_f <- min(lambda[selected])
  lp <- as.vector(XS %*% coef)
  p.pred  <- pnorm(lp)
  result <- list(coef=coef, lp = lp, pred=p.pred, lambda = lambda_f, formula = formula,
                 yout= YS, crit=crits, penalty = penalty, allowParallel = allowParallel)

  class(result) <- "ProbitLasso"
  return(result)
}




#####################################################################
# Predict function for ProbitLasso object
##########################################################################

predict.ProbitLasso <- function(object, newdata=NULL,...){
  if (missing(newdata)){
    pred <- object$pred
  } else {

    formula <- object$formula
    coef <- object$coef

    mf <- model.frame(formula, newdata)
    YS <- model.response(mf, "numeric")
    XS <- model.matrix(formula, newdata)

    lp <- as.vector((XS %*% coef))
    pred <-  pnorm(lp)

  }
  pred
}





#######################################################################
## Bootstrap internal validation technique to correct for overoptimism
# in predictions - mboot is the number of bootstrap samples.
# We fix BIC here for the probit regression
########################################################################

#' Title Bootstrap validation for probit regression with Lasso penalty
#'
#' @param object a fitted object of class inheriting from "ProbitLasso"
#' @param data data matrix containing both the outcome and selection variables
#' @param mboot number of bootstrap samples. 50 bootstrap samples are
#' used if not specified
#' @param seed integer value set for reproducibility
#' @param ...
#'
#' @return coefficients: same as in coef
#'
#' nonconvergence: the number of samples that did not converge
#' in the bootstrap sample
#'
#' lambda: optimal shrinkage parameter
#'
#' resu: data matrix containing apparent, bootstrap, test and optimism
#' corrected performance measures
#'
#' Vectors of performance measures for each bootstrap sample are also returned
#'
#' @export
#'
#' @examples
#' pp2 <- ProbitLasso(formula, data=data, allowParallel = TRUE, penalty="ALASSO", crit="bic")
#' boot_ProbitLasso(pp2, data=data, mboot=2, seed=1)
#'


boot_ProbitLasso <- function(object, data, mboot=50, seed,...){
  if (match("PRROC",.packages(),0)==0) require(PRROC)
  set.seed(seed)
  formula <- object$formula;  penalty <- object$penalty
  crit <- object$crit; allowParallel <- object$allowParallel

  auc <- pauc <- brier   <- ece <- mce <- c()
  auctest <- pauctest <- briertest <- ecetest <- mcetest <- c()

  # Looping to generate m bootstrap data
  data1 <- data
  bootdata <- list()
  for(i in 1:mboot)
  {
    iboot <- sample(1:nrow(data1), replace=TRUE)
    bootdata[[i]] <- data1[iboot,]

    options(warn=2)
    res <- try(ProbitLasso(formula, data=bootdata[[i]],allowParallel, penalty=penalty, crit=crit), silent=TRUE)
    if (any(class(res)=="try-error"))
    {
      auc[i] <- pauc[i] <-  brier[i] <- ece[i] <- mce[i] <-  NA
      auctest[i] <- pauctest[i] <- briertest[i] <- ecetest[i] <- mcetest[i] <- NA
    } else{

      p <- res$pred
      yboot <- res$yout
      fg <- p[yboot == 1]
      bg <- p[yboot == 0]
      roc <- roc.curve(scores.class0 = fg, scores.class1 = bg)
      pr <- pr.curve(scores.class0 = fg, scores.class1 = bg)
      auc[i] <- roc$auc
      pauc[i] <- pr$auc.integral
      brier[i] <- mean((p-yboot)^2)
      b_boot <- ece_mcefun(yboot, p, g=10)
      ece[i] <- b_boot$ece
      mce[i] <- b_boot$mce

      # predicting test data back in the original data

      yy <- object$yout
      prtest <- predict(res, data)
      fgt <- prtest[yy == 1]
      bgt <- prtest[yy == 0]

      roct <- roc.curve(scores.class0 = fgt, scores.class1 = bgt)
      prt <- pr.curve(scores.class0 = fgt, scores.class1 = bgt)
      auctest[i] <- roct$auc
      pauctest[i] <- prt$auc.integral
      briertest[i] <- mean((prtest-yy)^2)
      b_test <- ece_mcefun(yy, prtest, g=10)
      ecetest[i] <- b_test$ece
      mcetest[i] <- b_test$mce
    }
  }

  auc_r <- auc; auctest_r <- auctest; pauc_r <- pauc; pauctest_r <- pauctest
  brier_r <- brier; briertest_r <- briertest; ece_r <- ece; ecetest_r <- ecetest
  mce_r <- mce; mcetest_r <- mcetest

  ck <- sum(is.na(auc_r))
  auc <- mean(auc, na.rm=T); pauc <- mean(pauc, na.rm=T); brier <- mean(brier, na.rm=T)
  ece <- mean(ece, na.rm=T); mce <- mean(mce, na.rm=T)

  auct <- mean(auctest, na.rm=T); pauct <- mean(pauctest, na.rm=T)
  briert <- mean(briertest, na.rm=T); ecet <- mean(ecetest, na.rm=T)
  mcet <- mean(mcetest, na.rm=T)

  # Index measures from the orginal data set
  options(warn=2)
  orig <- try(object, silent=TRUE)

  if (any(class(orig)=="try-error"))
  {
    Oauc <- Opauc <-  Obrier <- Oece <- Omce <- NA
  } else{
    prorig <- orig$pred
    yy <- orig$yout
    fgo <- prorig[yy == 1]
    bgo <- prorig[yy == 0]
    oroc <- roc.curve(scores.class0 = fgo, scores.class1 = bgo)
    opr <- pr.curve(scores.class0 = fgo, scores.class1 = bgo)
    Oauc <- oroc$auc
    Opauc <- opr$auc.integral
    Obrier <- mean((prorig-yy)^2)
    b_orig <- ece_mcefun(yy, prorig, g=10)
    Oece <- b_orig$ece
    Omce <- b_orig$mce
  }
  selected <- orig$coef
  lambda <- orig$lambda
  index.orig <- c(Oauc,Opauc,Obrier, Oece, Omce)
  training <- c(auc,pauc,brier, ece, mce)
  test <- c(auct,pauct,briert, ecet, mcet)
  data.a <- data.frame(index.orig,training,test)
  data.a$optimism <- training-test
  data.a$index.corrected <- index.orig-data.a$optimism
  data.a$n <- rep(mboot-ck, length(training))
  data.a <- round(data.a,digits=4)
  xx <- c("auc","pauc","brier", "ece", "mce")
  rownames(data.a) <- xx
  boot <- list(coef = selected, lambda=lambda, nonconvergence=ck,resu=data.a,auc_r =auc_r,auctest_r =auctest_r,
               pauc_r =pauc_r,pauctest_r=pauctest_r,brier_r=brier_r,briertest_r=briertest_r,
               ece_r = ece_r, ecetest_r = ecetest_r, mce_r=mce_r, mcetest_r=mcetest_r)

  return(boot)
}



