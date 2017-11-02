################################################################################
################################## Functions ###################################
################################################################################


funcCrossValRF = function(dtIn, vDependentVars,vPredictorVars,nFolds, bOversample = T) {
  ######################### Cross-validate random forest #########################
  # Splits the data into multiple samples for cross-validation. At each fold, a
  # BayesNet is built. If there are too many independent variables for
  # the number of observations supplied, tikhonov regularization can be used to
  # reduce the dimensionality of the dataset.
  #
  # Inputs:
  # dt: data.table of the data
  # vDependentVars: list of column names in dt that are dependent variables
  # vPredictorsVars: list of column names in dt that may explain dependent
  #                     variables. Can be a list of list for multiple levels.
  #                     List has to be in reverse temporal order. The first
  #                     element are potential parents of the dependent variables
  #                         vPredictorVars[[n]][1...k]
  #                                   ...
  #                         vPredictorVars[[2]][1...k]
  #                         vPredictorVars[[1]][1...k]
  #                               vDependentVars
  # nFolds: num of subsamples
  # valAlpha: type 1 error threshold
  # bPenalize: boolean used to do ridge regression
  # bFitLocal: fit local neighborhood around each node
  #
  # Outputs:
  # bnOut: bn object which includes the network
  ################################################################################
  
  ### initialize empty vectors/matrices for later use
  nObs <- nrow(dtIn)
  
  # initialize intermediate variables
  # initialize outputs
  causal <- rep(0, times = nObs)
  probs <- rep(0, times = nObs)
  vTest <- rep(0, times = nObs)
  vRF.predicted <- vector("list",nFolds)
  
  # shuffle the data to get unbiased splits
  vTestInd <- createFolds(dtIn[,get(vDependentVars)], k = nFolds, list = TRUE)
  indStart <- 1
  indEnd <- 0
  
  cForm <- as.formula(paste(vDependentVars," ~ ",paste(unique(unlist(vPredictorVars)), collapse = "+")))
  
  for (kFold in seq(nFolds)) {
    print(c("KFOLD: ",kFold))
    ### shuffle the data to get unbiased splits
    indSamp <- vTestInd[[kFold]]
    indStart <- indEnd + 1
    indEnd <- indEnd + length(indSamp)
    
    ### split training and test
    dt.training = dtIn[-indSamp,]
    if (bOversample) {
      nObsTrain <- as.integer(nrow(dt.training)/2)
      perOver <- nObsTrain/table(dt.training[,get(vDependentVars)])[2]*100
      #perUnder <- table(dt.training[,get(vDependentVars)])[1]/nObsTrain*100
      cFormula <- as.formula(paste(vDependentVars, " ~ ",
                                   paste(unlist(vPredictorVars),collapse ="+")))
      dt.training <- SMOTE(form = cForm,
                           data = dt.training,
                           perc.over=perOver,
                           perc.under=100)
    }
    dt.test = dtIn[indSamp,]
    ### fit the model on the training data
    rf.model <- ranger(cForm, data = dt.training, num.trees=300,
                       importance="impurity", probability = T)
    
    ###predict each dependent variable, given all the parents
    cPredicted <- predict(rf.model, data=dt.test)#
    vPredicted <- rep(colnames(cPredicted$predictions)[1],nrow(dt.test))
    vPredicted[cPredicted$predictions[,1]<cPredicted$predictions[,2]] <- "yes"
    
    vRF.predicted[[kFold]] <- rf.model
    #vTestInd[[kFold]] <- dt.test.this[,]
    causal[indStart:indEnd] <- vPredicted
    probs[indStart:indEnd] <- cPredicted$predictions[,2]
    vTest[indStart:indEnd] <- dt.test[,get(vDependentVars)]
    
  }
  
  return(list(predicted = causal,
              probabilities = probs,
              #predcor = predcor,
              models = vRF.predicted,
              trueValues = vTest))
}
############################ End of funcCrossValRF #############################



funcCrossValGBM = function(dtIn,vDependentVars,vPredictorVars,nFolds, bOversample = T) {
  ############## Cross-validate gradient boosted regression models ###############
  # Splits the data into multiple samples for cross-validation. At each fold, a
  # GBM is built. If there are too many independent variables for
  # the number of observations supplied, tikhonov regularization can be used to
  # reduce the dimensionality of the dataset.
  #
  # Inputs:
  # dt: data.table of the data
  # vDependentVars: list of column names in dt that are dependent variables
  # vPredictorsVars: list of column names in dt that may explain dependent
  #                     variables. Can be a list of list for multiple levels.
  #                     List has to be in reverse temporal order. The first
  #                     element are potential parents of the dependent variables
  #                         vPredictorVars[[n]][1...k]
  #                                   ...
  #                         vPredictorVars[[2]][1...k]
  #                         vPredictorVars[[1]][1...k]
  #                               vDependentVars
  # nFolds: num of subsamples
  # bFitLocal: fit local neighborhood around each node
  #
  # Outputs:
  # bnOut: bn object which includes the network
  ################################################################################
  
  ### initialize empty vectors/matrices for later use
  nObs <- nrow(dtIn)
  
  # initialize intermediate variables
  # initialize outputs
  causal <- rep(0, nrow = nObs)
  probs <- rep(0, nrow = nObs)
  vTest <- rep(0, nrow = nObs)
  vGBM.predicted <- vector("list",nFolds)
  
  # shuffle the data to get unbiased splits
  vTestInd <- createFolds(dtIn[,get(vDependentVars)], k = nFolds, list = TRUE)
  indStart <- 1
  indEnd <- 0
  
  cForm <- as.formula(paste(vDependentVars," ~ ",paste(unique(unlist(vPredictorVars)), collapse = "+")))
  
  #if (nlevels(dtIn[,get(vDependentVars)])==2) {
  #  dtIn[,(vDependentVars):=as.integer(dtIn[,get(vDependentVars)])]
  #  set(dtIn,which(dtIn[,get(vDependentVars)]==1),vDependentVars,0)
  #  set(dtIn,which(dtIn[,get(vDependentVars)]==2),vDependentVars,1)
  #}
  
  for (kFold in seq(nFolds)) {
    print(c("KFOLD: ",kFold))
    ### shuffle the data to get unbiased splits
    indSamp <- vTestInd[[kFold]]
    indStart <- indEnd + 1
    indEnd <- indEnd + length(indSamp)
    # create a matrix to store the predicted values
    pred <- rep(0, nrow = length(indSamp))
    trueV <- rep(0, nrow = length(indSamp))
    
    ### split training and test
    dt.training = dtIn[-indSamp,]
    if (bOversample) {
      nObsTrain <- as.integer(nrow(dt.training)/2)
      perOver <- nObsTrain/table(dt.training[,get(vDependentVars)])[2]*100
      #perUnder <- table(dt.training[,get(vDependentVars)])[1]/nObsTrain*100
      dt.training <- SMOTE(form = cForm,
                           data = dt.training,
                           perc.over=perOver,
                           perc.under=100)
    }
    dt.test = dtIn[indSamp,]
    ### fit the model on the training data
    gbm.model <- gbm(cForm, data = dt.training,
                     distribution="gaussian",n.trees=1000,
                     shrinkage=0.05,interaction.depth=3,
                     bag.fraction = 0.5,train.fraction=1.0,
                     n.minobsinnode=10,cv.folds=0,
                     keep.data=F,verbose=F,n.cores=1)
    
    ###predict each dependent variable, given all the parents
    ind.Best <- gbm.perf(gbm.model,method="OOB")
    
    vGBM.predicted[[kFold]] <- gbm.model
    causal[indStart:indEnd] <- predict(gbm.model, dt.test,ind.Best)
    vTest[indStart:indEnd] <- dt.test[,get(vDependentVars)]
    
  }
  
  return(list(predicted = causal,
              models = vGBM.predicted,
              trueValues = vTest))
}
############################ End of funcCrossValRF #############################



funcCrossValSVM = function(dtIn,vDependentVars,vPredictorVars,nFolds, bOversample = T) {
  ################# Cross-validate support vector machine models #################
  # Splits the data into multiple samples for cross-validation. At each fold, a
  # SVM is built. If there are too many independent variables for
  # the number of observations supplied, tikhonov regularization can be used to
  # reduce the dimensionality of the dataset.
  #
  # Inputs:
  # dt: data.table of the data
  # vDependentVars: list of column names in dt that are dependent variables
  # vPredictorsVars: list of column names in dt that may explain dependent
  #                     variables. Can be a list of list for multiple levels.
  #                     List has to be in reverse temporal order. The first
  #                     element are potential parents of the dependent variables
  #                         vPredictorVars[[n]][1...k]
  #                                   ...
  #                         vPredictorVars[[2]][1...k]
  #                         vPredictorVars[[1]][1...k]
  #                               vDependentVars
  # nFolds: num of subsamples
  # bFitLocal: fit local neighborhood around each node
  #
  # Outputs:
  # bnOut: bn object which includes the network
  ################################################################################
  
  ### initialize empty vectors/matrices for later use
  nObs <- nrow(dtIn)
  
  # initialize intermediate variables
  # initialize outputs
  causal <- rep(0, nrow = nObs)
  probs <- rep(0, nrow = nObs)
  vTest <- rep(0, nrow = nObs)
  vModels <- vector("list",nFolds)
  
  # shuffle the data to get unbiased splits
  vTestInd <- createFolds(dtIn[,get(vDependentVars)], k = nFolds, list = TRUE)
  indStart <- 1
  indEnd <- 0
  
  cForm <- as.formula(paste(vDependentVars," ~ ",paste(unique(unlist(vPredictorVars)), collapse = "+")))
  
  #if (nlevels(dtIn[,get(vDependentVars)])==2) {
  #  dtIn[,(vDependentVars):=as.integer(dtIn[,get(vDependentVars)])]
  #  set(dtIn,which(dtIn[,get(vDependentVars)]==1),vDependentVars,0)
  #  set(dtIn,which(dtIn[,get(vDependentVars)]==2),vDependentVars,1)
  #}
  
  for (kFold in seq(nFolds)) {
    print(c("KFOLD: ",kFold))
    ### shuffle the data to get unbiased splits
    indSamp <- vTestInd[[kFold]]
    indStart <- indEnd + 1
    indEnd <- indEnd + length(indSamp)
    # create a matrix to store the predicted values
    pred <- rep(0, nrow = length(indSamp))
    trueV <- rep(0, nrow = length(indSamp))
    
    ### split training and test
    dt.training = dtIn[-indSamp,]
    if (bOversample) {
      nObsTrain <- as.integer(nrow(dt.training)/2)
      perOver <- nObsTrain/table(dt.training[,get(vDependentVars)])[2]*100
      #perUnder <- table(dt.training[,get(vDependentVars)])[1]/nObsTrain*100
      dt.training <- SMOTE(form = cForm,
                           data = dt.training,
                           perc.over=perOver,
                           perc.under=100)
    }
    dt.test = dtIn[indSamp,]
    ### fit the model on the training data
    svm_model <- svm(cForm, data=dt.training,
                     method="C-classification",
                     kernel="linear",
                     probability=T)
    
    ###predict each dependent variable, given all the parents
    vModels[[kFold]] <- svm_model
    causal[indStart:indEnd] <- as.vector(predict(svm_model, dt.test,probability=T))
    probs[indStart:indEnd] <- as.vector(attr(predict(svm_model, dt.test,probability=T),'probabilities')[,2])
    vTest[indStart:indEnd] <- dt.test[,get(vDependentVars)]
    
  }
  
  return(list(predicted = causal,
              probabilities = probs,
              models = vModels,
              trueValues = vTest))
}
############################ End of funcCrossValSVM ############################



funcCrossValReg = function(dtIn, vDependentVars, vPredictorVars,
                           nFolds = 10, bOversample = T) {
  ########################### Cross-validate bayes net ###########################
  # Splits the data into multiple samples for cross-validation. At each fold, a
  # BayesNet is built. If there are too many independent variables for
  # the number of observations supplied, tikhonov regularization can be used to
  # reduce the dimensionality of the dataset.
  #
  # Inputs:
  # dt: data.table of the data
  # vDependentVars: list of column names in dt that are dependent variables
  # vPredictorsVars: list of column names in dt that may explain dependent
  #                     variables. Can be a list of list for multiple levels.
  #                     List has to be in reverse temporal order. The first
  #                     element are potential parents of the dependent variables
  #                         vPredictorVars[[n]][1...k]
  #                                   ...
  #                         vPredictorVars[[2]][1...k]
  #                         vPredictorVars[[1]][1...k]
  #                               vDependentVars
  # nFolds: num of subsamples
  # valAlpha: type 1 error threshold
  # bPenalize: boolean used to do ridge regression
  # bFitLocal: fit local neighborhood around each node
  #
  # Outputs:
  # bnOut: bn object which includes the network
  ################################################################################
  
  
  ### initialize empty vectors/matrices for later use
  nObs <- nrow(dtIn)
  vAllVars <- c(vDependentVars, unlist(vPredictorVars))
  vPredictorVars <- unique(unlist(vPredictorVars))
  
  #dtIn[,(vDependentVars):=as.factor(dtIn[,get(vDependentVars)])]
  
  # initialize outputs
  causal <- rep(0, times = nObs)
  probs <- rep(0, times = nObs)
  vTest <- rep(0, times = nObs)
  vModels.predicted <- vector("list",nFolds)
  
  
  # shuffle the data to get unbiased splits
  vTestInd <- createFolds(dtIn[,get(vDependentVars)], k = nFolds, list = TRUE)
  indStart <- 1
  indEnd <- 0
  
  ############################## cross validation ##############################
  for (kFold in seq_len(nFolds)) {
    print(c("KFOLD: ",kFold))
    ### shuffle the data to get unbiased splits
    indSamp <- vTestInd[[kFold]]
    indStart <- indEnd + 1
    indEnd <- indEnd + length(indSamp)
    
    ### split training and test
    dt.training = dtIn[-indSamp,]
    if (length(vPredictorVars)==1) {
      cForm <- as.formula(paste(vDependentVars,' ~ ', vPredictorVars))
    } else {
      cForm <- as.formula(paste(vDependentVars,' ~ ', paste(unique(unlist(vPredictorVars)),collapse = " + ")))
    }
    if (bOversample) {
      nObsTrain <- as.integer(nrow(dt.training)/2)
      perOver <- nObsTrain/table(dt.training[,get(vDependentVars)])[2]*100
      #perUnder <- table(dt.training[,get(vDependentVars)])[1]/nObsTrain*100
      dt.training <- SMOTE(form = cForm,
                           data = dt.training,
                           perc.over=perOver,
                           perc.under=100)
    }
    dt.test = dtIn[indSamp,]
    ### fit the model on the training data,
    possibleError <- tryCatch({
      cModel <- glm(formula = cForm, data = dt.training, family = "gaussian")
    }, error = function(e) {
      print("Error glm")
    })
    
    if (inherits(possibleError, "error")) {
      causal[indStart:indEnd,] <- numeric(nrow(dt.test))
      vTest[indStart:indEnd,] <- dt.test[,get(vDependentVars)]
      cModel <- NULL
      next
    }
    
    
    ###predict each dependent variable, given all the parents
    vModels.predicted[[kFold]] <- cModel
    #vTestInd[[kFold]] <- dt.test.this[,]
    causal[indStart:indEnd] <- predict(cModel, newdata=dt.test)
    vTest[indStart:indEnd] <- dt.test[,get(vDependentVars)]
  } #end of kFold crossval
  
  ### overall cross-validated correlations
  #for (kDepVar in vDependentVars) {
  #  if (class(dt[,get(kDepVar)])=="factor") {
  #    dfTmp <- data.frame(a1=causal[,kDepVar], a2=dt[vIndObs, get(kDepVar)])
  #    tab <- xtabs(~a1+a2, data = dfTmp)
  #    predcor[kDepVar] <- vcd::assocstats(tab)$cramer
  #  } else {
  #    predcor[kDepVar] <- cor(causal[,kDepVar], dt[vIndObs, get(kDepVar)])
  #  }
  #}
  
  return(list(predicted = causal,
              models = vModels.predicted,
              trueValues = vTest))
} #end of funcCrossVal



funcSingleLogReg = function(dtIn, vDependentVars, vPredictorVars) {
  ##################### Single Variable Logistic Regression ######################
  #
  ################################################################################
  vModelsOut <- vector("list",length(vPredictorVars))
  names(vModelsOut) <- vPredictorVars
  
  for (kVar in vPredictorVars) {
    print(kVar)
    
    cForm <- as.formula(paste(vDependentVars,' ~ ', kVar))
    ### fit the model on the training data,
    possibleError <- tryCatch({
      vModelsOut[[kVar]] <- cModel <- glm(formula = cForm, data = dtIn, family = "binomial")
    }, error = function(e) {
      print("Error glm")
      vModelsOut[[kVar]] <- NA
    })
    
  } #end of kFold crossval
  
  return(vModelsOut)
}



funcStatOutput <- function(listModel, valPositive=T, valNegative=F,
                           bNumeric=T, threshType=3, valThresh=NULL) {
  ############################# Outputs Statistics  ##############################
  # Outputs "PPV","NPV","Sensitivity","accuracy","specificity","AUC" from model
  #
  # Inputs:
  # positive_label: label that is the positive outcome of scores
  # scores: vector predicted scores from model
  # truth: vector of truth
  # threshold: threshold of scores that predicts positive
  #
  # Outputs:
  # vOut: a vector of length (nrows(m))
  ################################################################################
  vProb <- listModel$predicted
  vPred <- rep(NA,times=length(vProb))
  vT <- listModel$trueValues
  
  if (threshType!='self') {
    vID <- as.character(seq(length(vProb)))
    vT2 <- rep(0,times=length(vProb))
    vT2[vT== valPositive] <- 1
    vProbNorm <- copy(vProb)
    cMin <- min(vProbNorm,na.rm=T)
    vProbNorm <- vProbNorm-cMin
    cMax <- max(vProbNorm,na.rm=T)
    vProbNorm <- vProbNorm/cMax
    
    dtPred <- data.table(id=vID, truth=vT2, prob=vProbNorm)
    dtPred <- dtPred[which(!is.nan(dtPred[,prob])),]
    
    valThresh <- optimal.thresholds(DATA=dtPred)[threshType,2]*cMax + cMin
  }
  
  
  if (bNumeric) {
    vPred[vProb>valThresh] <- valPositive
    vPred[!(vProb>valThresh)] <- valNegative
  } else {
    vPred <- vProb
  }
  
  
  nTP <- sum((vPred == valPositive) & (vT == valPositive),na.rm=T)
  nFN <- sum((vPred == valNegative) & (vT == valPositive),na.rm=T)
  nFP <- sum((vPred == valPositive) & (vT == valNegative),na.rm=T)
  nTN <- sum((vPred == valNegative) & (vT == valNegative),na.rm=T)
  
  valPPV <- nTP/(nTP + nFP)
  valNPV <- nTN/(nTN + nFN)
  valSensitivity <- nTP/sum(vT==valPositive)
  valAccuracy <- (nTP + nTN) / (length(vPred))
  valSpecificity <- nTN/sum(vT==valNegative)
  valF <- 2/(1/valSensitivity + 1/valPPV)
  valAUC <- colAUC(vProb,vT, plotROC=F)
  
  
  vOut <- c(valSensitivity, valSpecificity,valPPV,valNPV,valAccuracy,valF,valAUC)
  names(vOut) <- c("Sensitivity","specificity","Positive.P.V","Negative.P.V","Accuracy","F-score","AUC")
  
  return(vOut)
  
}


############################### End of Functions ###############################
