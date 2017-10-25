rm(list=ls(all=TRUE))
library(data.table)
#library(survival)
#library(DMwR)
library(glmnet)
library(caret)
#library(ranger)
#library(survivalROC)
#library(gbm)
#library(caTools)
#library(bnlearn)
#library(e1071)
library(PresenceAbsence)
library(FSelector)
library(corrgram)

#source("H:\\projects\\data_challenge_2017\\code\\.R")
source("H:\\Work\\R\\DataWranglingFunctions.R")


dt1 <- fread("H:\\projects\\data_challenge_2017\\data\\data.csv")


################################################################################
################################## Variables ###################################
################################################################################
vDependentVars <- "Y" 
vPredictorVars <- c("A","W1","W2","W3","W4","W5","W6","W7","W8","W9","W10","W11",
                    "W12","W13","W14","W15","W16","W17","W18","W19","W20","W21",
                    "W22","W23","W24","W25")
vFac <- c("A","W1","W2","W3","W4","W5","W6")
vNum <- c("W7","W8","W9","W10","W11","W12","W13","W14","W15","W16","W17","W18",
          "W19","W20","W21","W22","W23","W24","W25")

### Factorize data
for (kVar in vFac) {
  dt1[,(kVar):=as.factor(dt1[,get(kVar)])]
}

### normalize and binarize
dt1N <- funcNormalize(dtIn=dt1,
                      vNum=vNum,
                      bHard=F)
dt1N <- funcBinaryFactorize(dtIn=dt1N,
                          vFac=vFac)


################################################################################
############################### Basic Statistics ###############################
################################################################################
funcOddsRatio(dt1[,3:29,with=F], dt1[,A]==1)
vGainRatios <- gain.ratio(as.formula(paste(vDependentVars,"~ .")),dt1)

vCoefEstimate <- numeric(length(vPredictorVars))
vPValue <- numeric(length(vPredictorVars))
k <-1
for (kVar in vPredictorVars) {
  cForm <- as.formula(paste(vDependentVars,"~",kVar))
  cModel <- glm(formula = cForm, family = gaussian, data = dt1)
  vCoefEstimate[k] <- coef(summary(cModel))[2,1]
  vPValue[k] <- coef(summary(cModel))[2,4]
  k <- k+1
}
a <- data.table(variable=vPredictorVars,coef=vCoefEstimate,significance=vPValue)
print(a)

### Plot some stuff
#vInfoGain <- information.gain(as.formula(paste(vDependentVars,"~ .")),dt1.imp)
#p <- ggplot(data=data.frame(names=names(dt2)[1:68],gain.ratio=vGainRatios[,"attr_importance"]), 
#            aes(x=names, y=gain.ratio)) +
#  geom_bar(stat="identity")
#p + theme(axis.text.x = element_text(angle=90, hjust=1))


### correlations of the continuous data
corrgram(dt1, order=TRUE, lower.panel=panel.ellipse,
         upper.panel=panel.pts, text.panel=panel.txt,
         diag.panel=panel.minmax, 
         main="Simulated Data in PC2/PC1 Order")



################################################################################
################################ Model Learning ################################
################################################################################
vModelReg <- funcCrossValReg(dtIn=dt1, 
                             vDependentVars=vDependentVars,
                             vPredictorVars=vPredictorVars,
                             nFolds=10,
                             bOversample=T)
#vModelBN <- funcCrossValBN(dtIn=dt1F.imp, 
#                           vDependentVars=vDependentVars,
#                           vPredictorVars=vPredictorVarsBN,
#                           nFolds=10,
#                           bOversample=T)
vModelSVM <- funcCrossValSVM(dtIn=dt1N, 
                             vDependentVars=vDependentVars,
                             vPredictorVars=vPredictorVarsSVM,
                             nFolds=10,
                             bOversample=T)
vModelRF <- funcCrossValRF(dtIn=dt1, 
                           vDependentVars=vDependentVars,
                           vPredictorVars=vPredictorVars,
                           nFolds=10,
                           bOversample=T)
vModelGBM <- funcCrossValGBM(dtIn=dt1, 
                             vDependentVars=vDependentVars,
                             vPredictorVars=vPredictorVars,
                             nFolds=10,
                             bOversample=T)



################################################################################
############################### Model Evaluation ###############################
################################################################################
vModelBN$predicted <- vModelBN$probabilities
vModelRF$predicted <- vModelRF$probabilities
vModelSVM$predicted <- vModelSVM$probabilities

vReg.Stat <- funcStatOutput(listModel=vModelReg,valPositive=2,valNegative=1,bNumeric=T,threshType=3)#,valThresh=0.748431)
vBN.Stat <- funcStatOutput(listModel=vModelBN,valPositive=2,valNegative=1,bNumeric=T,threshType=3)
vSVM.Stat <- funcStatOutput(listModel=vModelSVM,valPositive=2,valNegative=1,bNumeric=T,threshType=3)#,valThresh=0.5209497)
vRF.Stat <- funcStatOutput(listModel=vModelRF,valPositive=2,valNegative=1,bNumeric=T,threshType=3)#,valThresh=0.5298393)
vGBM.Stat <- funcStatOutput(listModel=vModelGBM,valPositive=2,valNegative=1,bNumeric=T,threshType=3)#,valThresh=1.513645)
