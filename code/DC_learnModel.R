rm(list=ls(all=TRUE))
library(data.table)
library(survival)
library(DMwR)
library(glmnet)
library(caret)
library(ranger)
library(survivalROC)
library(gbm)
library(caTools)
#library(bnlearn)
library(e1071)
library(PresenceAbsence)
library(FSelector)
library(corrgram)

#source("H:\\projects\\data_challenge_2017\\code\\DC_learnFunctions.R")
#source("H:\\Work\\R\\DataWranglingFunctions.R")
source("~/projects/data_challenge_2017/code/DC_learnFunctions.R")
source("~/work/MyWork/R/DataWranglingFunctions.R")



#dt1 <- fread("H:\\projects\\data_challenge_2017\\data\\data.csv")
dt1 <- fread("~/projects/data_challenge_2017/data/data.csv")
dt1<-dt1[,-c("ID","V1")]

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
vPredictorVarsSVM <- c("A","W5","W6","W7","W8","W9","W10","W11","W12","W13",
                       "W14","W15","W16","W17","W18","W19","W20","W21","W22",
                       "W23","W24","W25","W1.A","W1.B","W1.C","W1.D","W2.A",
                       "W2.B","W2.C","W2.D","W3.A","W3.B","W3.C","W4.A",
                       "W4.B","W4.C")

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
                             bOversample=F)
#vModelBN <- funcCrossValBN(dtIn=dt1F.imp, 
#                           vDependentVars=vDependentVars,
#                           vPredictorVars=vPredictorVarsBN,
#                           nFolds=10,
#                           bOversample=T)
vModelSVM <- funcCrossValSVM(dtIn=dt1N, 
                             vDependentVars=vDependentVars,
                             vPredictorVars=vPredictorVarsSVM,
                             nFolds=10,
                             bOversample=F)
vModelRF <- funcCrossValRF(dtIn=dt1, 
                           vDependentVars=vDependentVars,
                           vPredictorVars=vPredictorVars,
                           nFolds=10,
                           bOversample=F)
vModelGBM <- funcCrossValGBM(dtIn=dt1, 
                             vDependentVars=vDependentVars,
                             vPredictorVars=vPredictorVars,
                             nFolds=10,
                             bOversample=F)



################################################################################
############################### Model Evaluation ###############################
################################################################################
plot(vModelReg$trueValues,vModelReg$predicted); cor.test(vModelReg$trueValues,vModelReg$predicted); sqrt(mean((vModelReg$trueValues-vModelReg$predicted)^2))
plot(vModelSVM$trueValues,vModelSVM$predicted); cor.test(vModelSVM$trueValues,vModelSVM$predicted); sqrt(mean((vModelSVM$trueValues-vModelSVM$predicted)^2))
plot(vModelRF$trueValues,vModelRF$predicted); cor.test(vModelRF$trueValues,vModelRF$predicted); sqrt(mean((vModelRF$trueValues-vModelRF$predicted)^2))
plot(vModelGBM$trueValues,vModelGBM$predicted); cor.test(vModelGBM$trueValues,vModelGBM$predicted); sqrt(mean((vModelGBM$trueValues-vModelGBM$predicted)^2))


