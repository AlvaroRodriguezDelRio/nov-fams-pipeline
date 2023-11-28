library(ggplot2)
library(tidyr)
library(dplyr)
library(FactoMineR)
library(data.table)
library(caret)
library(reshape2)
library(vegan)
library(patchwork)
library(ecodist)
library(MASS)
library(glmnet)
library(pROC)


# 70 -30 train test ml function (logistic regression) 
train_test = function(data, samples){
  
  data_ml = merge(data, samples, by = "sample")
  data_ml$Group = as.factor(data_ml$Group)
  data_ml$sample = NULL
  data_ml$population = NULL
  
  aucs = list()
  accs = list()
  n = 0
  for (j in 1:50) {
    set.seed(n)
    
    # get train and test data
    validation_index = createDataPartition(data_ml$Group, p=0.30, list=FALSE)
    train <- data_ml[-validation_index,]
    test <- data_ml[validation_index,]
    
    yt <- as.factor(train$Group)
    train$Group = NULL
    xt <- train
    
    ytest <- as.factor(test$Group)
    test$Group = NULL
    xtest <- test
    
    # fit model
    cvglmfit  <- cv.glmnet(as.matrix(xt), yt,family = "binomial")
    s0 = cvglmfit$lambda.min
    fit = glmnet(as.matrix(xt), yt,family = "binomial")
    pred = predict(fit,newx = as.matrix(xtest),type = "class", s = s0)
    
    # get AUC
    auc = assess.glmnet(fit, newx = as.matrix(xtest), newy = ytest,s = s0)
    
    # get accuracy
    cm = confusionMatrix(table(pred,ytest))  
    acc = cm$overall[1]
    
    # store data
    n = n + 1
    aucs[n] = auc$auc
    accs[n] = acc
    
    
  }
  r = as.data.frame(cbind(aucs,accs))
  names(r) = c("auc","acc")
  r$auc = as.numeric(r$auc)
  r$acc = as.numeric(r$acc)
  return(r)
}

# train test 70 - 30%, with random forest
train_test_rf = function(data,samples){
  
  data_ml = merge(data, samples, by = "sample")
  data_ml$Group = as.factor(data_ml$Group)
  data_ml$sample = NULL
  data_ml$population = NULL
  aucs = list()
  accs = list()
  n = 0
  for (j in 1:20) {
    
    set.seed(n)
    
    # random sample 1000 kos / novel fams 
    data_ml_sub = data_ml[,sample(ncol(data_ml), 1000)]
    data_ml_sub$Group = data_ml$Group
    data_ml = data_ml_sub
    
    # get train and test data
    validation_index = createDataPartition(data_ml$Group, p=0.30, list=FALSE)
    train <- data_ml[validation_index,]
    test <- data_ml[-validation_index,]
    
    xt <- train
    xtest <- test
    
    # train model
    control <- trainControl(method="cv", number=10)
    metric <- "Accuracy"
    fit.rf <- train(Group~., data=xt, method="rf", metric=metric, trControl=control)
    
    # calculate accuracy
    predictions <- predict(fit.rf, xtest)
    results = confusionMatrix(predictions, as.factor(xtest$Group))
    acc = results$overall[1]
    
    # calculate AUC
    predictions <- predict(fit.rf, xtest,  type="prob")
    auc = roc(xtest$Group, predictions$CTR)
    
    # combine in table
    n = n + 1
    aucs[n] = auc$auc
    accs[n] = acc
    
  }
  r = as.data.frame(cbind(aucs,accs))
  names(r) = c("auc","acc")
  r$auc = as.numeric(r$auc)
  r$acc = as.numeric(r$acc)
  return(r)
}

####
# load data
####

# load abundance data
data = read.csv("Abs_per_family.csv",sep = ',',header = T)


# load metadata
samples = read.table("metadata.tsv",header = T)


#####
# logistic regression 70% train - 30% test
#####


# run model
ttnfam = train_test(nfam_data,samples)

#####
# Random forest 30% train - 70% test, 1000 random nfams
#####


# run model
ttnfam = train_test_rf(nfam_data,samples)


