
rm(list = ls())
dev.off()
library("FFTrees")
library("rpart")
library("rpart.plot")
library(tidyverse)
library(caret)
library(dplyr)
library(ggplot2)
library(PRROC)
library(ROCit)

#obtains mean cue used of a cart decision tree
mcu_cart = function(tree){
  node <- as.numeric(row.names(tree$frame))
  depth <- floor(log(node, 2))
  
  t = table(tree$where)
  d = as.data.frame(t)
  d$Var1 = as.integer(as.character(d$Var1))
  cues = c()
  for (i in c(1:nrow(d)) ){
    cues  = c(cues, rep(depth[d$Var1[i]],d$Freq[i]) )
  }
  return(mean(cues))
}

# functions for signal detection analysis
sens = function(table){
  sens = table[4]/(table[4]+table[2])
  return(sens)}

spec = function(table){
  spec = table[1]/(table[1]+table[3])
  return(spec)}


dprime <- function(table_mat) {
  hit = table_mat[4]/ (table_mat[4] + table_mat[2])
  fa = table_mat[3]/ (table_mat[3] + table_mat[1])
    return(qnorm(hit) - qnorm(fa))
}

beta <- function(table_mat) {
  hit = table_mat[4]/ (table_mat[4] + table_mat[2])
  fa = table_mat[3]/ (table_mat[3] + table_mat[1])
  zhr <- qnorm(hit)
  zfar <- qnorm(fa)
  return(exp(-zhr*zhr/2+zfar*zfar/2))
}

bacc = function(table){
  bacc = 0.5*sens(table) + 0.5 * spec(table) 
  return(bacc)}

# some example datasets
datasets = c("titanic", "cmc", "yeast", "wine", "breastcancer")

# creating blank dataframe to store results
result = data.frame(fit = character(), 
                     tr = integer(), 
                    data  = character(),
                    bacc = integer(),
                    mcu= integer(),
                    dprime = integer(),
                    beta = integer())

# running 10 simulations on each dataset
# performs FFT, CART, and logistic regression on each dataset

for (sim in c(1:10)){
  for (dInd in c(1:length(datasets))){
    
    # reads and parses data
    dataname = paste(datasets[dInd], ".csv", sep = "")
    data = read.csv(dataname, sep = ",")
    data = subset(data, select = -c(X))
    
    # dividing into training and testing datsets
    for (tInd in c(1:5)){
      tr = c(30, 50, 100, 200, 400)[tInd]
      
      dataTrain = sample(1:nrow(data), nrow(data))
      dataTrain = dataTrain[1:tr]
      
      dataTest = c(1:nrow(data))[-dataTrain]
      
  # creates FFT
      tree = FFTrees(formula = pred ~.,
                     data = data[dataTrain,],
                     data.test = data[dataTest,])
  
      #plot(tree, data = "test")
      #summary(tree)
      #plot(tree, what = "cues")
      predictions = predict(tree, data[dataTest,])
      table_mat = table(data$pred[dataTest], predictions)
      
      summ = summary(tree)
      mcu = summ$train[6]
      
      add = data.frame(fit = "fft", tr = tr, data = dataname, bacc = bacc(table_mat),
                       mcu = mcu, dprime = dprime(table_mat), beta = beta(table_mat))
      result = rbind(result, add)
      
      # creates CART tree
      
      tree2 = rpart(pred~., 
                    data = data[dataTrain,],
                    method = "class")
      #rpart.plot(tree2)
      #summary(tree2)
      #plotcp(tree2)
      #printcp(tree2)
      #plot(tree2, uniform=TRUE)
      #text(tree2, use.n=TRUE, all=TRUE, cex=.8)
  
      predict_unseen = predict(tree2, data[dataTest,], type = "class")
      table_mat = table(data$pred[dataTest], predict_unseen)
      
      add = data.frame(fit = "cart", tr = tr, data = dataname, bacc = bacc(table_mat),
                       mcu = mcu_cart(tree2), dprime = dprime(table_mat), beta = beta(table_mat))
      result = rbind(result, add)
      
  
      # performs logistic regression
      
      prediction =   logfit %>% predict(newdata = data[dataTest,], type = "response")
      predicted.classes <- ifelse(prediction > 0.5, TRUE, FALSE)
      table_mat = table(data$pred[dataTest], predicted.classes)

      
      add = data.frame(fit = "lr", tr = tr, data = dataname, bacc = bacc(table_mat),
                       mcu = 0, dprime = dprime(table_mat), beta = beta(table_mat))
      result = rbind(result, add)
      
    }
  }
}

# averages across simulations
 sum_result = group_by(result, fit, tr, data) %>%
  summarise(bacc = mean(bacc), dprime = mean(dprime), beta = mean(beta), mcu = mean(mcu))
