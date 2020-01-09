rm(list=ls())
library(survival)
library(gbm)
library(glmnet)
library(randomForestSRC)
library(CoxBoost)
library(mboost)
library(ggplot2)
#########################################

#Convex combination of exsiting algorithms
convex.cox = function(Y, X, delta, iter.max = 10){
  iter = 1
  eps = 0.000001
  censor.indx = which(delta == 0)
  failure.indx = which(delta == 1)
  censor.X = X[censor.indx, , drop = FALSE]
  failure.X = X[failure.indx, , drop = FALSE]
  n = length(Y)  
  q = ncol(X)
  t = sort(unique(Y[failure.indx]))  # unique failure times in increasing order
  m = length(unique(Y[failure.indx])) 
  D = lapply(1:m, function(i) which(Y == t[i] & delta==1))
  d = sapply(D, length)
  R = lapply(1:m, function(i) which(Y >= t[i])) #risk set
  C = lapply(1:n, function(i) which(t < Y[i]))
  l = function(beta){
    (sum(sapply(1:m, function(i) sum(beta * X[D[[i]], ]) - d[i]*
                  log(sum(exp(beta*X[R[[i]], , drop=FALSE]))))))
  }
  
  beta.tilda = rep(0, q)
  eta.tilda = beta.tilda %*% t(X)
  
  while(iter < iter.max){
    beta_old = beta.tilda
    w = -sapply(1:n, function(i) 
      if(length(C[[i]])>0){
        sum(sapply(C[[i]], function(k) {
          (exp(eta.tilda[i])*sum(exp(eta.tilda[R[[k]]])) - 
             exp(eta.tilda[i])^2)*d[k] / sum(exp(eta.tilda[R[[k]]]))^2
        }))
      } else{0}
    )
    
    z = sapply(1:n, function(i) 
      if(length(C[[i]])>0){
        eta.tilda[i] - 1/w[i]*(delta[i] - 
                                 exp(eta.tilda[i])*sum(sapply(C[[i]], function(k){
                                   d[k]/sum(exp(eta.tilda[R[[k]]]))
                                 })) )
      } else {0} 
      
    )      
    
    beta.hat = rep(0, q) #initialize beta.hat
    for (k in 1:q){
      if (k == 1){
        beta = beta.tilda
      } else{
        beta = c(beta.hat[1:(k-1)], beta.tilda[k:q])
      }
      numerator = as.vector((w*(z-rowSums(sweep(X[,-k],2, beta[-k], FUN='*'))))%*%
                              X[ , k,drop=FALSE])
      denominator = w %*% X[,k, drop=FALSE]^2
      if(is.na(numerator)| is.na(denominator)){
        beta.hat[k]=0.001
      }
      else{
        b= numerator/denominator
        beta.hat[k] = ifelse(b < 0, 0.001, b)
      }
    }
    beta.tilda = beta.hat/sum(beta.hat)
    eta.tilda = beta.tilda %*% t(X)
    relative_diff = (abs(l(beta.tilda)-l(beta_old))/(abs(l(beta_old))+eps))
    if (relative_diff<0.000001)
      break
    iter = iter + 1
    if (iter == iter.max )
      message(paste('WARNING: algorithm did not converge in', iter.max, 'iterations.'))
  }
  beta.tilda = ifelse(beta.tilda<0.001,0, beta.tilda)
  return (list (eta = eta.tilda, iter=iter, init =beta.hat,coef =beta.tilda, coef_old =beta_old, diff =relative_diff))
}

############################################################
convex.cox_pair = function(Y, X, delta, iter.max = 10){
  eps = 0.000001
  Q = ncol(X)
  cor_mat = abs(cor(X))
  selected = as.numeric(sort(which(cor_mat == min(cor_mat), arr.ind = TRUE)[1,]))
  Orig.X = X
  X.bar = X[, -selected]
  X = X[, selected]
  censor.indx = which(delta == 0)
  failure.indx = which(delta == 1)
  failure.X = X[failure.indx, , drop = FALSE]
  n = length(Y)  
  q = ncol(X)
  t = sort(unique(Y[failure.indx]))  # unique failure times in increasing order
  m = length(unique(Y[failure.indx])) 
  D = lapply(1:m, function(i) which(Y == t[i] & delta==1))
  d = sapply(D, length)
  R = lapply(1:m, function(i) which(Y >= t[i]))
  C = lapply(1:n, function(i) which(t < Y[i]))
  l = function(beta){
    (sum(sapply(1:m, function(i) sum(beta * X[D[[i]], ]) - d[i]*
                  log(sum(exp(beta*X[R[[i]], , drop=FALSE]))))))
  }
  
  # iinitialize beta.tilda and eta.tilda
  beta.tilda = c(0,1)
  eta.tilda = beta.tilda %*% t(X)
  
  beta.sequence = c()
  selected.sequence = selected
  
  for (ii in seq(Q-1)){
    iter = 1
    while(iter < iter.max){
      beta_old =beta.tilda
      w = -sapply(1:n, function(i) 
        if(length(C[[i]])>0){
          sum(sapply(C[[i]], function(k) {
            (exp(eta.tilda[i])*sum(exp(eta.tilda[R[[k]]])) - 
               exp(eta.tilda[i])^2)*d[k] / sum(exp(eta.tilda[R[[k]]]))^2
          }))
        } else{0}
      )
      
      z = sapply(1:n, function(i) 
        
        if(length(C[[i]])>0){
          eta.tilda[i] - 1/w[i]*(delta[i] - 
                                   exp(eta.tilda[i])*sum(sapply(C[[i]], function(k){
                                     d[k]/sum(exp(eta.tilda[R[[k]]]))
                                   })) )
        } else {0} 
        
      )      
      
      if (isTRUE(all.equal(X[,1],X[,2]))){
        beta1.hat =0.001
      }else{
        numerator = sum(w*(X[, 1] - X[, 2])*(z-X[,2]))
        denominator = sum(w*(X[, 1] - X[, 2])^2)
        beta1.hat = numerator/denominator
      }
      
      if (beta1.hat >0 & beta1.hat<1)
        beta1.tilda = beta1.hat
      
      if (beta1.hat < 0)
        beta1.tilda = 0.001
      
      if(beta1.hat>1)
        beta1.tilda = 0.999
      
      
      beta2.tilda = 1 - beta1.tilda
      beta.tilda = c(beta1.tilda, beta2.tilda)
      eta.tilda = X %*% beta.tilda
      
      relative_diff = (abs(l(beta.tilda)-l(beta_old))/(abs(l(beta_old))+eps))
      if (relative_diff<0.000001)
        break
      
      iter = iter + 1
      if (iter == iter.max )
        message(paste('WARNING: algorithm did not converge in', 
                      iter.max, 'iterations.'))
    }
    
    ii = ii + 1
    selected = which.min(abs(cor(eta.tilda, X.bar)))
    selected.sequence = c(selected.sequence, selected)
    beta.sequence = c(beta.sequence, beta.tilda)
    
    X = matrix(c(eta.tilda, X.bar[,selected]), ncol=2)
    X.bar = X.bar[, -selected, drop=FALSE]
  }
  
  index = seq(Q)
  index2 = selected.sequence[1:2]
  index = index[-index2]
  for (x in selected.sequence[3:Q]){
    index2 = c(index2, index[x])
    index = index[-x]
  }
  
  beta = rep(0,Q)
  beta[1] = beta.sequence[1]*prod(beta.sequence[seq(3, 2*Q-3,2)])
  beta[Q] = beta.sequence[length(beta.sequence)]
  for (i in 2:(Q-1)){
    beta[i] = beta.sequence[2*(i-1)]*prod(beta.sequence[seq(2*i-1, 2*(Q-1),2)])
  }
  
  Match = match(seq(Q), index2)
  beta.tilda = beta[Match]  
  beta.tilda_new = ifelse(beta.tilda <0.0005, 0, beta.tilda)
  eta.tilda = Orig.X %*% beta.tilda
  eta.tilda_new = Orig.X %*% beta.tilda_new
  
  return (list (coef =beta.tilda, coef_new = beta.tilda_new, eta = eta.tilda, eta_new = eta.tilda_new, iter=iter) )
}
###############################################
convex.cox.pred = function (Y, X, delta, beta, X_new, times){
  censor.indx = which(delta == 0)
  failure.indx = which(delta == 1)
  censor.X = X[censor.indx, , drop = FALSE]
  failure.X = X[failure.indx, , drop = FALSE]
  n = length(Y)  
  q = ncol(X)
  t = sort(unique(Y[failure.indx])) 
  m = length(unique(Y[failure.indx])) 
  D = lapply(1:m, function(i) which(Y == t[i] & delta==1))
  d = sapply(D, length)
  R = lapply(1:m, function(i) which(Y >= t[i])) #risk set
  C = lapply(1:n, function(i) which(t < Y[i]))
  times_sorted = sort(unique(times)) 
  s = length(unique(times)) 
  p = lapply(1:n, function(i) which(times_sorted < Y[i]))
  #calculating baseline hazard at the predicted times
  basehaz = sapply(1:m, function(i){
    d[i]/sum(exp(beta*X[R[[i]], , drop=FALSE]))
    
  })
  pred = sapply(times_sorted, function(i){
    cumbasehaz = sum(basehaz[which(t<=i)]) 
    cumhaz = cumbasehaz*exp(X_new%*%beta)
    pred_surv = exp(-cumhaz)
    
  }) 
  
  return(pred)
}
################################################
c_indexCv_combined = function(data,k){
  y_dat = Surv(data$obs.time,data$status)
  set.seed(1)
  folds = sample(rep(1:k, length.out = nrow(data)))
  prediction_cox = c()
  prediction_lasso = c()
  prediction_ridge = c()
  prediction_net = c()
  prediction_coxboost =c()
  prediction_mboost =c()
  prediction_gbm = c()
  index =c()
  for (j in 1:k){
    idx = which(folds==j)
    train = data[-idx,]
    test = data[idx,]
    y_train = Surv(train$obs.time, train$status)
    y_test = Surv(test$obs.time,test$status)
    x_train = model.matrix(~., train[,-c(1,2)])
    x_test = model.matrix(~., test[,-c(1,2)])
    fit_cox = coxph(y_train~., data=train[,-c(1,2)])
    pred_cox = predict(fit_cox, test[,-c(1,2)], type="lp")  
    fit_lasso = glmnet(x_train,y_train, family="cox", alpha=1)
    cvFit_lasso = cv.glmnet(x_train,y_train, family="cox", alpha=1)
    pred_lasso = predict(fit_lasso,x_test, s=cvFit_lasso$lambda.min, type="link") 
    fit_ridge = glmnet(x_train,y_train, family="cox", alpha=0)
    cvFit_ridge = cv.glmnet(x_train,y_train, family="cox", alpha=0)
    pred_ridge = predict(fit_ridge,x_test, s=cvFit_ridge$lambda.min, type="link") 
    fit_net = glmnet(x_train,y_train, family="cox", alpha=0.5)
    cvFit_net = cv.glmnet(x_train,y_train, family="cox", alpha=0.5)
    pred_net = predict(fit_net,x_test, s=cvFit_net$lambda.min, type="link") 
    optimal.penalty = optimCoxBoostPenalty(time=train$obs.time, 
                                           status=train$status, x_train,        
                                           minstepno=50,maxstepno=200,start.penalty=9*sum(train$status==1),iter.max=10,
                                           upper.margin=0.05)
    fit_coxboost = CoxBoost(time= train$obs.time,status=train$status, x_train,
                            standardize=TRUE,subset=1:length(train$obs.time), 
                            penalty=optimal.penalty$penalty, 
                            stepno=optimal.penalty$cv.res$optimal.step) 
    pred_coxboost = predict(fit_coxboost, x_test , type="lp")  
    
    fit_mboost = glmboost(y_train~., data = train[,-c(1,2)],family=
                            CoxPH(),control=boost_control(mstop = 500))
    pred_mboost = predict(fit_mboost,test[,-c(1,2)], type="link")  
    
    fit_gbm = gbm(y_train~., data=train[,-c(1,2)], distribution="coxph",
                  n.trees=2000, shrinkage=0.01, interaction.depth=5 ,
                  bag.fraction = 0.5, cv.folds=5)
    iter = gbm.perf(fit_gbm,method="OOB")
    pred_gbm = predict(fit_gbm,test[,-c(1,2)],n.trees=iter, type="link") 
    index = c(index,idx)
    prediction_cox = c(prediction_cox, pred_cox)
    prediction_lasso = c(prediction_lasso, pred_lasso)
    prediction_ridge = c(prediction_ridge, pred_ridge)
    prediction_net = c(prediction_net, pred_net)
    prediction_coxboost = c(prediction_coxboost, pred_coxboost)
    prediction_mboost = c(prediction_mboost, pred_mboost)
    prediction_gbm = c(prediction_gbm, pred_gbm)
  }
  
  Match = match(seq(nrow(data)), index)
  prediction_cox = prediction_cox[Match]
  prediction_lasso = prediction_lasso[Match]
  prediction_ridge = prediction_ridge[Match]
  prediction_net = prediction_net[Match]
  prediction_coxboost = prediction_coxboost[Match]
  prediction_mboost = prediction_mboost[Match]
  prediction_gbm = prediction_gbm[Match]
  c_cox = survConcordance(y_dat~prediction_cox)$concordance
  c_lasso = survConcordance(y_dat~prediction_lasso)$concordance
  c_ridge = survConcordance(y_dat~prediction_ridge)$concordance
  c_net = survConcordance(y_dat~prediction_net)$concordance
  c_coxboost = survConcordance(y_dat~prediction_coxboost)$concordance
  c_mboost = survConcordance(y_dat~prediction_mboost)$concordance
  c_gbm = survConcordance(y_dat~prediction_gbm)$concordance
  final_pred = cbind(prediction_cox, prediction_lasso,prediction_ridge, prediction_net, 
                     prediction_coxboost,prediction_mboost, prediction_gbm)
  return(list(pred = final_pred, c_index=c(c_cox, c_lasso, c_ridge, c_net,c_coxboost, c_mboost, c_gbm)))
}

#######
#prediction of all algorithm for a new patient
new_patient_pred = function(train, test){
  pred_new = c()
  y_train = Surv(train$obs.time, train$status)
  x_train = model.matrix(~., train[,-c(1,2)])
  x_test = model.matrix(~., test)
  fit_cox = coxph(y_train~., data=train[,-c(1,2)])
  pred_cox = predict(fit_cox, newdata=test, type="lp")  
  fit_lasso = glmnet(x_train,y_train, family="cox", alpha=1)
  cvFit_lasso = cv.glmnet(x_train,y_train, family="cox", alpha=1)
  pred_lasso = predict(fit_lasso, x_test, s=cvFit_lasso$lambda.min, type="link") 
  fit_ridge = glmnet(x_train,y_train, family="cox", alpha=0)
  cvFit_ridge = cv.glmnet(x_train,y_train, family="cox", alpha=0)
  pred_ridge = predict(fit_ridge,x_test, s=cvFit_ridge$lambda.min, type="link") 
  fit_net = glmnet(x_train,y_train, family="cox", alpha=0.5)
  cvFit_net = cv.glmnet(x_train,y_train, family="cox", alpha=0.5)
  pred_net = predict(fit_net,x_test, s=cvFit_net$lambda.min, type="link") 
  optimal.penalty = optimCoxBoostPenalty(time=train$obs.time, 
                                         status=train$status, x_train,        
                                         minstepno=50,maxstepno=200,start.penalty=9*sum(train$status==1),iter.max=10,
                                         upper.margin=0.05)
  fit_coxboost = CoxBoost(time= train$obs.time,status=train$status, x_train,
                          standardize=TRUE,subset=1:length(train$obs.time), 
                          penalty=optimal.penalty$penalty, 
                          stepno=optimal.penalty$cv.res$optimal.step) 
  pred_coxboost = predict(fit_coxboost, x_test , type="lp")  
  
  fit_mboost = glmboost(y_train~., data = train[,-c(1,2)],family=
                          CoxPH(),control=boost_control(mstop = 500))
  pred_mboost = predict(fit_mboost,test, type="link")  
  
  fit_gbm = gbm(y_train~., data=train[,-c(1,2)], distribution="coxph",
                n.trees=2000, shrinkage=0.01, interaction.depth=5 ,
                bag.fraction = 0.5, cv.folds=5)
  iter = gbm.perf(fit_gbm,method="OOB")
  pred_gbm = predict(fit_gbm,test,n.trees=iter, type="link") 
  
  pred_new = c(pred_cox, pred_lasso,pred_ridge, pred_net, 
               pred_coxboost,pred_mboost, pred_gbm)
  return(pred_new)
}


#lung cancer data from survival package
names(lung)
naCount = sapply(lung, function(i) sum(is.na(i)))
lung= subset(lung, !is.na(ph.ecog) & !is.na(ph.karno) &!is.na(pat.karno)& !is.na(meal.cal) &!is.na(wt.loss) &!is.na(inst))
class_lung = sapply(lung, function(i) class(i))
head(lung)
lung$status = lung$status-1
options(na.action="na.pass")
lung1 = lung[,c(2, 3,1, 4, 5, 6, 7, 8, 9, 10)]
colnames(lung1)[1] ="obs.time"
dim(lung1)
table(lung1$status)
lung_c_idx = c_indexCv_combined(lung1,10) 
new_lung = cbind(lung1$obs.time,lung1$status, lung_c_idx$pred)
new_lung = as.data.frame(new_lung)
names(new_lung)= c("obs.time", "status","cox", "lasso", "ridge", "net", 
                   "coxboost", "mboost", "gbm")
X= as.matrix(new_lung[,3:9])
Y= as.matrix(new_lung[,1] ) 
delta = as.matrix(new_lung[,2])
beta1 = convex.cox(Y, X, delta, iter.max = 100)$coef
head(new_lung)
####################
#calculate the mode
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

##################################
#test = cbind(inst=1, age = 75, sex=1, ph.ecog =1, ph.karno=80, pat.karno=80, meal.cal=1500, wt.loss =10)
test = cbind(inst=1, age = 75, sex=1, ph.ecog =2, ph.karno=60, pat.karno=70, meal.cal=1500, wt.loss =10)
test = as.data.frame(test)
X_new = new_patient_pred (lung1, test)

#needs to be change prediction from each method seperately
#times = c(364, 728, 1092, 1456, 1820)
times = seq(from= 60, to=720, by=60)

survP= convex.cox.pred (Y, X, delta, beta1, X_new, times) #this work perfectly 

time = seq(from=0, to=24, by=2)

Dat = cbind(time, c(1,survP))
Dat = as.data.frame(Dat)

ggplot(data=Dat, aes(x=time, y=V2)) +
  geom_line()+scale_y_continuous(breaks = seq(0, 1, 0.1))+scale_x_continuous(breaks = seq(0, 24, 2))+
  xlab("Time in months") + ylab("Predicted survival probability")+
  ggtitle("24-month survival probability for the hypothetical patient")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.margin = unit(c(1,1,1,1), "cm"))
  

#this is for testing
# time_surv= Surv(new_lung$obs.time, new_lung$status)
# fit_cox = coxph(time_surv~., data=new_lung)
# pred_cox = predict(fit_cox, X_new, type="lp") 
# 




