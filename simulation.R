#Plotting survival curve
library(survival)
library(glmnet)
library(randomForestSRC)
library(CoxBoost)
library(mboost)
library(gbm)
library(ggplot2)

#function to plot survival curve using ggplot2
ggsurv = function(s, CI = 'def', plot.cens = T, surv.col = 'gg.def',
                  cens.col = 'red', lty.est = 1, lty.ci = 2,
                  cens.shape = 3, back.white = F, xlab = 'Time',
                  ylab = 'Survival', main = ''){
  
  library(ggplot2)
  strata = ifelse(is.null(s$strata) ==T, 1, length(s$strata))
  stopifnot(length(surv.col) == 1 | length(surv.col) == strata)
  stopifnot(length(lty.est) == 1 | length(lty.est) == strata)
  
  ggsurv.s = function(s, CI = 'def', plot.cens = T, surv.col = 'gg.def',
                      cens.col = 'red', lty.est = 1, lty.ci = 2,
                      cens.shape = 3, back.white = F, xlab = 'Time',
                      ylab = 'Survival', main = ''){
    
    dat = data.frame(time = c(0, s$time),
                     surv = c(1, s$surv),
                     up = c(1, s$upper),
                     low = c(1, s$lower),
                     cens = c(0, s$n.censor))
    dat.cens = subset(dat, cens != 0)
    
    col = ifelse(surv.col == 'gg.def', 'black', surv.col)
    
    pl = ggplot(dat, aes(x = time, y = surv)) +
      xlab(xlab) + ylab(ylab) + ggtitle(main) +
      geom_step(col = col, lty = lty.est)
    
    pl = if(CI == T | CI == 'def') {
      pl + geom_step(aes(y = up), color = col, lty = lty.ci) +
        geom_step(aes(y = low), color = col, lty = lty.ci)
    } else (pl)
    
    pl = if(plot.cens == T & length(dat.cens) > 0){
      pl + geom_point(data = dat.cens, aes(y = surv), shape = cens.shape,
                      col = cens.col)
    } else if (plot.cens == T & length(dat.cens) == 0){
      stop ('There are no censored observations')
    } else(pl)
    
    pl = if(back.white == T) {pl + theme_bw()
    } else (pl)
    pl
  }
  
  ggsurv.m = function(s, CI = 'def', plot.cens = T, surv.col = 'gg.def',
                      cens.col = 'red', lty.est = 1, lty.ci = 2,
                      cens.shape = 3, back.white = F, xlab = 'Time',
                      ylab = 'Survival', main = '') {
    n = s$strata
    
    groups = factor(unlist(strsplit(names
                                    (s$strata), '='))[seq(2, 2*strata, by = 2)])
    gr.name =  unlist(strsplit(names(s$strata), '='))[1]
    gr.df = vector('list', strata)
    ind = vector('list', strata)
    n.ind = c(0,n); n.ind = cumsum(n.ind)
    for(i in 1:strata) ind[[i]] = (n.ind[i]+1):n.ind[i+1]
    
    for(i in 1:strata){
      gr.df[[i]] = data.frame(
        time = c(0, s$time[ ind[[i]] ]),
        surv = c(1, s$surv[ ind[[i]] ]),
        up = c(1, s$upper[ ind[[i]] ]),
        low = c(1, s$lower[ ind[[i]] ]),
        cens = c(0, s$n.censor[ ind[[i]] ]),
        group = rep(groups[i], n[i] + 1))
    }
    
    dat = do.call(rbind, gr.df)
    dat.cens = subset(dat, cens != 0)
    
    pl = ggplot(dat, aes(x = time, y = surv, group = group)) +
      xlab(xlab) + ylab(ylab) + ggtitle(main) +
      geom_step(aes(col = group, lty = group))
    
    col = if(length(surv.col == 1)){
      scale_colour_manual(name = gr.name, values = rep(surv.col, strata))
    } else{
      scale_colour_manual(name = gr.name, values = surv.col)
    }
    
    pl = if(surv.col[1] != 'gg.def'){
      pl + col
    } else {pl + scale_colour_discrete(name = gr.name)}
    
    line = if(length(lty.est) == 1){
      scale_linetype_manual(name = gr.name, values = rep(lty.est, strata))
    } else {scale_linetype_manual(name = gr.name, values = lty.est)}
    
    pl = pl + line
    
    pl = if(CI == T) {
      if(length(surv.col) > 1 && length(lty.est) > 1){
        stop('Either surv.col or lty.est should be of length 1 in order
             to plot 95% CI with multiple strata')
      }else if((length(surv.col) > 1 | surv.col == 'gg.def')[1]){
        pl + geom_step(aes(y = up, color = group), lty = lty.ci) +
          geom_step(aes(y = low, color = group), lty = lty.ci)
      } else{pl +  geom_step(aes(y = up, lty = group), col = surv.col) +
          geom_step(aes(y = low,lty = group), col = surv.col)}
    } else {pl}
    
    
    pl = if(plot.cens == T & length(dat.cens) > 0){
      pl + geom_point(data = dat.cens, aes(y = surv), shape = cens.shape,
                      col = cens.col)
    } else if (plot.cens == T & length(dat.cens) == 0){
      stop ('There are no censored observations')
    } else(pl)
    
    pl = if(back.white == T) {pl + theme_bw()
    } else (pl)
    pl
  }
  pl = if(strata == 1) {ggsurv.s(s, CI , plot.cens, surv.col ,
                                 cens.col, lty.est, lty.ci,
                                 cens.shape, back.white, xlab,
                                 ylab, main)
  } else {ggsurv.m(s, CI, plot.cens, surv.col ,
                   cens.col, lty.est, lty.ci,
                   cens.shape, back.white, xlab,
                   ylab, main)}
  pl
}

#C-index for exsiting models using cross validation (oracle and rsf are not included)
#Thay are calculated seperately since the data for oracle and the structure of rsf 
#are different 
c_indexCv_combined = function(data,k){
  y_dat = Surv(data$obs.time,data$status)
  set.seed(369)
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
                  bag.fraction = 0.5, cv.folds=10)
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
  return(list(pred = final_pred, c_index =c(c_cox, c_lasso, c_ridge, c_net,c_coxboost, c_mboost, c_gbm)))
}
```


```{r eval=FALSE, message=FALSE}
#c_index for RSF
c_indexCv_rsf = function(data,k,ntree=500){
  set.seed(369)
  folds = sample(rep(1:k, length.out = nrow(data)))
  c_idx_cv = vector(length=k)
  for (j in 1:k){
    idx = which(folds==j)
    train = data[-idx,]
    test = data[idx,]
    v.obj = rfsrc(Surv(obs.time,status)~.,data=train,ntree=ntree,tree.err=T,
                  splitrule="logrankscore", nodesize=5, err.block = 1)
    pred = predict(v.obj, test, err.block=1)
    c_idx_cv[j] = 1-mean(pred$err.rate)
  }
  return(mean(c_idx_cv, na.rm=T))
}

#c_index for oracle model
c_indexCv_oracle =function(data,k){
  set.seed(369)
  y_dat = Surv(data$obs.time,data$status)
  folds = sample(rep(1:k, length.out = nrow(data)))
  prediction =c()
  index =c()
  for (j in 1:k){
    idx = which(folds==j)
    train = data[-idx,]
    test = data[idx,]
    y_train = Surv(train$obs.time, train$status)
    y_test = Surv(test$obs.time,test$status)
    fit_cox = coxph(y_train~., data=train[,-c(1,2)])
    pred = predict(fit_cox, test[,-c(1,2)], type="lp")
    prediction =c(prediction, pred)
    index = c(index, idx)
  }
  Match = match(seq(nrow(data)), index)
  prediction = prediction[Match]
  c_index = survConcordance(y_dat~prediction)$concordance
  return(c_index)
}

#winsorization function
winsor =function (x, multiple=2)
{
  med = median(x)
  y = x - med
  sc = mad(y, center=0) * multiple
  y[ y > sc ] = sc
  y[ y < -sc ] = -sc
  r.time = y + med
}

#Generating low dimensional time to event data
set.seed(456)
n = 500
p = 50
beta = c(rep(0.3,10),rep(-0.5,10),rep(0,p-20)) 
x = matrix(rnorm(n*p),n,p)
A = -(log(runif(n)))
real.time = A/(0.2*exp(drop(x %*% beta)))
real.time = winsor(real.time,2)
cens.time = rexp(n,rate=0.18)
cens.time = winsor(cens.time,2)
status = ifelse(real.time <= cens.time,1,0)
obs.time = ifelse(real.time <= cens.time,real.time,cens.time)
z = as.data.frame(x)
survData = cbind(obs.time,status,z)
surv_dat= survfit(Surv(obs.time,status) ~ 1, data = survData)
pl = ggsurv(surv_dat)
pl

#Generating low dimensional time to event data only with informative covariates
beta1 = c(rep(0.3,10),rep(-0.5,10))
x1 = x[,1:20]
real.time = A/(0.2*exp(drop(x1 %*% beta1)))
real.time = winsor(real.time,2)
status = ifelse(real.time <= cens.time,1,0)
obs.time = ifelse(real.time <= cens.time,real.time,cens.time)
z = as.data.frame(x1)
survData1 = cbind(obs.time,status,z)

#c_indices from different models when linear assumption holds

c_index1_linear =c_indexCv_combined(survData,10)

c_index2_linear = c_indexCv_rsf (survData,10,ntree=1000)

c_index3_linear = c_indexCv_oracle(survData1,10)

newCovs = c_indexCv_combined(survData,10)
newdata = cbind( survData$obs.time,survData$status, newCovs$pred)
newdata = as.data.frame(newdata)
names(newdata)= c("obs.time", "status","cox", "lasso", "ridge", "net", 
                  "coxboost", "mboost", "gbm")
X= as.matrix(newdata[,3:9])
Y= as.matrix(newdata[,1] ) 
delta = as.matrix(newdata[,2])

#generating low dimensional time to event data with some quadratic effects
set.seed(456)
n = 500
p = 50
beta = c(rep(0.2,10),rep(-0.3,10),rep(0,p-20), rep(0.01,3), rep(-0.01,2))
x = matrix(rnorm(n*p),n,p)
s= sample(1:20, size=5)
y = apply(as.matrix(x[,s], drop=F), 2, function(i) i^2)
z = cbind(x,y)
A=-(log(runif(n))) 
real.time = A/(0.2*exp(drop(z %*% beta)))
real.time = winsor(real.time,2)
cens.time = rexp(n,rate=0.18)
cens.time = winsor(cens.time,2)
status = ifelse(real.time <= cens.time,1,0)
obs.time = ifelse(real.time <= cens.time,real.time,cens.time)
z = as.data.frame(z)
survData = cbind(obs.time,status,z)
surv_dat= survfit(Surv(obs.time,status) ~ 1, data = survData)
pl = ggsurv(surv_dat)
pl

#Generating low dimensional data with quadratic effcts only with informative covariates
beta1 = c(rep(0.2,10),rep(-0.3,10),rep(0.01,3), rep(-0.01,2))
x1 = cbind(x[,1:20],y)
real.time = A/(0.2*exp(drop(x1 %*% beta1)))
real.time = winsor(real.time,2)
status = ifelse(real.time <= cens.time,1,0)
obs.time = ifelse(real.time <= cens.time,real.time,cens.time)
z1 = as.data.frame(x1)
survData1 = cbind(obs.time,status,z1)

#c_indices from different models when quadratic terms are present
c_index1_quad =c_indexCv_combined(survData,10)

c_index2_quad = c_indexCv_rsf(survData,10,ntree=500)

c_index3_quad = c_indexCv_oracle(survData1,10)

newCovs1 = c_indexCv_combined(survData,10)
newdata1 = cbind( survData$obs.time,survData$status, newCovs1$pred)
newdata1 = as.data.frame(newdata1)
names(newdata1)= c("obs.time", "status","cox", "lasso", "ridge", "net", "coxboost", "mboost", "gbm")


X= as.matrix(newdata1[,3:9])
Y= as.matrix(newdata1[,1] ) 
delta = as.matrix(newdata1[,2])

#Gererating low dimensional time to event data with some interaction effects
set.seed(456)
n = 500
p = 50
beta = c(rep(0.2,10),rep(-0.3,10),rep(0,p-20), rep(0.01,3), rep(-0.01,2))
x = matrix(rnorm(n*p),n,p)
s=combn(1:20,2)
y1=sample(complex(real=s[1,],imaginary=s[2,]),5,replace = F)
y2=rbind(Re(y1),Im(y1))
y = matrix(c(x[,y2[1,1]]*x[,y2[2,1]],x[,y2[1,2]]*x[,y2[2,2]],x[,y2[1,3]]*x[,y2[2,3]],
             x[,y2[1,4]]*x[,y2[2,4]],x[,y2[1,5]]*x[,y2[2,5]]), nrow=500, ncol=5)
z = cbind(x,y)
A= -(log(runif(n)))
real.time = A/(0.2*exp(drop(z %*% beta)))
real.time = winsor(real.time,2)
cens.time = rexp(n,rate=0.18)
cens.time = winsor(cens.time,2)
status = ifelse(real.time <= cens.time,1,0)
obs.time = ifelse(real.time <= cens.time,real.time,cens.time)
z = as.data.frame(z)
survData = cbind(obs.time,status,z)
surv_dat= survfit(Surv(obs.time,status) ~ 1, data = survData)
pl = ggsurv(surv_dat)
pl

#Generating low dimensional time to event data with some interaction effects only with informative covariates
beta1 = c(rep(0.2,10),rep(-0.3,10),rep(0.01,3), rep(-0.01,2))
x1 = cbind(x[,1:20],y)
real.time = A/(0.2*exp(drop(x1 %*% beta1)))
real.time = winsor(real.time,2)
status = ifelse(real.time <= cens.time,1,0)
obs.time = ifelse(real.time <= cens.time,real.time,cens.time)
z = as.data.frame(x1)
survData1 = cbind(obs.time,status,z)

#c_indices from different models when interaction terms are present

c_index1_inter =c_indexCv_combined(survData,10)

c_index2_inter = c_indexCv_rsf(survData,10,ntree=500)

c_index3_inter = c_indexCv_oracle(survData1,10)

newCovs2 = c_indexCv_combined(survData,10)
newdata2 = cbind( survData$obs.time,survData$status, newCovs2$pred)
newdata2 = as.data.frame(newdata2)
names(newdata2)= c("obs.time", "status","cox", "lasso", "ridge", "net", 
                   "coxboost", "mboost", "gbm")

X= as.matrix(newdata2[,3:9])
Y= as.matrix(newdata2[,1] ) 
delta = as.matrix(newdata2[,2])

#Generating high dimensional time to event data
set.seed(456)
n = 800
p = 1000
beta = c(rep(0.2,10),rep(-0.3,10),rep(0,p-20))
x = matrix(rnorm(n*p),n,p)
A= -(log(runif(n)))
real.time = A/(0.2*exp(drop(x %*% beta)))
real.time = winsor(real.time,2)
cens.time = rexp(n,rate=0.18)
cens.time = winsor(cens.time,2)
status = ifelse(real.time <= cens.time,1,0)
obs.time = ifelse(real.time <= cens.time,real.time,cens.time)
z = as.data.frame(x)
survData = cbind(obs.time,status,z)
surv_dat= survfit(Surv(obs.time,status) ~ 1, data = survData)
pl = ggsurv(surv_dat)

#Generating high dimensional time to event data with informative covariates only
beta1 =  c(rep(0.2,10),rep(-0.3,10))
x1 = x[,1:20]
real.time = A/(0.2*exp(drop(x1 %*% beta1)))
real.time = winsor(real.time,2)
status = ifelse(real.time <= cens.time,1,0)
obs.time = ifelse(real.time <= cens.time,real.time,cens.time)
z = as.data.frame(x1)
survData1 = cbind(obs.time,status,z)

#c-indices for high dimensional linear cases
c_index11_linear =c_indexCv_combined(survData,10)

c_index22_linear = c_indexCv_rsf(survData,10,ntree=500)

c_index33_linear = c_indexCv_oracle(survData1,10)

newCovs3 = c_indexCv_combined(survData,10)
newdata3 = cbind( survData$obs.time,survData$status, newCovs3$pred)
newdata3 = as.data.frame(newdata3)
names(newdata3)= c("obs.time", "status","cox", "lasso", "ridge", "net", 
                   "coxboost", "mboost", "gbm")


X= as.matrix(newdata3[,3:9])
Y= as.matrix(newdata3[,1] ) 
delta = as.matrix(newdata3[,2])

#Generating high dimensional time to event data with some quadratic effects
set.seed(456)
n = 800
p = 1000
beta = c(rep(0.2,10),rep(-0.3,10),rep(0,p-20), rep(-0.01,3), rep(0.01,2))
x = matrix(rnorm(n*p),n,p)
s= sample(1:20, size=5)
y = apply(as.matrix(x[,s], drop=F), 2, function(i) i^2)
z = cbind(x,y)
A= -(log(runif(n)))
real.time = -(log(runif(n)))/(0.2*exp(drop(z %*% beta)))
real.time = winsor(real.time,2)
cens.time = rexp(n,rate=0.18)
cens.time = winsor(cens.time,2)
status = ifelse(real.time <= cens.time,1,0)
obs.time = ifelse(real.time <= cens.time,real.time,cens.time)
z = as.data.frame(z)
survData = cbind(obs.time,status,z)
surv_dat= survfit(Surv(obs.time,status) ~ 1, data = survData)
pl = ggsurv(surv_dat)

#Generating high dimensional time to event data with some quadratic effects and informative covariates only
beta1 = c(rep(0.2,10),rep(-0.3,10),rep(-0.01,3), rep(0.01,2))
x1 = cbind(x[,1:20],y)
real.time = A/(0.2*exp(drop(x1 %*% beta1)))
real.time = winsor(real.time,2)
status = ifelse(real.time <= cens.time,1,0)
obs.time = ifelse(real.time <= cens.time,real.time,cens.time)
z = as.data.frame(x1)
survData1 = cbind(obs.time,status,z)

#c-indices for high dimensional cases with quadratic effects
c_index11_quad =c_indexCv_combined(survData,10)

c_index22_quad = c_indexCv_rsf(survData,10,ntree=500)

c_index33_quad = c_indexCv_oracle(survData1,10)

newCovs4 = c_indexCv_combined(survData,10)
newdata4 = cbind( survData$obs.time,survData$status, newCovs4$pred)
newdata4 = as.data.frame(newdata4)
names(newdata4)= c("obs.time", "status","cox", "lasso", "ridge", "net", 
                   "coxboost", "mboost", "gbm")


X= as.matrix(newdata4[,3:9])
Y= as.matrix(newdata4[,1] ) 
delta = as.matrix(newdata4[,2])

#Generating high dimensional data with some interaction effects
set.seed(456)
n = 800
p = 1000
beta = c(rep(0.2,10),rep(-0.3,10),rep(0,p-20), rep(0.01,3),rep(-0.01,2))
x = matrix(rnorm(n*p),n,p)
s=combn(1:20,2)
y1=sample(complex(real=s[1,],imaginary=s[2,]),5,replace = F)
y2=rbind(Re(y1),Im(y1))
y = matrix(c(x[,y2[1,1]]*x[,y2[2,1]],x[,y2[1,2]]*x[,y2[2,2]],x[,y2[1,3]]*x[,y2[2,3]],
             x[,y2[1,4]]*x[,y2[2,4]],x[,y2[1,5]]*x[,y2[2,5]]), nrow=800, ncol=5)
z = cbind(x,y)
A=-(log(runif(n)))
real.time = -(log(runif(n)))/(0.2*exp(drop(z %*% beta)))
real.time = winsor(real.time,2)
cens.time = rexp(n,rate=0.18)
cens.time = winsor(cens.time,2)
status = ifelse(real.time <= cens.time,1,0)
obs.time = ifelse(real.time <= cens.time,real.time,cens.time)
z = as.data.frame(z)
survData = cbind(obs.time,status,z)
surv_dat= survfit(Surv(obs.time,status) ~ 1, data = survData)
pl = ggsurv(surv_dat)

#Generating high dimensional data with some interaction effects and with informative covaraites only
beta1 = c(rep(0.2,10),rep(-0.3,10),rep(0.01,3),rep(-0.01,2))
x1 = cbind(x[,1:20],y)
real.time = A/(0.2*exp(drop(x1 %*% beta1)))
real.time = winsor(real.time,2)
status = ifelse(real.time <= cens.time,1,0)
obs.time = ifelse(real.time <= cens.time,real.time,cens.time)
z = as.data.frame(x1)
survData1 = cbind(obs.time,status,z)

#c_indices for high dimensional cases with interaction effects
c_index11_inter =c_indexCv_combined(survData,10)

c_index22_inter = c_indexCv_rsf(survData,10,ntree=500)

c_index33_inter = c_indexCv_oracle(survData1,10)

newCovs6 = c_indexCv_combined(survData,10)
#Creating a low dimensional data set with linear predictors as covarites
newdata6 = cbind( survData$obs.time,survData$status, newCovs6$pred)
newdata6 = as.data.frame(newdata6)
names(newdata6)= c("obs.time", "status","cox", "lasso", "ridge", "net", "coxboost", "mboost", "gbm")


X= as.matrix(newdata6[,3:9])
Y= as.matrix(newdata6[,1] ) 
delta = as.matrix(newdata6[,2])

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
  R = lapply(1:m, function(i) which(Y >= t[i]))
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

c_indexCv_convex = function(data,k,niter){
  y_dat = Surv(data$obs.time,data$status)
  set.seed(369)
  folds = sample(rep(1:k, length.out = nrow(data)))
  prediction = c()
  index = c()
  for (j in 1:k){
    idx = which(folds==j)
    train = data[-idx,]
    test = data[idx,]
    Y_train= as.matrix(train[,1]) 
    X_train= as.matrix(train[,3:ncol(train)])
    delta_train = as.matrix(train[,2])
    X_test= as.matrix(test[,3:ncol(test)])
    y_test = Surv(test$obs.time,test$status)
    beta = convex.cox(Y_train,X_train,delta_train,iter.max=niter)$coef 
    pred = X_test%*%as.matrix(beta) 
    prediction = c(prediction, pred)
    index = c(index, idx)
  } 
  
  Match = match(seq(nrow(data)), index)
  prediction = prediction[Match]
  c_index = survConcordance(y_dat~prediction)$concordance
  return(c_index)
}

#sequentional convex combinations of existing algorithms (unsing min correlation)
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

c_indexCv_convex_pair = function(data,k,niter){
  y_dat = Surv(data$obs.time,data$status)
  set.seed(369)
  folds = sample(rep(1:k, length.out = nrow(data)))
  prediction = c()
  index = c()
  for (j in 1:k){
    idx = which(folds==j)
    train = data[-idx,]
    test = data[idx,]
    Y_train= as.matrix(train[,1] ) 
    X_train= as.matrix(train[,3:ncol(train)])
    delta_train = as.matrix(train[,2])
    X_test= as.matrix(test[,3:ncol(test)])
    y_test = Surv(test$obs.time,test$status)
    beta = convex.cox_pair(Y_train,X_train,delta_train,iter.max=niter)$coef_new
    pred = X_test%*%as.matrix(beta)
    prediction = c(prediction, pred)
    index = c(index, idx)
  } 
  
  Match = match(seq(nrow(data)), index)
  prediction = prediction[Match]
  c_index = survConcordance(y_dat~prediction)$concordance
  return(c_index)
}

p_new0= c_indexCv_convex(newdata,10,1000)
p_new00 = c_indexCv_convex_pair(newdata,10,1000) 
p_new1= c_indexCv_convex(newdata1,10,1000)
p_new11 = c_indexCv_convex_pair(newdata1,10,1000) 
p_new2= c_indexCv_convex(newdata2,10,1000)
p_new22 = c_indexCv_convex_pair(newdata2,10,1000) 
p_new3= c_indexCv_convex(newdata3,10,1000)
p_new33 = c_indexCv_convex_pair(newdata3,10,1000) 
p_new4= c_indexCv_convex(newdata4,10,1000)
p_new44 = c_indexCv_convex_pair(newdata4,10,1000) 
p_new5= c_indexCv_convex(newdata5,10,1000)
p_new55 = c_indexCv_convex_pair(newdata5,10,1000) 
p_new6= c_indexCv_convex(newdata5,10,1000)
p_new66 = c_indexCv_convex_pair(newdata5,10,1000) 