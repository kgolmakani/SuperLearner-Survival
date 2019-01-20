library("survival")
library(glmnet)
library(randomForestSRC)
library(CoxBoost)
library(mboost)
library(gbm)
library(ggplot2)
library(RColorBrewer)
#function for computing k-fold cross-validated concordance index based on different algorithms
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

#1. lung cancer data from survival package
names(lung)
naCount = sapply(lung, function(i) sum(is.na(i)))
lung= subset(lung, !is.na(ph.ecog) & !is.na(ph.karno) &!is.na(pat.karno)& !is.na(meal.cal) &!is.na(wt.loss))
class_lung = sapply(lung, function(i) class(i))
head(lung)
lung$status = lung$status-1
options(na.action="na.pass")
lung = lung[,-1]
colnames(lung)[1] ="obs.time"
dim(lung)
table(lung$status)
lung_c_idx = c_indexCv_combined(lung,10) 
lung_c_idx$c_index
new_lung = cbind(lung$obs.time,lung$status, lung_c_idx$pred)
new_lung = as.data.frame(new_lung)
names(new_lung)= c("obs.time", "status","cox", "lasso", "ridge", "net", 
                   "coxboost", "mboost", "gbm")
X= as.matrix(new_lung[,3:9])
Y= as.matrix(new_lung[,1] ) 
delta = as.matrix(new_lung[,2])
p_convex_lung= c_indexCv_convex(new_lung,10,1000)
p_sconvex_lung = c_indexCv_convex_pair(new_lung,10,1000) 
##############################################################################
#2.rats data from survival package
names(rats)
naCount = sapply(rats, function(i) sum(is.na(i)))
class_rats = sapply(rats, function(i) class(i))
colnames(rats)[3] ="obs.time"
rats = subset(rats, select=c(3,4,1,2,5))
names(rats)
rats$sex = as.factor(rats$sex)
table(rats$status)
dim(rats)
rats_c_idx = c_indexCv_combined(rats,10) 
rats_c_idx$c_index
new_rats = cbind(rats$obs.time,rats$status, rats_c_idx$pred)
new_rats = as.data.frame(new_rats)
names(new_rats)= c("obs.time", "status","cox", "lasso", "ridge", "net", 
                   "coxboost", "mboost", "gbm")
X= as.matrix(new_rats[,3:9])
Y= as.matrix(new_rats[,1] ) 
delta = as.matrix(new_rats[,2])
p_convex_rats= c_indexCv_convex(new_rats,10,1000)
p_sconvex_rats = c_indexCv_convex_pair(new_rats,10,1000) 
#############################################################################
#3.colon data from survival package
names(colon)
colon = colon[,-c(1,2)]
naCount = sapply(colon, function(i) sum(is.na(i)))
colon = subset(colon, !is.na(nodes) & !is.na(differ))
class_colon = sapply(colon, function(i) class(i))
colnames(colon)[13] ="obs.time"
names(colon)
colon = subset(colon, select=c(13,8,1,2,3,4,5,6,7,9,10,11,12,14))
dim(colon)
table(colon$status)
colon_c_idx = c_indexCv_combined(colon,10) 

new_colon = cbind(colon$obs.time,colon$status, colon_c_idx$pred)
new_colon = as.data.frame(new_colon)
names(new_colon) = c("obs.time", "status","cox", "lasso", "ridge", "net", 
                    "coxboost", "mboost", "gbm")
X= as.matrix(new_colon[,3:9])
Y= as.matrix(new_colon[,1] ) 
delta = as.matrix(new_colon[,2])
p_convex_colon = c_indexCv_convex(new_colon,10,2)
p_sconvex_colon = c_indexCv_convex_pair(new_colon,10,2) 

############################################################################
#4.veteran data from survival package
names(veteran)
naCount = sapply(veteran, function(i) sum(is.na(i)))
class_veteran = sapply(veteran, function(i) class(i))
colnames(veteran)[3] ="obs.time"
names(veteran)
veteran = subset(veteran, select=c(3,4,1,2,5,6,7,8))
names(veteran)
dim(veteran)
table(veteran$status)

veteran_c_idx = c_indexCv_combined(veteran,10) 
veteran_c_idx$c_index
new_veteran = cbind(veteran$obs.time,veteran$status, veteran_c_idx$pred)
new_veteran = as.data.frame(new_veteran)
names(new_veteran)= c("obs.time", "status","cox", "lasso", "ridge", "net", 
                      "coxboost", "mboost", "gbm")
X= as.matrix(new_veteran[,3:9])
Y= as.matrix(new_veteran[,1] ) 
delta = as.matrix(new_veteran[,2])
p_convex_veteran = c_indexCv_convex(new_veteran,10,1000)
p_sconvex_veteran = c_indexCv_convex_pair(new_veteran,10,1000) 
############################################################################
#5.kidney data d from survival package
names(kidney)
naCount = sapply(kidney, function(i) sum(is.na(i)))
class_veteran = sapply(kidney, function(i) class(i))
kidney=kidney[,-1]
colnames(kidney)[1] ="obs.time"
dim(kidney)
table(kidney$status)
kidney_c_idx = c_indexCv_combined(kidney,10) 
kidney_c_idx$c_index
new_kidney = cbind(kidney$obs.time,kidney$status, kidney_c_idx$pred)
new_kidney = as.data.frame(new_kidney)
names(new_kidney)= c("obs.time", "status","cox", "lasso", "ridge", "net", 
                     "coxboost", "mboost", "gbm")
X= as.matrix(new_kidney[,3:9])
Y= as.matrix(new_kidney[,1] ) 
delta = as.matrix(new_kidney[,2])
p_convex_kidney = c_indexCv_convex(new_kidney,10,1000)
p_sconvex_kidney = c_indexCv_convex_pair(new_kidney,10,1000) 
#######################################################
#6.retinopathy data from survival package
names(retinopathy)
naCount = sapply(retinopathy, function(i) sum(is.na(i)))
class_retinopathy = sapply(retinopathy, function(i) class(i))

retinopathy=retinopathy[,-1]
retinopathy = subset(retinopathy, select=c(6,7,1,2,3,4,5,8))
colnames(retinopathy)[1] ="obs.time"
dim(retinopathy)
table(retinopathy$status)
retinopathy_c_idx = c_indexCv_combined(retinopathy,10) 
retinopathy_c_idx$c_index

new_retinopathy = cbind(retinopathy$obs.time,retinopathy$status, retinopathy_c_idx$pred)
new_retinopathy = as.data.frame(new_retinopathy)
names(new_retinopathy)= c("obs.time", "status","cox", "lasso", "ridge", "net", 
                          "coxboost", "mboost", "gbm")
X= as.matrix(new_retinopathy[,3:9])
Y= as.matrix(new_retinopathy[,1] ) 
delta = as.matrix(new_retinopathy[,2])
p_convex_retinopathy = c_indexCv_convex(new_retinopathy,10,1000)
p_sconvex_retinopathy = c_indexCv_convex_pair(new_retinopathy,10,1000) 
###########################################
#boxplots
methods = read.table("~/Desktop/methods.csv", header=TRUE, sep=",")
methods1 = subset(methods, Data!="lung" & Data!="colon"& Data!="retinopathy")
ggplot(methods1, aes(x=Data, y=c.index,color=Method)) +
  geom_point()+ scale_color_brewer(palette="Set3")+
  ggtitle("Performance across different data sets")+
  xlab("Data sets") +
  ylab("10-folds cross-validated Harrell’s C-index (%)")+
  theme(plot.margin=unit(c(1,0.5,0.5,1),"cm"))+
  theme(plot.title=element_text(hjust=0.5, size=14))+theme(text=element_text(size=16))

methods2 = subset(methods, Data=="lung" | Data=="colon"| Data=="retinopathy")

ggplot(methods2, aes(x = Data, y = c.index, fill = Method)) +
  stat_boxplot(geom ='errorbar') + geom_boxplot() +
  scale_fill_brewer(palette="Set3")+
  ggtitle("Performance across different data sets")+
  xlab("Data sets") +
  ylab("10-folds cross-validated Harrell’s C-index (%)")+
  theme(plot.margin=unit(c(1,0.5,0.5,1),"cm"))+
  theme(plot.title=element_text(hjust=0.5, size=14))+theme(text=element_text(size=16))
############################################################
source("https://bioconductor.org/biocLite.R")
biocLite("RTCGAToolbox")
library(RTCGAToolbox)
# Low grade glioma
lggData = getFirehoseData(dataset = 'LGG', forceDownload = FALSE, Clinic = TRUE, RNAseq2_Gene_Norm = TRUE, runDate = "20160128")
clin = getData(lggData, 'Clinical')
dim(clin)
naCount = sapply(clin, function(i) sum(is.na(i)))
clin= clin[,-c(1,6,8,9,11)]
clin$obs.time = ifelse(is.na(clin$days_to_death),clin$days_to_last_followup, clin$days_to_death)
clin = subset(clin, !is.na(clin$obs.time) & clin$obs.time>0)
clin = clin[,-c(3,4)]
clin = subset(clin, select=c(8,2,1,3,4,5,6,7))
colnames(clin)[c(2,3)]=c("status","age")
clin = subset(clin, !is.na(race) & !is.na(ethnicity) & !is.na(radiation_therapy) &!is.na(age))
naCount = sapply(clin, function(i) sum(is.na(i)))
dim(clin)
table(clin$status)
gene = getData(lggData, 'RNASeq2GeneNorm')
dim(gene)
gene = as.data.frame(gene)
stcs = rownames(clin)
stcs = gsub(pattern="\\.",replacement="-",stcs)
stcs = toupper(stcs)
samplesDat = data.frame(matrix(nrow=length(stcs),ncol=3))
rownames(samplesDat) = stcs
for(j in 1:length(stcs))
{
  tmpRow = unlist(strsplit(stcs[j],split="-"))
  samplesDat[stcs[j],] = tmpRow
}

stcs = paste(samplesDat[,1],samplesDat[,2],samplesDat[,3],sep="-")
rownames(clin) = stcs

sampleIDs1 = colnames(gene)
sampleIDs1 = gsub(pattern="\\.",replacement="-",sampleIDs1)
samplesDat = data.frame(matrix(nrow=length(sampleIDs1),ncol=7))
rownames(samplesDat) = sampleIDs1
for(j in 1:length(sampleIDs1))
{
  tmpRow = unlist(strsplit(sampleIDs1[j],split="-"))
  samplesDat[sampleIDs1[j],] = tmpRow
}

sampleIDs1 = as.character(samplesDat[,4])
sampleIDs1 <- substr(sampleIDs1,1,nchar(sampleIDs1)-1)
sampleIDs1 <- as.numeric(sampleIDs1)
samplesDat[,4] <- sampleIDs1
sampleIDs1 <- paste(samplesDat[,1],samplesDat[,2],samplesDat[,3],samplesDat[,4],sep="-")
sampleIDs11 <- paste(samplesDat[,1],samplesDat[,2],samplesDat[,3],sep="-")
colnames(gene) <- sampleIDs1
gene = gene[,!duplicated(sampleIDs11)]
colnames(gene) = sampleIDs11[!duplicated(sampleIDs11)]
commonSamples <- intersect(rownames(clin),names(gene))
length(commonSamples)
class(commonSamples)
sd_gene = apply(t(gene),2, sd)
qt = quantile(sd_gene , probs =0.7)
gene1= t(gene)[,sd_gene>qt]
lgg = merge(clin,gene1,by="row.names")
lgg = lgg[,-1]
row.has.na = apply(lgg, 1, function(x){any(is.na(x))})
sum(row.has.na)
lgg = lgg[!row.has.na,]
dim(lgg)
table(lgg$status)
lgg$obs.time= as.numeric(lgg$obs.time)
lgg$status = as.numeric(lgg$status)
lgg$age = as.numeric(lgg$age)
lgg$gender = as.factor(lgg$gender)
lgg$radiation_therapy = as.factor(lgg$radiation_therapy)
lgg$histological_type = as.factor(lgg$histological_type)
lgg$race = as.factor(lgg$race)
lgg$ethnicity = as.factor(lgg$ethnicity)
#lgg data set is ready now

#function to compute k-fold cross-validated concordance index for Lasso-Cox, Ridge-Cox, EN-Cox
  c_indexCv_combined1 = function(data,k){
  y_dat = Surv(data$obs.time,data$status)
  set.seed(1)
  folds = sample(rep(1:k, length.out = nrow(data)))
  prediction_lasso = c()
  prediction_ridge = c()
  prediction_net = c()
  index =c()
  for (j in 1:k){
    idx = which(folds==j)
    train = data[-idx,]
    test = data[idx,]
    y_train = Surv(train$obs.time, train$status)
    y_test = Surv(test$obs.time,test$status)
    x = model.matrix(~., data[,-c(1,2)])
    fit_lasso = glmnet(x[-idx,],y_train, family="cox", alpha=1)
    cvFit_lasso = cv.glmnet(x[-idx,],y_train, family="cox", alpha=1)
    pred_lasso = predict(fit_lasso,x[idx,], s=cvFit_lasso$lambda.min, type="link") 
    fit_ridge = glmnet(x[-idx,],y_train, family="cox", alpha=0)
    cvFit_ridge = cv.glmnet(x[-idx,],y_train, family="cox", alpha=0)
    pred_ridge = predict(fit_ridge,x[idx,], s=cvFit_ridge$lambda.min, type="link") 
    fit_net = glmnet(x[-idx,],y_train, family="cox", alpha=0.5)
    cvFit_net = cv.glmnet(x[-idx,],y_train, family="cox", alpha=0.5)
    pred_net = predict(fit_net,x[idx,], s=cvFit_net$lambda.min, type="link") 
    index = c(index,idx)
    prediction_lasso = c(prediction_lasso, pred_lasso)
    prediction_ridge = c(prediction_ridge, pred_ridge)
    prediction_net = c(prediction_net, pred_net)
  }
  
  Match = match(seq(nrow(data)), index)
  prediction_lasso = prediction_lasso[Match]
  prediction_ridge = prediction_ridge[Match]
  prediction_net = prediction_net[Match]
  c_lasso = survConcordance(y_dat~prediction_lasso)$concordance
  c_ridge = survConcordance(y_dat~prediction_ridge)$concordance
  c_net = survConcordance(y_dat~prediction_net)$concordance
  final_pred = cbind(prediction_lasso, prediction_ridge, prediction_net)
  return(list(pred = final_pred, c_index=c(c_lasso, c_ridge, c_net)))
}

lgg_c_idx1 = c_indexCv_combined1(lgg,5) 
lgg_c_idx1$c_index


#function to compute k-fold cross-validated concordance index for Coxboost, mboost, gbm
c_indexCv_combined2 = function(data,k){
  y_dat = Surv(data$obs.time,data$status)
  set.seed(1)
  folds = sample(rep(1:k, length.out = nrow(data)))
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
    x = model.matrix(~., data[,-c(1,2)])
    optimal.penalty = optimCoxBoostPenalty(time=train$obs.time, 
    status=train$status, x[-idx,],        
    minstepno=50,maxstepno=200,start.penalty=9*sum(train$status==1),
    iter.max=10, upper.margin=0.05)
    fit_coxboost = CoxBoost(time= train$obs.time,status=train$status, x[-idx,],
    standardize=TRUE,subset=1:length(train$obs.time), 
    penalty=optimal.penalty$penalty, 
    stepno=optimal.penalty$cv.res$optimal.step) 
    pred_coxboost = predict(fit_coxboost, x[idx,] , type="lp")  
    
    fit_mboost = glmboost(y_train~., data = train[,-c(1,2)],family=
                            CoxPH(),control=boost_control(mstop = 500))
    pred_mboost = predict(fit_mboost,test[,-c(1,2)], type="link")  
    
    fit_gbm = gbm(y_train~., data=train[,-c(1,2)], distribution="coxph",
                  n.trees=2000, shrinkage=0.01, interaction.depth=5 ,
                  bag.fraction = 0.5, cv.folds=5)
    iter = gbm.perf(fit_gbm,method="OOB")
    pred_gbm = predict(fit_gbm,test[,-c(1,2)],n.trees=iter, type="link") 
    index = c(index,idx)
    prediction_coxboost = c(prediction_coxboost, pred_coxboost)
    prediction_mboost = c(prediction_mboost, pred_mboost)
    prediction_gbm = c(prediction_gbm, pred_gbm)
  }
  
  Match = match(seq(nrow(data)), index)
  prediction_coxboost = prediction_coxboost[Match]
  prediction_mboost = prediction_mboost[Match]
  prediction_gbm = prediction_gbm[Match]
  c_coxboost = survConcordance(y_dat~prediction_coxboost)$concordance
  c_mboost = survConcordance(y_dat~prediction_mboost)$concordance
  c_gbm = survConcordance(y_dat~prediction_gbm)$concordance
  final_pred = cbind(prediction_coxboost,prediction_mboost, prediction_gbm)
  return(list(pred = final_pred, c_index=c(c_coxboost, c_mboost, c_gbm)))
}

lgg_c_idx2 = c_indexCv_combined2(lgg,5) 
lgg_c_idx2$c_index
