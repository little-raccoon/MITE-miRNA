setwd("D:/Data/gzl/")
library(xgboost)
library(magrittr)
library(ggplot2)
library(reshape2)
set.seed(19191)
sample_nums=100 #times down-samplings
sample_scale = 20000 # pos samples everytime 

#load and construct input data
array<-read.table("RNA_array.tsv",header = T)

array_label<-array$type
array_mat<-xgb.DMatrix(data = as.matrix(array[,-1]),label=array_label)
feat_names<-colnames(array_mat)
xgb<-xgboost(array_mat,
             max_depth = 2,
             eta = 1, nthread = 2, nrounds = 50,
             objective = "binary:logistic",eval_metric = "error") #0.975
importance_matrix <- xgb.importance(model = xgb,feature_names = feat_names)
imp_melt<-melt(importance_matrix,measure.vars = 2:4,id.vars = 1)
rm(array_mat)

#Create results matrix
result_freq<-matrix(rep(0,length(feat_names)*sample_nums),
                    ncol = length(feat_names))
colnames(result_freq)<-feat_names
result_gain<-matrix(rep(0,length(feat_names)*sample_nums),
                    ncol = length(feat_names))
colnames(result_gain)<-feat_names
result_cover<-matrix(rep(0,length(feat_names)*sample_nums),
                    ncol = length(feat_names))
colnames(result_cover)<-feat_names

array_pos<-as.matrix(array[array$type == 1,-1])
array_neg<-as.matrix(array[array$type == 0,-1])


#Run xgboost to constrcut models and evaluate the feature importance metrics
for(i in 1:sample_nums){
  array_mat<-xgb.DMatrix(data = rbind(array_pos[sample(x = 1:nrow(array_pos),
                                                       size = sample_scale,
                                                       replace = F),],
                                      array_neg[sample(x = 1:nrow(array_neg),
                                                       size = sample_scale,
                                                       replace = F),]),
                         label=c(rep(1,sample_scale),rep(0,sample_scale)))
  
  xgb<-xgboost(array_mat[sample(1:nrow(array_mat),
                                size = sample_scale,
                                replace = F),],
               max_depth = 2,
               eta = 1, nthread = 2, nrounds = 500,
               objective = "binary:logistic")
  imp_mat<-xgb.importance(model = xgb,feature_names = feat_names) %>% as.data.frame()
  rownames(imp_mat)<-imp_mat$Feature
  result_freq[i,]<-imp_mat[feat_names,"Frequency"]
  result_cover[i,]<-imp_mat[feat_names,"Cover"]
  result_gain[i,]<-imp_mat[feat_names,"Gain"]
}

save(result_cover,result_freq,result_gain,file="Down_sampling.RData")
save(importance_matrix,file="all_Data_xgboost_importance.RData")
