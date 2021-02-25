### Scalability analysis of BRM and other methods for different missing value patterns and number of rows

library(mice)
library(randomForest)
library(Hmisc)
library(gbm)
library(data.table)
library(Metrics)
library(gtools)
library(MASS)
library(caret)
library(rpart)
library(tree)
library(DescTools)
library(gtools)
library(dplyr)
library(iml)
library(DALEX)
library(partykit)

gen_missing_bike <- function(bike,propMissing)
{  
  set.seed(1234)
  miss_1 <- sample(nrow(bike),propMissing*nrow(bike),replace=FALSE )
  set.seed(123)
  miss_2 <- sample(nrow(bike),propMissing*nrow(bike),replace=FALSE )
  
  bike_mis <- bike  
  
  ### Generate missing values
  
  bike_mis[miss_1,names(bike_mis) %in% c('windspeed','hr','ToD','windsp2')] <- NA
  bike_mis[miss_2,names(bike_mis) %in% c('temp','hum','weathersit','temp2','hum2')] <- NA
  
  ## Add random noise (only little = 5% ) ## Preserve the noise dataset (important)
  
  miss_3 <- bike_mis[,names(bike_mis)%in% c('ToD','hr','weathersit','windspeed','temp','hum','windsp2','temp2','hum2')]
  miss_4 <- bike_mis[,!names(bike_mis)%in% c('ToD','hr','weathersit','windspeed','temp','hum','windsp2','temp2','hum2')]
  
  for(col in 1:ncol(miss_3))
  {
    m <- 12 + col  ##For reproducibility
    set.seed(m)
    miss_temp <- sample(nrow(miss_3),nrow(miss_3)*0.05,replace=FALSE )
    miss_3[miss_temp,col] <- NA 
  }
  
  bike_mis_n <- cbind(miss_3,miss_4)
  return(bike_mis_n)
}

## Into train and test 
tr_test_st <- function(data_in)
{
  smp_size <- floor(0.75 * nrow(data_in))
  set.seed(12345)
  train_ind <- sample(seq_len(nrow(data_in)), size = smp_size)
  train <- data_in[train_ind, ]
  test <- data_in[-train_ind, ]
  return (list(train,test))
}

## Model errors of subsets
maec <- function(actual,predicted)
{ 
  mean(as.matrix(abs(actual-predicted)), na.rm=T) 
}

rmsec <- function(actual,predicted)
{ 
  sqrt(mean((actual-predicted)^2, na.rm=T))  
}

mapec <- function(actual,predicted)
{ 
  mean(as.matrix(abs(actual-predicted)/actual), na.rm=T)*100 
}

smapec <- function(actual,predicted)
{ 
  mean(as.matrix(2*abs(actual-predicted)/( abs(actual)+abs(predicted) )) , na.rm=T)*100 
}

impute_manual <- function(data_in)
{
  for(j in 1:ncol(data_in))
  {
    if(is.factor(data_in[,j]))
    {
      data_in[is.na(data_in[,j]),j] <- names(sort(-table(data_in[,j])))[1] 
    }
    else
    {
      data_in[is.na(data_in[,j]),j]  <- mean(data_in[,j],na.rm=T)
    }  
  }
  return(data_in)
}

impute_test <- function(data_tr,data_test)
{
  for(j in 1:(ncol(data_test)-2) )  ## Two extra columns added to data_test as compared to data_train: (?)
  {
    if(is.factor(data_test[,j]))
    {
      data_test[is.na(data_test[,j]),j] <- names(sort(-table(data_tr[,j])))[1]
    }
    else
    {
      data_test[is.na(data_test[,j]),j]  <- mean(data_tr[,j],na.rm=T)
    }  
  }
  return(data_test)
}

## Missing proportion in training data:
miss_prop <- function(data_in)
{
  miss_prop_cols <- sum(sapply(data_in, function(y) sum(length(which(is.na(y))))))
  miss_proportion <- miss_prop_cols/(dim(data_in)[1]*dim(data_in)[2])
  return(miss_proportion)
}

#miss_prop(data_train) 
#miss_prop(data_test)

brm_num_blocks <- function(data_train_x,low_threshold = 0.05)
{
  data_present <- data.frame( matrix(as.integer(!is.na(data_train_x)),nrow = nrow(data_train_x), ncol = ncol(data_train_x))) 
  names(data_present) <- names(data_train_x) 
  
  missing_prop <- function(MVI,k, low_threshold)
  {
    determine_best_cluster <- function(data_in,clus_number)
    {
      n <- 10
      kmeans_o <- lapply(1:n,function(x) { set.seed(x); a<-kmeans(data_in,clus_number); return(a) } ) 
      opt_kmeans <- which.max(lapply(kmeans_o,function(x) {100-x$tot.withinss/x$totss*100})) 
      best_cluster <- kmeans_o[[opt_kmeans]]
      return(best_cluster)
    }	
    
    member_clusters <- determine_best_cluster(MVI,k)
    members <- member_clusters$cluster
    data_subs_1 <- split(MVI,members)  
    data_subs_2 <- list(length(data_subs_1)) ## This is for identifying the columns
    missingness <- vector(length = length(data_subs_1))
    for(list in 1:length(data_subs_1))
    {
      completeness <- colSums(data_subs_1[[list]])
      data_subs_2[[list]] <- data_subs_1[[list]][,which(completeness>low_threshold*nrow(data_subs_1[[list]])) ]  
      ## data_subs_2[[list]] <- data_subs_1[[list]][,which(completeness>(1-low_threshold)*nrow(data_subs_1[[list]])) ]
      missing_vals <- nrow(data_subs_2[[list]]) - colSums(data_subs_2[[list]])
      missingness[list] <- sum(missing_vals)
    }
    return(sum(missingness)/(dim(MVI)[1]*dim(MVI)[2]) )
  }
  
  M <- (ncol(data_train_x))  ## Let the maximum number of clusters be equal to the number of columns - 1. 
  
  missing_prop_val <- vector(length=M)
  
  for(k in 1:M)
  {
    missing_prop_val[k] <- missing_prop(data_present,k, low_threshold)
  }
  return(missing_prop_val)
}

brm_ov_subsets <- function(num_blocks, X_in, Y_in,low_threshold = 0.05)
{
  data_xy <- cbind(X_in, Y_in)
  data_present <- data.frame( matrix(as.integer(!is.na(X_in)),nrow = nrow(X_in), ncol = ncol(X_in))) 
  names(data_present) <- names(X_in) 
  
  columns_in_subsets <- function(data_xy, MVI,num_blocks,low_threshold )
  {
    columns <- list()
    determine_best_cluster <- function(data_in,clus_number)
    {
      n <- 10
      kmeans_o <- lapply(1:n,function(x) { set.seed(x); a<-kmeans(data_in,clus_number); return(a) } ) 
      opt_kmeans <- which.max(lapply(kmeans_o,function(x) {100-x$tot.withinss/x$totss*100})) 
      best_cluster <- kmeans_o[[opt_kmeans]]
      return(best_cluster)
    }	
    member_clusters <- determine_best_cluster(MVI,num_blocks)
    members <- member_clusters$cluster
    mcc <- as.data.frame(member_clusters$centers)
    mcc_rounded <- round(mcc,0)
    
    data_subs_1 <- split(MVI,members)  
    data_subs_out <- list(length(data_subs_1))
    data_values <- split(data_xy, members)
    for(list in 1:length(data_subs_1))
    {
      completeness <- colSums(data_subs_1[[list]])
      inputs_in_subset <- which(completeness>low_threshold*nrow(data_subs_1[[list]])) 
      ## inputs_in_subset <- which(completeness>(1-low_threshold)*nrow(data_subs_1[[list]]))
      outcome_index <- which(names(data_values[[list]]) ==  names(Y_in) )
      data_subs_out[[list]] <- data_values[[list]][, c(inputs_in_subset,outcome_index) ]
      columns[[list]] <-names(data_subs_out[[list]])
    }
    return(list(mcc_rounded,columns,data_subs_out ))
  }
  
  n_ov_subsetting <- columns_in_subsets(data_xy,data_present,num_blocks,low_threshold)
  cluster_centers_subsets <- n_ov_subsetting[[1]]
  columns_n_ov_subsets <- n_ov_subsetting[[2]]
  n_ov_subsets <- n_ov_subsetting[[3]]
  ## 
  ### Logic for Set theory determination of overlapping subsets
  
  ov_subsets <- list(length(n_ov_subsets))
  
  for(list1 in 1:length(n_ov_subsets))
  {
    ov_subsets[[list1]] <- n_ov_subsets[[list1]]
    for(list2 in 1:length(n_ov_subsets))
    {
      if( all(columns_n_ov_subsets[[list1]] %in% columns_n_ov_subsets[[list2]]) && (list2 != list1 ) )
      {
        subset_cols <- columns_n_ov_subsets[[list1]]
        ov_subsets[[list1]] <- rbind(ov_subsets[[list1]], n_ov_subsets[[list2]][,subset_cols] )
      }
    }
    ov_subsets[[list1]] <- impute_manual(ov_subsets[[list1]])
    n_ov_subsets[[list1]] <- impute_manual(n_ov_subsets[[list1]])
  }
  return(list(ov_subsets,cluster_centers_subsets,n_ov_subsets, columns_n_ov_subsets ))
}

### Run regression models for each subset
stepreg <- function(data_in)
{
  X <- data_in[,names(data_in) %in% names(data_train_x) ]
  Y <- data_in[,names(data_in) %in% names(data_train_y) ]
  datam <- data.table(X,Y)
  fmla <- as.formula(paste("Y ~ ", paste(names(X), collapse= "+")))
  # reg <- stepAIC(glm( fmla, data=data_in,family="poisson"),direction="both",trace=FALSE)
  reg <- stepAIC(glm.nb( fmla, data=data_in),direction="both",trace=FALSE)
  # reg <- stepAIC(lm( fmla, data=data_in),direction="both",trace=FALSE)
  # reg <- lm( fmla, data=datam)
  # reg <- stepAIC(glm.nb( fmla, data=data_in),direction="both",trace=FALSE)
  return(reg)
}

gbm_train <- function(data_in)
{
  X <- data_in[,names(data_in) %in% names(data_train_x) ]
  Y <- data_in[,names(data_in) %in% names(data_train_y) ]
  datam <- data.table(X,Y) 
  fmla <- as.formula(paste("Y ~ ", paste(names(X), collapse= "+")))
  #reg <- gbm( fmla, data=datam ,verbose=FALSE,n.trees = 500,shrinkage=0.1,interaction.depth = 3)
  reg <- gbm( fmla, data=datam ,distribution ="poisson" ,verbose=FALSE,n.trees = 500,shrinkage=0.1,interaction.depth = 3)
  return(reg)
}

pred_reg <- function(models, data_test,X_in, Y_in,cluster_centers_subsets)
{
  data_train <- cbind(X_in, Y_in)
  ## Check test data, to see which is the best matched model subset
  t_data_present <- data.frame( matrix(as.integer(!is.na(data_test)),nrow = nrow(data_test), ncol = ncol(data_test))) 
  names(t_data_present) <- names(data_test) 
  t_data_present <- t_data_present[,!names(t_data_present) %in% names(Y_in)]
  
  nearest_clus <- apply(t_data_present, 1, function(d)
  {
    temp <- apply(cluster_centers_subsets,1,function(x) dist(rbind(x,d)))
    return(which.min(temp))
  }
  )
  
  data_test_in <- data_test
  data_test_in$nearest_clus <- as.factor(nearest_clus)
  data_test_in$rownum <- 1:nrow(data_test_in)
  data_groups <- split(data_test_in,nearest_clus)
  get_back_row_order <- do.call("rbind",data_groups)
  
  data_igroups <- lapply(data_groups,function(x) impute_test(data_train,x))
  reg_outcomes <- list()
  
  for(i in 1: length(data_igroups))
  {   
    data_in <- data_igroups[[i]]
    X <- data_in[,names(data_in) %in% names(X_in)]
    Y <- data_in[,names(data_in) %in% names(Y_in)]
    datam <- data.table(X,Y)
    if(class(models[[i]])[1] == "negbin")
    {
      reg_outcomes[[i]] <- predict(models[[i]],newdata = datam,type="response")
    }
    else if(class(models[[i]]) == "gbm"){
      reg_outcomes[[i]] <- predict(models[[i]],newdata = datam,n.trees=models[[i]]$n.trees,type="response")
    }
    # reg_outcomes[[i]] <- predict(reg_m[[i]],newdata = datam,type="response")
  }
  
  reg_out_a <- unlist(reg_outcomes)
  reg_ord <- reg_out_a[order(get_back_row_order[,"rownum"])]
  return(reg_ord)
}

model_fit_comparisons_reg <- function(models, ov_subsets_)
{
  R_sq_l <- list()
  
  for(l in 1:length(models))
  {
    pred = predict(models[[l]],ov_subsets_[[l]],  type= "response")
    act = ov_subsets_[[l]][,"cnt"]
    R_sq_l[[l]] = 1 - (sum((act - pred)^2)/(sum((act - mean(act))^2)))
  }  
  return(R_sq_l)
  #return(lapply(models, function(x) summary(x)$r.squared))
}

## Get the global coefficient scores 
get_global_wts <- function(models)
{
  global_reg_par <- function(model)
  {
    nrow <- length(model$residuals)
    coef <- summary(model)$coefficients[,1]
    se <- summary(model)$coefficients[,2]
    coef*nrow/se
  }
  
  global_reg_par_wt <- function(model)
  {
    nrow <- length(model$residuals)
    se <- summary(model)$coefficients[,2]
    nrow/se
  }
  
  global_wt_num <- lapply(models,global_reg_par)
  global_wt_denom <- lapply(models,global_reg_par_wt)
  
  global_wt_denom2 <- lapply(global_wt_denom,function(x) data.frame(t(x)))
  global_wt_denom_mat <- do.call(smartbind,global_wt_denom2)
  global_wt_denom_mat_avg <- apply(global_wt_denom_mat,2,function(x) sum(x,na.rm=T))
  
  global_wt_num2 <- lapply(global_wt_num,function(x) data.frame(t(x)))
  global_wt_num_mat <- do.call(smartbind,global_wt_num2)
  global_wt_num_mat_avg <- apply(global_wt_num_mat,2,function(x) sum(x,na.rm=T))
  
  global_wts <- global_wt_num_mat_avg/global_wt_denom_mat_avg
  return(global_wts)  
}

## Get the global coefficient scores 
get_global_wts_gbm <- function(models)
{
  global_gbm_par_wt <- function(model){
    return(sqrt(length(model$var.names)*length(model$fit))/mse_gbm(model))
  }
  
  wt2 <- do.call(sum,lapply(models,global_gbm_par_wt))
  
  global_gbm_par <- function(model){
    varimp <- summary(model)*sqrt(length(model$var.names)*length(model$fit))/(mse_gbm(model)*wt2)
    # variable = data.frame(summary(model)[,1]
    return(varimp)
  }  ## Wont work as it is ordered for importance
  
  gbm_wt_parameters <- lapply(models,global_gbm_par)
  
  weighted_gbm <- do.call(smartbind,lapply(gbm_wt_parameters, function(x) t(x)))
  #weigthed_gbm2 <- weighted_gbm[!is.na(weighted_gbm$ToD),]
  weighted_gbm_sum <- apply(weighted_gbm,2, function(x) sum(x,na.rm=T))
  
  return(weighted_gbm_sum)  
}


feature_drop_loss <- function(subsets,X_in,Y_in, model, models,cluster_centers_subsets)
{
  data_train <- cbind(X_in, Y_in)
  pred_insample <- pred_reg(models, data_train,X_in,Y_in,cluster_centers_subsets)
  pred_insample_rmse <- rmsec(unlist(Y_in),pred_insample)
  feature_imp <- vector(length = ncol(X_in))
  
  for(feature in 1:ncol(X_in))
  {
    subsets_modified <- lapply(subsets, function(sub) 
    {
      sub[,names(sub) == names(X_in)[feature]] <- runif(n= nrow(sub))
      return(sub)
    }
    )
    models_modified <- lapply(subsets_modified,model)
    pred_insample_modified <- pred_reg(models_modified, data_train,X_in,Y_in,cluster_centers_subsets)
    pred_insample_rmse_modified <- rmsec(unlist(Y_in),pred_insample_modified)
    feature_imp[feature] =  pred_insample_rmse_modified - pred_insample_rmse
  }
  feature_imp_matrix <- data.frame(variable = names(X_in),delta_rmse = feature_imp)
  return(feature_imp_matrix)
}

#mse_reg <- function(model) { sum((model$residuals)^2)/length(model) }   ### THIS IS WRONG! 
mse_reg <- function(model) { mean((model$residuals)^2) } ##NN
mse_gbm <- function(model) {model$train.error[length(model$train.error)]} ##NN

##### main{} #####

## Get the data and do some pre-processing
data_all <- read.csv("data/bike_all.csv")
bike <- data_all[,!names(data_all) %in% c("instant","dteday","atemp","yr","workingday","holiday")] 
bike$weekday <- as.factor(bike$weekday)
bike$weathersit[bike$weathersit==4] <- 3
bike$weathersit <- as.factor(bike$weathersit)
bike$season <- as.factor(bike$season)
#bike$yr <- as.factor(bike$yr)
bike$ToD <- cut(bike$hr, breaks = c(-0.01,5,10,16,19,24))
#bike$windsp2 <- bike$windspeed^2
#bike$temp2 <- bike$temp^2
#bike$hum2 <- bike$hum^2

## Renormalize
bike$temp <- bike$temp*(39+8) -8 

dat_list <- list()

dat_list[[1]] <- gen_missing_bike(bike,0) 
dat_list[[2]] <- gen_missing_bike(bike,0.05)
dat_list[[3]] <- gen_missing_bike(bike,0.10) 
dat_list[[4]] <- gen_missing_bike(bike,0.15)
dat_list[[5]] <- gen_missing_bike(bike,0.20) 
dat_list[[6]] <- gen_missing_bike(bike,0.25)
dat_list[[7]] <- gen_missing_bike(bike,0.5)
dat_list[[8]] <- gen_missing_bike(bike,0.75) 

### RE-ORDER TO:
## n_ov_subset[[1]] has weathersit, temp, hum missing
## n_ov_subset[[2]] has hr, Windspeed, ToD as missing
## n_ov_subset[[3]] has no missing
## n_ov_subset[[4]] has weather, temp, hum, hr, Windspeed, ToD missing


## Instead of using the matrix estimating functions, we do an analogous scalar weighted execution of MBI  
MBI_model <- function(n_ov_subsets)
{ 
  # n_ov_subsets_x <- list(n_ov_subsets_ordered[[3]],n_ov_subsets_ordered[[4]],n_ov_subsets_ordered[[2]],n_ov_subsets_ordered[[1]])
  
  n_ov_subsets_ordered <- n_ov_subsets
  
  columns_n_ov_subsets <- lapply(n_ov_subsets, function(x) names(x))
  
  for(list in 1:length(n_ov_subsets))
  {
    if(all(columns_n_ov_subsets[[list]] == c("hr",        "windspeed", "ToD",       "season",    "mnth",      "weekday",   "cnt")) )
    { n_ov_subsets_ordered[[1]] <- n_ov_subsets[[list]] }
    if(all(columns_n_ov_subsets[[list]] == c("weathersit", "temp",       "hum" ,       "season",    "mnth",      "weekday",   "cnt")) )
    {  n_ov_subsets_ordered[[2]] <- n_ov_subsets[[list]] }
    if(all(columns_n_ov_subsets[[list]] == c("hr",  "weathersit","temp", "hum" ,     "windspeed", "ToD",       "season",    "mnth",      "weekday",   "cnt")) )
    {  n_ov_subsets_ordered[[3]] <- n_ov_subsets[[list]] }
    if(all(columns_n_ov_subsets[[list]] == c(  "season",    "mnth",      "weekday",   "cnt")) )
    {  n_ov_subsets_ordered[[4]] <- n_ov_subsets[[list]] }
  }
  
  n_ov_subsets_x = n_ov_subsets_ordered
  
  for(list in 1:length(n_ov_subsets_ordered))
  {
    n_ov_subsets_x[[list]] <- n_ov_subsets_x[[list]][,!names(n_ov_subsets_x[[list]]) %in% c("cnt")]
  }
  
  all_imputed_subsets <- n_ov_subsets_x
  subs <- n_ov_subsets_x 
  
  ### SUBSET 1:
  
  temp_pred <- list()
  temp_m <- list()
  fmla_temp <- list()
  cols_temp <- list()
  cols_temp[[1]] <- names(subs[[2]])[names(subs[[2]]) %in% names(subs[[1]])]
  cols_temp[[2]] <- names(subs[[3]])[names(subs[[3]]) %in% names(subs[[1]])]
  fmla_temp[[1]] <- as.formula(paste("temp ~ ", paste(cols_temp[[1]], collapse= "+")))
  fmla_temp[[2]] <- as.formula(paste("temp ~ ", paste(cols_temp[[2]], collapse= "+")))
  temp_m[[1]] <- lm(fmla_temp[[1]],data=subs[[2]])
  temp_m[[2]] <- lm(fmla_temp[[2]],data=subs[[3]])
  temp_pred[[1]] <- predict(temp_m[[1]],newdata = subs[[1]]) 
  temp_pred[[2]] <- predict(temp_m[[2]],newdata = subs[[1]]) 
  temp_pred_m <- apply(do.call(cbind, temp_pred),1, mean)
  all_imputed_subsets[[1]][,"temp"] <- temp_pred_m
  
  hum_pred <- list()
  hum_m <- list()
  fmla_hum <- list()
  cols_hum <- list()
  cols_hum[[1]] <- names(subs[[2]])[names(subs[[2]]) %in% names(subs[[1]])]
  cols_hum[[2]] <- names(subs[[3]])[names(subs[[3]]) %in% names(subs[[1]])]
  fmla_hum[[1]] <- as.formula(paste("hum ~ ", paste(cols_hum[[1]], collapse= "+")))
  fmla_hum[[2]] <- as.formula(paste("hum ~ ", paste(cols_hum[[2]], collapse= "+")))
  hum_m[[1]] <- lm(fmla_hum[[1]],data=subs[[2]])
  hum_m[[2]] <- lm(fmla_hum[[2]],data=subs[[3]])
  hum_pred[[1]] <- predict(hum_m[[1]],newdata = subs[[1]]) 
  hum_pred[[2]] <- predict(hum_m[[2]],newdata = subs[[1]]) 
  hum_pred_m <- apply(do.call(cbind, hum_pred),1, mean)
  all_imputed_subsets[[1]][,"hum"] <- hum_pred_m
  
  weathersit_pred <- list()
  weathersit_m <- list()
  fmla_weathersit <- list()
  cols_weathersit <- list()
  cols_weathersit[[1]] <- names(subs[[2]])[names(subs[[2]]) %in% names(subs[[1]])]
  cols_weathersit[[2]] <- names(subs[[3]])[names(subs[[3]]) %in% names(subs[[1]])]
  fmla_weathersit[[1]] <- as.formula(paste("weathersit ~ ", paste(cols_weathersit[[1]], collapse= "+")))
  fmla_weathersit[[2]] <- as.formula(paste("weathersit ~ ", paste(cols_weathersit[[2]], collapse= "+")))
  weathersit_m[[1]] <- ctree(fmla_weathersit[[1]],data=subs[[2]])
  weathersit_m[[2]] <- ctree(fmla_weathersit[[2]],data=subs[[3]])
  weathersit_pred[[1]] <- predict(weathersit_m[[1]],newdata = subs[[1]]) 
  weathersit_pred[[2]] <- predict(weathersit_m[[2]],newdata = subs[[1]]) 
  weathersit_pred_m <- apply(data.frame(weathersit_pred[[1]],weathersit_pred[[2]]),1, function(x) names(sort(-table(x)))[1]  )
  all_imputed_subsets[[1]][,"weathersit"] <- weathersit_pred_m
  
  ### SUBSET 2:
  
  hr_pred <- list()
  hr_m <- list()
  fmla_hr <- list()
  cols_hr <- list()
  cols_hr[[1]] <- names(subs[[1]])[names(subs[[1]]) %in% names(subs[[2]])]
  cols_hr[[2]] <- names(subs[[3]])[names(subs[[3]]) %in% names(subs[[2]])]
  fmla_hr[[1]] <- as.formula(paste("hr ~ ", paste(cols_hr[[1]], collapse= "+")))
  fmla_hr[[2]] <- as.formula(paste("hr ~ ", paste(cols_hr[[2]], collapse= "+")))
  hr_m[[1]] <- lm(fmla_hr[[1]],data=subs[[1]])
  hr_m[[2]] <- lm(fmla_hr[[2]],data=subs[[3]])
  hr_pred[[1]] <- predict(hr_m[[1]],newdata = subs[[2]]) 
  hr_pred[[2]] <- predict(hr_m[[2]],newdata = subs[[2]]) 
  hr_pred_m <- apply(do.call(cbind, hr_pred),1, mean)
  all_imputed_subsets[[2]][,"hr"] <- hr_pred_m
  
  windspeed_pred <- list()
  windspeed_m <- list()
  fmla_windspeed <- list()
  cols_windspeed <- list()
  cols_windspeed[[1]] <- names(subs[[1]])[names(subs[[1]]) %in% names(subs[[2]])]
  cols_windspeed[[2]] <- names(subs[[3]])[names(subs[[3]]) %in% names(subs[[2]])]
  fmla_windspeed[[1]] <- as.formula(paste("windspeed ~ ", paste(cols_windspeed[[1]], collapse= "+")))
  fmla_windspeed[[2]] <- as.formula(paste("windspeed ~ ", paste(cols_windspeed[[2]], collapse= "+")))
  windspeed_m[[1]] <- lm(fmla_windspeed[[1]],data=subs[[1]])
  windspeed_m[[2]] <- lm(fmla_windspeed[[2]],data=subs[[3]])
  windspeed_pred[[1]] <- predict(windspeed_m[[1]],newdata = subs[[2]]) 
  windspeed_pred[[2]] <- predict(windspeed_m[[2]],newdata = subs[[2]]) 
  windspeed_pred_m <- apply(do.call(cbind, windspeed_pred),1, mean)
  all_imputed_subsets[[2]][,"windspeed"] <- windspeed_pred_m 
  
  ToD_pred <- list()
  ToD_m <- list()
  fmla_ToD <- list()
  cols_ToD <- list()
  cols_ToD[[1]] <- names(subs[[1]])[names(subs[[1]]) %in% names(subs[[2]])]
  cols_ToD[[2]] <- names(subs[[3]])[names(subs[[3]]) %in% names(subs[[2]])]
  fmla_ToD[[1]] <- as.formula(paste("ToD ~ ", paste(cols_ToD[[1]], collapse= "+")))
  fmla_ToD[[2]] <- as.formula(paste("ToD ~ ", paste(cols_ToD[[2]], collapse= "+")))
  ToD_m[[1]] <- ctree(fmla_ToD[[1]],data=subs[[1]])
  ToD_m[[2]] <- ctree(fmla_ToD[[2]],data=subs[[3]])
  ToD_pred[[1]] <- predict(ToD_m[[1]],newdata = subs[[2]]) 
  ToD_pred[[2]] <- predict(ToD_m[[2]],newdata = subs[[2]]) 
  ToD_pred_m <- apply(data.frame(ToD_pred[[1]],ToD_pred[[2]]),1, function(x) names(sort(-table(x)))[1] )
  all_imputed_subsets[[2]][,"ToD"] <- ToD_pred_m
  
  ### SUBSET 4:
  
  hr_pred <- list()
  hr_m <- list()
  fmla_hr <- list()
  cols_hr <- list()
  cols_hr[[1]] <- names(subs[[1]])[names(subs[[1]]) %in% names(subs[[4]])]
  cols_hr[[2]] <- names(subs[[3]])[names(subs[[3]]) %in% names(subs[[4]])]
  fmla_hr[[1]] <- as.formula(paste("hr ~ ", paste(cols_hr[[1]], collapse= "+")))
  fmla_hr[[2]] <- as.formula(paste("hr ~ ", paste(cols_hr[[2]], collapse= "+")))
  hr_m[[1]] <- lm(fmla_hr[[1]],data=subs[[1]])
  hr_m[[2]] <- lm(fmla_hr[[2]],data=subs[[3]])
  hr_pred[[1]] <- predict(hr_m[[1]],newdata = subs[[4]]) 
  hr_pred[[2]] <- predict(hr_m[[2]],newdata = subs[[4]]) 
  hr_pred_m <- apply(do.call(cbind, hr_pred),1, mean)
  all_imputed_subsets[[4]][,"hr"] <- hr_pred_m
  
  windspeed_pred <- list()
  windspeed_m <- list()
  fmla_windspeed <- list()
  cols_windspeed <- list()
  cols_windspeed[[1]] <- names(subs[[1]])[names(subs[[1]]) %in% names(subs[[4]])]
  cols_windspeed[[2]] <- names(subs[[3]])[names(subs[[3]]) %in% names(subs[[4]])]
  fmla_windspeed[[1]] <- as.formula(paste("windspeed ~ ", paste(cols_windspeed[[1]], collapse= "+")))
  fmla_windspeed[[2]] <- as.formula(paste("windspeed ~ ", paste(cols_windspeed[[2]], collapse= "+")))
  windspeed_m[[1]] <- lm(fmla_windspeed[[1]],data=subs[[1]])
  windspeed_m[[2]] <- lm(fmla_windspeed[[2]],data=subs[[3]])
  windspeed_pred[[1]] <- predict(windspeed_m[[1]],newdata = subs[[4]]) 
  windspeed_pred[[2]] <- predict(windspeed_m[[2]],newdata = subs[[4]]) 
  windspeed_pred_m <- apply(do.call(cbind, windspeed_pred),1, mean)
  all_imputed_subsets[[4]][,"windspeed"] <- windspeed_pred_m 
  
  ToD_pred <- list()
  ToD_m <- list()
  fmla_ToD <- list()
  cols_ToD <- list()
  cols_ToD[[1]] <- names(subs[[1]])[names(subs[[1]]) %in% names(subs[[4]])]
  cols_ToD[[2]] <- names(subs[[3]])[names(subs[[3]]) %in% names(subs[[4]])]
  fmla_ToD[[1]] <- as.formula(paste("ToD ~ ", paste(cols_ToD[[1]], collapse= "+")))
  fmla_ToD[[2]] <- as.formula(paste("ToD ~ ", paste(cols_ToD[[2]], collapse= "+")))
  ToD_m[[1]] <- ctree(fmla_ToD[[1]],data=subs[[1]])
  ToD_m[[2]] <- ctree(fmla_ToD[[2]],data=subs[[3]])
  ToD_pred[[1]] <- predict(ToD_m[[1]],newdata = subs[[4]]) 
  ToD_pred[[2]] <- predict(ToD_m[[2]],newdata = subs[[4]]) 
  ToD_pred_m <- apply(data.frame(ToD_pred[[1]],ToD_pred[[2]]),1, function(x) names(sort(-table(x)))[1]     ) 
  all_imputed_subsets[[4]][,"ToD"] <- ToD_pred_m
  
  temp_pred <- list()
  temp_m <- list()
  fmla_temp <- list()
  cols_temp <- list()
  cols_temp[[1]] <- names(subs[[2]])[names(subs[[2]]) %in% names(subs[[4]])]
  cols_temp[[2]] <- names(subs[[3]])[names(subs[[3]]) %in% names(subs[[4]])]
  fmla_temp[[1]] <- as.formula(paste("temp ~ ", paste(cols_temp[[1]], collapse= "+")))
  fmla_temp[[2]] <- as.formula(paste("temp ~ ", paste(cols_temp[[2]], collapse= "+")))
  temp_m[[1]] <- lm(fmla_temp[[1]],data=subs[[2]])
  temp_m[[2]] <- lm(fmla_temp[[2]],data=subs[[3]])
  temp_pred[[1]] <- predict(temp_m[[1]],newdata = subs[[4]]) 
  temp_pred[[2]] <- predict(temp_m[[2]],newdata = subs[[4]]) 
  temp_pred_m <- apply(do.call(cbind, temp_pred),1, mean)
  all_imputed_subsets[[4]][,"temp"] <- temp_pred_m
  
  hum_pred <- list()
  hum_m <- list()
  fmla_hum <- list()
  cols_hum <- list()
  cols_hum[[1]] <- names(subs[[2]])[names(subs[[2]]) %in% names(subs[[4]])]
  cols_hum[[2]] <- names(subs[[3]])[names(subs[[3]]) %in% names(subs[[4]])]
  fmla_hum[[1]] <- as.formula(paste("hum ~ ", paste(cols_hum[[1]], collapse= "+")))
  fmla_hum[[2]] <- as.formula(paste("hum ~ ", paste(cols_hum[[2]], collapse= "+")))
  hum_m[[1]] <- lm(fmla_hum[[1]],data=subs[[2]])
  hum_m[[2]] <- lm(fmla_hum[[2]],data=subs[[3]])
  hum_pred[[1]] <- predict(hum_m[[1]],newdata = subs[[4]]) 
  hum_pred[[2]] <- predict(hum_m[[2]],newdata = subs[[4]]) 
  hum_pred_m <- apply(do.call(cbind, hum_pred),1, mean)
  all_imputed_subsets[[4]][,"hum"] <- hum_pred_m
  
  weathersit_pred <- list()
  weathersit_m <- list()
  fmla_weathersit <- list()
  cols_weathersit <- list()
  cols_weathersit[[1]] <- names(subs[[2]])[names(subs[[2]]) %in% names(subs[[4]])]
  cols_weathersit[[2]] <- names(subs[[3]])[names(subs[[3]]) %in% names(subs[[4]])]
  fmla_weathersit[[1]] <- as.formula(paste("weathersit ~ ", paste(cols_weathersit[[1]], collapse= "+")))
  fmla_weathersit[[2]] <- as.formula(paste("weathersit ~ ", paste(cols_weathersit[[2]], collapse= "+")))
  weathersit_m[[1]] <- ctree(fmla_weathersit[[1]],data=subs[[2]])
  weathersit_m[[2]] <- ctree(fmla_weathersit[[2]],data=subs[[3]])
  weathersit_pred[[1]] <- predict(weathersit_m[[1]],newdata = subs[[4]]) 
  weathersit_pred[[2]] <- predict(weathersit_m[[2]],newdata = subs[[4]]) 
  weathersit_pred_m <- apply(data.frame(weathersit_pred[[1]],weathersit_pred[[2]]),1, function(x) names(sort(-table(x)))[1] )
  all_imputed_subsets[[4]][,"weathersit"] <- weathersit_pred_m
  
  for(list in 1:length(n_ov_subsets_ordered))
  {
    all_imputed_subsets[[list]][,"cnt"] <- n_ov_subsets_ordered[[list]][,"cnt"]
  }
  
  imputed_data <- do.call(smartbind,all_imputed_subsets)
  imputed_data$ToD <- factor(imputed_data$ToD)
  imputed_data$weathersit <- factor(imputed_data$weathersit)
  
  return(imputed_data)
}
## Can build whatever model is required over this and predict for imputed test dataset

# mbi_imputed_data <- MBI_model(ov_subsets)
mbi_imputed_data <- MBI_model(n_ov_subsets)

#reg_mbi <- stepAIC(glm( cnt~., family= "poisson", data=mbi_imputed_data))
#reg_mbi <- stepAIC(glm.nb( cnt~., data=mbi_imputed_data))
reg_mbi <- stepAIC(lm( cnt~., data=mbi_imputed_data), direction='both',trace=FALSE)

# gbm_mbi <- gbm(cnt~., data=mbi_imputed_data, distribution ="poisson", verbose=FALSE,n.trees = 500,shrinkage=0.1,interaction.depth = 3)


mbi_predictions <- function(data_train, data_test)
{
  # data_xy <- cbind(X_in, Y_in)
  data_xy <- data_train
  X_in <- data_train[,!names(data_train) %in% c("cnt")]
  Y_in <- data.frame(cnt = data_train[,"cnt"])
  data_present <- data.frame( matrix(as.integer(!is.na(X_in)),nrow = nrow(X_in), ncol = ncol(X_in))) 
  names(data_present) <- names(X_in) 
  
  columns_in_subsets <- function(data_xy, MVI,num_blocks,low_threshold = 0.05)
  {
    columns <- list()
    determine_best_cluster <- function(data_in,clus_number)
    {
      n <- 10
      kmeans_o <- lapply(1:n,function(x) { set.seed(x); a<-kmeans(data_in,clus_number); return(a) } )
      opt_kmeans <- which.max(lapply(kmeans_o,function(x) {100-x$tot.withinss/x$totss*100}))
      best_cluster <- kmeans_o[[opt_kmeans]]
      return(best_cluster)
    }
    member_clusters <- determine_best_cluster(MVI,num_blocks)
   #  member_clusters <- kmeans(MVI,num_blocks)
    members <- member_clusters$cluster
    mcc <- as.data.frame(member_clusters$centers)
    mcc_rounded <- round(mcc,0)
    
    data_subs_1 <- split(MVI,members)  
    data_subs_out <- list(length(data_subs_1))
    data_values <- split(data_xy, members)
    for(list in 1:length(data_subs_1))
    {
      completeness <- colSums(data_subs_1[[list]])
      ## inputs_in_subset <- which(completeness>low_threshold*nrow(data_subs_1[[list]]))
      inputs_in_subset <- which(completeness>(low_threshold)*nrow(data_subs_1[[list]]))
      outcome_index <- which(names(data_values[[list]]) ==  names(Y_in) )
      data_subs_out[[list]] <- data_values[[list]][, c(inputs_in_subset,outcome_index) ]
      columns[[list]] <-names(data_subs_out[[list]])
    }
    return(list(mcc_rounded,columns,data_subs_out ))
  }
  
  n_ov_subsetting <- columns_in_subsets(data_xy,data_present,num_blocks,low_threshold = 0.05)
  cluster_centers_subsets <- n_ov_subsetting[[1]]
  columns_n_ov_subsets <- n_ov_subsetting[[2]]
  n_ov_subsets <- n_ov_subsetting[[3]]
  
  n_ov_subsets <- lapply(n_ov_subsets, impute_manual) 
  
  ### Hard-code the order of subsets here:
  ## n_ov_subset[[1]] has weathersit, temp, hum missing
  ## n_ov_subset[[2]] has hr, Windspeed, ToD as missing
  ## n_ov_subset[[3]] has no missing
  ## n_ov_subset[[4]] has weather, temp, hum, hr, Windspeed, ToD missing
  
  n_ov_subsets_ordered <- n_ov_subsets
  
  for(list in 1:length(n_ov_subsets))
  {
    if(all(columns_n_ov_subsets[[list]] == c("hr",        "windspeed", "ToD",       "season",    "mnth",      "weekday",   "cnt")) )
    { n_ov_subsets_ordered[[1]] <- n_ov_subsets[[list]] }
    if(all(columns_n_ov_subsets[[list]] == c("weathersit", "temp",       "hum" ,       "season",    "mnth",      "weekday",   "cnt")) )
    {  n_ov_subsets_ordered[[2]] <- n_ov_subsets[[list]] }
    if(all(columns_n_ov_subsets[[list]] == c("hr",  "weathersit","temp", "hum" ,     "windspeed", "ToD",       "season",    "mnth",      "weekday",   "cnt")) )
    {  n_ov_subsets_ordered[[3]] <- n_ov_subsets[[list]] }
    if(all(columns_n_ov_subsets[[list]] == c(  "season",    "mnth",      "weekday",   "cnt")) )
    {  n_ov_subsets_ordered[[4]] <- n_ov_subsets[[list]] }
  }
  
  mbi_imputed_data <- MBI_model(n_ov_subsets_ordered)
  
  reg_mbi <- stepAIC(lm( cnt~., data=mbi_imputed_data), direction="both", trace=FALSE)
  # reg_mbi <- glm( cnt~., family= "poisson", data=mbi_imputed_data)
  # reg_mbi <- glm.nb( cnt~., data=mbi_imputed_data)
  #gbm_mbi <- gbm(cnt~., data=mbi_imputed_data, verbose=FALSE,n.trees = 500,shrinkage=0.1,interaction.depth = 3)
  
  # reg_pred_mbi <- predict(reg_mbi, data_test_imputed,type="response")
  reg_pred_mbi <- predict(reg_mbi, data_test_imputed)
  #gbm_pred_mbi <- predict(gbm_mbi, data_test_imputed,n.trees=gbm_mbi$n.trees)
  # return(list(reg_pred_mbi, gbm_pred_mbi))
  return(reg_pred_mbi)
}



bike_mis_n <- dat_list[[7]]
na_count <-sapply(bike_mis_n, function(y) sum(length(which(is.na(y)))))

### Here, we take different data sizes into consideration ##

sc_data_list <- list()
set.seed(12345)
sc_data_list[[1]] <- bike_mis_n[sample(nrow(bike_mis_n),5000,replace=TRUE ),]
set.seed(12345)
sc_data_list[[2]] <- bike_mis_n[sample(nrow(bike_mis_n),10000,replace=TRUE ),]
set.seed(12345)
sc_data_list[[3]] <- bike_mis_n[sample(nrow(bike_mis_n),15000,replace=TRUE ),]
set.seed(12345)
sc_data_list[[4]] <- bike_mis_n[sample(nrow(bike_mis_n),50000,replace=TRUE ),]
set.seed(12345)
sc_data_list[[5]] <- bike_mis_n[sample(nrow(bike_mis_n),100000,replace=TRUE ),]
set.seed(12345)
sc_data_list[[6]] <- bike_mis_n[sample(nrow(bike_mis_n),500000,replace=TRUE ),]
set.seed(12345)
sc_data_list[[7]] <- bike_mis_n[sample(nrow(bike_mis_n),1000000,replace=TRUE ),]


impute_manual <- function(data_in)
{
  for(j in 1:ncol(data_in))
  {
    if(is.factor(data_in[,j]))
    {
      data_in[is.na(data_in[,j]),j] <- names(sort(-table(data_in[,j])))[1] 
    }
    else
    {
      data_in[is.na(data_in[,j]),j]  <- mean(data_in[,j],na.rm=T)
    }  
  }
  return(data_in)
}


impute_test <- function(data_tr,data_test)
{
  for(j in 1:(ncol(data_test)-2) )  ## Two extra columns added to data_test as compared to data_train
  {
    if(is.factor(data_test[,j]))
    {
          data_test[is.na(data_test[,j]),j] <- names(sort(-table(data_tr[,j])))[1]
    }
    else
    {
      data_test[is.na(data_test[,j]),j]  <- mean(data_tr[,j],na.rm=T)
    }  
  }
  return(data_test)
}


red_model <- function(data_train, data_test)
{
  ## Clustering
  
  data_present <- data.frame( matrix(as.integer(!is.na(data_train)),nrow = nrow(data_train), ncol = ncol(data_train))) 
  names(data_present) <- names(data_train) 
  
  ## This is a very simplified version where we know k as 4.
  member_clusters <- kmeans(data_present,4)
  members <- member_clusters$cluster
  
  ## Once the cluster memberships are determined
  data_subs <- split(data_train,members)  
  
  col_names <- list()
  data_sets <- list()
  
  mcc <- as.data.frame(member_clusters$centers)
  mcc_rounded <- round(mcc,0)
  for(row in 1:nrow(mcc))
  {
    temp <- vector()
    for( col in 1:ncol(mcc)) 
    {
      if (mcc[row,col] > 0.1)
      {
        temp <- c(temp, colnames(mcc)[col] )  
      }
    }
    col_names[[row]] <- temp
  } 
  
  dat_inl <- list() 
  
  ## Impute for each cluster using HMisc (mean impute):
  
  ### This is an important step!
    for (l in 1:length(data_subs))
    {
      temp <- data_subs[[l]][,names(data_subs[[l]]) %in% col_names[[l]]] 
      temp2 <- impute_manual(temp)
      dat_inl[[l]] <- temp2
    }

  pars_te <- do.call(smartbind, dat_inl)
  
  pars2 <- lapply(col_names, function(x) {na.omit(pars_te[,names(pars_te)%in% x])}) 
  
  match_levels <- function(data_original, data_subset)
  {
    for(col in 1:ncol(data_subset))
    {
      if(is.factor(data_subset[,col]))
      {
        levels(data_subset[,col]) <- levels(data_original[,names(data_subset)[col]])  
      }
    }
    return(data_subset)
  }
  
  par <- lapply(pars2, function(x) {match_levels(data_train,x)}) 
  
  reg_models <- function(data_in)
  {
    X <- data_in[,!names(data_in) %in% c('cnt','instant','dteday','atemp')]
    Y <- data_in[,c('cnt')]
    datam <- data.table(X,Y) 
    datam <- match_levels(data_train,datam)
    fmla <- as.formula(paste("Y ~ ", paste(names(X), collapse= "+")))
    reg <- stepAIC(lm(fmla , data=datam),trace=F)
    return(reg)
  }
  
  reg_mas <- lapply(par,reg_models)
  
  mse_reg <- function(model) { sum((model$residuals)^2)/length(model) } 
  
  ## Check test data, to see which is the best matched model subset
  
  ## Clustering
  t_data_present <- data.frame( matrix(as.integer(!is.na(data_test)),nrow = nrow(data_test), ncol = ncol(data_test))) 

  names(t_data_present) <- names(data_test) 
  
  system.time(
  nearest_clus <- apply(t_data_present, 1, function(d)
  {
    temp <- apply(mcc_rounded,1,function(x) dist(rbind(x,d)))
    return(which.min(temp))
  }
  )
  )
  
  data_test_in <- data_test
  data_test_in$nearest_clus <- as.factor(nearest_clus)
  data_test_in$rownum <- 1:nrow(data_test_in)
  data_groups <- split(data_test_in,nearest_clus)
  get_back_row_order <- do.call("rbind",data_groups)
  
  data_igroups <- lapply(data_groups,function(x) impute_test(data_train,x))
  
  ## Write separate for OLS and GBM as syntax is different
  data_i1groups <- lapply(data_igroups,function(x) x[,!names(x) %in% c("nearest_clus")]) 
  data_i2groups <- lapply(data_i1groups,function(x) match_levels(data_train,x))
  
  reg_outcomes <- list()
  
  for(i in 1: length(data_i2groups))
  {   
    data_in <- data_i2groups[[i]]
    X <- data_in[,!names(data_in) %in% c('cnt','instant','dteday','atemp')]
    Y <- data_in[,c('cnt')]
    datam <- data.table(X,Y)
    datam <- match_levels(data_train,datam)
    reg_outcomes[[i]] <- predict(reg_mas[[i]],newdata=datam)
  }
  
  reg_out_a <- unlist(reg_outcomes)
  reg_ord <- reg_out_a[order(get_back_row_order[,"rownum"])]
  
  return(reg_ord)
}


#red_ord15K <- red_model(data_train, data_test)

ens_model <- function(data_train, data_test)
{
  ## Clustering
  data_present <- data.frame( matrix(as.integer(!is.na(data_train)),nrow = nrow(data_train), ncol = ncol(data_train))) 
  names(data_present) <- names(data_train) 
  
  ## This is a very simplified version where we know k as 4.
  
  determine_best_cluster <- function(data_in,clus_number)
  {
    n <- 10
    kmeans_o <- lapply(1:n,function(x) { set.seed(x); a<-kmeans(data_in,clus_number); return(a) } ) 
    opt_kmeans <- which.max(lapply(kmeans_o,function(x) {100-x$tot.withinss/x$totss*100})) 
    best_cluster <- kmeans_o[[opt_kmeans]]
    return(best_cluster)
  }	
  
  member_clusters <- determine_best_cluster(data_present,4)
  members <- member_clusters$cluster
  
  ## Once the cluster memberships are determined
  data_subs <- split(data_train,members)  
  
  col_names <- list()
  data_sets <- list()
  
  mcc <- as.data.frame(member_clusters$centers)
  mcc_rounded <- round(mcc,0)
  for(row in 1:nrow(mcc))
  {
    temp <- vector()
    for( col in 1:ncol(mcc)) 
    {
      if (mcc[row,col] > 0.1)
      {
        temp <- c(temp, colnames(mcc)[col] )  
      }
    }
    col_names[[row]] <- temp
  } 
  
  dat_inl <- list() 
  
  match_levels <- function(data_original, data_subset)
  {
    for(col in 1:ncol(data_subset))
    {
      if(is.factor(data_subset[,col]))
      {
        levels(data_subset[,col]) <- levels(data_original[,names(data_subset)[col]])  
      }
    }
    return(data_subset)
  }
  
  ## Impute for each cluster using HMisc (mean impute):
  
  ### This is an important step!
  for (l in 1:length(data_subs))
  {
    temp <- data_subs[[l]][,names(data_subs[[l]]) %in% col_names[[l]]] 
    temp2 <- impute_manual(temp)
    dat_inl[[l]] <- temp2
  }
  
  par <- lapply(dat_inl, function(x) {match_levels(data_train,x)}) 
  
  reg_models <- function(data_in)
  {
    X <- data_in[,!names(data_in) %in% c('cnt','instant','dteday','atemp')]
    Y <- data_in[,c('cnt')]
    datam <- data.table(X,Y) 
    datam <- match_levels(data_train,datam)
    fmla <- as.formula(paste("Y ~ ", paste(names(X), collapse= "+")))
    reg <- stepAIC(lm(fmla , data=datam),trace=F)
    reg$xlevels[['weathersit']] <- levels(data_train$weathersit)
    return(reg)
  }
  
  
  reg_mas <- lapply(par,reg_models)
  
  mse_reg <- function(model) { sum((model$residuals)^2)/length(model$residuals) } 
  
  data_test_imputed <- impute_test(data_train,data_test)
  data_test_imputed <- match_levels(data_train,data_test_imputed)
  
  ensemble_reg <- function(model_list,data_train,test_data)
  {
    #test_data <- match_levels(data_train,test_data)
    reg_ensemble_pred_1 <- 0
    reg_ensemble_pred_2 <- 0
    mse_list <- list()
    predict_model_list <- list()
    for(k in 1:length(model_list) ) ## Last one is baseline
    {
      
      m <- model_list[[k]]
      reg_ensemble_pred_1 <- reg_ensemble_pred_1 + (predict(m,newdata=test_data)*((length(m$residuals)*length(m$coefficients))^0.5/mse_reg(m)) )
      reg_ensemble_pred_2 <- reg_ensemble_pred_2 + ((length(m$residuals)*length(m$coefficients))^0.5/mse_reg(m))
      
      mse_list[[k]] <- mse_reg(m)
      predict_model_list[[k]] <-  predict(m,newdata=test_data)
    }                     
    
    reg_ensemble_pred <- reg_ensemble_pred_1/reg_ensemble_pred_2 
    reg_ensemble_best <- predict_model_list[[which.min(mse_list)]]
    return(list(reg_ensemble_pred,reg_ensemble_best))
  }
  
  reg_ens_out <- ensemble_reg(reg_mas[1:length(par)],data_train,data_test_imputed)	
  return(reg_ens_out)
}

## CART 

cart_model <- function(data_train,data_test)
{  
  #X <- data_train[,!names(data_train) %in% c('cnt','instant','dteday')]
  #Y <- data_train[,c('cnt')]
  #datam <- data.table(X,Y) 
  #train_control <- trainControl(method="cv",number=10, search = "random")  # 10 fold CV
  #rpart.grid <- expand.grid(cp=c(0.001,0.01,0.1,0.3)) 
  #model.cart <- train(Y~.,data=datam ,trControl=train_control, method="rpart"  ,tuneGrid=rpart.grid )
  #cart.model <- rpart(Y ~ .,          method="anova", data=datam) 
  cart.model <- rpart(cnt ~ .,          method="anova", data=data_train)  
  cart_pred <- predict(cart.model, data=data_test)
  return(cart_pred)
}

svi_model <- function(data_train,data_test)
{
  data_test_imputed <-  impute_test(data_train,data_test)  
  svi.data <- rfImpute(data_train[,!names(data_train)%in% c('cnt')],data_train[,names(data_train)%in% c('cnt')],ntree=50,iter=3)  
  names(svi.data)[1] <- "cnt"
  svi.model <- stepAIC(lm(cnt~.,data=svi.data),trace=F)
  svi.pred <- predict(svi.model, newdata=data_test_imputed)
  return(svi.pred)
}

mvi_model <- function(data_train,data_test)
{
  data_test_imputed <-  impute_test(data_train,data_test)  
  data_mi <- mice(data_train,m=5,maxit=3,meth='pmm',seed=1234) 
  completedData <- complete(data_mi,1)
  model_mi <- list()
  for(list in 1:data_mi$m)
  {
    model_mi[[list]] <- stepAIC(lm(cnt ~ .,data=complete(data_mi,list)))
  }
  model_meta_var <- data.frame(matrix(NA, nrow = nrow(completedData), ncol = data_mi$m )) 
  for(i in 1:data_mi$m)
  {
    model_meta_var[,i] <- predict(model_mi)
  }
  
  Y <- data_train$cnt
  X <- model_meta_var
  datam <- data.table(X,Y)
  rf.model_meta <- randomForest(Y~.,data=datam,importance =FALSE,ntree=100)
  
  ## Prediction is same procedure
  meta_var_pred <- data.frame(matrix(NA, nrow = nrow(data_test_imputed), ncol = data_mi$m)) 
  for(i in 1:data_mi$m)
  {
    meta_var_pred[,i] <- predict(model_mi,data_test_imputed)
  }
  
  rf_meta_pred <- predict(rf.model_meta,meta_var_pred)
  return(rf_meta_pred)
}

maec <- function(actual,predicted)
{ 
  mean(as.matrix(abs(actual-predicted))) ## Use your own function ##
}

rmsec <- function(actual,predicted)
{ 
  sqrt(mean((actual-predicted)^2))  ## Use your own function ##
}

red_ord <- list()
ens_ord <- list()
cart_ord <- list()
svi_ord <- list()
mvi_ord <- list()
mbi_ord <- list()

Time_red <- c(
  system.time(red_ord[[1]] <- red_model(tr_test_st(sc_data_list[[1]])[[1]], tr_test_st(sc_data_list[[1]])[[2]]))[3],
  system.time(red_ord[[2]] <- red_model(tr_test_st(sc_data_list[[2]])[[1]], tr_test_st(sc_data_list[[2]])[[2]]))[3],
  system.time(red_ord[[3]] <- red_model(tr_test_st(sc_data_list[[3]])[[1]], tr_test_st(sc_data_list[[3]])[[2]]))[3],
  system.time(red_ord[[4]] <- red_model(tr_test_st(sc_data_list[[4]])[[1]], tr_test_st(sc_data_list[[4]])[[2]]))[3],
  system.time(red_ord[[5]] <- red_model(tr_test_st(sc_data_list[[5]])[[1]], tr_test_st(sc_data_list[[5]])[[2]]))[3],
  system.time(red_ord[[6]] <- red_model(tr_test_st(sc_data_list[[6]])[[1]], tr_test_st(sc_data_list[[6]])[[2]]))[3],
  system.time(red_ord[[7]] <- red_model(tr_test_st(sc_data_list[[7]])[[1]], tr_test_st(sc_data_list[[7]])[[2]]))[3]
)

write.csv(Time_red, "BRM_model_times.csv")

Time_ens <- c(
  system.time(ens_ord[[1]] <- ens_model(tr_test_st(sc_data_list[[1]])[[1]], tr_test_st(sc_data_list[[1]])[[2]]))[3],
  system.time(ens_ord[[2]] <- ens_model(tr_test_st(sc_data_list[[2]])[[1]], tr_test_st(sc_data_list[[2]])[[2]]))[3],
  system.time(ens_ord[[3]] <- ens_model(tr_test_st(sc_data_list[[3]])[[1]], tr_test_st(sc_data_list[[3]])[[2]]))[3],
  system.time(ens_ord[[4]] <- ens_model(tr_test_st(sc_data_list[[4]])[[1]], tr_test_st(sc_data_list[[4]])[[2]]))[3],
  system.time(ens_ord[[5]] <- ens_model(tr_test_st(sc_data_list[[5]])[[1]], tr_test_st(sc_data_list[[5]])[[2]]))[3],
  system.time(ens_ord[[6]] <- ens_model(tr_test_st(sc_data_list[[6]])[[1]], tr_test_st(sc_data_list[[6]])[[2]]))[3],
  system.time(ens_ord[[7]] <- ens_model(tr_test_st(sc_data_list[[7]])[[1]], tr_test_st(sc_data_list[[7]])[[2]]))[3]
)

write.csv(Time_ens, "Ensembled_model_times.csv")

Time_cart <- c(
  system.time(cart_ord[[1]] <- cart_model(tr_test_st(sc_data_list[[1]])[[1]], tr_test_st(sc_data_list[[1]])[[2]]))[3],
  system.time(cart_ord[[2]] <- cart_model(tr_test_st(sc_data_list[[2]])[[1]], tr_test_st(sc_data_list[[2]])[[2]]))[3],
  system.time(cart_ord[[3]] <- cart_model(tr_test_st(sc_data_list[[3]])[[1]], tr_test_st(sc_data_list[[3]])[[2]]))[3],
  system.time(cart_ord[[4]] <- cart_model(tr_test_st(sc_data_list[[4]])[[1]], tr_test_st(sc_data_list[[4]])[[2]]))[3],
  system.time(cart_ord[[5]] <- cart_model(tr_test_st(sc_data_list[[5]])[[1]], tr_test_st(sc_data_list[[5]])[[2]]))[3],
  system.time(cart_ord[[6]] <- cart_model(tr_test_st(sc_data_list[[6]])[[1]], tr_test_st(sc_data_list[[6]])[[2]]))[3],
  system.time(cart_ord[[7]] <- cart_model(tr_test_st(sc_data_list[[7]])[[1]], tr_test_st(sc_data_list[[7]])[[2]]))[3]
)

write.csv(Time_cart, "CART_times.csv")

Time_mbi <- c(
  system.time(mbi_ord[[1]] <- mbi_predictions(tr_test_st(sc_data_list[[1]])[[1]], tr_test_st(sc_data_list[[1]])[[2]]))[3],
  system.time(mbi_ord[[2]] <- mbi_predictions(tr_test_st(sc_data_list[[2]])[[1]], tr_test_st(sc_data_list[[2]])[[2]]))[3],
  system.time(mbi_ord[[3]] <- mbi_predictions(tr_test_st(sc_data_list[[3]])[[1]], tr_test_st(sc_data_list[[3]])[[2]]))[3],
  system.time(mbi_ord[[4]] <- mbi_predictions(tr_test_st(sc_data_list[[4]])[[1]], tr_test_st(sc_data_list[[4]])[[2]]))[3],
  system.time(mbi_ord[[5]] <- mbi_predictions(tr_test_st(sc_data_list[[5]])[[1]], tr_test_st(sc_data_list[[5]])[[2]]))[3],
  system.time(mbi_ord[[6]] <- mbi_predictions(tr_test_st(sc_data_list[[6]])[[1]], tr_test_st(sc_data_list[[6]])[[2]]))[3],
  system.time(mbi_ord[[7]] <- mbi_predictions(tr_test_st(sc_data_list[[7]])[[1]], tr_test_st(sc_data_list[[7]])[[2]]))[3]
)
### Sometimes all of them do not run, but separately, we can run for each mbi object separately.
## In one run, we got c(1.25, 2.3, 3.15, 11.22, 24.49, 129.64, 262.53 )

Time_svi <- vector()

for(i in 1:7)
{
  try(temp <- system.time(svi_ord[[i]] <- svi_model(tr_test_st(sc_data_list[[i]])[[1]], tr_test_st(sc_data_list[[i]])[[2]]))[3])
Time_svi <- c(Time_svi,temp)
}  

Time_mvi <- vector()

i <- 1
temp <- try(temp <- system.time(mvi_ord[[i]] <- mvi_model(tr_test_st(sc_data_list[[i]])[[1]], tr_test_st(sc_data_list[[i]])[[2]]))[3]) 
temp
Time_mvi <- c(Time_mvi,temp)

i <- i + 1
temp <- try(temp <- system.time(mvi_ord[[i]] <- mvi_model(tr_test_st(sc_data_list[[i]])[[1]], tr_test_st(sc_data_list[[i]])[[2]]))[3]) 
temp
Time_mvi <- c(Time_mvi,temp)
# 58.42
i <- i + 1
temp <- try(temp <- system.time(mvi_ord[[i]] <- mvi_model(tr_test_st(sc_data_list[[i]])[[1]], tr_test_st(sc_data_list[[i]])[[2]]))[3]) 
temp
Time_mvi <- c(Time_mvi,temp)

i <- i + 1
temp <- try(temp <- system.time(mvi_ord[[i]] <- mvi_model(tr_test_st(sc_data_list[[i]])[[1]], tr_test_st(sc_data_list[[i]])[[2]]))[3]) 
temp
Time_mvi <- c(Time_mvi,temp)


