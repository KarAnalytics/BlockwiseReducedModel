#### Evaluating the blockwise reduced modeling (BRM) method on UCI bike sharing dataset

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
  }  
  
  gbm_wt_parameters <- lapply(models,global_gbm_par)
  
  weighted_gbm <- do.call(smartbind,lapply(gbm_wt_parameters, function(x) t(x)))
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

# sapply(bike, function(y) sum(length(which(is.na(y)))))
# sapply(dat_list[[1]], function(y) sum(length(which(is.na(y)))))
# sapply(dat_list[[2]], function(y) sum(length(which(is.na(y)))))
# sapply(dat_list[[3]], function(y) sum(length(which(is.na(y)))))
# sapply(dat_list[[4]], function(y) sum(length(which(is.na(y)))))

data_train_no_miss <- tr_test_st(bike)[[1]]
data_test_no_miss <- tr_test_st(bike)[[2]]

# reg_all <- stepAIC(lm( cnt~., data=data_train_no_miss),direction="both",trace=FALSE)
reg_ps <- stepAIC(glm( cnt~., family= "poisson", data=data_train_no_miss),direction="both",trace=FALSE)
summary(reg_ps)
# exp(coef(reg_ps))

#hist(data_train_no_miss$cnt)
#mean (data_train_no_miss$cnt)
#var(data_train_no_miss$cnt)

reg_all<- glm.nb( cnt~., data=data_train_no_miss)
summary(reg_all)

reg_all_pred <- predict(reg_all, type="response")

R_sq_all <- 1 - (sum((data_train_no_miss[,"cnt"] - reg_all_pred)^2)/(sum((data_train_no_miss[,"cnt"] - mean(data_train_no_miss[,"cnt"]))^2)))


chiq = 2 * (logLik(reg_p) - logLik(reg_all))
pchisq( chiq, df = 1, lower.tail = FALSE)

temp_insample <- predict(reg_p,type="response")

gbm_all <- gbm(cnt~., data=data_train_no_miss, distribution ="poisson", verbose=FALSE,n.trees = 500,shrinkage=0.1,interaction.depth = 3)

temp <- predict(gbm_all,newdata = data_train_no_miss,n.trees=gbm_all$n.trees,type="response")

head(data_train_no_miss[,"cnt"])
head(temp)
### 10% missing:
# data_train <- tr_test_st(dat_list[[3]])[[1]]
# data_test <- tr_test_st(dat_list[[3]])[[2]]

### 50% missing:
data_train <- tr_test_st(dat_list[[7]])[[1]]
data_test <- tr_test_st(dat_list[[7]])[[2]]


data_test_imputed <- impute_test(data_train,data_test)
data_train_imputed <- impute_manual(data_train)

data_train_y <- data.frame(data_train[,"cnt"])
data_train_x <- data_train[,!names(data_train) %in% c("cnt")]
names(data_train_y) = "cnt" 

missing_prop_val <- brm_num_blocks(data_train_x,low_threshold = 0.05)  

### We assume four blocks identified here!
plot(missing_prop_val)
imputation_thres <- 0.05

num_blocks <- min(which(missing_prop_val <= imputation_thres))
#num_blocks <- 4

### TO DO:
## When to use BRM and when not to? With a quantitative measure for the dataset. 
## Number of blocks identified automatically using the plot (derivative or something) 

## We decide how many blocks to consider based on plot

brm_model <- brm_ov_subsets(num_blocks,data_train_x, data_train_y,low_threshold = 0.05)

ov_subsets <- brm_model[[1]]
cluster_centers_subsets <- brm_model[[2]]
n_ov_subsets <- brm_model[[3]]
subsets_columns <- brm_model[[4]]

## From here standardized code for R stops (for now)!

reg_m <- lapply(ov_subsets,stepreg)
gbm_m <- lapply(ov_subsets,gbm_train)

rsquares_list <- unlist(model_fit_comparisons_reg(reg_m, ov_subsets ))

pred_reg_outcomes <- pred_reg(reg_m, data_test,data_train_x, data_train_y,cluster_centers_subsets)
pred_reg_insample <- pred_reg(reg_m, data_train,data_train_x, data_train_y,cluster_centers_subsets)
pred_per_reg <- c(rmsec(data_test[,"cnt"],pred_reg_outcomes), maec(data_test[,"cnt"],pred_reg_outcomes), smapec(data_test[,"cnt"],pred_reg_outcomes))
pred_per_reg_insample <- c(rmsec(data_train[,"cnt"],pred_reg_insample), maec(data_train[,"cnt"],pred_reg_insample), smapec(data_train[,"cnt"],pred_reg_insample))
global_wts_reg <- get_global_wts(reg_m)

pred_gbm_outcomes <- pred_reg(gbm_m, data_test,data_train_x, data_train_y,cluster_centers_subsets)
pred_gbm_insample <- pred_reg(gbm_m, data_train,data_train_x, data_train_y,cluster_centers_subsets)
pred_per_gbm <- c(rmsec(data_test[,"cnt"],pred_gbm_outcomes), maec(data_test[,"cnt"],pred_gbm_outcomes), smapec(data_test[,"cnt"],pred_gbm_outcomes))
pred_per_gbm_insample <- c(rmsec(data_train[,"cnt"],pred_gbm_insample), maec(data_train[,"cnt"],pred_gbm_insample), smapec(data_train[,"cnt"],pred_gbm_insample))

## BRM execution codes in one place as a function:

brm_predictions <- function(num_blocks,X_in, Y_in, data_test,low_threshold = 0.05)
{
  brm_model <- brm_ov_subsets(num_blocks,X_in, Y_in,low_threshold)
  ov_subsets <- brm_model[[1]]
  cluster_centers_subsets <- brm_model[[2]]
  n_ov_subsets <- brm_model[[3]]
  subsets_columns <- brm_model[[4]]
  reg_m <- lapply(ov_subsets,stepreg)
  gbm_m <- lapply(ov_subsets,gbm_train)
  pred_reg_outcomes <- pred_reg(reg_m, data_test,data_train_x, data_train_y,cluster_centers_subsets)
  pred_gbm_outcomes <- pred_reg(gbm_m, data_test,data_train_x, data_train_y,cluster_centers_subsets)
  return(list(pred_reg_outcomes,pred_gbm_outcomes))
}

## Compare this with the our own weighted method:
feature_imp_mat <- feature_drop_loss(ov_subsets,data_train_x, data_train_y, gbm_train, gbm_m, cluster_centers_subsets)
global_wts_gbm <- get_global_wts_gbm(gbm_m)

### SE and p-value using semi-parametric bootstrapping

boot_n <- 100
resid <- data_train[,"cnt"] - pred_reg_insample
var_list <- list()
system.time(
for(s in 1:boot_n)
          {
              set.seed(s*10)
              newY <- pred_reg_insample + sample(resid, size = nrow(data_train), replace = TRUE)
              Y_dat <- data.frame(round(newY))
              names(Y_dat) <- names(data_train_y) 
              Y_dat[Y_dat<= round(min(pred_reg_insample))] <- round(min(pred_reg_insample))
              brm_model_bs <- brm_ov_subsets(num_blocks,data_train_x, Y_dat, low_threshold = 0.05)
              ov_subsets_bs <- brm_model_bs[[1]]
              cluster_centers_subsets_bs <- brm_model_bs[[2]]
              reg_m_bs <- lapply(ov_subsets_bs,stepreg)
              var_list[[s]] <- data.frame(t(get_global_wts(reg_m_bs))) 
         }
)
### 100 seconds for bootstrapp = 100

summary(reg_m[[1]])
summary(reg_m[[2]])

var_list_t0 <- lapply(var_list, function(x) data.frame(t(x)))
# var_list_t <- do.call(rbind,var_list)
#var_list_t <- smartbind(var_list_t0[[1]],var_list_t0[[2]] )
var_list_t <- do.call(smartbind,var_list_t0)
var_means <- apply(var_list_t, 2, function(x) mean(x, na.rm=T) )
var_se <- apply(var_list_t, 2, function(x) sd(x, na.rm=T))
var_t <- (var_means/var_se)
pvals = 2*pt(-abs(var_t), df=nrow(var_list_t)-1)

pvals_dat = data.frame(var = names(pvals), pval= pvals)

global_wts_df = data.frame(var = names(get_global_wts(reg_m)), coef = get_global_wts(reg_m))

global_wts_df = left_join(global_wts_df, pvals_dat, by =c("var"="var"))

reg_all$coefficients  ## Quite close!

## R-squared: https://internal.ncl.ac.uk/ask/numeracy-maths-statistics/statistics/regression-and-correlation/coefficient-of-determination-r-squared.html
R_sq = 1 - (sum((data_train[,"cnt"] - pred_reg_insample)^2)/(sum((data_train[,"cnt"] - mean(data_train[,"cnt"]))^2)))

#####
### Benchmarks: MBI, Ensemble, SVI, MVI, CART for bike sharing:
#####

## MBI: Steps - Find the blocks (using our method, then the different J(j,k) and then do the weighted imputation)

## Subset by Subset:
                
## n_ov_subsets[[1]] has "hr"         "weathersit" "temp"       "hum"        "windspeed"  "ToD" missing
## n_ov_subsets[[2]] has no missing
## n_ov_subsets[[3]] has "weathersit" "temp"       "hum" missing
## n_ov_subsets[[4]] has "hr"        "windspeed" "ToD" missing

names(data_train)[! names(data_train) %in% names(n_ov_subsets[[1]]) ]
names(data_train)[! names(data_train) %in% names(n_ov_subsets[[2]]) ]
names(data_train)[! names(data_train) %in% names(n_ov_subsets[[3]]) ]
names(data_train)[! names(data_train) %in% names(n_ov_subsets[[4]]) ]

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

# mbi_imputed_data <- MBI_model(ov_subsets)
mbi_imputed_data <- MBI_model(n_ov_subsets)

#reg_mbi <- stepAIC(glm( cnt~., family= "poisson", data=mbi_imputed_data))
reg_mbi <- stepAIC(glm.nb( cnt~., data=mbi_imputed_data))

gbm_mbi <- gbm(cnt~., data=mbi_imputed_data, distribution ="poisson", verbose=FALSE,n.trees = 500,shrinkage=0.1,interaction.depth = 3)

mbi_predictions <- function(X_in, Y_in, data_test)
{
  data_xy <- cbind(X_in, Y_in)
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
  
  #reg_mbi <- stepAIC(lm( cnt~., data=mbi_imputed_data))
  # reg_mbi <- glm( cnt~., family= "poisson", data=mbi_imputed_data)
  # reg_mbi <- glm.nb( cnt~., data=mbi_imputed_data)
  reg_mbi <- glm.nb( cnt~., data=mbi_imputed_data)
  #gbm_mbi <- gbm(cnt~., data=mbi_imputed_data, verbose=FALSE,n.trees = 500,shrinkage=0.1,interaction.depth = 3)
  
  reg_pred_mbi <- predict(reg_mbi, data_test_imputed,type="response")
  #gbm_pred_mbi <- predict(gbm_mbi, data_test_imputed,n.trees=gbm_mbi$n.trees)
  # return(list(reg_pred_mbi, gbm_pred_mbi))
  return(reg_pred_mbi)
}

#pred_per_mbi_reg <- c(rmsec(data_test_imputed[,"cnt"],predict(reg_mbi, data_test_imputed, type= "response")), maec(data_test_imputed[,"cnt"],predict(reg_mbi, data_test_imputed)), smapec(data_test_imputed[,"cnt"],predict(reg_mbi, data_test_imputed)))
# pred_per_mbi_gbm <- c(rmsec(data_test_imputed[,"cnt"],predict(gbm_mbi, data_test_imputed,n.trees=gbm_mbi$n.trees)), maec(data_test_imputed[,"cnt"],predict(gbm_mbi, data_test_imputed,n.trees=gbm_mbi$n.trees)), smapec(data_test_imputed[,"cnt"],predict(gbm_mbi, data_test_imputed,n.trees=gbm_mbi$n.trees)))

## No missing:
#pred_per_reg_all <- c(rmsec(data_test_no_miss[,"cnt"],predict(reg_all, data_test_no_miss)), maec(data_test_no_miss[,"cnt"],predict(reg_all, data_test_no_miss)), smapec(data_test_no_miss[,"cnt"],predict(reg_all, data_test_no_miss)))
#pred_per_gbm_all <- c(rmsec(data_test_no_miss[,"cnt"],predict(gbm_all, data_test_no_miss,n.trees=gbm_all$n.trees)), maec(data_test_no_miss[,"cnt"],predict(gbm_all, data_test_no_miss,n.trees=gbm_all$n.trees)), smapec(data_test_no_miss[,"cnt"],predict(gbm_all, data_test_no_miss,n.trees=gbm_all$n.trees)))

## List-wise deletion

data_train_listwise <- na.omit(data_train)
# reg_listwise <- stepAIC(glm( cnt~., family= "poisson", data=data_train_listwise))
reg_listwise <- stepAIC(glm.nb( cnt~., data=data_train_listwise))
gbm_listwise <- gbm(cnt~., data=data_train_listwise, distribution ="poisson", verbose=FALSE,n.trees = 500,shrinkage=0.1,interaction.depth = 3)

listwise_model <- function(data_train, data_test)
{
  data_train_listwise <- na.omit(data_train)
  # reg_listwise <- stepAIC(lm( cnt~., data=data_train_listwise))
  # reg_listwise <- stepAIC(glm( cnt~., family = "poisson",data=data_train_listwise))
  reg_listwise <- stepAIC(glm.nb( cnt~.,data=data_train_listwise))
  gbm_listwise <- gbm(cnt~., data=data_train_listwise, distribution ="poisson", verbose=FALSE,n.trees = 500,shrinkage=0.1,interaction.depth = 3)
  #data_test_listwise <- na.omit(data_test)
  data_test_imputed <- impute_test(data_train_listwise,data_test)
  listwise_pred_reg <- predict(reg_listwise, data_test_imputed, type = "response")
  listwise_pred_gbm <- predict(gbm_listwise, data_test_imputed, type = "response" ,n.trees=gbm_listwise$n.trees)
  return(list(listwise_pred_reg, listwise_pred_gbm))
}

#pred_per_reg_listwise <- c(rmsec(data_test_imputed[,"cnt"],predict(reg_listwise, data_test_imputed,type="response")), maec(data_test_imputed[,"cnt"],predict(reg_listwise, data_test_imputed)), smapec(data_test_imputed[,"cnt"],predict(reg_listwise, data_test_imputed)))
#pred_per_gbm_listwise <- c(rmsec(data_test_imputed[,"cnt"],predict(gbm_listwise, data_test_imputed,type="response",n.trees=gbm_listwise$n.trees)), maec(data_test_imputed[,"cnt"],predict(gbm_listwise, data_test_imputed,n.trees=gbm_listwise$n.trees)), smapec(data_test_imputed[,"cnt"],predict(gbm_listwise, data_test_imputed,n.trees=gbm_listwise$n.trees)))

## Ensemble:
ens_model <- function(num_blocks, X_in, Y_in,test_data, low_threshold = 0.05)
{
  
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
  
    data_xy <- cbind(X_in, Y_in)
    data_test_imputed <- impute_test(data_xy,test_data)
    data_test_imputed <- match_levels(data_xy,data_test_imputed)
    
    data_xy <- match_levels(data_test_imputed,data_xy)
    
    data_present <- data.frame( matrix(as.integer(!is.na(X_in)),nrow = nrow(X_in), ncol = ncol(X_in))) 
    names(data_present) <- names(X_in) 
    
    columns_in_subsets <- function(data_xy, MVI,num_blocks, low_threshold )
    {
      columns <- list()
      determine_best_cluster <- function(data_in,clus_number)
      {
        n <- 10
        kmeans_o <- lapply(1:n,function(x) { set.seed(x); a <- kmeans(data_in,clus_number); return(a) } ) 
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
    
    n_ov_subsetting <- columns_in_subsets(data_xy,data_present,num_blocks, low_threshold)
    cluster_centers_subsets <- n_ov_subsetting[[1]]
    columns_n_ov_subsets <- n_ov_subsetting[[2]]
    n_ov_subsets <- n_ov_subsetting[[3]]
    
    n_ov_subsets <- lapply(n_ov_subsets, impute_manual) 
    
    stepreg_ens <- function(data_in)
    {
      X <- data_in[,names(data_in) %in% names(X_in) ]
      Y <- data_in[,names(data_in) %in% names(Y_in) ]
      datam <- data.table(X,Y)
      fmla <- as.formula(paste("Y ~ ", paste(names(X), collapse= "+")))
      # reg <- stepAIC(glm( fmla, data=data_in, family= "poisson"),direction="both",trace=FALSE)
      reg <- stepAIC(glm.nb( fmla, data=data_in),direction="both",trace=FALSE)
      return(reg)
    }
    
    gbm_train_ens <- function(data_in)
    {
      X <- data_in[,names(data_in) %in% names(X_in) ]
      Y <- data_in[,names(data_in) %in% names(Y_in) ]
      datam <- data.table(X,Y) 
      fmla <- as.formula(paste("Y ~ ", paste(names(X), collapse= "+")))
      reg <- gbm( fmla, data=datam ,verbose=FALSE, distribution = "poisson", n.trees = 500,shrinkage=0.1,interaction.depth = 3)
      return(reg)
    }
    
  reg_ens <- lapply(n_ov_subsets,stepreg_ens)
  gbm_ens <- lapply(n_ov_subsets,gbm_train_ens)
  
  #mse_reg <- function(model) { sum((model$residuals)^2)/length(model$residuals) } 
  mse_reg <- function(model) { mean((model$residuals)^2) }
  mse_gbm <- function(model) {model$train.error[length(model$train.error)]}
  
  ensemble_reg <- function(model_list,test_data)
  {
    reg_ensemble_pred_1 <- 0
    reg_ensemble_pred_2 <- 0
    mse_list <- list()
    predict_model_list <- list()
    for(k in 1:length(model_list) ) ## Last one is baseline
    {
      m <- model_list[[k]]
      
      if(class(m)[1] == "negbin")
      {
        try(reg_ensemble_pred_1 <- reg_ensemble_pred_1 + (predict(m,newdata=test_data,type="response")*((length(m$residuals)*length(m$coefficients))^0.5/mse_reg(m)) ) ) 
        reg_ensemble_pred_2 <- reg_ensemble_pred_2 + ((length(m$residuals)*length(m$coefficients))^0.5/mse_reg(m))
        mse_list[[k]] <- mse_reg(m)
        predict_model_list[[k]] <-  predict(m,newdata=test_data,type="response")
      }else if(class(m) == "gbm"){
        try(reg_ensemble_pred_1 <- reg_ensemble_pred_1 + (predict(m,test_data,n.trees=m$n.trees,type="response")*((length(m$var.names)*length(m$fit))^0.5/(mse_reg(m))) ) )
        reg_ensemble_pred_2 <- reg_ensemble_pred_2 + ((length(m$var.names)*length(m$fit))^0.5/(mse_gbm(m)))
        mse_list[[k]] <- mse_gbm(m)
        predict_model_list[[k]] <-  predict(m,test_data,n.trees=m$n.trees,type="response")
      }
      
    }                     
    
    if(class(model_list[[1]])[1] == "negbin")
    {
      return(reg_ensemble_pred_1/reg_ensemble_pred_2)
    }
    else if(class(model_list[[1]]) == "gbm")
    {
      return(predict_model_list[[which.min(mse_list)]])
    }
    # reg_ensemble_best <- predict_model_list[[which.min(mse_list)]]
    # reg_ensemble_pred <- reg_ensemble_pred_1/reg_ensemble_pred_2 
     # return(reg_ensemble_best)
  }
  
  reg_ens_out <- ensemble_reg(reg_ens,data_test_imputed)
  gbm_ens_out <- ensemble_reg(gbm_ens,data_test_imputed)
  return(list(reg_ens_out, gbm_ens_out))
}

pred_ens <- ens_model(num_blocks, data_train_x, data_train_y,data_test,low_threshold = 0.05)
ens_reg <- pred_ens[[1]]
ens_gbm <- pred_ens[[2]]

pred_ens_reg <- c(rmsec(data_test_imputed[,"cnt"],ens_reg), maec(data_test_imputed[,"cnt"],ens_reg), smapec(data_test_imputed[,"cnt"],ens_reg))
pred_per_gbm <- c(rmsec(data_test_imputed[,"cnt"],ens_gbm), maec(data_test_imputed[,"cnt"],ens_gbm), smapec(data_test_imputed[,"cnt"],ens_gbm))

### MI:

mi_model <- function(data_train,data_test)
{
  data_test_imputed <-  impute_test(data_train,data_test)  
  data_mi <- mice(data_train,m=5,maxit=3,meth='pmm',seed=1234) 
  completedData <- complete(data_mi,1)
  model_mi_gbm <- list()
  model_mi_reg <- list()
  for(list in 1:data_mi$m)
  {
    model_mi_gbm[[list]] <- gbm( cnt~., distribution = "poisson", data=complete(data_mi,list) ,verbose=FALSE,n.trees = 500,shrinkage=0.1,interaction.depth = 3)
    # model_mi_reg[[list]] <- glm( cnt~., family="poisson", data=complete(data_mi,list) )
    model_mi_reg[[list]] <- glm.nb( cnt~.,  data=complete(data_mi,list) ) 
  }
  model_meta_var_gbm <- data.frame(matrix(NA, nrow = nrow(completedData), ncol = data_mi$m ))
  model_meta_var_reg <- data.frame(matrix(NA, nrow = nrow(completedData), ncol = data_mi$m ))
  for(i in 1:data_mi$m)
  {
    model_meta_var_gbm[,i] <- predict(model_mi_gbm[[i]],n.trees = model_mi_gbm[[i]]$n.trees,type="response")
    model_meta_var_reg[,i] <- predict(model_mi_reg[[i]],type="response")
  }
  
  Y <- data_train$cnt
  Xgbm <- model_meta_var_gbm
  datam1 <- data.table(Xgbm,Y)
  rf.model_meta_gbm <- randomForest(Y~.,data=datam1,importance =FALSE,ntree=50)

  Xreg <- model_meta_var_reg
  datam2 <- data.table(Xreg,Y)
  rf.model_meta_reg <- randomForest(Y~.,data=datam2,importance =FALSE,ntree=50)
  
  ## Prediction is same procedure
  meta_var_pred_gbm <- data.frame(matrix(NA, nrow = nrow(data_test_imputed), ncol = data_mi$m)) 
  meta_var_pred_reg <- data.frame(matrix(NA, nrow = nrow(data_test_imputed), ncol = data_mi$m)) 
  for(i in 1:data_mi$m)
  {
    meta_var_pred_gbm[,i] <- predict(model_mi_gbm[[i]],data_test_imputed,type="response",n.trees = model_mi_gbm[[i]]$n.trees)
    meta_var_pred_reg[,i] <- predict(model_mi_reg[[i]],data_test_imputed,type="response")
  }
  
  rf_meta_pred_gbm <- predict(rf.model_meta_gbm,meta_var_pred_gbm)
  rf_meta_pred_reg <- predict(rf.model_meta_reg,meta_var_pred_reg)
  return(list(rf_meta_pred_reg,rf_meta_pred_gbm))
}

### SVI:
svi_model <- function(data_train,data_test)
{
  data_test_imputed <-  impute_test(data_train,data_test)  
  svi.data <- rfImpute(data_train[,!names(data_train)%in% c('cnt')],data_train[,names(data_train)%in% c('cnt')],ntree=50,iter=3)  
  names(svi.data)[1] <- "cnt"
  svi.model_gbm <- gbm( cnt~., data=svi.data , distribution = "poisson",verbose=FALSE,n.trees = 500,shrinkage=0.1,interaction.depth = 3) 
  svi.pred_gbm <- predict(svi.model_gbm, type="response",  newdata=data_test_imputed,n.trees = svi.model_gbm$n.trees)
  # svi.model_reg <- glm( cnt~., family="poisson", data=svi.data) 
  svi.model_reg <- glm.nb( cnt~.,  data=svi.data) 
  svi.pred_reg <- predict(svi.model_reg,type="response", newdata=data_test_imputed)
  return(list(svi.pred_reg, svi.pred_gbm))
}

##Evaluation:
### CART:

cart_model <- function(data_train, data_test)
{
  cart_m <- rpart(cnt ~ .,          method="anova", data=data_train )
  cart_prediction <- predict(cart_m, data_test)
  return(cart_prediction)
}

system.time( pred_mi <- mi_model(data_train,data_test) )
system.time( pred_svi <- svi_model(data_train,data_test) ) 
system.time( pred_cart <- cart_model(data_train,data_test) ) 

mi_reg <- pred_mi[[1]]
mi_gbm <- pred_mi[[2]]

pred_mi_reg <- c(rmsec(data_test_imputed[,"cnt"],mi_reg), maec(data_test_imputed[,"cnt"],mi_reg), smapec(data_test_imputed[,"cnt"],mi_reg))
pred_per_gbm <- c(rmsec(data_test_imputed[,"cnt"],mi_gbm), maec(data_test_imputed[,"cnt"],mi_gbm), smapec(data_test_imputed[,"cnt"],mi_gbm))

svi_reg <- pred_svi[[1]]
svi_gbm <- pred_svi[[2]]

pred_svi_reg <- c(rmsec(data_test_imputed[,"cnt"],svi_reg), maec(data_test_imputed[,"cnt"],svi_reg), smapec(data_test_imputed[,"cnt"],svi_reg))
pred_per_gbm <- c(rmsec(data_test_imputed[,"cnt"],svi_gbm), maec(data_test_imputed[,"cnt"],svi_gbm), smapec(data_test_imputed[,"cnt"],svi_gbm))

pred_per_cart <- c(rmsec(data_test[,"cnt"],pred_cart), maec(data_test[,"cnt"],pred_cart), smapec(data_test[,"cnt"],pred_cart))

#### All predictions in one place

predperf <- function(predictions, test_dat)
{
  predictions <- round(predictions)  ## Since they are count
  perf <- c(rmsec(test_dat[,"cnt"],predictions), maec(test_dat[,"cnt"],predictions), smapec(test_dat[,"cnt"],predictions), mapec(test_dat[,"cnt"],predictions))
  return(perf)    
}

p_brm_reg <- matrix(NA, nrow = 8, ncol = 4)
p_brm_gbm <- matrix(NA, nrow = 8, ncol = 4)
p_mbi_reg <- matrix(NA, nrow = 8, ncol = 4)
p_ens_reg <- matrix(NA, nrow = 8, ncol = 4)
p_ens_gbm <- matrix(NA, nrow = 8, ncol = 4)
p_listwise_reg <- matrix(NA, nrow = 8, ncol = 4)
p_listwise_gbm <- matrix(NA, nrow = 8, ncol = 4)
p_svi_reg <- matrix(NA, nrow = 8, ncol = 4)
p_svi_gbm <- matrix(NA, nrow = 8, ncol = 4)
p_mi_reg <- matrix(NA, nrow = 8, ncol = 4)
p_mi_gbm <- matrix(NA, nrow = 8, ncol = 4)
p_cart <- matrix(NA, nrow = 8, ncol = 4)

brm_pred <- list()
ens_pred <- list()
mbi_pred <- list()
listwise_pred <- list()
mi_pred <- list()
svi_pred <- list()
cart_pred <- list()

data_train <- tr_test_st(dat_list[[1]])[[1]]
data_test <- tr_test_st(dat_list[[1]])[[2]]
data_train_y <- data.frame(data_train[,"cnt"])
data_train_x <- data_train[,!names(data_train) %in% c("cnt")]
names(data_train_y) = "cnt" 
data_test_imputed <- impute_test(data_train,data_test)
data_train_imputed <- impute_manual(data_train)


brm_pred[[1]] <- brm_predictions(num_blocks,data_train_x, data_train_y, data_test)
p_brm_reg[1,] <- predperf(brm_pred[[1]][[1]], data_test)
p_brm_gbm[1,] <- predperf(brm_pred[[1]][[2]], data_test)

ens_pred[[1]] <- ens_model(num_blocks,data_train_x, data_train_y, data_test,low_threshold = 0.05)
p_ens_reg[1,] <- predperf(ens_pred[[1]][[1]], data_test_imputed)
p_ens_gbm[1,] <- predperf(ens_pred[[1]][[2]], data_test_imputed)

mbi_pred[[1]] <- mbi_predictions(data_train_x, data_train_y, data_test)
p_mbi_reg[1,] <- predperf(mbi_pred[[1]], data_test_imputed)

#data_test_listwise <- na.omit(data_test)
data_test_listwise <- data_test
listwise_pred[[1]] <- listwise_model(data_train, data_test_listwise)
p_listwise_reg[1,] <- predperf(listwise_pred[[1]][[1]], data_test_listwise)
p_listwise_gbm[1,] <- predperf(listwise_pred[[1]][[2]], data_test_listwise)

mi_pred[[1]] <- mi_model(data_train, data_test)
p_mi_reg[1,] <- predperf(mi_pred[[1]][[1]], data_test_imputed)
p_mi_gbm[1,] <- predperf(mi_pred[[1]][[2]], data_test_imputed)

svi_pred[[1]] <- svi_model(data_train, data_test)
p_svi_reg[1,] <- predperf(svi_pred[[1]][[1]], data_test_imputed)
p_svi_gbm[1,] <- predperf(svi_pred[[1]][[2]], data_test_imputed)

cart_pred[[1]] <- cart_model(data_train, data_test)
p_cart[1,] <- predperf(cart_pred[[1]], data_test_imputed)

#### 2 ####

data_train <- tr_test_st(dat_list[[2]])[[1]]
data_test <- tr_test_st(dat_list[[2]])[[2]]
data_train_y <- data.frame(data_train[,"cnt"])
data_train_x <- data_train[,!names(data_train) %in% c("cnt")]
names(data_train_y) = "cnt" 
data_test_imputed <- impute_test(data_train,data_test)
data_train_imputed <- impute_manual(data_train)

brm_pred[[2]] <- brm_predictions(num_blocks,data_train_x, data_train_y, data_test)
p_brm_reg[2,] <- predperf(brm_pred[[2]][[1]], data_test)
p_brm_gbm[2,] <- predperf(brm_pred[[2]][[2]], data_test)

ens_pred[[2]] <- ens_model(num_blocks,data_train_x, data_train_y, data_test,low_threshold = 0.05)
p_ens_reg[2,] <- predperf(ens_pred[[2]][[1]], data_test_imputed)
p_ens_gbm[2,] <- predperf(ens_pred[[2]][[2]], data_test_imputed)

mbi_pred[[2]] <- mbi_predictions(data_train_x, data_train_y, data_test)
p_mbi_reg[2,] <- predperf(mbi_pred[[2]], data_test_imputed)

#data_test_listwise <- na.omit(data_test)
data_test_listwise <- data_test
listwise_pred[[2]] <- listwise_model(data_train, data_test_listwise)
p_listwise_reg[2,] <- predperf(listwise_pred[[2]][[1]], data_test_listwise)
p_listwise_gbm[2,] <- predperf(listwise_pred[[2]][[2]], data_test_listwise)

mi_pred[[2]] <- mi_model(data_train, data_test)
p_mi_reg[2,] <- predperf(mi_pred[[2]][[1]], data_test_imputed)
p_mi_gbm[2,] <- predperf(mi_pred[[2]][[2]], data_test_imputed)

svi_pred[[2]] <- svi_model(data_train, data_test)
p_svi_reg[2,] <- predperf(svi_pred[[2]][[1]], data_test_imputed)
p_svi_gbm[2,] <- predperf(svi_pred[[2]][[2]], data_test_imputed)

cart_pred[[2]] <- cart_model(data_train, data_test)
p_cart[2,] <- predperf(cart_pred[[2]], data_test_imputed)

#### subset 3

data_train <- tr_test_st(dat_list[[3]])[[1]]
data_test <- tr_test_st(dat_list[[3]])[[2]]
data_train_y <- data.frame(data_train[,"cnt"])
data_train_x <- data_train[,!names(data_train) %in% c("cnt")]
names(data_train_y) = "cnt" 
data_test_imputed <- impute_test(data_train,data_test)
data_train_imputed <- impute_manual(data_train)

brm_pred[[3]] <- brm_predictions(num_blocks,data_train_x, data_train_y, data_test)
p_brm_reg[3,] <- predperf(brm_pred[[3]][[1]], data_test)
p_brm_gbm[3,] <- predperf(brm_pred[[3]][[2]], data_test)

ens_pred[[3]] <- ens_model(num_blocks,data_train_x, data_train_y, data_test,low_threshold = 0.05)
p_ens_reg[3,] <- predperf(ens_pred[[3]][[1]], data_test_imputed)
p_ens_gbm[3,] <- predperf(ens_pred[[3]][[2]], data_test_imputed)

mbi_pred[[3]] <- mbi_predictions(data_train_x, data_train_y, data_test)
p_mbi_reg[3,] <- predperf(mbi_pred[[3]], data_test_imputed)

#data_test_listwise <- na.omit(data_test)
data_test_listwise <- data_test
listwise_pred[[3]] <- listwise_model(data_train, data_test_listwise)
p_listwise_reg[3,] <- predperf(listwise_pred[[3]][[1]], data_test_listwise)
p_listwise_gbm[3,] <- predperf(listwise_pred[[3]][[2]], data_test_listwise)

mi_pred[[3]] <- mi_model(data_train, data_test)
p_mi_reg[3,] <- predperf(mi_pred[[3]][[1]], data_test_imputed)
p_mi_gbm[3,] <- predperf(mi_pred[[3]][[2]], data_test_imputed)

svi_pred[[3]] <- svi_model(data_train, data_test)
p_svi_reg[3,] <- predperf(svi_pred[[3]][[1]], data_test_imputed)
p_svi_gbm[3,] <- predperf(svi_pred[[3]][[2]], data_test_imputed)

cart_pred[[3]] <- cart_model(data_train, data_test)
p_cart[3,] <- predperf(cart_pred[[3]], data_test_imputed)

#### subset 4

data_train <- tr_test_st(dat_list[[4]])[[1]]
data_test <- tr_test_st(dat_list[[4]])[[2]]
data_train_y <- data.frame(data_train[,"cnt"])
data_train_x <- data_train[,!names(data_train) %in% c("cnt")]
names(data_train_y) = "cnt" 
data_test_imputed <- impute_test(data_train,data_test)
data_train_imputed <- impute_manual(data_train)

brm_pred[[4]] <- brm_predictions(num_blocks,data_train_x, data_train_y, data_test)
p_brm_reg[4,] <- predperf(brm_pred[[4]][[1]], data_test)
p_brm_gbm[4,] <- predperf(brm_pred[[4]][[2]], data_test)

ens_pred[[4]] <- ens_model(num_blocks,data_train_x, data_train_y, data_test,low_threshold = 0.05)
p_ens_reg[4,] <- predperf(ens_pred[[4]][[1]], data_test_imputed)
p_ens_gbm[4,] <- predperf(ens_pred[[4]][[2]], data_test_imputed)

mbi_pred[[4]] <- mbi_predictions(data_train_x, data_train_y, data_test)
p_mbi_reg[4,] <- predperf(mbi_pred[[4]], data_test_imputed)

#data_test_listwise <- na.omit(data_test)
data_test_listwise <- data_test
listwise_pred[[4]] <- listwise_model(data_train, data_test_listwise)
p_listwise_reg[4,] <- predperf(listwise_pred[[4]][[1]], data_test_listwise)
p_listwise_gbm[4,] <- predperf(listwise_pred[[4]][[2]], data_test_listwise)

mi_pred[[4]] <- mi_model(data_train, data_test)
p_mi_reg[4,] <- predperf(mi_pred[[4]][[1]], data_test_imputed)
p_mi_gbm[4,] <- predperf(mi_pred[[4]][[2]], data_test_imputed)

svi_pred[[4]] <- svi_model(data_train, data_test)
p_svi_reg[4,] <- predperf(svi_pred[[4]][[1]], data_test_imputed)
p_svi_gbm[4,] <- predperf(svi_pred[[4]][[2]], data_test_imputed)

cart_pred[[4]] <- cart_model(data_train, data_test)
p_cart[4,] <- predperf(cart_pred[[4]], data_test_imputed)

##
#### subset 5

data_train <- tr_test_st(dat_list[[5]])[[1]]
data_test <- tr_test_st(dat_list[[5]])[[2]]
data_train_y <- data.frame(data_train[,"cnt"])
data_train_x <- data_train[,!names(data_train) %in% c("cnt")]
names(data_train_y) = "cnt" 
data_test_imputed <- impute_test(data_train,data_test)
data_train_imputed <- impute_manual(data_train)

brm_pred[[5]] <- brm_predictions(num_blocks,data_train_x, data_train_y, data_test)
p_brm_reg[5,] <- predperf(brm_pred[[5]][[1]], data_test)
p_brm_gbm[5,] <- predperf(brm_pred[[5]][[2]], data_test)

ens_pred[[5]] <- ens_model(num_blocks,data_train_x, data_train_y, data_test,low_threshold = 0.05)
p_ens_reg[5,] <- predperf(ens_pred[[5]][[1]], data_test_imputed)
p_ens_gbm[5,] <- predperf(ens_pred[[5]][[2]], data_test_imputed)

mbi_pred[[5]] <- mbi_predictions(data_train_x, data_train_y, data_test)
p_mbi_reg[5,] <- predperf(mbi_pred[[5]], data_test_imputed)

#data_test_listwise <- na.omit(data_test)
data_test_listwise <- data_test
listwise_pred[[5]] <- listwise_model(data_train, data_test_listwise)
p_listwise_reg[5,] <- predperf(listwise_pred[[5]][[1]], data_test_listwise)
p_listwise_gbm[5,] <- predperf(listwise_pred[[5]][[2]], data_test_listwise)

mi_pred[[5]] <- mi_model(data_train, data_test)
p_mi_reg[5,] <- predperf(mi_pred[[5]][[1]], data_test_imputed)
p_mi_gbm[5,] <- predperf(mi_pred[[5]][[2]], data_test_imputed)

svi_pred[[5]] <- svi_model(data_train, data_test)
p_svi_reg[5,] <- predperf(svi_pred[[5]][[1]], data_test_imputed)
p_svi_gbm[5,] <- predperf(svi_pred[[5]][[2]], data_test_imputed)

cart_pred[[5]] <- cart_model(data_train, data_test)
p_cart[5,] <- predperf(cart_pred[[5]], data_test_imputed)

### 
#### subset 6

data_train <- tr_test_st(dat_list[[6]])[[1]]
data_test <- tr_test_st(dat_list[[6]])[[2]]
data_train_y <- data.frame(data_train[,"cnt"])
data_train_x <- data_train[,!names(data_train) %in% c("cnt")]
names(data_train_y) = "cnt" 
data_test_imputed <- impute_test(data_train,data_test)
data_train_imputed <- impute_manual(data_train)

brm_pred[[6]] <- brm_predictions(num_blocks,data_train_x, data_train_y, data_test)
p_brm_reg[6,] <- predperf(brm_pred[[6]][[1]], data_test)
p_brm_gbm[6,] <- predperf(brm_pred[[6]][[2]], data_test)

ens_pred[[6]] <- ens_model(num_blocks,data_train_x, data_train_y, data_test,low_threshold = 0.05)
p_ens_reg[6,] <- predperf(ens_pred[[6]][[1]], data_test_imputed)
p_ens_gbm[6,] <- predperf(ens_pred[[6]][[2]], data_test_imputed)

mbi_pred[[6]] <- mbi_predictions(data_train_x, data_train_y, data_test)
p_mbi_reg[6,] <- predperf(mbi_pred[[6]], data_test_imputed)

#data_test_listwise <- na.omit(data_test)
data_test_listwise <- data_test
listwise_pred[[6]] <- listwise_model(data_train, data_test_listwise)
p_listwise_reg[6,] <- predperf(listwise_pred[[6]][[1]], data_test_listwise)
p_listwise_gbm[6,] <- predperf(listwise_pred[[6]][[2]], data_test_listwise)

mi_pred[[6]] <- mi_model(data_train, data_test)
p_mi_reg[6,] <- predperf(mi_pred[[6]][[1]], data_test_imputed)
p_mi_gbm[6,] <- predperf(mi_pred[[6]][[2]], data_test_imputed)

svi_pred[[6]] <- svi_model(data_train, data_test)
p_svi_reg[6,] <- predperf(svi_pred[[6]][[1]], data_test_imputed)
p_svi_gbm[6,] <- predperf(svi_pred[[6]][[2]], data_test_imputed)

cart_pred[[6]] <- cart_model(data_train, data_test)
p_cart[6,] <- predperf(cart_pred[[6]], data_test_imputed)

#### subset 7

data_train <- tr_test_st(dat_list[[7]])[[1]]
data_test <- tr_test_st(dat_list[[7]])[[2]]
data_train_y <- data.frame(data_train[,"cnt"])
data_train_x <- data_train[,!names(data_train) %in% c("cnt")]
names(data_train_y) = "cnt" 
data_test_imputed <- impute_test(data_train,data_test)
data_train_imputed <- impute_manual(data_train)

brm_pred[[7]] <- brm_predictions(num_blocks,data_train_x, data_train_y, data_test)
p_brm_reg[7,] <- predperf(brm_pred[[7]][[1]], data_test)
p_brm_gbm[7,] <- predperf(brm_pred[[7]][[2]], data_test)

ens_pred[[7]] <- ens_model(num_blocks,data_train_x, data_train_y, data_test,low_threshold = 0.05)
p_ens_reg[7,] <- predperf(ens_pred[[7]][[1]], data_test_imputed)
p_ens_gbm[7,] <- predperf(ens_pred[[7]][[2]], data_test_imputed)

mbi_pred[[7]] <- mbi_predictions(data_train_x, data_train_y, data_test)
p_mbi_reg[7,] <- predperf(mbi_pred[[7]], data_test_imputed)

#data_test_listwise <- na.omit(data_test)
data_test_listwise <- data_test
listwise_pred[[7]] <- listwise_model(data_train, data_test_listwise)
p_listwise_reg[7,] <- predperf(listwise_pred[[7]][[1]], data_test_listwise)
p_listwise_gbm[7,] <- predperf(listwise_pred[[7]][[2]], data_test_listwise)

mi_pred[[7]] <- mi_model(data_train, data_test)
p_mi_reg[7,] <- predperf(mi_pred[[7]][[1]], data_test_imputed)
p_mi_gbm[7,] <- predperf(mi_pred[[7]][[2]], data_test_imputed)

svi_pred[[7]] <- svi_model(data_train, data_test)
p_svi_reg[7,] <- predperf(svi_pred[[7]][[1]], data_test_imputed)
p_svi_gbm[7,] <- predperf(svi_pred[[7]][[2]], data_test_imputed)

cart_pred[[7]] <- cart_model(data_train, data_test)
p_cart[7,] <- predperf(cart_pred[[7]], data_test_imputed)

#### subset 8

data_train <- tr_test_st(dat_list[[8]])[[1]]
data_test <- tr_test_st(dat_list[[8]])[[2]]
data_train_y <- data.frame(data_train[,"cnt"])
data_train_x <- data_train[,!names(data_train) %in% c("cnt")]
names(data_train_y) = "cnt" 
data_test_imputed <- impute_test(data_train,data_test)
data_train_imputed <- impute_manual(data_train)

brm_pred[[8]] <- brm_predictions(num_blocks,data_train_x, data_train_y, data_test)
p_brm_reg[8,] <- predperf(brm_pred[[8]][[1]], data_test)
p_brm_gbm[8,] <- predperf(brm_pred[[8]][[2]], data_test)

ens_pred[[8]] <- ens_model(num_blocks,data_train_x, data_train_y, data_test,low_threshold = 0.05)
p_ens_reg[8,] <- predperf(ens_pred[[8]][[1]], data_test_imputed)
p_ens_gbm[8,] <- predperf(ens_pred[[8]][[2]], data_test_imputed)

mbi_pred[[8]] <- mbi_predictions(data_train_x, data_train_y, data_test)
p_mbi_reg[8,] <- predperf(mbi_pred[[8]], data_test_imputed)

#data_test_listwise <- na.omit(data_test)
data_test_listwise <- data_test
listwise_pred[[8]] <- listwise_model(data_train, data_test_listwise)
p_listwise_reg[8,] <- predperf(listwise_pred[[8]][[1]], data_test_listwise)
p_listwise_gbm[8,] <- predperf(listwise_pred[[8]][[2]], data_test_listwise)

mi_pred[[8]] <- mi_model(data_train, data_test)
p_mi_reg[8,] <- predperf(mi_pred[[8]][[1]], data_test_imputed)
p_mi_gbm[8,] <- predperf(mi_pred[[8]][[2]], data_test_imputed)

svi_pred[[8]] <- svi_model(data_train, data_test)
p_svi_reg[8,] <- predperf(svi_pred[[8]][[1]], data_test_imputed)
p_svi_gbm[8,] <- predperf(svi_pred[[8]][[2]], data_test_imputed)

cart_pred[[8]] <- cart_model(data_train, data_test)
p_cart[8,] <- predperf(cart_pred[[8]], data_test_imputed)


### Print outcomes:

pred_all_reg <- predperf(predict(reg_all, data_test_no_miss,type="response"), data_test_no_miss)
pred_all_reg_p <- predperf(predict(reg_p, data_test_no_miss,type="response"), data_test_no_miss)
pred_all_gbm <- predperf(predict(gbm_all, data_test_no_miss,n.trees=gbm_all$n.trees, type="response"), data_test_no_miss)

p_brm_reg 
p_brm_gbm 
p_mbi_reg 
p_ens_reg 
p_ens_gbm 
p_listwise_reg 
p_listwise_gbm 
p_svi_reg 
p_svi_gbm 
p_mi_reg 
p_mi_gbm 
p_cart 

## Arrange these for plots:

rmse_reg_mat <- t(cbind(p_listwise_reg[,1],p_svi_reg[,1],p_mi_reg[,1],p_ens_reg[,1],p_mbi_reg[,1],p_brm_reg[,1]))
mae_reg_mat <- t(cbind(p_listwise_reg[,2],p_svi_reg[,2],p_mi_reg[,2],p_ens_reg[,2],p_mbi_reg[,2],p_brm_reg[,2]))
smape_reg_mat <- t(cbind(p_listwise_reg[,3],p_svi_reg[,3],p_mi_reg[,3],p_ens_reg[,3],p_mbi_reg[,3],p_brm_reg[,3]))

rmse_gbm_mat <- t(cbind(p_listwise_gbm[,1],p_svi_gbm[,1],p_mi_gbm[,1],p_ens_gbm[,1],p_cart[,1],p_brm_gbm[,1]))
mae_gbm_mat <- t(cbind(p_listwise_gbm[,2],p_svi_gbm[,2],p_mi_gbm[,2],p_ens_gbm[,2], p_cart[,2],p_brm_gbm[,2]))
smape_gbm_mat <- t(cbind(p_listwise_gbm[,3],p_svi_gbm[,3],p_mi_gbm[,3],p_ens_gbm[,3], p_cart[,3],p_brm_gbm[,3]))

out_reg = data.frame(brm1 = rmse_reg_mat[6,c(2,6,7,8)], 
                     mi1 = rmse_reg_mat[3,c(2,6,7,8)],
                     ensemble1 = rmse_reg_mat[4,c(2,6,7,8)],
                     svi1 = rmse_reg_mat[2,c(2,6,7,8)],
                     ifr1 = rmse_reg_mat[5,c(2,6,7,8)],
                     list1 = rmse_reg_mat[1,c(2,6,7,8)], 
                     brm2 = mae_reg_mat[6,c(2,6,7,8)],
                     mi2 = mae_reg_mat[3,c(2,6,7,8)],
                     ensemble2 = mae_reg_mat[4,c(2,6,7,8)],
                     svi2 = mae_reg_mat[2,c(2,6,7,8)],
                     ifr2 = mae_reg_mat[5,c(2,6,7,8)],
                     list2 = mae_reg_mat[1,c(2,6,7,8)],
                     brm3 = smape_reg_mat[6,c(2,6,7,8)],
                     mi3 = smape_reg_mat[3,c(2,6,7,8)],
                     ensemble3 = smape_reg_mat[4,c(2,6,7,8)],
                     svi3 = smape_reg_mat[2,c(2,6,7,8)],
                     ifr3 = smape_reg_mat[5,c(2,6,7,8)],
                     list3 = smape_reg_mat[1,c(2,6,7,8)]                     
                     )

out_gbm = data.frame(brm1 = rmse_gbm_mat[6,c(2,6,7,8)], 
                     mi1 = rmse_gbm_mat[3,c(2,6,7,8)],
                     ensemble1 = rmse_gbm_mat[4,c(2,6,7,8)],
                     svi1 = rmse_gbm_mat[2,c(2,6,7,8)],
                     cart1 = rmse_gbm_mat[5,c(2,6,7,8)],
                     list1 = rmse_gbm_mat[1,c(2,6,7,8)], 
                     brm2 = mae_gbm_mat[6,c(2,6,7,8)],
                     mi2 = mae_gbm_mat[3,c(2,6,7,8)],
                     ensemble2 = mae_gbm_mat[4,c(2,6,7,8)],
                     svi2 = mae_gbm_mat[2,c(2,6,7,8)],
                     cart2 = mae_gbm_mat[5,c(2,6,7,8)],
                     list2 = mae_gbm_mat[1,c(2,6,7,8)],
                     brm3 = smape_gbm_mat[6,c(2,6,7,8)],
                     mi3 = smape_gbm_mat[3,c(2,6,7,8)],
                     ensemble3 = smape_gbm_mat[4,c(2,6,7,8)],
                     svi3 = smape_gbm_mat[2,c(2,6,7,8)],
                     cart3 = smape_gbm_mat[5,c(2,6,7,8)],
                     list3 = smape_gbm_mat[1,c(2,6,7,8)]
)

write.csv(cbind(out_reg,out_gbm), "Prediction_performance_comparison.csv",row.names=F)
