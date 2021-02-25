### Prepare the dataset - Input and outcomes
### Dec has wave 5 Nov 24 till Dec 19, and Wave 6 onwards
## Outcome is y_test = Covid_positive tested in last 14 days


library(ROCR)

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

library(dplyr)
library(iml)
library(DALEX)
library(partykit)

library(tidyverse)
library(qdapTools)

library(janitor)
library(pryr)


Dec_d0 = read.csv("symtrack/data/2020-12.csv")
dim(Dec_d0)
names(Dec_d)

# Dec_d = Dec_d0 %>% filter( B10a == 1 | B10a == 2 )

Dec_dp = Dec_d0 %>% filter( B10a == 1 )

Dec_dn = Dec_d0 %>% filter( B10a == 2 )

set.seed(12345)
ind = sample(nrow(Dec_dn),nrow(Dec_dp), replace = FALSE) 

Dec_d = rbind(Dec_dp, Dec_dn[ind,])  

## Randomize the positive and negative cases
set.seed(12345)
random_ind <- sample(nrow(Dec_d))
Dec_d = Dec_d[random_ind,]

# Dec_d = subset(Dec_d0, B10a == 1 | B10a == 2)

## Symptoms = Dec_d[,c("A1_1" , "A1_2",  "A1_3"  ,    "A1_4",          "A1_5")]

## ppl_house_sick = Dec_d[,c("A2")]
## People in house
ppl_house = Dec_d[,c("A5_1" , "A5_2",  "A5_3")]

## Only integers, positive numbers and <=20 allowed. 
ppl_house$A5_1[ppl_house$A5_1 > 20] = NA
ppl_house$A5_2[ppl_house$A5_2 > 20] = NA
ppl_house$A5_3[ppl_house$A5_3 > 20] = NA

ppl_house$A5_1[ppl_house$A5_1 < 0] = NA
ppl_house$A5_2[ppl_house$A5_2 < 0] = NA
ppl_house$A5_3[ppl_house$A5_3 < 0] = NA

ppl_house$A5_1[round(ppl_house$A5_1) != ppl_house$A5_1 ] = NA
ppl_house$A5_2[round(ppl_house$A5_2) != ppl_house$A5_2] = NA
ppl_house$A5_3[round(ppl_house$A5_3) != ppl_house$A5_3 ] = NA

table(ppl_house$A5_3)

ppl_community_sick = Dec_d[,c("A4")]
ppl_community_sick[round(ppl_community_sick) != ppl_community_sick ] = NA
ppl_community_sick[ppl_community_sick < 0] = NA

ppl_community_sick[ppl_community_sick >20 ] = 20

### Need to clean this variable!
pre_conditions = as.data.frame(as.matrix(Dec_d[,c("C1")]), stringsAsFactors = FALSE)

#options(stringsAsFactors=FALSE)
#temp = pre_conditions[1:10,1]
#mtabulate(strsplit(temp, ","))

pre_cond_mat = mtabulate(strsplit(pre_conditions[,1], ","))
index_na_cond <- rowSums(pre_cond_mat)
pre_cond_mat[index_na_cond == 0,] = NA


## Gender
demo_gender = Dec_d[,c("D1")]
demo_gender[demo_gender > 2] = NA
demo_gender = factor(demo_gender)

table(demo_gender)

### Age:
demo_age = Dec_d[,c("D2")]
demo_age = factor(demo_age)

### Origin and Race missing: 

demo_edu = Dec_d[,c("D8")]
demo_edu = factor(demo_edu)

demo_workOutside = Dec_d[,c("D10")]

work_cat = Dec_d[,c("Q64")]
work_cat = factor(work_cat)


work_other = Dec_d[,c("Q80")]
work_other = factor(work_other)

social_distancing = Dec_d[,c("C7")]

wear_mask_freq = Dec_d[,c("C14")]
wear_mask_freq[wear_mask_freq == 6] <- NA

others_wear_mask_freq = Dec_d[,c("C16")]
others_wear_mask_freq[others_wear_mask_freq == 6] = NA

### Children in household (can be an important factor):
summary(as.factor(Dec_d$E1_1))

children_kinder = Dec_d$E1_1
children_kinder[children_kinder == 5] = NA

children_grade1_5 = Dec_d$E1_2
children_grade1_5[children_grade1_5 == 5] = NA

children_grade6_8 = Dec_d$E1_3
children_grade6_8[children_grade6_8 == 5] = NA

children_grade9_12 = Dec_d$E1_4
children_grade9_12[children_grade9_12 == 5] = NA

children_fulltime = Dec_d$E2_1
children_fulltime[children_fulltime == 4] = NA

children_parttime = Dec_d$E2_2
children_parttime[children_parttime == 4] = NA

how_many_people = Dec_d[,c("C10_1_1","C10_2_1","C10_3_1","C10_4_1")]

how_many_people$C10_1_1[how_many_people$C10_1_1 > 100] = NA
how_many_people$C10_2_1[how_many_people$C10_2_1 > 100] = NA
how_many_people$C10_3_1[how_many_people$C10_3_1 > 100] = NA
how_many_people$C10_4_1[how_many_people$C10_4_1 > 100] = NA

how_many_people$C10_1_1[how_many_people$C10_1_1 < 0] = NA
how_many_people$C10_2_1[how_many_people$C10_2_1 < 0] = NA
how_many_people$C10_3_1[how_many_people$C10_3_1 < 0] = NA
how_many_people$C10_4_1[how_many_people$C10_4_1 < 0] = NA

how_many_people$C10_1_1[round(how_many_people$C10_1_1) != how_many_people$C10_1_1   ] = NA
how_many_people$C10_2_1[round(how_many_people$C10_2_1) != how_many_people$C10_2_1   ] = NA
how_many_people$C10_3_1[round(how_many_people$C10_3_1) != how_many_people$C10_3_1   ] = NA
how_many_people$C10_4_1[round(how_many_people$C10_4_1) != how_many_people$C10_4_1   ] = NA

flu_vaccine = Dec_d[,c("C17")]
flu_vaccine[flu_vaccine == 2] = NA

test = factor(ifelse(Dec_d$B10a == 1, 1, 0)) 

# all_data = data.frame(test = test, symptoms = Symptoms,ppl_house_sick= ppl_house_sick,ppl_house = ppl_house, ppl_community_sick = ppl_community_sick,
#                       pre_cond_mat = pre_cond_mat,demo_gender = demo_gender, demo_age = demo_age, demo_edu = demo_edu,
#                       demo_workOutside = demo_workOutside, work_cat = work_cat,work_other =  work_other, 
#                       social_distancing = social_distancing, wear_mask_freq = wear_mask_freq, how_many_people = how_many_people, 
#                       flu_vaccine = flu_vaccine)

all_data = data.frame(test = test, ppl_house = ppl_house, 
                      pre_cond_mat = pre_cond_mat,demo_gender = demo_gender, demo_age = demo_age, demo_edu = demo_edu,
                      demo_workOutside = demo_workOutside, work_cat = work_cat,work_other =  work_other, 
                      social_distancing = social_distancing, wear_mask_freq = wear_mask_freq,
                      others_wear_mask_freq = others_wear_mask_freq, children_kinder = children_kinder,
                      children_grade1_5 = children_grade1_5,
                      children_grade6_8 = children_grade6_8, children_grade9_12 = children_grade9_12,
                      children_fulltime = children_fulltime, children_parttime = children_parttime,
                      how_many_people = how_many_people, 
                      flu_vaccine = flu_vaccine)

sapply(all_data, function(y) sum(length(which(is.na(y)))))

data.frame(sapply(all_data, function(y) round(sum(length(which(is.na(y))))/length(y)*100,2) ))

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

brm_num_blocks <- function(data_train_x,low_threshold = 0.10)
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
    # data_subs_2 <- list(length(data_subs_1)) ## This is for identifying the columns
    data_subs_2 <- vector("list", length = length(data_subs_1))
    missingness <- vector(length = length(data_subs_1))
    block_info <- vector(length = length(data_subs_1))
    for(list in 1:length(data_subs_1))
    {
      completeness <- colSums(data_subs_1[[list]])
      
      ## data_subs_2[[list]] <- data.frame(data_subs_1[[list]][,which(completeness>(1-low_threshold)*nrow(data_subs_1[[list]])) ] )
      data_subs_2[[list]] <- data.frame(data_subs_1[[list]][,which(completeness>(low_threshold*nrow(data_subs_1[[list]]))) ] )
      
      
      missing_vals <- nrow(data_subs_2[[list]]) - colSums(data_subs_2[[list]])
      missingness[list] <- sum(missing_vals)
      block_info[list] <- sum(colSums(data_subs_2[[list]]))
    }
    ## return(list(  sum(block_info), (sum(missingness)/(dim(MVI)[1]*dim(MVI)[2]))   ) )
    return(sum(missingness)/(dim(MVI)[1]*dim(MVI)[2]))
  }
  
  M <- (ncol(data_train_x))  ## Let the maximum number of clusters be equal to the number of columns - 1. 
  
  # total_info = sum(colSums(data_present))
  # missing_prop_list <- vector("list", length = M) 
  ## missing_prop_score <- vector(length = M) 
  missing_prop_val <- vector(length = M) 

  for(k in 1:M)
  {
    # missing_prop_list[[k]] <- missing_prop(data_present,k, low_threshold)
    ### The score is the sum of info missed out + total missing proportion across blocks
    #missing_prop_score[k] = (total_info - missing_prop_list[[k]][[1]])/total_info +  missing_prop_list[[k]][[2]]
    missing_prop_val[k] = missing_prop(data_present,k, low_threshold)
  }
  
  return(missing_prop_val)
}

brm_ov_subsets <- function(num_blocks, X_in, Y_in,low_threshold = 0.10)
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
    data_subs_out <- vector("list", length = length(data_subs_1))
    data_values <- split(data_xy, members)
    for(list in 1:length(data_subs_1))
    {
      completeness <- colSums(data_subs_1[[list]])
      # inputs_in_subset <- which(completeness>(1-low_threshold)*nrow(data_subs_1[[list]]))
      inputs_in_subset <- which(completeness>(low_threshold*nrow(data_subs_1[[list]])))
      outcome_index <- which(names(data_values[[list]]) ==  names(Y_in) )
      
      #data_subs_out[[list]] <- data_values[[list]][, c(inputs_in_subset,outcome_index) ]
      data_subs_out[[list]] <- data.frame(data_values[[list]][, c(inputs_in_subset,outcome_index) ])
      
      ## 10 low_threshold reduces missing values, but also truncates information. 
      ## Since info more important than reducing missingness for each block, we capture max info per block
      ## A more complex formulation of optimization of information and missingness 
      ###   for different values of alpha can be done, but that is out of scope for this paper.   
      
      columns[[list]] <- names(data_subs_out[[list]])
    }
    return(list(mcc_rounded,columns,data_subs_out ))
  }
  
  n_ov_subsetting <- columns_in_subsets(data_xy,data_present,num_blocks, low_threshold)
  cluster_centers_subsets <- n_ov_subsetting[[1]]
  columns_n_ov_subsets <- n_ov_subsetting[[2]]
  n_ov_subsets <- n_ov_subsetting[[3]]
  
  ## 
  ### Logic for Set theory determination of overlapping subsets
  
  ov_subsets <- vector("list", length(n_ov_subsets))
  
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

### UPDATE: Adding stepreg to this instead of plain logistic to avoid an error of 
### computational signularity while evaluating MI

logitreg <- function(data_in)
{
  X <- data_in[,names(data_in) %in% names(data_train_x) ]
  Y <- data_in[,names(data_in) %in% names(data_train_y) ]
  
  X <- remove_constant(X, na.rm= FALSE, quiet = TRUE)
  
  datam <- data.table(X,Y)
  fmla <- as.formula(paste("Y ~ ", paste(names(X), collapse= "+")))
  reg <- glm( fmla, data=data_in,family="binomial")
  # reg <- stepAIC(glm( fmla, data=data_in,family="binomial"), direction = "both", trace = FALSE)
  return(reg)
}

logitregstepwise <- function(data_in)
{
  X <- data_in[,names(data_in) %in% names(data_train_x) ]
  Y <- data_in[,names(data_in) %in% names(data_train_y) ]
  
  X <- remove_constant(X, na.rm= FALSE, quiet = TRUE)
  
  datam <- data.table(X,Y)
  fmla <- as.formula(paste("Y ~ ", paste(names(X), collapse= "+")))
  # reg <- glm( fmla, data=data_in,family="binomial")
  reg <- stepAIC(glm( fmla, data=data_in,family="binomial"), direction = "both", trace = FALSE)
  return(reg)
}


gbm_train <- function(data_in)
{
  X <- data_in[,names(data_in) %in% names(data_train_x) ]
  Y <- data_in[,names(data_in) %in% names(data_train_y) ]
  Y = as.logical(as.integer(Y)-1)
  X <- remove_constant(X, na.rm= FALSE, quiet = TRUE)
  datam <- data.table(X,Y) 
  fmla <- as.formula(paste("Y ~ ", paste(names(X), collapse= "+")))
  reg <- gbm( fmla, data=datam ,verbose=FALSE,n.trees = 500,shrinkage=0.05,interaction.depth = 3)
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
    if(class(models[[i]]) == "lm")
    {
      reg_outcomes[[i]] <- predict(models[[i]],newdata = datam)
    }
    else if(class(models[[i]]) == "glm"){
      reg_outcomes[[i]] <- predict(models[[i]],newdata = datam,type = "response")
    }
    else if(class(models[[i]]) == "gbm"){
      reg_outcomes[[i]] <- predict(models[[i]],newdata = datam,n.trees=models[[i]]$n.trees)
    }
    # reg_outcomes[[i]] <- predict(reg_m[[i]],newdata = datam,type="response")
  }
  
  reg_out_a <- unlist(reg_outcomes)
  reg_ord <- reg_out_a[order(get_back_row_order[,"rownum"])]
  return(reg_ord)
}

model_fit_comparisons_reg <- function(models)
{
  return(lapply(models, function(x) summary(x)$r.squared))
}

# https://thestatsgeek.com/2014/02/08/r-squared-in-logistic-regression/
model_fit_comparisons_logitreg <- function(models)
{
  return(lapply(models, function(x) 
    {
    y_ = x$y
    nullmod <- glm(y_ ~ 1, family="binomial")
    denom <- logLik(nullmod)
    1-logLik(x)/denom
    }
    ))
}


## Get the global coefficient scores 
get_global_wts_logit <- function(models)
{
  
  global_reg_par_wt <- function(model)
  {
    nrow <- length(model$residuals)
    #se <- summary(model)$coefficients[,2]
    #nrow/se
    coef <- summary(model)$coefficients[,1]
    coef_ind = coef
    coef_ind <- coef_ind/coef_ind
    e_reg_ <- e_logit(model)  
    return(coef_ind*nrow/e_reg_)    
    # return(coef_ind/e_reg_) 
  }
  
  
  global_reg_par <- function(model)
  {
    nrow <- length(model$residuals)
    coef <- summary(model)$coefficients[,1]
    # se <- summary(model)$coefficients[,2]
    # coef*nrow/se
    e_reg_ <- e_logit(model) 
    return(coef*nrow/e_reg_)
    # return(coef/e_reg_)
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
    mse_gbm_ <- mse_gbm(model)
    nrow <- length(model$fit)
    varimp <- summary(model)[2]
    varimp_ind <- varimp
    varimp_ind <- varimp_ind/varimp_ind
    # return(sqrt(length(model$var.names)*length(model$fit))/mse_gbm(model))
    return(varimp_ind*nrow/mse_gbm_)
    # return(varimp_ind/mse_gbm_)
  }
  
  
  global_gbm_par <- function(model){
    # varimp <- summary(model)*sqrt(length(model$var.names)*length(model$fit))/(mse_gbm(model)*wt2)
    varimp <- summary(model)[2]
    nrow <- length(model$fit)
    mse_gbm_ <- mse_gbm(model)
    return(varimp*nrow/mse_gbm_)
    # return(varimp/mse_gbm_)
  }  
  
  global_wt_num <- lapply(models,global_gbm_par)
  global_wt_denom <- lapply(models,global_gbm_par_wt)
  
  global_wt_denom2 <- lapply(global_wt_denom,function(x) data.frame(t(x)))
  global_wt_denom_mat <- do.call(smartbind,global_wt_denom2)
  global_wt_denom_mat_avg <- apply(global_wt_denom_mat,2,function(x) sum(x,na.rm=T))
  
  global_wt_num2 <- lapply(global_wt_num,function(x) data.frame(t(x)))
  global_wt_num_mat <- do.call(smartbind,global_wt_num2)
  global_wt_num_mat_avg <- apply(global_wt_num_mat,2,function(x) sum(x,na.rm=T))
  
  global_wts <- global_wt_num_mat_avg/global_wt_denom_mat_avg
  
  return(global_wts)  
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

#mse_reg <- function(model) { mean((model$residuals)^2) } 
mse_gbm <- function(model) {model$train.error[length(model$train.error)]} 

e_logit <- function(model) {summary(model)$deviance }

## BRM execution codes in one place as a function:
brm_predictions <- function(num_blocks,X_in, Y_in, data_test)
{
  brm_model <- brm_ov_subsets(num_blocks,X_in, Y_in, low_threshold = 0.10)
  ov_subsets <- brm_model[[1]]
  cluster_centers_subsets <- brm_model[[2]]
  n_ov_subsets <- brm_model[[3]]
  subsets_columns <- brm_model[[4]]
  reg_m <- lapply(ov_subsets,logitreg)
  gbm_m <- lapply(ov_subsets,gbm_train)
  pred_reg_outcomes_sc <- pred_reg(reg_m, data_test,data_train_x, data_train_y,cluster_centers_subsets)
  pred_gbm_outcomes_sc <- pred_reg(gbm_m, data_test,data_train_x, data_train_y,cluster_centers_subsets)
  
  pred_reg_outcomes <- factor(ifelse(pred_reg_outcomes_sc >= 0.5, 1, 0))
  pred_reg_outcomes[is.na(pred_reg_outcomes)] = 0
  
  pred_gbm_outcomes = factor(ifelse(pred_gbm_outcomes_sc > 0, 1, 0)) 
  pred_gbm_outcomes[is.na(pred_gbm_outcomes)] = 0
  
  return(list(pred_reg_outcomes,pred_gbm_outcomes))
}




acc_out <- function(y_actual, y_predict)
{ 
  ### Use own formula, works best!!
  
  # y_actual = data_test_imputed[,"test"]
  # y_predict = listwise_pred_reg
  # y_predict = pred_reg_outcomes
  ## ref = levels(y_actual)[2]
  acc_ = sum(as.integer(y_actual == y_predict))/length(y_actual) 
    ## accuracy(y_actual,y_predict)
  auc_ = auc(y_actual,y_predict)
  recall_ = sum(as.integer(y_actual == 1 & y_predict == 1))/sum(as.integer(y_actual ==1)) 
    ## recall(data = y_predict, reference = y_actual, relevant = ref)  ## From caret
  precision_ =  sum(as.integer(y_actual == 1 & y_predict == 1))/sum(as.integer(y_predict ==1)) 
    ## precision(data = y_predict, reference = y_actual, relevant = ref)
  f1_ =  2*(recall_ * precision_)/(recall_ + precision_ )
    ## F_meas(data = y_predict, reference = y_actual, relevant = ref)
  pred_ <- c(round(acc_,4), round(auc_,4) , round(recall_,4) , round(precision_,4) ,round(f1_,4) )
  return(pred_)
}


##### main{} #####

### Since each input has some missing values, we have a dummy variable that is constant 
## but has no missing values, for data processing convenience
# all_data$dummy = 1

data_train <- tr_test_st(all_data)[[1]]
data_test <- tr_test_st(all_data)[[2]]

data_test_imputed <- impute_test(data_train,data_test)
data_train_imputed <- impute_manual(data_train)

round(sapply(data_train, function(y) sum(length(which(is.na(y)))))/dim(data_train)[1],2)

# sapply(data_test, function(y) sum(length(which(is.na(y)))))

data_train_y <- data.frame(data_train[,"test"])
data_train_x <- data_train[,!names(data_train) %in% c("test")]
names(data_train_y) = "test" 

system.time( missing_prop_val <- brm_num_blocks(data_train_x,low_threshold = 0.10)  )
### Takes around 10 minutes for 150K data, so take sample for now

### We assume four blocks identified here!
plot(missing_prop_val)
#imputation_thres <- 0.05
imputation_thres <- 0.10

num_blocks <- min(which(missing_prop_val <= imputation_thres))
#num_blocks <- 10

## When to use BRM and when not to? With a quantitative measure for the dataset. 
## Number of blocks identified automatically using the plot (derivative or something) 
## We decide how many blocks to consider based on plot

### low_threshold is the proportion of data up to which you are ready to impute. If it is lower, 
### only a few columns which have full data within block will be included. Since clustering is stochastic, 
### We set values between 0.05 to 0.20. For higher values such as 0.95, almost all columns will get selected 
### if clustering quality is bad. For MAR, we will observe the blocks are missing, and we need to use higher alpha. 

brm_model <- brm_ov_subsets(num_blocks,data_train_x, data_train_y,low_threshold = 0.10) 

ov_subsets <- brm_model[[1]]
cluster_centers_subsets <- brm_model[[2]]
n_ov_subsets <- brm_model[[3]]
subsets_columns <- brm_model[[4]]

lapply(ov_subsets, nrow)

## Binary classification problem

##system.time( reg_m <- lapply(ov_subsets,stepreg)  ) 
## 3 mins only logit, 15 minutes for stepwise logit
system.time( logit_m <- lapply(ov_subsets,logitreg) )
summary(logit_m[[3]])
summary(logit_m[[10]])

gbm_m <- lapply(ov_subsets,gbm_train)
### Converted the outcomes to logical using as.logical(as.integer(y)-1) in gbm_train
summary(gbm_m[[3]])
summary(gbm_m[[7]])

# rsquares_list <- unlist(model_fit_comparisons_reg(reg_m))
rsquares_list <- unlist(model_fit_comparisons_logitreg(logit_m))
round(rsquares_list,2)

pred_reg_outcomes_sc <- pred_reg(logit_m, data_test,data_train_x, data_train_y,cluster_centers_subsets)
pred_reg_insample_sc <- pred_reg(logit_m, data_train,data_train_x, data_train_y,cluster_centers_subsets)

pred_reg_outcomes = pred_reg_outcomes_sc
pred_reg_insample = pred_reg_insample_sc

pred_reg_outcomes[is.na(pred_reg_outcomes)] = 0
pred_reg_insample[is.na(pred_reg_insample)] = 0

pred_reg_outcomes <- factor(ifelse(pred_reg_outcomes >= 0.5, 1, 0))
pred_reg_insample <- factor(ifelse(pred_reg_insample >= 0.5, 1, 0))

global_wts_logit <- get_global_wts_logit(logit_m)

pred_gbm_outcomes_sc <- pred_reg(gbm_m, data_test,data_train_x, data_train_y,cluster_centers_subsets)
pred_gbm_insample_sc <- pred_reg(gbm_m, data_train,data_train_x, data_train_y,cluster_centers_subsets)

pred_gbm_outcomes = factor(ifelse(pred_gbm_outcomes_sc > 0, 1, 0)) 
pred_gbm_insample = factor(ifelse(pred_gbm_insample_sc > 0, 1, 0))

pred_gbm_outcomes[is.na(pred_gbm_outcomes)] = 0
pred_gbm_insample[is.na(pred_gbm_insample)] = 0

#pred_gbm_outcomes_sc <- pred_reg(gbm_m, data_test,data_train_x, data_train_y,cluster_centers_subsets)
#pred_gbm_insample_sc <- pred_reg(gbm_m, data_train,data_train_x, data_train_y,cluster_centers_subsets)

# ## Convert [-1,1] to [0,1]
# gbm_normal_2_risk <- function(vec_in)
# {
#   vec_out = (vec_in - min(vec_in)) / ( max(vec_in) - min(vec_in))
#   vec_out[is.na(vec_out)] = 0
#   return( vec_out )
# }

# pred_gbm_outcomes = gbm_normal_2_risk(pred_gbm_outcomes_sc)
# pred_gbm_insample = gbm_normal_2_risk(pred_gbm_insample_sc)

## accuracy, auc, recall, precision, f1 

brm_acc_reg = acc_out(data_test_imputed[,"test"],pred_reg_outcomes)
brm_acc_gbm = acc_out(data_test_imputed[,"test"],pred_gbm_outcomes)

brm_acc_reg_in = acc_out(data_train_imputed[,"test"],pred_reg_insample)
brm_acc_gbm_in = acc_out(data_train_imputed[,"test"],pred_gbm_insample)

## List-wise deletion

sapply(data_train, function(y) sum(length(which(is.na(y)))))

### Remove columns with >50% missing
#data_train_columnwise = data_train[,c(1:25,29,30,35)]

# data_train_listwise <- na.omit(data_train_columnwise)
data_train_listwise <- na.omit(data_train)

sapply(data_train_listwise, function(y) sum(length(which(is.na(y)))))

system.time( reg_listwise <- logitreg(data_train_listwise)  )

# system.time( reg_listwise <- logitregstepwise(data_train_listwise)  )
### 1 min for 2000 observations

gbm_listwise = gbm_train(data_train_listwise)

listwise_pred_reg_sc <- predict(reg_listwise, data_test_imputed)
listwise_pred_gbm_sc <- predict(gbm_listwise, data_test_imputed,n.trees=gbm_listwise$n.trees)

# listwise_pred_reg = listwise_pred_reg_sc
# listwise_pred_reg[is.na(listwise_pred_reg)] = 0
# listwise_pred_reg[listwise_pred_reg < 0 ] = 0
# listwise_pred_reg[listwise_pred_reg > 0 ] = 1
# listwise_pred_gbm = gbm_normal_2_risk(listwise_pred_gbm_sc)

listwise_pred_reg <- factor(ifelse(listwise_pred_reg_sc >= 0.5, 1, 0))
listwise_pred_reg[is.na(listwise_pred_reg)] = 0
listwise_pred_gbm = factor(ifelse(listwise_pred_gbm_sc > 0, 1, 0))
listwise_pred_gbm[is.na(listwise_pred_gbm)] = 0

listwise_acc_reg = acc_out(data_test_imputed[,"test"], listwise_pred_reg)
listwise_acc_gbm = acc_out(data_test_imputed[,"test"], listwise_pred_gbm)




### MI:

mi_model <- function(data_train_mi,data_test_mi)
{
  # data_test_imputed_mi <-  impute_test(data_train_mi,data_test_mi)  
  data_mi <- mice(data_train_mi,m=5,maxit=3,meth='pmm',seed=1234) 
  completedData <- complete(data_mi,1)
  model_mi_gbm <- list()
  model_mi_reg <- list()
  #test_gbm = as.logical(as.integer(data_train$test)-1)
  
  for(list in 1:data_mi$m)
  {
    # model_mi_gbm[[list]] <- gbm( test_gbm ~., data=complete(data_mi,list) ,verbose=FALSE,n.trees = 500,shrinkage=0.1,interaction.depth = 3)
    try( model_mi_gbm[[list]] <- gbm( as.logical(as.integer(test)-1) ~., data=complete(data_mi,list) ,verbose=FALSE,n.trees = 50,shrinkage=0.1,interaction.depth = 3) ) 
    
    try( model_mi_reg[[list]] <- glm( test ~., data=complete(data_mi,list), family="binomial"  ) )
  }
  model_meta_var_gbm <- data.frame(matrix(NA, nrow = nrow(completedData), ncol = data_mi$m ))
  model_meta_var_reg <- data.frame(matrix(NA, nrow = nrow(completedData), ncol = data_mi$m ))
  for(i in 1:data_mi$m)
  {
    try( model_meta_var_gbm_sc <- predict(model_mi_gbm[[i]],n.trees = model_mi_gbm[[i]]$n.trees) ) 
    try(model_meta_var_reg_sc <- predict(model_mi_reg[[i]], type = "response" )) 
    
    model_meta_var_gbm[,i] <- factor(ifelse(model_meta_var_gbm_sc >= 0.5, 1, 0))
    model_meta_var_gbm[,i][is.na(model_meta_var_gbm[,i])] = 0
    
    model_meta_var_reg[,i] = factor(ifelse(model_meta_var_reg_sc > 0, 1, 0)) 
    model_meta_var_reg[,i][is.na(model_meta_var_reg[,i])] = 0
  }
  
  # Y <- data_train$test
  # Y_gbm = test_gbm
  Y = data_train_mi$test
  
  Xgbm <- model_meta_var_gbm
  datam1 <- data.table(Xgbm,Y)
  
  rf.model_meta_gbm <- randomForest(Y ~.,data=datam1,importance =FALSE,ntree=50)
  
  Xreg <- model_meta_var_reg
  datam2 <- data.table(Xreg,Y)
  rf.model_meta_reg <- randomForest(Y~.,data=datam2,importance =FALSE,ntree=50)
  
  ## Prediction is same procedure
  meta_var_pred_gbm <- data.frame(matrix(NA, nrow = nrow(data_test_mi), ncol = data_mi$m)) 
  meta_var_pred_reg <- data.frame(matrix(NA, nrow = nrow(data_test_mi), ncol = data_mi$m)) 
  for(i in 1:data_mi$m)
  {
    try( meta_var_pred_gbm[,i] <- predict(model_mi_gbm[[i]],data_test_mi,n.trees = model_mi_gbm[[i]]$n.trees)  )
    meta_var_pred_reg[,i] <- predict(model_mi_reg[[i]],data_test_mi)
  }
  
  rf_meta_pred_gbm <- predict(rf.model_meta_gbm,meta_var_pred_gbm)
  rf_meta_pred_reg <- predict(rf.model_meta_reg,meta_var_pred_reg)
  return(list(rf_meta_pred_reg,rf_meta_pred_gbm))
}

### SVI:
svi_model <- function(data_train_svi,data_test_svi)
{
  # data_test_imputed <-  impute_test(data_train,data_test)  
  svi.data <- rfImpute(data_train_svi[,!names(data_train_svi)%in% c('test')],data_train_svi[,names(data_train_svi)%in% c('test')],ntree=100,iter=1)  
  names(svi.data)[1] <- "test"
  # svi.data$test_gbm = as.logical(as.integer(svi.data$test)-1)
  # svi.model_gbm <- gbm( test~., data=svi.data ,verbose=FALSE,n.trees = 500,shrinkage=0.1,interaction.depth = 1) 
  svi.model_gbm <- gbm( as.logical(as.integer(test)-1) ~., data=svi.data ,verbose=FALSE,n.trees = 100,shrinkage=0.1,interaction.depth = 1) 
  svi.pred_gbm_sc <- predict(svi.model_gbm, newdata=data_test_svi,n.trees = svi.model_gbm$n.trees)
  svi.model_reg <- glm( as.logical(as.integer(test)-1) ~., data=svi.data, family = "binomial") 
  svi.pred_reg_sc <- predict(svi.model_reg, newdata=data_test_svi)
  
  svi.pred_gbm <- factor(ifelse(svi.pred_gbm_sc >= 0.5, 1, 0))
  svi.pred_gbm[is.na(svi.pred_gbm)] = 0
  
  svi.pred_reg = factor(ifelse(svi.pred_reg_sc > 0, 1, 0)) 
  svi.pred_reg[is.na(svi.pred_reg)] = 0
  
  return(list(svi.pred_reg, svi.pred_gbm))
}

##Evaluation:
### CART:

cart_model <- function(data_train, data_test_cart)
{
  cart_m <- rpart(test ~ ., method="class", data=data_train )
  cart_prediction <- predict(cart_m, data_test_cart)
  return(cart_prediction)
}

#system.time( pred_mi <- mi_model(data_train[1:5000,],data_test_imputed) )
system.time( pred_svi <- svi_model(data_train[1:10000,],data_test_imputed) )  ## Data was too big for SVI (~168GB)
system.time( pred_cart <- cart_model(data_train,data_test_imputed) ) 

pred_cart_b = factor(ifelse(pred_cart[,2] > 0.5, 1, 0)) 

#mi_reg <- pred_mi[[1]]
#mi_gbm <- pred_mi[[2]]

svi_reg <- pred_svi[[1]]
svi_gbm <- pred_svi[[2]]

cart_acc = acc_out(data_test_imputed[,"test"], pred_cart_b)

svi_acc_reg = acc_out(data_test_imputed[,"test"], svi_reg)
svi_acc_gbm = acc_out(data_test_imputed[,"test"], svi_gbm)

### If predictions are only 1s:

acc_out(data_test_imputed[,"test"], factor(c(0,rep(1,length(data_test_imputed[,"test"])-1))))

## IF only zeros:
acc_out(data_test_imputed[,"test"], factor(c(1,rep(0,length(data_test_imputed[,"test"])-1))))


## Ensemble:
ens_model <- function(num_blocks, X_in, Y_in,test_data, low_threshold = 0.10)
{
  X_in = data_train_x
  Y_in = data_train_y
  test_data = data_test

    
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
  
  n_ov_subsetting <- columns_in_subsets(data_xy,data_present,num_blocks, low_threshold = 0.10)
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
    reg <- glm( fmla, data=data_in, family= "binomial")
    return(reg)
  }
  
  gbm_train_ens <- function(data_in)
  {
    X <- data_in[,names(data_in) %in% names(X_in) ]
    Y <- data_in[,names(data_in) %in% names(Y_in) ]
    Y_logit <- as.logical(as.integer(Y)-1)
    datam <- data.table(X,Y_logit) 
    fmla <- as.formula(paste("Y_logit ~ ", paste(names(X), collapse= "+")))
    reg <- gbm( fmla, data=datam ,verbose=FALSE,n.trees = 500,shrinkage=0.1,interaction.depth = 3)
    return(reg)
  }
  
  reg_ens <- lapply(n_ov_subsets,stepreg_ens)
  gbm_ens <- lapply(n_ov_subsets,gbm_train_ens)
  
  mse_reg <- function(model) {summary(model)$deviance }  ## Being lazy here
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
      
      if(class(model_list[[k]]) == "glm")
      {
        try(reg_ensemble_pred_1 <- reg_ensemble_pred_1 + (predict(m,newdata=test_data)*((length(m$residuals)*length(m$coefficients))^0.5/mse_reg(m)) ) ) 
        reg_ensemble_pred_2 <- reg_ensemble_pred_2 + ((length(m$residuals)*length(m$coefficients))^0.5/mse_reg(m))
        mse_list[[k]] <- mse_reg(m)
        predict_model_list[[k]] <-  predict(m,newdata=test_data)
      }else if(class(model_list[[k]]) == "gbm"){
        try(reg_ensemble_pred_1 <- reg_ensemble_pred_1 + (predict(m,test_data,n.trees=m$n.trees)*((length(m$var.names)*length(m$fit))^0.5/(mse_reg(m))) ) )
        reg_ensemble_pred_2 <- reg_ensemble_pred_2 + ((length(m$var.names)*length(m$fit))^0.5/(mse_gbm(m)))
        mse_list[[k]] <- mse_gbm(m)
        predict_model_list[[k]] <-  predict(m,test_data,n.trees=m$n.trees)
      }
      
    }                     
    
    reg_ensemble_best <- predict_model_list[[which.min(mse_list)]]
    # reg_ensemble_pred <- reg_ensemble_pred_1/reg_ensemble_pred_2 
    return(reg_ensemble_best)
  }
  
  reg_ens_out <- ensemble_reg(reg_ens,data_test_imputed)
  gbm_ens_out <- ensemble_reg(gbm_ens,data_test_imputed)
  return(list(reg_ens_out, gbm_ens_out))
}

# pred_ens <- ens_model(num_blocks, data_train_x, data_train_y,data_test,low_threshold = 0.10)
#ens_reg <- pred_ens[[1]]
#ens_gbm <- pred_ens[[2]]

ens_reg <- reg_ens_out
ens_gbm <- gbm_ens_out

ens_reg <- factor(ifelse(ens_reg >= 0.5, 1, 0))
ens_reg[is.na(ens_reg)] = 0
ens_gbm = factor(ifelse(ens_gbm > 0, 1, 0))
ens_gbm[is.na(ens_gbm)] = 0

ense_acc_reg = acc_out(data_test_imputed[,"test"], ens_reg)
ens_acc_gbm = acc_out(data_test_imputed[,"test"], ens_gbm)


### Do this later for binary classification!
## Two ways of getting coefficients: 
feature_imp_reg <- feature_drop_loss(ov_subsets,data_train_x, data_train_y, gbm_train, gbm_m, cluster_centers_subsets)
global_wts_logit <- get_global_wts_logit(logit_m)

data.frame(global_wts_logit)

feature_imp_gbm <- feature_drop_loss(ov_subsets,data_train_x, data_train_y, gbm_train, gbm_m, cluster_centers_subsets)
global_wts_gbm <- get_global_wts_gbm(gbm_m)

### SEMI_PARAMETRIC bootstrapp, as non-parametric may not be ideal.

boot_n <- 100
err_p <- sum((data_train[,"test"] != pred_reg_insample))/length(data_train[,"test"])  ## Error
brm_acc_reg_in

err_p2 <- 0.2
#err_p2 <- 0.1
### We sample a binomial dist with the err_p and flip the predicted responses in each 
### bootstrapp

set.seed(20)
t = rbinom(length(pred_reg_insample),size = 1,prob = err_p)
#t[1:20]
pred_reg_insample_int = as.integer(pred_reg_insample)
pred_reg_insample_inv = pred_reg_insample_int
pred_reg_insample_inv[pred_reg_insample_inv == 2 ] <- 0
pred_reg_insample_ = pred_reg_insample_int - 1
pred_reg_insample_inv <- factor(pred_reg_insample_inv)

var_list2 <- list()
system.time(
  for(s in 1:boot_n)
  {
    set.seed(s*10)
    ind = as.logical(rbinom(length(pred_reg_insample),size = 1,prob = err_p2))
    newY <- pred_reg_insample
    newY[ind] <- pred_reg_insample_inv[ind] ### Over-writes some samples
    Y_dat <- data.frame(newY)
    names(Y_dat) <- names(data_train_y) 
    brm_model_bs <- brm_ov_subsets(num_blocks,data_train_x, Y_dat, low_threshold = 0.10)
    ov_subsets_bs <- brm_model_bs[[1]]
    cluster_centers_subsets_bs <- brm_model_bs[[2]]
    reg_m_bs <- lapply(ov_subsets_bs,logitreg)
    var_list2[[s]] <- get_global_wts(reg_m_bs)
  }
)
### 100 seconds for bootstrap = 100

# summary(reg_m[[1]])
# summary(reg_m[[2]])

var_list_t2 <- do.call(rbind,var_list2)
var_means2 <- apply(var_list_t2, 2, mean)
var_se2 <- apply(var_list_t2, 2, sd)
var_t2 <- (var_means2/var_se2)
pvals2 = 2*pt(-abs(var_t2), df=nrow(var_list_t2)-1)
round(pvals2,2)
data.frame(round(pvals2,4))

#ep_0.4 <- list(var_list_t, var_means, var_se, var_t, pvals)
ep_0.2 <- list(var_list_t2, var_means2, var_se2, var_t2, pvals2)

#ep_0.1 <- list(var_list_t2, var_means2, var_se2, var_t2, pvals2)

tail(pred_reg_insample)

### Smart join of all the coeff 

t = data.frame(summary(logit_m[[1]])$coef)
rownames(t)
subsett <- data.frame(var = rownames(t), coef = t[,1], pval = t[,4] )

global_reg <- data.frame(var = names(global_wts_logit), coef = global_wts_logit, pval = ep_0.2[[5]])
global_reg$coef[global_reg$pval > 0.05] <- NA

global_reg2 <- global_reg

for(m in 1:length(logit_m))
{
  t = data.frame(summary(logit_m[[m]])$coef)
  subsett <- data.frame(var = rownames(t), coef = t[,1], pval = t[,4] )
  subsett$coef[subsett$pval > 0.05] <- NA
  global_reg2 <- left_join(global_reg2, subsett[,c(1,2)], by =c("var"="var"))
  global_reg <- left_join(global_reg, subsett, by =c("var"="var"))
}
    

write.csv(global_reg2, "Covid19_coef_brm.csv",na="", row.names= FALSE)

accuracy_mat = rbind(cart_acc, listwise_acc_reg,svi_acc_reg, ense_acc_reg, brm_acc_reg, 
                     listwise_acc_gbm,svi_acc_gbm,ens_acc_gbm, brm_acc_gbm     )

write.csv(accuracy_mat, "Covid19_brm_predictions.csv",na="", row.names= FALSE)

