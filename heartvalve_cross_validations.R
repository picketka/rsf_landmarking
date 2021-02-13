rm(list=ls())

# Load packages -------------------------------------------------------------------------------
library("MASS") #mvnorm
library(survival) #survfit
library(ggplot2) #ggplot
library(JM)
#install.packages("randomForestSRC") #-- update package since there has been some changes
library(randomForestSRC)
library(timeROC)
library(tidyverse)
library(KMsurv)
library(km.ci)
library(pec)
library(prodlim)
library(lme4)
library(pROC)
library(joineRML)
library(dynpred)
library("JMbayes")
library(boot)
library(JM)
library(joineR)
#set global option that warnings be treated as errors for loops
#options(warn = 2)
library(gridExtra)
library(RColorBrewer)


# Set main directory (can replace with a string of the directory location)
dir0<-'C:/Users/'
source('C:/Users/heartvalve_prediction_functions.R')
#source('C:/Users/sureshk/Downloads/Kaci_App/2019_11_11_heartvalve_prediction_functions.R')

dat <- heart.valve
dat$age_10 <- dat$age/10
dat$prenyha <- factor(dat$prenyha)
dat$lv <- factor(dat$lv)
dat$emergenc <- factor(dat$emergenc)
dat$hc <- factor(dat$hc)
dat$sten.reg.mix <- factor(dat$sten.reg.mix)
dat$hs <- factor(dat$hs)

dat <- dat %>% mutate %>%
  group_by(num)%>%
  mutate(
    base_ef = first(ef, order_by = time),
    base_log.grad = first(log.grad, order_by = time),
    base_log.lvmi = first(log.lvmi, order_by = time),
    first_time = first(time),
    n = n()
  )

dat.id<-dat[!duplicated(dat$num),] 




#split the data into 5 individual datasets:
names_dat <- names(dat.id)
n_folds <- 5
#randomly order individuals
set.seed(1234)
dat.id    <- dat.id[order(runif(nrow(dat.id))), ]
bins  <- rep(1:n_folds, nrow(dat.id) / n_folds)
#split_ids <- split(tac.id$idn, sample(1:n_folds, nrow(tac.id), replace = T))
split_ids <- split(dat.id$num, c(bins,1))

####################################start of cross validation loop################################### 
i = 1
for(i in 1:n_folds){
  # j <- 2 #for testing individual runs
  # Parameters to change ------------------------------------------------------------------------
  ## 1. Simulate baseline covariates (how many? Continuous? Binary?) 
    
    #pull training data from i'th split dataset
    set <- as.data.frame(split_ids[i])
    names(set) <- c("num")
    train_data <-subset(dat,!dat$num%in%set$num)
    test_data <-subset(dat,dat$num%in%set$num) # longitudinal data set for Test
    train_data.id<-train_data[!duplicated(train_data$num),]
    test_data.id<-test_data[!duplicated(test_data$num),]
    
    #for loop to loop through all prediction times and save predictions for models of interest
    #need these variables later for many things so define them outside the function
    w_predict<-3
    LMx<-c(0.5,1,1.5,2,2.5,3) #c(2.5,3,4,5)#, 7)
    
    
    # Joint Model ---------------------------------------------------------------------------------
    ## Fit model from which data was simulated and compare estimated coefficients to true values
    ## Where are all of the parameters you used found in the fit model
    joint_model_time<- system.time({
      lmeFit<-try(lme(log.lvmi~time+ sex+age_10+ lv + con.cabg  ,data=train_data, #factor(lvh) + size +
                      random=~time|num)) #+creat, +sex*time
      survFit<-try(coxph(Surv(fuyrs,status)~sex+age_10+ lv + con.cabg ,data=train_data.id,x=TRUE)) #+creat
      jointFit<-try(jointModel(lmeFit,survFit,timeVar="time"))
    })
    
    ############################################Cox Supermodel based on stacked data sets###########################################################
    ############################################create stacked dataset and modify longitudinal dataset
    #functions of time
    dat$Time_noadmin <- dat$fuyrs
    dat$event_noadmin <- dat$status
   
    i <- 1
    g1 <- function(t) (t)
    g2 <- function(t) (t)^2
    LMdata <- NULL
    for (i in 1:length(LMx)){
      LMdat1<- cutLM(data=dat, outcome = list(time = "fuyrs", status = "status"),
                     LM=LMx[i],horizon=LMx[i]+w_predict,covs=list(fixed=c("sex","bsa","lvh" , "prenyha" ,
                                                                          "redo","size"  , "con.cabg" , "creat" ,  "dm"  , 
                                                                          "acei" ,"lv"   , "emergenc", "hc"  , "sten.reg.mix" ,
                                                                          "hs" , "age_10", "base_ef", "base_log.grad","base_log.lvmi", "Time_noadmin", "event_noadmin"), varying = c("log.lvmi")),
                     format = "long", id = "num", rtime = "time", right = F) #prior to 7/12/20 this was train_data
      LMdat2<- cutLM(data=dat, outcome = list(time = "fuyrs", status = "status"),
                     LM=LMx[i],horizon=LMx[i]+w_predict,covs=list(varying = c("log.grad")),
                     format = "long", id = "num", rtime = "time", right = F) 
      LMdat2 <- LMdat2[,c("num", "LM", "log.grad")]
      LMdat3<- cutLM(data=dat, outcome = list(time = "fuyrs", status = "status"),
                     LM=LMx[i],horizon=LMx[i]+w_predict,covs=list(varying = c("ef")),
                     format = "long", id = "num", rtime = "time", right = F) 
      LMdat3 <- LMdat3[,c("num", "LM", "ef")]
      
      all_data <- Reduce(function(x,y) merge(x = x, y = y, by = c("num", "LM")), 
                         list(LMdat1,LMdat2, LMdat3))
      
      LMdata <- rbind(LMdata, all_data)
    }
    
    
    LMdata$y_tau<-LMdata$log.lvmi*LMdata$LM
    LMdata$y_tau2<-LMdata$log.lvmi*LMdata$LM^2
    LMdata$LM1 <- g1(LMdata$LM)
    LMdata$LM2 <- g2(LMdata$LM)
    
    #create administrative censoring variables for longi dataset (keep same functions of LM and y for easier predictions)
    #create a "LM" variable that matches observation time
    train_data$LM <- train_data$time
    train_data$y_tau<-train_data$log.lvmi*train_data$LM
    train_data$y_tau2<-train_data$log.lvmi*train_data$LM^2
    train_data$LM1 <- g1(train_data$LM)
    train_data$LM2 <- g2(train_data$LM)
    
    
    train_data$event_admin <- ifelse(train_data$fuyrs <= train_data$LM+w_predict, train_data$status, 0) #do we need data in counting process format (should the event have not occured for longi measures before Event?)
    train_data$Time_admin <- ifelse(train_data$fuyrs <= train_data$LM+w_predict, train_data$fuyrs, train_data$LM+w_predict)
    
  
   LMdata_train <- subset(LMdata,!LMdata$num%in%set$num)
    ##############################################################################################################################################################################
    ###########################################################start of for loop for 'output predictions' for loop###########################################################
    for(t0 in LMx){
      print(t0)
      
      #w_predict<-3 #could eventually put this on the outside fo the loop too?
      summ_window <- 2 #amount of time to summarize longi measure prior to LM
      
      # 
      trunc_full.lm <- subset(LMdata, LM==t0)
      
      train_dat_sub.lm <-subset(trunc_full.lm,!trunc_full.lm$num%in%set$num) #longitudinal data set for Training 
      test_dat_sub.lm <-subset(trunc_full.lm,trunc_full.lm$num%in%set$num)
      
      
      out<-data.frame("num"=test_dat_sub.lm$num)
      out$pred_jm<-rep(NA,nrow(test_dat_sub.lm))
      
      #########################predictions for Joint model -- can only predict if joint model is not NA
      jm_predtime <- system.time({
        if(inherits(jointFit, "try-error")|nrow(test_dat_sub.lm) == 0) {
          out$pred_jm <- NA
        }else{
          x.new <- JM::survfitJM(jointFit,newdata=test_data[test_data$fuyrs > t0,],idVar="num",survTimes=t0+w_predict,last.time=rep(t0,nrow(test_dat_sub.lm)),simulate=FALSE)
          y.new <- t(matrix(unlist(x.new$summaries),nrow=2))
          out$pred_jm<-1-y.new[,2]
        }
      })#end of jm_pred_time
      
      ###########model with tuning 
      rsf_tuning_mod_time <- system.time({
        set.seed(j)
        model_form1 <- Surv(fuyrs, status) ~ sex+age_10+ lv + con.cabg+log.lvmi +log.grad +ef # + median + cv + min + avg + change
        rsf1_test <- try(tune.rfsrc(model_form1, data = train_dat_sub.lm, mtryStart = 1, nodesizeTry = c(1, 5, seq(10, 100, by = 5)), ntreeTry = 50,na.action= "na.impute",
                                    sampsize = 0.623*(nrow(train_dat_sub.lm)),
                                    nsplit = 10, stepFactor = 1.25, improve = 1e-2, strikeout = 3, maxIter = 25,
                                    trace = FALSE, doBest = TRUE)) #mtry is going to go from 1 to p = total number of parameters
        nodesize1 <- ifelse(inherits(rsf1_test, "try-error"), 15, rsf1_test$optimal[1])#default is 15
        mtry1 <- ifelse(inherits(rsf1_test, "try-error"), 5, rsf1_test$optimal[2]) #the recommended is p/3 = 4 but highly correlated vars so increase
        rsf1 <- try(rfsrc(model_form1, data = train_dat_sub.lm, ntree = 1000,nodesize = nodesize1, mtry = mtry1, na.action= "na.impute",
                          block.size = 1, statistics = TRUE, forest = TRUE, importance=TRUE, nsplit = 10))
        
        
        rsf_imp_tun<- as.data.frame(t(as.data.frame(rsf1$importance)))
        rsf_imp_tun$pred_time<- t0
        rsf_importance_tun <- rbind(rsf_importance_tun, rsf_imp_tun)
        #keep predictions  
        for (k in 1:nrow(test_dat_sub.lm)) {
          # Solve for cox ph probability at specified time
          preds_rsf_tun <- try(predictSurvProb.rsf_imp(rsf1, newdata = test_dat_sub.lm[k,], times = t0+w_predict))
          out$pred_rsf_tun[k] <- ifelse(inherits(preds_rsf_tun, "try-error"), NA, 1- preds_rsf_tun) #[k]
        }
      }) #end of tuning_mod_time
      
      #RSF with no tuning parameters
      rsf_mod_time<- system.time({
        set.seed(j)
        rsf2 <- try(rfsrc(model_form1, data = train_dat_sub.lm, ntree = 1000,nodesize = 15, mtry = 5, na.action= "na.impute",
                          block.size = 1, statistics = TRUE, forest = TRUE, importance=TRUE,nsplit = 10)) #nodesize and mtry optimized in 1 simulation then kept standard across sims and LM times
        
        rsf_imp <- as.data.frame(t(as.data.frame(rsf2$importance)))
        rsf_imp$pred_time<- t0
        rsf_importance <- rbind(rsf_importance, rsf_imp)  
      })#end of rsf_mod_time
       
       #keep predictions -- determine if we should predict if rsf_tac runs an error for no deaths in train data.... 
      rsf_pred_time_imp<- system.time({   
        for (k in 1:nrow(test_dat_sub.lm)) {
          # Solve for cox ph probability at specified time
          preds_rsf1 <- try(predictSurvProb.rsf_imp(rsf2, newdata = test_dat_sub.lm[k,], times = t0+w_predict))
          out$pred_rsf1[k] <- ifelse(inherits(preds_rsf1, "try-error"), NA, 1-preds_rsf1)#[k]
        }
      })#end of rsf_mod_time
      
      
      #do not impute on predictions BUT imput in model building
      rsf_pred_time<- system.time({   
        for (k in 1:nrow(test_dat_sub.lm)) {
          # Solve for cox ph probability at specified time
          preds_rsf_noimp1 <- try(predictSurvProb.rsf(rsf2, newdata = test_dat_sub.lm[k,], times = t0+w_predict))
          out$pred_rsf_noimp1[k] <- ifelse(inherits(preds_rsf_noimp1, "try-error"), NA, 1-preds_rsf_noimp1)#[k]
        }
      })#end of rsf_mod_time
      
      
      #run an rsf model with no imputation (versus above where we predict outcomes with no imputation)
      rsf_mod_time_noimp<- system.time({
        set.seed(j)
        rsf2b <- try(rfsrc(model_form1, data = train_dat_sub.lm, ntree = 1000,nodesize = 15, mtry = 5, na.action= "na.omit",
                          block.size = 1, statistics = TRUE, forest = TRUE, importance=TRUE,nsplit = 10)) #nodesize and mtry optimized in 1 simulation then kept standard across sims and LM times
        
        rsf_imp_noimp <- as.data.frame(t(as.data.frame(rsf2b$importance)))
        rsf_imp_noimp$pred_time<- t0
        rsf_importance_noimp <- rbind(rsf_importance_noimp, rsf_imp_noimp) 
        for (k in 1:nrow(test_dat_sub.lm)) {
          # Solve for cox ph probability at specified time
          preds_rsf_noimp2 <- try(predictSurvProb.rsf(rsf2b, newdata = test_dat_sub.lm[k,], times = t0+w_predict))
          out$pred_rsf_noimp2[k] <- ifelse(inherits(preds_rsf_noimp2, "try-error"), NA, 1-preds_rsf_noimp2)#[k]
        }
      })#end of rsf_mod_time
      
      
      
      #without administrative censoring -- using Time_noadmin instead of Time
      set.seed(j)
      model_form3 <- Surv(Time_noadmin, event_noadmin) ~ sex+age_10+ lv + con.cabg+log.lvmi +log.grad +ef
      rsf3 <- try(rfsrc(model_form3, data = train_dat_sub.lm, ntree = 1000,nodesize = 15, mtry = 5, na.action= "na.impute",
                        block.size = 1, statistics = TRUE, forest = TRUE, importance=TRUE,nsplit = 10))
      
      rsf_imp_noadmin <- as.data.frame(t(as.data.frame(rsf3$importance)))
      rsf_imp_noadmin $pred_time<- t0
      rsf_importance_noadmin <- rbind(rsf_importance_noadmin, rsf_imp_noadmin )
      
      #keep predictions -- determine if we should predict if rsf_tac runs an error for no deaths in train data.... 
      
      for (k in 1:nrow(test_dat_sub.lm)) {
        # Solve for cox ph probability at specified time
        preds3 <- try(predictSurvProb.rsf_imp(rsf3, newdata = test_dat_sub.lm[k,], times = t0+w_predict)) #impute for rn
        out$pred_rsf_noadmin[k] <- ifelse(inherits(preds3, "try-error"), NA, 1-preds3) #[k]
      }
      ###with many more additional covariates
      rsf_vars_mod_time<- system.time({
        set.seed(j)
        model_form4 <- Surv(fuyrs, status) ~ sex + age_10 + bsa + lvh + prenyha + redo + size + con.cabg+creat+dm + acei + 
          lv + emergenc+hc+sten.reg.mix + hs + log.lvmi  + ef + log.grad  #most likely too much missingness for these two 
        rsf4 <- try(rfsrc(model_form4, data = train_dat_sub.lm, ntree = 1000,nodesize = 15, mtry = 5,na.action= "na.impute",  
                          block.size = 1, statistics = TRUE, forest = TRUE, importance=TRUE,nsplit = 10))
        
        rsf_imp_vars <- as.data.frame(t(as.data.frame(rsf4$importance)))
        rsf_imp_vars$pred_time<- t0
        rsf_importance_vars <- rbind(rsf_importance_vars, rsf_imp_vars )
        
        #keep predictions -- determine if we should predict if rsf_tac runs an error for no deaths in train data.... 
        
        for (k in 1:nrow(test_dat_sub.lm)) {
          preds4<- try(predictSurvProb.rsf_imp(rsf4, newdata = test_dat_sub.lm[k,], times = t0+w_predict))
          # Solve for cox ph probability at specified time
          out$pred_rsf_vars[k] <- ifelse(inherits(preds4, "try-error"), NA, 1-preds4)#[k]
        }
      }) #end of rsf_vars_mod_time 
      
      ###run a non-administratively censored all variable rsf
      set.seed(j)
      model_form5 <- Surv(Time_noadmin, event_noadmin) ~ sex + age_10 + bsa + lvh + prenyha + redo + size + con.cabg+creat+dm + acei + 
        lv + emergenc+hc+sten.reg.mix + hs + log.lvmi  + ef + log.grad  #most likely too much missingness for these two 
      rsf_nv <- try(rfsrc(model_form5, data = train_dat_sub.lm, ntree = 1000,nodesize = 15, mtry = 5,na.action= "na.impute",  
                        block.size = 1, statistics = TRUE, forest = TRUE, importance=TRUE,nsplit = 10))
      #keep predictions
      for (k in 1:nrow(test_dat_sub.lm)) {
        preds_nv<- try(predictSurvProb.rsf_imp(rsf_nv, newdata = test_dat_sub.lm[k,], times = t0+w_predict))
        # Solve for cox ph probability at specified time
        out$pred_rsf_noadmin_vars[k] <- ifelse(inherits(preds_nv, "try-error"), NA, 1-preds_nv)#[k]
      }
      
      
      
      ###########model with tuning 
      rsf_tuning_vars_mod_time <- system.time({
        set.seed(j)
        rsf_vars_test <- try(tune.rfsrc(model_form4, data = train_dat_sub.lm, mtryStart = 1, nodesizeTry = c(1, 5, seq(10, 100, by = 5)), ntreeTry = 50,na.action= "na.impute",
                                    sampsize = 0.623*(nrow(train_dat_sub.lm)),
                                    nsplit = 10, stepFactor = 1.25, improve = 1e-2, strikeout = 3, maxIter = 25,
                                    trace = FALSE, doBest = TRUE)) #mtry is going to go from 1 to p = total number of parameters
        nodesize_vars <- ifelse(inherits(rsf_vars_test, "try-error"), 15, rsf_vars_test$optimal[1])#default is 15
        mtry_vars <- ifelse(inherits(rsf_vars_test, "try-error"), 5, rsf_vars_test$optimal[2]) #the recommended is p/3 = 4 but highly correlated vars so increase
        rsf_vars_tun <- try(rfsrc(model_form4, data = train_dat_sub.lm, ntree = 1000,nodesize = nodesize_vars, mtry = mtry_vars,na.action= "na.impute", 
                          block.size = 1, statistics = TRUE, forest = TRUE, importance=TRUE, nsplit = 10))
        
        rsf_imp_tun_vars<- as.data.frame(t(as.data.frame(rsf_vars_tun$importance)))
        rsf_imp_tun_vars$pred_time<- t0
        rsf_importance_tun_vars <- rbind(rsf_importance_tun_vars, rsf_imp_tun_vars)
        #keep predictions  
        for (k in 1:nrow(test_dat_sub.lm)) {
          # Solve for cox ph probability at specified time
          preds_rsf_tun_vars <- try(predictSurvProb.rsf_imp(rsf_vars_tun, newdata = test_dat_sub.lm[k,], times = t0+w_predict))
          out$pred_rsf_tun_vars[k] <- ifelse(inherits(preds_rsf_tun_vars, "try-error"), NA, 1- preds_rsf_tun_vars) #[k]
        }
      }) #end of tuning_vars_mod_time
      
     
      ############################################Cox Model at Landmark##########################################
      #correctly specified Cox model (w/o noise not with interactions)
      cox_time <-system.time({
        LMcox <- coxph( Surv(fuyrs, status) ~ sex+age_10+ lv+ con.cabg + log.lvmi , #
                        data=train_dat_sub.lm, method="breslow") #make sure the order matches cox_prob function
        
        cox_sm<-data.frame(Independent_variable=unlist(c("sex", "age_10", "lv2","lv3", "con.cabg", "log.lvmi")),
                        Estimate=paste0(round(coef(LMcox),2)," (",
                                        round(confint(LMcox)[,1],2),", ",
                                        round(confint(LMcox)[,2],2),")"),
                        Pvalue=round(summary(LMcox)$coef[,5],4),
                        row.names=NULL)
        cox_sm$pred_time <- t0
        cox_sm$ppl_train <- nrow(train_dat_sub.lm)
        cox_sm$ppl_test <- nrow(test_dat_sub.lm)
        
        cox_sm_mods <- rbind(cox_sm_mods, cox_sm)
        
        pred_cox <- NULL
        for (k in 1:nrow(test_dat_sub.lm)) {
          # Solve for cox ph probability at specified time
          preds_c <- try(cox_prob(LMcox,test_dat_sub.lm, t0, w_predict, k))
          pred_cox[k] <- ifelse(inherits(preds_c, "try-error"), NA, 1-preds_c)#function listed above
          out$pred_cox[k] <- ifelse(is.na(LMcox$coefficients[1]), NA, pred_cox[k])
        }
      })#end of cox_time
      
      
      LMcox_lg <- coxph( Surv(fuyrs, status) ~sex + age_10 + bsa + lvh + prenyha + redo + size + con.cabg+creat+dm + acei + 
                        lv + emergenc + hc + sten.reg.mix + hs + log.lvmi , 
                      data=train_dat_sub.lm, method="breslow") #make sure the order matches cox_prob function
      
      out$pred_time <- t0
      output_predictions <- rbind(output_predictions, out)
      out <- out[0,]
      
       print(t0)
    } #end of prediction function
    
} # end of full function   