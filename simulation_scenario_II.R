N <- 1000 #number of simulated individuals 
K <- 30 #number of planned repeated measurements per subject, per outcome
cens_horiz <- 15 #maximum follow-up time
insp.rate <- 2 #rate at which inspection times are made (increase to make inspections more frequent; mean time bw insp= 1/insp.rate)
## 2. Set true values for longitudinal model
## What values do you need? Fixed effects (intercept, slope); Random effects (intercept, slope)
## Measurement errror for longitudinal measurements 
sigma.y <- 0.6 #measurement error standard deviation
betas <- c("(Intercept)" = 5.6, #try: 5.6
           "time"=-0.45, #try:  -0.45
           "X1" = -0.25, #sex try: -0.25
           "X2" = -0.11, #"age" try: -0.01
           "X3" = -0.3) #, #BSA Try: -0.3
#            "X1:time" = -0.15,#does sex get more important over time? #try: -0.15
#            "X1:X2"  = -0.1, #try: -0.1
#            "X1:X3" = 0.1 #try 0.09
# ) #betas for longitudinal model
#modify to include additional covariates, slopes, interactions, etc. 
D <- matrix(c(1, 0.5, 
              0.5, 1), 2, 2) #covariance matrix for random effects 
D <- (D + t(D)) / 2 
#square matrix, dimension should be same as number of random effects

## 3. Set true values for survival model 
## What model are you using for the hazard (Weibull?)
## What parameters does that model need (Shape? Scale?)
gammas <- c("(Intercept)" = -8, #try: -5.75, -6.75, -7.75, -8.75, -10.75, - 11.75
            "X1" = 1.90, #try: 0.15, 0.75, 0.08, -0.75, 1.75
            "X2" = 2.15, #age effect on survival is strong try: 1.04, 2.04
            "X3" = 2.65)  #try: -0.65, - 0.35, 0.65(all events way too early), 0.35, 0.25, 0.15
#gammas: coefficients for Weibull hazard 
phi <- 2 #shape of Weibull hazard
#modify based on achieving desired shape/scale 

## 4. Association Parameters 
## How will the longitudinal model and the survival model be associated
## What kind of association (Positive? Negative?)
alpha0 <- 0.07 #association parameter #try: -0.5, -0.25, -0.75


######################################################################################
#######################################Functions Necessary############################
######################################################################################
###########change invS, cox_prob, and probability based on covariates included########
######################################################################################
######################################################################################
## Write a function of the hazard for each individual as a function of "s" and "u", a uniform
## random variable 
#!!!!!!!!!!!!!!!!needs to be adusted if additional covariates added!!!!!!!!!!!!!!!!
## Integrate this function w.r.t "s" over the interval (0,t), where you want to solve for "t" 
invS <- function (t, u, i) {
  #inverse function used to obtain survival time 
  #t: time 
  #u: draw from a uniform[0,1] distribution
  #i: patient id indicator 
  h <- function (s) {
    #hazard function (based on survival and longitudinal submodels)
    X1 <- X1[i]
    X2 <- X2[i]
    X3 <- X3[i]
    XX <- cbind(1, s, X1, X2, X3) #, X1*s, X1*X2, X1*X3
    ZZ <- cbind(1, s)
    f1 <- as.vector(XX %*% betas + rowSums(ZZ * b[rep(i, nrow(ZZ)), ]))
    haz<- exp(log(phi) + (phi - 1) * log(s) + eta.t[i] +  f1^2 * alpha0*log(1+s))
    return(haz)
  }
  integrate(h, lower = 0, upper = t)$value + log(1-u) 
}

#for each individual calculate the probability of conditional survival from a LM super model 
#!!!!!!!!!!!!!!!!needs to be adusted if additional covariates added!!!!!!!!!!!!!!!!
predict_supercox_3cov<-function(bh,bet,t0,xdata,w_predict, i)
{
  Fw <- NULL
  tti<-t0
  sfi<-bh
  sfi$haz0<-diff(c(0,sfi$hazard))
  xdata <- xdata[i,]
  sfi$hazard<-sfi$haz0*exp(bet[1]*xdata$X1 + bet[2]*xdata$X2 +
                             bet[3]*xdata$X3 +bet[4]*xdata$y + bet[5]*xdata$y*g1(tti) + bet[6]*xdata$y*g2(tti)+
                             bet[7]*g1(tti) + bet[8]*g2(tti)) 
  sfi$Haz <- cumsum(sfi$hazard)
  tmp <- evalstep(sfi$time,sfi$Haz,c(tti,tti+w_predict),subst=0)
  Fw <- c(Fw,exp(-(tmp[2]-tmp[1])))
  return(Fw)
}
#prediction for mispecified/noise Cox model(s?)
predict_supercox_10cov<-function(bh,bet,t0,xdata,w_predict, i)
{
  Fw <- NULL
  tti<-t0
  sfi<-bh
  sfi$haz0<-diff(c(0,sfi$hazard))
  xdata <- xdata[i,]
  sfi$hazard<-sfi$haz0*exp(bet[1]*xdata$X1 + bet[2]*xdata$X2 +
                             bet[3]*xdata$X3 +bet[4]*xdata$X4 + bet[5]*xdata$X5 +
                             bet[6]*xdata$X6 +bet[7]*xdata$X7 + bet[8]*xdata$X8 +bet[9]*xdata$X9 + 
                             bet[10]*xdata$X10 +bet[11]*xdata$y + bet[12]*xdata$y*g1(tti) + bet[13]*xdata$y*g2(tti)+
                             bet[14]*g1(tti) + bet[15]*g2(tti)) 
  sfi$Haz <- cumsum(sfi$hazard)
  tmp <- evalstep(sfi$time,sfi$Haz,c(tti,tti+w_predict),subst=0)
  Fw <- c(Fw,exp(-(tmp[2]-tmp[1])))
  return(Fw)
}
#for each individual calculate the probability of conditional survival from a coxph() model
cox_prob <- function(model,dat, t0, w, i){
  sf <- survfit(model, newdata = dat[i,]) #get the predicted cummulative hazard for each individual based on fit model
  sf2 <- data.frame(time = sf$time, surv = sf$surv, Haz = -log(sf$surv))
  tmp <- evalstep(sf2$time,sf2$Haz,c(t0,t0+w),subst=0) #landmark time and event horizon defined here to evaluate CHF
  Fw <- exp(-(tmp[2]-tmp[1])) #this is calculating probability of event?
  Fw
}

#function only works if betas, eta.t, b, alpha0 and phi are already named above it
#!!!!!!!!!!!!!!!!needs to be adusted if additional covariates added!!!!!!!!!!!!!!!!
probability <- function (t, w, i) {
  #inverse function used to obtain survival time 
  #t: time 
  #u: draw from a uniform[0,1] distribution
  #i: patient id indicator 
  h <- function (s) {
    X1<- X1[i] #i is defined by the function call--- we made it equal Test_data_sub.id$id
    X2 <- X2[i]
    X3 <- X3[i]
    XX <- cbind(1, s, X1, X2, X3) #, X1*s, X1*X2, X1*X3
    ZZ <- cbind(1, s)
    f1 <- as.vector(XX %*% betas + rowSums(ZZ * b[rep(i, nrow(ZZ)), ])) 
    #f1 <- as.vector(XX %*% betas + rowSums(ZZ * b[i, ])) #put in their true values of b here
    haz<- exp(log(phi) + (phi - 1) * log(s) + eta.t[i] +  f1^2 * alpha0*log(1+s)) 
    return(haz)
  }
  exp(-integrate(h, lower = t, upper = t+w, subdivisions = 10000)$value) #-- integrate from lower t0 to upper t0+w_predict exp(-)
} #end of probability function

###prediction function based on Morgensen PEC paper
predictSurvProb.rsf <- function (object, newdata, times, ...) {
  N <- NROW(newdata)
  class(object) <- c("rfsrc", "grow")
  S <- predict.rfsrc(object, newdata = newdata)$survival
  if(N == 1) S <- matrix(S, nrow = 1)
  Time <- object$time.interest
  p <- cbind(1, S)[, 1 + sindex(Time, times),drop = FALSE]
  if(NROW(p) != NROW(newdata) || NCOL(p) != length(times))
    stop("Prediction failed") #prediction fails if predictions cannot be made for every subject
  p
}

############################Calculate time dependent Brier Score and AUC#################################
#custom made BS calculation function (edited from Blanche paper)
BS <- function(timepoints,times,status,pred,cause=1){ 
  n <- length(times)
  n_times <- length(timepoints)
  timepoints <- timepoints[order(timepoints)]
  times_names <- paste("t=",timepoints,sep="")
  # output initialisation 
  BS <- rep(NA,n_times)
  CumInci <- rep(NA,n_times)
  surv <- rep(NA,n_times)
  Stats <- matrix(NA,nrow=n_times,ncol=4)
  hit1_all <- matrix(NA,nrow=n,ncol=n_times)
  hit2_all <- matrix(NA,nrow=n,ncol=n_times)
  epsilon_i <- matrix(NA,nrow=n,ncol=n_times)
  #adds name to outputs
  names(BS) <- times_names
  names(CumInci) <- times_names
  names(surv) <- times_names
  colnames(Stats) <- c("Cases","survivor at t","Other events at t","Censored at t")
  rownames(Stats) <- times_names
  colnames(epsilon_i) <- times_names
  colnames(hit1_all) <-  times_names
  colnames(hit2_all)  <- times_names 
  # we need to order to use the ipcw() function of the pec package
  #browser()
  order_T <- order(times)
  times <-  times[order_T]
  delta  <-  status[order_T]
  pred <-  pred[order_T,,drop=FALSE]
  #compute KM weights
  weights <- ipcw(Surv(failure_time,status)~1,
                  data=data.frame(failure_time=times,status=as.numeric(delta!=0)),
                  method="marginal",times=timepoints,subjectTimes=times,subjectTimesLag=1)
  Mat_data <- cbind(times,delta,as.numeric(delta==0))
  colnames(Mat_data) <- c("T","delta","indic_Cens")
  # computate weights of cases
  Weights_cases_all <- 1/(weights$IPCW.subjectTimes*n)
  # loop on all time points
  for(t in 1:n_times){
    Cases <- (Mat_data[,"T"]<= timepoints[t] &  Mat_data[,"delta"]==cause)
    Controls_1 <- (Mat_data[,"T"]> timepoints[t] )
    Controls_2 <- (Mat_data[,"T"]<= timepoints[t] &  Mat_data[,"delta"]!=cause & Mat_data[,"delta"]!=0)  
    # compute weights
    Weights_controls_1 <- rep(1/(weights$IPCW.times[t]*n),times=n)
    Weights_cases <- Weights_cases_all
    Weights_controls_2 <- Weights_cases_all
    Weights_cases[!Cases] <- 0
    Weights_controls_1[!Controls_1] <- 0
    Weights_controls_2[!Controls_2] <- 0   
    #compute outputs
    CumInci[t] <- c(sum(Weights_cases))
    surv[t] <- c(sum(Weights_controls_1))
    Stats[t,] <- c(sum(Cases),sum(Controls_1),sum(Controls_2),n-sum(Cases)-sum(Controls_1)-sum(Controls_2)) 
    hit1_all[,t] <- (Weights_controls_1*((pred[,t])^2))*n
    hit2_all[,t] <- (Weights_cases*((1-pred[,t])^2) + Weights_controls_2*((pred[,t])^2))*n
    BS[t] <- (sum(hit1_all[,t], na.rm = T) +sum(hit2_all[,t], na.rm = T))/n #added na.rm = T to this to fix BS calc for RSF...
  } 
  
  out <- list(BS=BS,res=(hit1_all+hit2_all),
              CumulativeIncidence=CumInci,survProb=surv,n=n,Stats=Stats,timepoints=timepoints
  )
  class(out) <- "ipcwEBS"
  out 
}

#caluclate AUC and BS based on predictions
PE<-function(model,LMx,w_predict,BS_dat)
{    
  AUCst.M1 <- rep(NA,length(LMx))
  BrierS.s.M1 <- rep(NA,length(LMx))
  for (s in LMx){
    # print(s)
    # Create landmark data set
    d.s<-BS_dat[,c("Time","event")]
    d.s$Pred.s.M1<-BS_dat[,model]
    d.s<-subset(d.s,BS_dat$LM==s) #subset(d.s,d.s$Time>s)
    d.s$time.s<-d.s$Time-s #subtract landmark time from Time of end event (censor or event)
    # AUC and BS for prediction based on M1
    # estimate ROC curve and AUC
    ROC.s.M1<-timeROC(T=d.s$time.s,
                      delta=d.s$event,
                      marker=d.s$Pred.s.M1,
                      cause=1,weighting="marginal",
                      times=c(w_predict),
                      iid=TRUE)
    # estimate expected Brier score
    BS.s.M1 <- BS(timepoints=c(w_predict),
                  times=d.s$time.s,
                  status=d.s$event,
                  pred=as.matrix(d.s$Pred.s.M1),
                  cause=1)
    # save useful results (estimates, s.e. and iid decompositions)
    BrierS.s.M1[which(s==LMx)] <- BS.s.M1$BS # BS estimate
    AUCst.M1[which(s==LMx)] <- ROC.s.M1$AUC[2] # AUC estimate
  }
  return(list(BrierS.s.M1,AUCst.M1))
}


######################################END NECESSARY FUNCTIONS BEGIN SIMULATION##############################
seed <- j
set.seed(seed)
X1 <- rbinom(N,size=1,prob=0.60) #tried 0.4 first. #simulate binary covariate "sex"(proportion = 60% male) #decreased from 0.71 male for better sim performance?
X2 <- rbinom(N,size=1,prob=0.50)#rnorm(N,mean=6.5,sd=1.2) #simulate continuous covariate "age/10" -- left skewed dist???
X3 <- rbinom(N,size=1,prob=0.55)#rnorm(N, mean = 1.86, sd = 0.23) #simulate continuous covariate "bsa"
X4 <- rbinom(N,size=1,prob=0.50)
X5 <- rbinom(N,size=1,prob=0.65)
X6 <- rbinom(N,size=1,prob=0.45)
X7 <- rbinom(N,size=1,prob=0.55)
X8 <- rbinom(N,size=1,prob=0.70)
X9 <- rbinom(N,size=1,prob=0.60)
X10 <- rbinom(N,size=1,prob=0.50)
## 3. Simulate random effects
## b: simulated random effects b using covariance matrix D 
b <- mvrnorm(N, rep(0, nrow(D)), D)
# Create data set with baseline covariate data ------------------------------------------------
##1. Create a data set with simulated baseline covariates 
##2. Assign a patient ID to each row
base_dat<-data.frame(id=1:N,X1, X2, X3, X4, X5, X6, X7, X8, X9, X10)

# Longitudinal submodel -----------------------------------------------------------------------
#yi(t)=mi(t)+ei(t)=Xi'(t)*beta+zi'(t)*bi+ei(t)
#bi~N(0,D); ei(t)~N(0,sigma.y^2)
## 1. What are the time points at which to generate longitudinal measurements 
## (Fixed?, Random?)
# Fixed measurement times: 
#times <- c(replicate(N, 0:(K-1)))
# Random measurement times: 
times <- c(replicate(N,cumsum(c(0,rexp(n=K-1,rate=insp.rate))))) #exponential 
t.max <- 15
#times <- c(replicate(N, c(0, sort(runif(K-1, 0, t.max))))) #uniform 
max.FUtime <- max(times) + 2 * IQR(times)

times_dat<-data.frame(id=rep(1:N,each=K),time=times)
DF<-merge(times_dat,base_dat,by="id",all.x=TRUE)
## 2. Create design matrices X and Z 
## X: design matrix for fixed effects (should have as many columns as elements in "betas")
## Z: design  matrix for random effects (should have as many columns as elements in "D")
X <- model.matrix(~ time + X1+X2 +X3 , data = DF) #+ X1:time+ X1:X2+ X1:X3 
Z <- model.matrix(~ time, data = DF)

id <- rep(1:N, each = K)
eta.y <- as.vector(X %*% betas + rowSums(Z * b[id, ])) #Xi'(t)*beta+zi'(t)*bi

# Survival submodel ---------------------------------------------------------------------------
## 1. Create design matrix for survival model 
## W: design matrix for survival (should have as many columns as "gammas")
W<- cbind("(Intercept)" = 1, "X1" = X1, "X2" = X2, "X3" = X3)


eta.t <- as.vector(W %*% gammas)

## 4. Simulate longitudinal responses
y <- rnorm(N * K, eta.y, sigma.y) #~Normal(mean=eta.y, sd=sigma.y) 

## 2. Simulate event times using function 
## 3. Solve for survival times 
u <- runif(N) #generate random U[0,1] random variable 
trueTimes <- numeric(N) #initiate vector of true survival times
for (i in 1:N) {
  # Solve for true survival time 
  Root <- try(uniroot(invS, interval = c(1e-05, max.FUtime), u = u[i],
                      i = i,extendInt="upX")$root, TRUE)
  # "Cure" patients have event time set to Infinity 
  trueTimes[i] <- ifelse(inherits(Root, "try-error"),Inf,Root)
}

## 4. Simulate censoring times 
## What kind of distribution do you want for this? (Uniform?, Exponential?)
Ctimes<-runif(N,0,cens_horiz)
#Ctimes <- cens_horiz
Time <- pmin(trueTimes, Ctimes,rep(cens_horiz,N))
event <- ifelse(trueTimes<=Time,1,0)

# Create simulated data set -------------------------------------------------------------------
## Create the datasets required for analysis (joint model, landmark, etc.)
## 1. Create longitudinal data set of simulated individuals (including baseline covariates)
Time_dat<-data.frame(id=rep(1:N,each=K),Time=rep(Time,each=K))
## Drop the longitudinal measurements that were taken after the observed event time for each subject
ind <- times <= Time_dat$Time
y <- y[ind]
X <- X[ind, , drop = FALSE]
Z <- Z[ind, , drop = FALSE]
id <- id[ind]
id <- match(id, unique(id))

dat <- DF[ind, ]
dat$id <- id
dat$y <- y #longitudinal marker value
dat$Time <- Time[id] #event time
dat$event <- event[id] #indication of event at "Time" 

#make baseline (first) values a covariate
dat_base <- dat[dat$time == 0, ] #baseline longi measure as predictor
dat_base$base_y <- dat_base$y
dat_base <- dat_base[,c('id', 'base_y')]
dat <- merge(dat, dat_base, by = 'id')

## 2. Create baseline data set of simulation individuals (num rows = N)
set <- sort(sample(unique(id), N/2))
train_data <-subset(dat,!dat$id%in%set) #longitudinal data set for Training 
test_data <-subset(dat,dat$id%in%set) # longitudinal data set for Test
train_data.id<-train_data[!duplicated(train_data$id),]
test_data.id<-test_data[!duplicated(test_data$id),]

#for loop to loop through all prediction times and save predictions for models of interest
#need these variables later for many things so define them outside the function
w_predict<-5
LMx<-c(0.5,1.5,2.5, 3.5, 4.5)     

  # Joint Model ---------------------------------------------------------------------------------
  ## Fit model from which data was simulated and compare estimated coefficients to true values
  ## Where are all of the parameters you used found in the fit model
  #correctly specified model like how simulated
  joint_model_time<- system.time({
    lmeFit<-try(lme(y~time+X1+ X2 + X3 ,data=train_data, random=~time|id)) 
    survFit<-try(coxph(Surv(Time,event)~X1+ X2 + X3 ,data = train_data.id,x=TRUE))
    jointFit<-try(jointModel(lmeFit,survFit,timeVar="time"))
  })
  
  
  ############################################JM Prediction + Cox and RSF Models (Same models fit in Scenario I & II)#############################
  for(t0 in LMx){
    trunc_full <- dat[dat$Time > t0,]
    
   #sort by id and time
    trunc_full <- trunc_full[order(trunc_full$id, trunc_full$time),]
    #calculate  slope from baseline
    trunc_full$slope <- (trunc_full$y-trunc_full$base_y)/trunc_full$time
    #create a non-administratve censored Time variable to see how RSF performs
    trunc_full$Time_noadmin <- trunc_full$Time
    trunc_full$event_noadmin <- trunc_full$event
    
    train_dat_sub <-subset(trunc_full,!trunc_full$id%in%set) #longitudinal data set for Training 
    test_dat_sub <-subset(trunc_full,trunc_full$id%in%set) # longitudinal data set for Test
    
    train_dat_sub.lm <- cutLM(data = train_dat_sub, outcome = list(time = "Time", status = "event"),
                              LM = t0, horizon = t0+w_predict, covs = list(fixed = c("X1", "X2", "X3", "X4","X5", "X6", "X7", "X8","X9", "X10",
                                                                                     "base_y","median", "cv", "min", "avg", "change","slope", "Time_noadmin", "event_noadmin"), varying = c("y")), #removed , "median", "cv", "min" from fixed
                              format = "long", id = "id", rtime = "time", right = F) #right = F makes the value of y come from time of landmark not last obs prior
    
    test_dat_sub.lm <- cutLM(data = test_dat_sub, outcome = list(time = "Time", status = "event"),
                             LM = t0, horizon = t0+w_predict, covs = list(fixed = c("X1", "X2", "X3", "X4","X5", "X6", "X7", "X8","X9", "X10",
                                                                                    "base_y","median", "cv", "min", "avg", "change","slope", "Time_noadmin", "event_noadmin"), varying = c("y")), #, "median", "cv", "min" in fixed
                             format = "long", id = "id", rtime = "time", right = F)
    
    #create variables to match long test dataset
    test_dat_sub.lm$y_tau<-test_dat_sub.lm$y*test_dat_sub.lm$LM
    test_dat_sub.lm$y_tau2<-test_dat_sub.lm$y*test_dat_sub.lm$LM^2
    
    out<-data.frame("id"=test_dat_sub.lm$id)
    
    
    #########################predictions for Joint model -- can only predict if joint model is not NA
    jm_predtime <- system.time({
      if(inherits(jointFit, "try-error")) {
        out$pred_jm <- NA
      }else{
        x.new <- JM::survfitJM(jointFit,newdata=test_dat_sub,idVar="id",survTimes=t0+w_predict,last.time=rep(t0,nrow(test_dat_sub.lm)),simulate=FALSE)
        y.new <- t(matrix(unlist(x.new$summaries),nrow=2))
        out$pred_jm<-1-y.new[,2]
      }
    })#end of jm_pred_time
    
    
    ##########################Random Survival Forest with Tuning of Parameters###############################
    model_form1 <- Surv(Time, event) ~ X1+X2+X3+ y 
    rsf1_test <- try(tune.rfsrc(model_form1, data = train_dat_sub.lm, mtryStart = 1, nodesizeTry = c(1, 5, seq(10, 100, by = 5)), ntreeTry = 50,
                                sampsize = 0.623*(nrow(train_dat_sub.lm)),
                                nsplit = 10, stepFactor = 1.25, improve = 1e-2, strikeout = 3, maxIter = 25,
                                trace = FALSE, doBest = TRUE)) #mtry is going to go from 1 to p = total number of parameters
    nodesize1 <- ifelse(inherits(rsf1_test, "try-error"), 5, rsf1_test$optimal[1])
    mtry1 <- ifelse(inherits(rsf1_test, "try-error"), 2, rsf1_test$optimal[2]) #the recommended is p/3 = 4 but highly correlated vars so increase
    rsf1 <- try(rfsrc(model_form1, data = train_dat_sub.lm, ntree = 1000,nodesize = nodesize1, mtry = mtry1, 
                      block.size = 1, statistics = TRUE, forest = TRUE, importance=TRUE, nsplit = 10))
    
    
    #keep predictions  
    preds_rsf_tun <- try(predictSurvProb.rsf(rsf1, newdata = test_dat_sub.lm, times = t0+w_predict))
    for (k in 1:nrow(test_dat_sub.lm)) {
      # Solve for cox ph probability at specified time
      out$pred_rsf_tun[k] <- ifelse(inherits(preds_rsf_tun, "try-error"), NA, 1- preds_rsf_tun[k])
    }
    
    ############################################Cox Model at each Landmark##########################################
    #correctly specified Cox model (w/o noise not with interactions)
      LMcox <- try(coxph( Surv(Time, event) ~ X1+X2+X3 + y  , #"mispecified but not +X1:X2 + X1:X3
                          data=train_dat_sub.lm, method="breslow")) #make sure the order matches cox_prob function
      pred_cox <- NULL
      for (k in 1:nrow(test_dat_sub.lm)) {
        # Solve for cox ph probability at specified time
        pred_cox <- NULL
        out$pred_cox[k] <- ifelse(inherits(LMcox, "try-error"), NA, 1-cox_prob(LMcox,test_dat_sub.lm, t0, w_predict, k))#function listed above
      }
    
  }#end of loop for each landmark
    