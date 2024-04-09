
dataset = 'PBC2' #PBC2/ADNI/ROSMAP
if (dataset == 'PBC2') {
  load("data/pbc2.RData")
  lmark = c(2.5, 3, 3.5)
  tmax = 8
}
if (dataset == 'ADNI') {
  load("data/adni.RData")
  lmark = 2:4
  tmax = 10
}
if (dataset == 'ROSMAP') {
  load("data/rosmap.RData")
  lmark = 2:6
  tmax = 15
}

# info for RCV
n.fold = ifelse(dataset == 'PBC2', 5, 10)
n.RCV = ifelse(dataset == 'PBC2', 20, 10)
SEED = 2023
seeds.RCV = SEED:(SEED + n.RCV-1)

comp.time = rep(NA, length(lmark))

source('00_functions.R')
library(pec)
library(riskRegression)
library(survcomp)
library(magrittr)
library(dplyr)
library(data.table)

for (l in lmark){
  matrix.C.index = matrix(NA,n.RCV,n.fold)
  # Prediction times
  times = c((l+1):tmax)
  
  # use original or transformed variables?
  long = long.original
  # covariate names
  fixed.cov.names = all.fixed.cov.names
  long.cov.names = all.long.cov.names
  if(dataset == 'ROSMAP') long.cov.names = names.1percMiss
  
  # LANDMARKING - LOCF dataset --------------------------------------------------
  #remove subjects who had event/were censored before landmark
  surv = subset(surv, time.to.event >= l)
  long = subset(long, id %in% surv$id)
  # remove measurements after landmark
  long = subset(long, time.fup <= l)
  
  surv.LOCF = surv[,c('id','time.to.event','event',fixed.cov.names)]
  
  # LOCF - last observation carried forward
  # [NB] I am assuming that long data are ordered according to time.fup
  i = 0
  for(id in surv$id){
    i = i+1
    for(q in long.cov.names){
      tmp = long[long$id == id,q]
      if(length(tail(tmp[!is.na(tmp)],1)) == 0)
        surv.LOCF[i,q] = NA
      else
        surv.LOCF[i,q] = tail(tmp[!is.na(tmp)],1)
    }
  }
  
  # remove subjects that have a NA values
  surv.LOCF = surv.LOCF[complete.cases(surv.LOCF),]
  
  # Repeated k-folds cross validation --------------------------------------------
  q = 0
  mean.perf.RCV = list()
  
  # Repeated CV loop
  start.time = Sys.time()
  for(seed in seeds.RCV){
    perf.list = list()

    # k-fold CV loop
    folds = Create_folds(surv.LOCF, n.fold, seed) #ids of surv.data for each fold
    for(k.fold in 1:n.fold){
      
      # Print loop information
      print(paste0('Lmark: ',l,' Repetition: ',(seed+1-SEED),' Fold: ',k.fold))
      
      #training data and test data
      id.test = folds[[k.fold]]$ids.test
      surv.LOCF.train = surv.LOCF[!surv.LOCF$id%in%id.test, ]
      surv.LOCF.test = surv.LOCF[surv.LOCF$id%in%id.test, ]
      
      # Cox model
      fit.cox = coxph(Surv(time.to.event, event) ~., data = surv.LOCF.train[,-1], 
                      model = TRUE, x = TRUE)
      
      # Predict survival probabilities
      surv.prob = predictSurvProb(fit.cox, surv.LOCF.test, times=times)
      fail.prob = as.data.frame(1 - surv.prob)
      colnames(fail.prob) = times
      
      # Brier Score
      X = Score(as.list(fail.prob), data=surv.LOCF.test, formula=Surv(time.to.event, event) ~1,
                metrics=c('brier','auc'),times=times, exact = FALSE,conf.int = FALSE,
                cens.model = "cox",splitMethod = "none", B = 0,
                verbose = FALSE)
      auc = subset(X$AUC$score, model == times)[,-1]
      brier = subset(X$Brier$score, model == times)[,-1]
      
      # Concordance index
      X = concordance.index(x=fail.prob[,1], surv.time=surv.LOCF.test$time,
                            surv.event=surv.LOCF.test$event, method='noether')
      CI = X$c.index
      
      # Store C index to compute standard deviation
      matrix.C.index[seed+1-SEED,k.fold] = CI
      
      # Performance indices Table
      table = cbind(rep(l, length(times)), times, auc$AUC, brier$Brier, rep(CI, length(times)))
      colnames(table) = c('lmark', 't', 'AUC', 'BS', 'CI')
      
      # Store performance table
      perf.list[[k.fold]] = as.data.frame(table)
    }
    # average over the CV run
    q = q+1
    mean.perf.RCV[[q]] = rbindlist(perf.list)[,lapply(.SD,mean,na.rm=TRUE),list(lmark,t)]
  }
  
  # average over the RCV
  PERF = rbindlist(mean.perf.RCV)[,lapply(.SD,mean,na.rm=TRUE),list(lmark,t)]
  print(PERF)
  
  ## Save PERFORMANCE table ------------------------------------------------------
  model = 'landmarking'
  file.name = paste0(dataset,'_',model,'_',l*10,".RData")
  file.path = paste0('results/',dataset,'/',file.name)
  save(PERF, mean.perf.RCV, file = file.path)
  
  # Computing Time
  end.time = Sys.time()
  comp.time[which(lmark == l)] = as.numeric(difftime(end.time, start.time, units = 'secs'))
}

save(comp.time, file = paste0('ctime/', dataset, '_', model, '.RData'))

