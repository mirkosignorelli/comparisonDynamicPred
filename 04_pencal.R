ncores = 10

dataset = 'ROSMAP' #PBC2/ADNI/ROSMAP
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
library(pencal)
library(survival)
library(riskRegression)
library(dplyr)
library(magrittr)
library(survcomp)
library(data.table)

for (l in lmark){
  matrix.C.index = matrix(NA,n.RCV,n.fold)
  # Prediction times
  times = c((l+1):tmax)

  # use original or transformed variables?
  long = long.transformed
  # covariate names
  fixed.cov.names = all.fixed.cov.names
  long.cov.names = all.long.cov.names
  if(dataset == 'ROSMAP') long.cov.names = names.1percMiss
  
  # Set time dependent variable used in lmms
  long$t.var.lmm = long$time.fup
  # Fixed and Random effects
  fixefs = ~ t.var.lmm
  ranefs = ~ t.var.lmm | id
  # Baseline covariates
  if (dataset == 'PBC2'){
    baseline.covs =  ~ age.baseline + sex + drug
    pfac.base.covs = c(0,1,1)
  }
  if (dataset == 'ADNI'){
    baseline.covs =  ~ age.baseline + PTGENDER + PTEDUCAT + status.bl + APOE4
    pfac.base.covs = c(0,1,1,1,1,1)
  }
  if (dataset == 'ROSMAP'){
    baseline.covs =  ~ age_bl + msex + educ + dx_bl + cancer_bl
    pfac.base.covs = c(0,1,1,1,1)
  }
  
  # Pencal needs the variable 'time.to.event' to be called 'time'
  surv$time = surv$time.to.event
  
  # Landmarking ------------------------------------------------------------------
  #remove subjects who had event/were censored before landmark
  surv = subset(surv, time.to.event >= l)
  long = subset(long, id %in% surv$id)
  
  # remove measurements after landmark
  long = subset(long, time.fup <= l)
  
  # number of subjects retained in the analysis:
  print(paste0('Landmarking at ',l,', number of subjects retained in the analysis: ',nrow(surv)))
  
  # number of events (1) and censored observations:
  print('Number of censored (0) and events (1):')
  table(surv$event)
  
  # number of repeated measurements per subject:
  print('Number of repeated measurements per subject:')
  table(table(long$id))
  
  # ------------------------------------------------------------------------------
  # Looking/removing subjects with all NAs for a given covariate (after landmarking)
  # Store into a list the ids with all NAs for each covariate
  j = 0
  ids = list()
  rep.meas.per.sub = aggregate(time.fup ~ id, long, function(x){ sum(!is.na(x)) }, na.action = NULL)
  for(i in long.cov.names){
    j = j + 1
    tmp.NAperSub = aggregate(long[,c(i)] ~ id, long, function(x){ sum(is.na(x)) }, na.action = NULL)
    tmp.ids = rep.meas.per.sub$id[which(tmp.NAperSub[,2] == rep.meas.per.sub[,2])]
    ids[[j]] = tmp.ids
  }
  
  names(ids) = long.cov.names
  
  # Remove all the ids found
  print(paste0('Number of subjects with all NAs for at least one covariate: ',length(unique(unlist(ids))),' (removed)'))
  ids.to.remove = unique(unlist(ids))
  surv = surv[!surv$id %in% ids.to.remove,]
  long = long[long$id %in% surv$id,]
  
  #-------------------------------------------------------------------------------
  # Repeated k-fold cross validation
  # List to store mean of performance measures for each repetition
  q = 0
  mean.perf.RCV = list()
  
  # Repeated CV loop
  start.time = Sys.time()
  for(seed in seeds.RCV){
    perf.list = list()

    # k-fold CV loop
    folds = Create_folds(surv, n.fold, seed) #ids of surv.data for each fold
    for(k.fold in 1:n.fold){
      # Print loop information
      print(paste0('Lmark: ',l,' Repetition: ',(seed+1-SEED),' Fold: ',k.fold))
      
      #training data and test data
      id.test = folds[[k.fold]]$ids.test
      long.train = long[!long$id%in%id.test, ]
      long.test = long[long$id%in%id.test, ]
      surv.train = surv[!surv$id%in%id.test, ]
      surv.test = surv[surv$id%in%id.test, ]
      
      #Estimation of the linear mixed models
      step1 = fit_lmms(y.names = long.cov.names, fixefs = fixefs, ranefs = ranefs, 
                       long.data = long.train, surv.data = surv.train, 
                       t.from.base = time.fup, n.boots = 0, n.cores = ncores,
                       max.ymissing = 0.2, seed = 123)
      
      #Summarize lmms
      step2 = summarize_lmms(object = step1, n.cores = ncores)
      
      #Fit penalized regression calibration
      step3 = fit_prclmm(object = step2, surv.data = surv.train,
                            baseline.covs =  baseline.covs,
                            penalty='ridge', standardize=TRUE,
                            pfac.base.covs=pfac.base.covs, n.cores = ncores)
      
      #Performance indeces
      surv.preds = survpred_prclmm(step1, step2, step3, times = times,
                                   new.longdata = long.test,
                                   new.basecovs = surv.test)
      
      fail.preds = 1 - surv.preds$predicted_survival[ , -1]
      colnames(fail.preds) = times
      
      # Brier score
      X = Score(as.list(fail.preds), data=surv.test, formula=Surv(time.to.event,event)~1,
                metric=c('brier','auc'), times=times,exact = FALSE,conf.int = FALSE,cens.model = "cox",
                splitMethod = "none", B = 0, verbose = FALSE)
      
      auc = subset(X$AUC$score, model == times)
      brier = subset(X$Brier$score, model == times)
      
      # Concordance index
      X = concordance.index(x=fail.preds[,1], surv.time=surv.test$time.to.event,
                            surv.event=surv.test$event, method='noether')
      CI = X$c.index
      
      # Store C index to compute standard deviation
      matrix.C.index[seed+1-SEED,k.fold] = CI
      
      # Performance indices table
      table = cbind(rep(l, length(times)), times, auc$AUC, brier$Brier, rep(CI,length(times)))
      colnames(table) = c('lmark', 't', 'AUC', 'BS', 'CI')
      
      # Store performance table
      perf.list[[k.fold]] = as.data.frame(table)
    }
    # average over the CV run
    q = q+1
    mean.perf.RCV[[q]] = rbindlist(perf.list)[,lapply(.SD,mean,na.rm=TRUE), list(lmark,t)]
  }
  
  # average over the RCV
  PERF = rbindlist(mean.perf.RCV)[,lapply(.SD,mean,na.rm=TRUE), list(lmark,t)]
  print(PERF)
  
  #-------------------------------------------------------------------------------
  ## Save PERFORMANCE table
  model = 'pencal'
  file.name = paste0(dataset,"_",model,"_",l*10,".RData")
  file.path = paste0('results/',dataset,'/',file.name)
  save(PERF, mean.perf.RCV, file = file.path)

  # Computing time
  end.time = Sys.time()
  comp.time[which(lmark == l)] = as.numeric(difftime(end.time, start.time, units='secs'))
}

save(comp.time, file = paste0('ctime/', dataset, '_', model, '.RData'))
