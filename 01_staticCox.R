
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

# use original or transformed variables?
long = long.original
# covariate names
fixed.cov.names = all.fixed.cov.names
long.cov.names = all.long.cov.names
if(dataset == 'ROSMAP') long.cov.names = names.1percMiss

# info for RCV
n.fold = ifelse(dataset == 'PBC2', 5, 10)
n.RCV = ifelse(dataset == 'PBC2', 20, 10)
SEED = 2023
seeds.RCV = SEED:(SEED + n.RCV-1)

comp.time = rep(NA, length(lmark))

source('00_functions.R')
library(magrittr)
library(dplyr)
library(survcomp)
library(pec)
library(riskRegression)
library(data.table)
library(survival)

for (l in lmark){
  matrix.C.index = matrix(NA,n.RCV,n.fold)
  # Prediction times
  times = c((l+1):tmax)
  
  # Landmark the data ------------------------------------------------------------------
  # remove subjects who had event/were censored before landmark
  surv = subset(surv, time.to.event >= l)
  long = subset(long, id %in% surv$id)
  # remove measurements after landmark
  long = subset(long, time.fup <= l)
  
  # ------------------------------------------------------------------------------
  # [DO THIS STEP TO HAVE THE SAME DATASET FOR ALL THE METHODS]
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
  
  #-----------------------------------------------------------------------------
  
  # Keep just values at baseline
  df.baseline = subset(long, time.fup == 0)
  df.baseline = cbind(surv[,c('id','time.to.event','event',fixed.cov.names)],
                      df.baseline[,long.cov.names])
    
  # Remove subject with NA values at baseline
  df.baseline = na.omit(df.baseline)
  print(paste0('Number of subjects removed because of NA values at baseline: ',
               nrow(surv) - nrow(df.baseline)))
  
  #-----------------------------------------------------------------------------
  # [TO SUM UP]
  # number of subjects retained in the analysis:
  print(paste0('Landmarking at ',l,', number of subjects retained in the analysis: ',nrow(df.baseline)))
  
  # number of events (1) and censored observations:
  print('Number of censored (0) and events (1):')
  table(df.baseline$event)
  
  #-------------------------------------------------------------------------------
  # Repeated k-fold cross validation
  q = 0
  mean.perf.RCV = list()
  
  # Repeated CV loop
  start.time = Sys.time()
  for(seed in seeds.RCV){
    perf.list = list()

    # k-fold CV loop
    folds = Create_folds(df.baseline, n.fold, seed) #ids of surv.data for each fold
    for(k.fold in 1:n.fold){
      
      # Print loop information
      print(paste0('Lmark: ',l,' Repetition: ',(seed+1-SEED),' Fold: ',k.fold))
      
      #training data and test data
      id.test = folds[[k.fold]]$ids.test
      df.train = df.baseline[!df.baseline$id%in%id.test, ]
      df.test = df.baseline[df.baseline$id%in%id.test, ]
      
      # Cox model
      fit.cox = coxph(Surv(time.to.event, event) ~., data = df.baseline[,-1], 
                      model = TRUE, x = TRUE)
      surv.prob = predictSurvProb(fit.cox, df.test, times=times)

      # Predict survival probabilities
      fail.prob = as.data.frame(1 - surv.prob)
      colnames(fail.prob) = times
      
      #-------------------------------------------------------------------------
      # PERFORMANCE INDICES
      # Brier score
      X = Score(as.list(fail.prob), data=df.test, formula=Surv(time.to.event,event)~1,
                metric=c('brier','auc'), times=times,exact = FALSE,conf.int = FALSE,cens.model = "cox",
                splitMethod = "none", B = 0, verbose = FALSE)
      
      auc = subset(X$AUC$score, model == times)$AUC
      brier = subset(X$Brier$score, model == times)$Brier
      length(auc) = length(times)
      length(brier) = length(times)
    
      # Concordance index
      X = concordance.index(x=fail.prob[,1], surv.time=df.test$time.to.event,
                            surv.event=df.test$event, method='noether')
      CI = X$c.index
      
      # Store C index to compute standard deviation
      matrix.C.index[seed+1-SEED,k.fold] = CI
      
      table = cbind(rep(l, length(times)), times, auc, brier, rep(CI,length(times)))
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
  
  #-------------------------------------------------------------------------------
  ## Save PERFORMANCE table
  model = 'staticCox'
  file.name = paste0(dataset,'_',model,'_',l*10,".RData")
  file.path = paste0('results/',dataset,'/',file.name)
  save(PERF, mean.perf.RCV, file = file.path)
  
  # Computing Time
  end.time = Sys.time()
  comp.time[which(lmark == l)] = as.numeric(difftime(end.time, start.time, units = 'secs'))
}

save(comp.time, file = paste0('ctime/', dataset, '_', model, '.RData'))

