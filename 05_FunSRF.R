
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
source('00_functions4MFPCA.R')
library(MASS)
library(MFPCA)
library(randomForestSRC)
library(pec)
library(survival)
library(abind)
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
  Q = length(long.cov.names)
  
  # mtry value: set to sqrt(number of covs), rounded up
  # number of covs: ADNI 26, PBC2 11, ROSMAP 35
  mtry = ifelse(dataset == 'PBC2', 4, 6)

  #-------------------------------------------------------------------------------
  # Landmarking
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
  # From long to wide format
  nPat = dim(surv)[1]
  patID = surv$id
  
  # transfer longitudinal outcomes from long to wide
  multivar = array(NA, c(nPat, length(obstime), Q))
  for(i in 1:nPat){
    visits = which(obstime %in% (long$time.regular[long$id == patID[i]]))
    k.cov = 0
    for(cov in long.cov.names){
      k.cov = k.cov+1
      multivar[i,visits,k.cov] = long[long$id == patID[i],cov]
    }
  }
  
  #-------------------------------------------------------------------------------
  # Repeated cross validation
  
  # List to store mean of performance measures for each repetition
  q = 0
  mean.perf.RCV = list()
  
  count.fails = 0
  # Repeated CV loop
  start.time = Sys.time()
  for(seed in seeds.RCV){
    
    # To store performance measures 
    perf.list = list()

    # k-fold CV loop
    folds = Create_folds(surv, n.fold, seed) #ids of surv.data for each fold
    for(k.fold in 1:n.fold){
      
      # Print loop information
      print(paste0('Lmark: ',l,' Repetition: ',(seed+1-SEED),' Fold: ',k.fold))
      
      # Training data and test data
      id.test = folds[[k.fold]]$ids.test
      multivar.train = multivar[!surv$id %in% id.test, , ]
      multivar.test = multivar[surv$id %in% id.test,,]
      long.train = long[!long$id%in%id.test, ]
      long.test = long[long$id%in%id.test, ]
      surv.train = surv[!surv$id%in%id.test, ]
      surv.test = surv[surv$id%in%id.test, ]
      
      # Train the model ----------------------------------------------------------
      # Univariate FPCA via PACE
      Xi.train = L = phi.train = meanFun.train =  NULL
      pve = 0.90   # proportion of variance explained
      npc = NULL     # number of principal components - NOT USED IN THIS CASE
      argvals =  obstime/max(obstime) # scale the time domain to [0,1]
      
      res = try(for(p in 1:Q){
        tmp.ufpca = uPACE(multivar.train[,,p], argvals, nbasis=3, pve = pve, npc = npc)
        Xi.train = cbind(Xi.train, tmp.ufpca$scores) # FPC scores
        L = c(L, dim(tmp.ufpca$scores)[2])
        phi.train[[p]] = t(tmp.ufpca$functions@X) # FPC eigenfunctions
        meanFun.train[[p]] = tmp.ufpca$mu@X # estimated mean functions
      })
      
      if (class(res)!='try-error'){
        # Multivariate FPCA
        I = dim(surv.train)[1] # this value is needed in mFPCA
        mFPCA.train = mFPCA(Xi=Xi.train, phi=phi.train, p=Q, L=L )
        rho.train = mFPCA.train$rho  #MFPC scores
        pve = mFPCA.train$pve
        psi = mFPCA.train$psi
        Cms = mFPCA.train$Cms
        
        # Survival model - RSF
        colnames(rho.train) = paste0('rho.',(1:ncol(rho.train)))
        surv.train.temp = cbind(surv.train[,c('time.to.event','event',fixed.cov.names)], rho.train)
        rsf.fit = rfsrc(Surv(time.to.event, event)~., data=surv.train.temp, 
                        ntree=1000, seed=1234, mtry = mtry)
        
        # Prediction on test set ---------------------------------------------------
        # Univariate FPC
        Xi.test = NULL
        for(p in 1:Q){
          tmp.ufpca = uPACE(multivar.train[,,p], argvals, multivar.test[,,p], nbasis=3, pve=0.9)
          Xi.test = cbind(Xi.test, tmp.ufpca$scores) # dynamic FPC scores for test subjects 
        }
        
        # estimate MFPC scores for test subjects
        rho.test = mfpca.score(Xi.test, Cms)
        colnames(rho.test) = paste0('rho.',(1:ncol(rho.test)))
        surv.test.rho = cbind(surv.test[,c('time.to.event','event',fixed.cov.names)],rho.test)
        
        # predict longitudinal trajectories 
        long.pred = mfpca.pred(rho.test, meanFun.train, psi)
        
        # prediction for different times
        ith = 0
        surv.prob = NULL
        for(dt in times){
          ith = ith + 1
          surv.prob[[ith]] = cond.prob.pec(rsf.fit, surv.test.rho, l, dt)  # risk prediction
        }  
        
        # Brier Score and tdAUC
        fail.prob = 1 - as.data.frame(do.call(cbind, surv.prob))
        names(fail.prob) = times
        X = Score(as.list(fail.prob), data=surv.test, formula=Surv(time.to.event, event) ~1,
                  metrics=c('brier','auc'),times=times, exact = FALSE,conf.int = FALSE,
                  cens.model = "cox",splitMethod = "none", B = 0,
                  verbose = FALSE)
        auc = subset(X$AUC$score, model == times)[,-1]
        brier = subset(X$Brier$score, model == times)[,-1]
        
        # Concordance index
        X = concordance.index(x=(1-surv.prob[[1]]), surv.time=surv.test$time.to.event,
                              surv.event=surv.test$event, method='noether')
        CI = X$c.index
        
        # Store C index to compute standard deviation
        matrix.C.index[seed+1-SEED,k.fold] = CI
        
        # Performance indeces Table
        table = cbind(rep(l, length(times)), times, auc$AUC, brier$Brier, rep(CI, length(times)))
        colnames(table) = c('lmark', 't', 'AUC', 'BS', 'CI')
        
        # Store performance table
        perf.list[[k.fold]] = as.data.frame(table)
      }
      
      if (class(res)=='try-error'){
        print(paste0('Model has more coefficient than data, go to the next fold :('))
        count.fails = count.fails + 1
      }
    }
    
    # average over the CV run
    q = q+1
    if(length(perf.list) != 0) mean.perf.RCV[[q]] = rbindlist(perf.list)[,lapply(.SD,mean,na.rm=TRUE),list(lmark,t)]
  }
  
  # average over the RCV
  if(length(mean.perf.RCV) != 0){
    PERF = rbindlist(mean.perf.RCV)[,lapply(.SD,mean,na.rm=TRUE),list(lmark,t)]
    print(PERF)
  }
  
  print(paste0('Number of folds for which it was not possible to train the model: ',
               count.fails))
  
  #-------------------------------------------------------------------------------
  if(length(mean.perf.RCV) != 0){
    ## Save PERFORMANCE table
    model = 'FunRSF'
    file.name = paste0(dataset,"_",model,"_",l*10,".RData")
    file.path = paste0('results/',dataset,'/',file.name)
    save(PERF, mean.perf.RCV, count.fails, file = file.path)
  }
  if(length(mean.perf.RCV) == 0) print(paste0('For landmark ',l,' model can not be trained for any of the folds!'))
  
  # Computing Time
  end.time = Sys.time()
  comp.time[which(lmark == l)] = as.numeric(difftime(end.time, start.time, units = 'secs'))
}

save(comp.time, file = paste0('ctime/', dataset, '_', model, '.RData'))

