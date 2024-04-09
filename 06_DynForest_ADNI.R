ncores = 10

dataset = 'ADNI' #PBC2/ADNI/ROSMAP
run = 1 # run separately for each landmark due to high computing time
if (dataset=='PBC2') {
  load("data/pbc2.RData")
  lmark = c(2.5, 3, 3.5)[run]
  tmax = 8
}
if (dataset=='ADNI') {
  load("data/adni.RData")
  lmark = c(2:4)[run]
  tmax = 10
}
if (dataset=='ROSMAP') {
  load("data/rosmap.RData")
  lmark = c(2:6)[run]
  tmax = 15
}

long = long.transformed
long$id = as.numeric(long$id)
surv$id = as.numeric(surv$id)

long$t.var.lmm = long$time.fup 
long$t.var.pred = long$time.fup
fixed.cov.names = all.fixed.cov.names

if (dataset == 'PBC2'){
  long.cov.names = c("serChol","serBilir","albumin","alkaline","SGOT","platelets","prothrombin","histologic")
  timeVarModel = list(serChol = list(fixed = serChol ~ t.var.lmm, random = ~ t.var.lmm),
                      serBilir = list(fixed = serBilir ~ t.var.lmm, random = ~ t.var.lmm),
                      albumin = list(fixed = albumin ~ t.var.lmm, random = ~ t.var.lmm),
                      alkaline = list(fixed = alkaline ~ t.var.lmm, random = ~ t.var.lmm),
                      SGOT = list(fixed = SGOT ~ t.var.lmm, random = ~ t.var.lmm),
                      platelets = list(fixed = platelets ~ t.var.lmm, random = ~ t.var.lmm),
                      prothrombin = list(fixed = prothrombin ~ t.var.lmm, random = ~ t.var.lmm),
                      histologic = list(fixed = histologic ~ t.var.lmm, random = ~ t.var.lmm))
}

if (dataset == 'ADNI'){
  long.cov.names = c("ADAS11","ADASQ4","MMSE","RAVLT.immediate","RAVLT.learning","RAVLT.forgetting",
                     "LDELTOTAL","mPACCdigit","mPACCtrailsB","ADAS13","TRABSCOR","ICV",
                     "WholeBrain","Ventricles","Hippocampus","Entorhinal","Fusiform","MidTemp")
  timeVarModel = list(ADAS11 = list(fixed = ADAS11 ~ t.var.lmm, random = ~ t.var.lmm),
                      ADASQ4 = list(fixed = ADASQ4 ~ t.var.lmm, random = ~ t.var.lmm),
                      MMSE = list(fixed = MMSE ~ t.var.lmm, random = ~ t.var.lmm),
                      RAVLT.immediate = list(fixed = RAVLT.immediate ~ t.var.lmm, random = ~ t.var.lmm),
                      RAVLT.learning = list(fixed = RAVLT.learning ~ t.var.lmm, random = ~ t.var.lmm),
                      RAVLT.forgetting = list(fixed = RAVLT.forgetting ~ t.var.lmm, random = ~ t.var.lmm),
                      LDELTOTAL = list(fixed = LDELTOTAL ~ t.var.lmm, random = ~ 1|id),
                      mPACCdigit = list(fixed = mPACCdigit ~ t.var.lmm, random = ~ t.var.lmm),
                      mPACCtrailsB = list(fixed = mPACCtrailsB ~ t.var.lmm, random = ~ t.var.lmm),
                      ADAS13 = list(fixed = ADAS13 ~ t.var.lmm, random = ~ t.var.lmm),
                      TRABSCOR = list(fixed = TRABSCOR ~ t.var.lmm, random = ~ t.var.lmm),
                      ICV = list(fixed = ICV ~ t.var.lmm, random = ~ t.var.lmm),
                      WholeBrain = list(fixed = WholeBrain ~ t.var.lmm, random = ~ t.var.lmm),
                      Ventricles = list(fixed = Ventricles ~ t.var.lmm, random = ~ t.var.lmm),
                      Hippocampus = list(fixed = Hippocampus ~ t.var.lmm, random = ~ t.var.lmm),
                      Entorhinal = list(fixed = Entorhinal ~ t.var.lmm, random = ~ t.var.lmm),
                      Fusiform = list(fixed = Fusiform ~ t.var.lmm, random = ~ t.var.lmm),
                      MidTemp = list(fixed = MidTemp ~ t.var.lmm, random = ~ t.var.lmm))
}

if (dataset == 'ROSMAP'){
  long.cov.names = names.1percMiss
  timeVarModel = list(vasc_risks_sum = list(fixed = vasc_risks_sum ~ t.var.lmm, random = ~ t.var.lmm),
                      cogn_wo = list(fixed = cogn_wo ~ t.var.lmm, random = ~ t.var.lmm),
                      cogn_global = list(fixed = cogn_global ~ t.var.lmm, random = ~ t.var.lmm),
                      r_depres = list(fixed = r_depres ~ t.var.lmm, random = ~ t.var.lmm),
                      katzsum = list(fixed = katzsum ~ t.var.lmm, random = ~ t.var.lmm),
                      iadlsum = list(fixed = iadlsum ~ t.var.lmm, random = ~ t.var.lmm),
                      rosbsum = list(fixed = rosbsum ~ t.var.lmm, random = ~ t.var.lmm|id),
                      phys5itemsum = list(fixed = phys5itemsum ~ t.var.lmm, random = ~ t.var.lmm),
                      cesdsum = list(fixed = cesdsum ~ t.var.lmm, random = ~ t.var.lmm),
                      cogn_ep = list(fixed = cogn_ep ~ t.var.lmm, random = ~ t.var.lmm),
                      cogn_se = list(fixed = cogn_se ~ t.var.lmm, random = ~ t.var.lmm),
                      motor_dexterity = list(fixed = motor_dexterity ~ t.var.lmm, random = ~ t.var.lmm),
                      cogn_po = list(fixed = cogn_po ~ t.var.lmm, random = ~ t.var.lmm),
                      cogn_ps = list(fixed = cogn_ps ~ t.var.lmm, random = ~ t.var.lmm),
                      visions = list(fixed = vision ~ t.var.lmm, random = ~ t.var.lmm),
                      r_stroke = list(fixed = r_stroke ~ t.var.lmm, random = ~ t.var.lmm),
                      dbp_avg = list(fixed = dbp_avg ~ t.var.lmm, random = ~ t.var.lmm),
                      sbp_avg = list(fixed = sbp_avg ~ t.var.lmm, random = ~ t.var.lmm),
                      r_pd = list(fixed = r_pd ~ t.var.lmm, random = ~ t.var.lmm),
                      bmi = list(fixed = bmi ~ t.var.lmm, random = ~ t.var.lmm),
                      soc_net = list(fixed = soc_net ~ t.var.lmm, random = ~ t.var.lmm),
                      motor_gait = list(fixed = motor_gait ~ t.var.lmm, random = ~ t.var.lmm),
                      parkinsonism_tri = list(fixed = parkinsonism_tri ~ t.var.lmm, random = ~ t.var.lmm),
                      bradysc = list(fixed = bradysc ~ t.var.lmm, random = ~ t.var.lmm),
                      rigidsc = list(fixed = rigidsc ~ t.var.lmm, random = ~ t.var.lmm),
                      tremsc = list(fixed = tremsc ~ t.var.lmm, random = ~ t.var.lmm),
                      gaitsc = list(fixed = gaitsc ~ t.var.lmm, random = ~ t.var.lmm),
                      motor10 = list(fixed = motor10 ~ t.var.lmm, random = ~ t.var.lmm),
                      d_frailty = list(fixed = d_frailty ~ t.var.lmm, random = ~ t.var.lmm),
                      motor_handstreng = list(fixed = motor_handstreng ~ t.var.lmm, random = ~ t.var.lmm))
}

#-------------------------------------------------------------------------------
source('00_functions.R')
library(riskRegression)
library(survcomp)
library(magrittr)
library(dplyr)
library(DynForest)
library(data.table)

n.fold = ifelse(dataset == 'PBC2', 5, 10)
n.RCV = ifelse(dataset == 'PBC2', 20, 10)
SEED = 2023
seeds.RCV = SEED:(SEED + n.RCV-1)

# Prediction times
times = c((lmark+1):tmax)

# Hyperparameters
if (dataset == 'PBC2'){
  mtry = 4
  nodesize1 = 4
  minsplit1 = 4
  nodesize2 = 6
  minsplit2 = 4
  nodesize3 = 8
  minsplit3 = 4
  nodesize4 = 10
  minsplit4 = 4
} 
if (dataset == 'ADNI'){
  mtry = 6
  nodesize1 = 5
  minsplit1 = 5
  nodesize2 = 10
  minsplit2 = 5
  nodesize3 = 15
  minsplit3 = 5
  nodesize4 = 20
  minsplit4 = 5
}
if (dataset == 'ROSMAP'){
  mtry = 6
  nodesize1 = 15
  minsplit1 = 5
  nodesize2 = 20
  minsplit2 = 5
  nodesize3 = 25
  minsplit3 = 5
  nodesize4 = 30
  minsplit4 = 5
}

#-------------------------------------------------------------------------------
# Landmarking
#remove subjects who had event/were censored before landmark
surv = subset(surv, time.to.event >= lmark)
long = subset(long, id %in% surv$id)

# remove measurements after landmark
long = subset(long, time.fup <= lmark)

# number of subjects retained in the analysis:
print(paste0('Landmarking at ',lmark,', number of subjects retained in the analysis: ',nrow(surv)))

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
# Repeated k-folds cross validation
# List to store mean of performance measures for each repetition
q = 0
mean.perf.RCV = list()

matrix.C.index = matrix(NA,n.RCV,n.fold)

# Repeated CV loop
t0 = start.time = Sys.time()
count.fails = 0 # to count how many times the funcion dynForest returns an error
for(seed in seeds.RCV){
  perf.list = list()

  # k-fold CV loop
  folds = Create_folds(surv, n.fold, seed) #ids of surv.data for each fold
  for(k.fold in 1:n.fold){
    print(paste0('Lmark: ',lmark,' Repetition: ',(seed+1-SEED),' Fold: ',k.fold))
    
    #training data and test data
    id.test = folds[[k.fold]]$ids.test
    long.train = long[!long$id %in% id.test,]
    long.test = long[long$id %in% id.test,]
    surv.train = surv[!surv$id %in% id.test,]
    surv.test = surv[surv$id %in% id.test,]
    
    # Longitudinal predictors & fixed-time predictors
    timeData_train = long.train[,c("id","t.var.lmm",long.cov.names)]
    fixedData_train = surv.train[,c("id",fixed.cov.names)]
    timeData_pred = long.test[,c("id",'t.var.lmm','t.var.pred',long.cov.names)]
    fixedData_pred = surv.test[,c("id",fixed.cov.names)]
    
    # Grow the random forest
    Y = list(type = "surv", Y = surv.train[,c("id","time.to.event","event")])
    res_dyn = try(DynForest(timeData = timeData_train, fixedData = fixedData_train,
                        timeVar = "t.var.lmm", idVar = "id", timeVarModel = timeVarModel,
                        Y = Y, ntree = 200, mtry = mtry, nodesize = nodesize1, minsplit = minsplit1, 
                        seed = 1234, ncores = ncores))

    # If the function dynForest returns an error, increment count and go to the next loop
    if (class(res_dyn)=='try-error'){
      print(paste0('Failed to train the model! Trying again with nodesize=',nodesize2,' and minsplit=',minsplit2))
      res_dyn = try(DynForest(timeData = timeData_train, fixedData = fixedData_train,
                              timeVar = "t.var.lmm", idVar = "id", timeVarModel = timeVarModel,
                              Y = Y, ntree = 200, mtry = mtry, nodesize = nodesize2, minsplit = minsplit2, 
                              seed = 1234, ncores = ncores))

      if (class(res_dyn)=='try-error'){
        print(paste0('Failed to train the model AGAIN! Trying again with nodesize=',nodesize3,' and minsplit=',minsplit3))
        res_dyn = try(DynForest(timeData = timeData_train, fixedData = fixedData_train,
                                timeVar = "t.var.lmm", idVar = "id", timeVarModel = timeVarModel,
                                Y = Y, ntree = 200, mtry = mtry, nodesize = nodesize3, minsplit = minsplit3, 
                                seed = 1234, ncores = ncores))
        
        if (class(res_dyn)=='try-error'){
          print(paste0('Failed to train the model AGAIN! Trying again with nodesize=',nodesize4,' and minsplit=',minsplit4))
          res_dyn = try(DynForest(timeData = timeData_train, fixedData = fixedData_train,
                                  timeVar = "t.var.lmm", idVar = "id", timeVarModel = timeVarModel,
                                  Y = Y, ntree = 200, mtry = mtry, nodesize = nodesize4, minsplit = minsplit4, 
                                  seed = 1234, ncores = ncores))
        }
      }
    }
    
    # If the function DynForest does not return an error, keep going with the computations
    if (class(res_dyn)!='try-error'){
      
      # Prediction
      start.time = Sys.time()
      pred_dyn <- predict(object = res_dyn, 
                          timeData = timeData_pred, fixedData = fixedData_pred,
                          idVar = 'id', timeVar = 't.var.pred', t0 = lmark)
      end.time = Sys.time()
      print(paste0('Computation Time for prediction: ',round(end.time - start.time,2)))
      
      #estimate failure prob at times
      fail.preds = matrix(NA, nrow=nrow(fixedData_pred), length(times))
      ith = 1
      for(t in times){
        nth = length(which((pred_dyn$times - t) <= 0)) #select the closest time point from left
        fail.preds[,ith] = pred_dyn$pred_indiv[,nth]
        ith = ith +1
      }
      fail.preds = as.data.frame(fail.preds)
      colnames(fail.preds) = times
      
      # Brier Score and tdAUC with riskRegression package
      X = Score(as.list(fail.preds), data=surv.test, formula=Surv(time.to.event, event) ~1,
                metrics=c('brier','auc'),times=times,exact = FALSE,conf.int = FALSE,
                cens.model = "cox", splitMethod = "none", B = 0, verbose = FALSE)
      
      auc = subset(X$AUC$score, model == times)
      brier = subset(X$Brier$score, model == times)
      
      # Concordance index
      X = concordance.index(x=fail.preds[,1], surv.time=surv.test$time.to.event, surv.event=
                              surv.test$event, method='noether')
      CI = X$c.index
      
      # Store C index to compute standard deviation
      matrix.C.index[seed+1-SEED,k.fold] = CI
      
      # Performance indeces Table
      table = cbind(rep(lmark, length(times)), times, auc$AUC, brier$Brier, rep(CI, length(times)))
      colnames(table) = c('lmark', 't', 'AUC', 'BS', 'CI')
      
      # Store performance table
      perf.list[[k.fold]] = as.data.frame(table)
    }
    
    if (class(res_dyn)=='try-error'){
      print(paste0('Failed to train the model, go to the next fold :('))
      count.fails = count.fails + 1
    }
  }
  
  # Mean performance indices
  q = q+1
  if(length(perf.list) != 0) mean.perf.RCV[[q]] = rbindlist(perf.list)[,lapply(.SD,mean,na.rm=TRUE),list(lmark,t)]
}

# Mean performance indices
PERF = rbindlist(mean.perf.RCV)[,lapply(.SD,mean,na.rm=TRUE),list(lmark,t)]

# Total computation time for training
print(paste0('Number of times the function DynForest failed: ', count.fails))

#-------------------------------------------------------------------------------
## Save PERFORMANCE table
model = 'dynForest'
file.name = paste0(dataset,"_",model,"_",lmark*10,".RData")
file.path = paste0('results/',dataset,'/',file.name)
save(PERF, mean.perf.RCV, count.fails, file = file.path)

t1 = Sys.time()
comp.time = as.numeric(difftime(t1, t0, units='secs'))
save(comp.time, file = paste0('ctime/', dataset, '_', model,
                              '_', lmark*10, '.RData'))
