# Function to create k folds for CV

Create_folds <- function(surv, n_fold, seed) {
  
  set.seed(seed)
  
  folds <- vector(mode = "list", length = n_fold)
  
  fold_size <- round(1 / n_fold * length(surv$id), 0)
  
  ids.remaining <- as.numeric(surv$id) # Initialize
  for (i in 1:n_fold) {
    j <- n_fold - (i - 1) # Count down
    if (j > 1) { # Not last fold
      tmp <- surv %>%
        filter(id %in% ids.remaining) %>% # Excluded subjects from other split
        splitstackshape::stratified(., group = "event", size = 1 / j) 
      # size is set to be proportional to the number of observations per group
      # Store ids
      folds[[i]]$ids.test <- as.numeric(tmp$id)
      # Update remaining subjects for selection
      ids.remaining <- setdiff(ids.remaining, tmp$id)
    } else { # Last fold
      folds[[i]]$ids.test <- ids.remaining
    }
    folds[[i]]$seed <- seed
  }
  
  return (folds)
}
