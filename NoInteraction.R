#######################################################
#######################################################
### R Code for 1-1 matching analyses, no interaction term
###
###
### Contents:
###   I.   Library, data sourcing, working directory
###   II.  Sourced functions
###   III. Example syntax to run code
###
#######################################################
#######################################################


#######################################################
#######################################################
###   I. Library and data sourcing

library("survival")

setwd("~/")

mydata = readRDS("yourdata.rds")

#######################################################
#######################################################

#######################################################
#######################################################
###   II. Functions
###      a. expit
###      b. int_df
###      c. int_df_no_cov
###      d. fc_df
###      e. fc_df_no_cov
###      f. twostage_df
###      g. main_x

expit = function(x){return(exp(x)/(exp(x)+1))}

int_df = function(data, X, S, H, G, W, Y, strata, nstud, nref, covar=NA){
  
  #### 1. Create matrix to store output
  output_int = matrix(NA, ncol=3, nrow=1)
  colnames(output_int) = c("Point estimate of betax","Estimated RR", "Estimated variance of betahatx")
  
  #### 2. Rename variables names given to function to ensure consistency
  data$X = data[[X]]
  data$S = data[[S]]
  data$H = data[[H]]
  data$G = data[[G]]
  data$W = data[[W]]
  data$Y = data[[Y]]
  data$strata = data[[strata]]
  
  #### 3. Compute other useful quantities
  nz = length(covar)
  np = dim(data)[1]/2 # total number of strata
  n1 = c(rep(NA,nstud)) # number of matched pairs in each study
  nc = nstud - nref # number of studies needing calibration
  for(k in 1:nstud){
    n1[k] = sum(data$S==k)/2
  }
  
  #### 4. Complete calibration studies and add appropriate ahat, bhat to data frame
  a_hat = c(rep(0,nref),rep(NA,nc))
  b_hat = c(rep(1,nref),rep(NA,nc))
  
  for(k in (nref+1):nstud){
    cal_data_s  = subset(data, S==k & H==1)
    fit         = lm(X~W, data=cal_data_s)
    a_hat[k]    = fit$coefficients[1]
    b_hat[k]    = fit$coefficients[2]
  }
  
  data$a_hat = a_hat[data$S] ## adding a_hat and b_hat to the dataframe
  data$b_hat = b_hat[data$S]
  
  #### 5. Create xhat_fc variable- use H==2 to indicate when using X ref lab
  data$xhat_fc  = ifelse(data$H==2, data$X, data$a_hat + data$b_hat*data$W)
  data$xhat_int = ifelse(data$H==0, data$xhat_fc, data$X) 
  
  
  #### 6. Obtain point estimate from standard logistic regression
  formula        = as.formula(paste("Y~xhat_int+strata(strata)",paste(covar,collapse ="+"),sep="+"))
  int_fit        = clogit(formula, data=data)
  beta_hat       = int_fit$coefficients[1] # betahat_x value
  betazhat       = int_fit$coefficients[-1] # extract betazhat values
  output_int[1]   = beta_hat
  output_int[2]   = exp(beta_hat)
  
  #### 6.b. Compute Xd for each observation, Zd, and BX (linear combo)
  cases         = subset(data, Y==1)
  controls      = subset(data, Y==0)
  xd_int        = cases$xhat_int - controls$xhat_int
  data$xd_int   = rep(xd_int, each=2) # put in data frame for later reference in variance computation
  wd            = cases$W - controls$W
  data$wd       = rep(wd,each=2)
  
  ## Z differences
  col_nums        = which(names(data) %in% covar)
  Zvec_cases      = as.matrix(cases[,col_nums])
  Zvec_controls   = as.matrix(controls[,col_nums])
  Zd              = Zvec_cases - Zvec_controls
  BZ              = Zd %*% betazhat
  data$BXd        = rep(beta_hat*xd_int + BZ, each=2)
  
  ### Zd columns; update naming, add to data frame
  covdiffs           = Zd[rep(1:nrow(Zd), each = 2), ]
  covard             = paste(covar,"d",sep = "")
  colnames(covdiffs) = covard
  data               = cbind(data,covdiffs)
  
  #### 7. Compute variance; prepare matrices 
  dim_sand  = 2*nstud+1+nz
  A         = matrix(0, ncol = dim_sand, nrow = dim_sand)
  B         = matrix(0, ncol = dim_sand, nrow = dim_sand)
  
  ###################### A MATRIX #######################################
  
  #### 7.a.i: The upper block diagonals of A
  for(k in (nref+1):nstud){
    cal_data_s               = subset(data, (H==1 & S==k)) # Specific cal data with only controls
    A[2*k-1, 2*k-1]          = (1/np)*sum((cal_data_s$X-a_hat[k]-b_hat[k]*cal_data_s$W)^2)
    A[2*k-1,(2*k-1)+1]       = (1/np)*sum(cal_data_s$W*(cal_data_s$X-a_hat[k]-b_hat[k]*cal_data_s$W)^2)
    A[(2*k-1)+1,(2*k-1)]     = A[2*k-1, (2*k-1)+1]
    A[(2*k-1)+1,(2*k-1)+1]   = (1/np)*sum((cal_data_s$W^2)*(cal_data_s$X-a_hat[k]-b_hat[k]*cal_data_s$W)^2)
  }
  
  #### 7.a.ii. Upper right entries of A (A12)
  A12    = matrix(NA, ncol=(nz+1), nrow=(2*nstud))
  
  for(k in 1:nstud){
    
    ### Betax and a,b ###
    data_s_h     = subset(data, S==k & H==1) # Calibration controls
    Xd_s_h       = data_s_h$xd_int
    BXd_s_h      = data_s_h$BXd
    psi1         = data_s_h$X-a_hat[k]-b_hat[k]*data_s_h$W # estimating equ 1
    psi2         = (data_s_h$X-a_hat[k]-b_hat[k]*data_s_h$W)*data_s_h$W # estimating equ 2
    psiR         = Xd_s_h/(exp(BXd_s_h)+1) # last estimating equation
    A12[(2*k-1),1]   = (1/np)*sum(psi1*psiR)
    A12[(2*k),1]     = (1/np)*sum(psi2*psiR)
    
    ### Betaz and a,b ###
    ## Grab the Z covariate vector of interest; do computations with that particular Z
    for(i in 1:nz){ 
      col_num  = which(names(data_s_h) == covard[i])
      ZZ       = data_s_h[,col_num]
      psiZ     = ZZ/(1+exp(BXd_s_h))
      
      A12[(2*k-1),(1+i)]  = (1/np)*sum(psi1*psiZ)
      A12[(2*k),(1+i)]    = (1/np)*sum(psi2*psiZ) 
    }
    
  }
  
  A[1:(2*nstud), (2*nstud+1):dim_sand] = A12 # Fill appropriate piece of A matrix
  A21                                  = t(A12) # Reflect entries for A21
  A[(2*nstud+1):dim_sand,1:(2*nstud)]  = A21
  
  
  ## Now compute A22 entries
  A22    = matrix(NA, ncol=(nz+1), nrow=(nz+1))
  Xd     = data$xd_int[seq(1, nrow(data), 2)] # vector of all Xd in study
  BXd    = data$BXd[seq(1, nrow(data), 2)] # vector of all BXd in study
  
  ## Betax^2 entry
  psiX      = Xd/(exp(BXd)+1) 
  A22[1,1]  = (1/np)*sum(psiX^2)
  
  ### Computations in A22 involving betax and betaz, betaz squared
  
  for(i in 1:nz){ 
    col_num  = which(names(data) == covard[i])
    ZZ       = data[,col_num]
    ZZ       = ZZ[seq(1, length(ZZ), 2)]
    psiZ     = ZZ/(1+exp(BXd))
    
    A22[(1),(1+i)]  = (1/np)*sum(psiX*psiZ)
    A22[(1+i),(1)]  = (1/np)*sum(psiX*psiZ)
    A22[(1+i),(1+i)]    = (1/np)*sum(psiZ*psiZ) 
  }
  
  ### Computations involving cross terms of Z
  
  ## Interact the covariate psi terms (Z1 and Z2...) for all pairs
  
  pairs = combn(seq(1,nz,1),m=2)
  for(j in 1:choose(nz,2)){ #for each pairwise combination
    pair = pairs[,j]
    
    # Obtain appropriate matrix column vector based on pairwise selection
    col_num1 = which(names(data) == covard[pair[1]])
    ZZ1 = data[,col_num1]
    ZZ1 = ZZ1[seq(1, length(ZZ1), 2)] # get unique Zd obs
    psi_betaz1    = ZZ1/(1+exp(BXd))
    
    col_num2 = which(names(data) == covar[pair[2]])
    ZZ2 = data[,col_num2]
    ZZ2 = ZZ2[seq(1, length(ZZ2), 2)] # get unique Zd obs
    psi_betaz2    = ZZ2/(1+exp(BXd))
    
    #Fill matrix elements
    A22[(1+pair[2]),(1+pair[1])] = (1/np)*sum(psi_betaz1*psi_betaz2)
    A22[(1+pair[1]),(1+pair[2])] = (1/np)*sum(psi_betaz1*psi_betaz2)
  }
  
  ## Now place A22 in the main A matrix
  A[(2*nstud+1):dim_sand, (2*nstud+1):dim_sand] = A22
  
  ## Eliminate extra entries from A not associated with calibration studies
  A = A[(2*nref+1):dim_sand,(2*nref+1):dim_sand]
  
  #################### B Matrix ############################  
  
  ### 8.a.i. B11 entries: derivatives associated with calibration studies
  for(k in (nref+1):nstud){
    cal_data_s               = subset(data, S==k & H==1) # Specific cal data with only controls
    B[2*k-1,2*k-1]           = -1
    B[2*k-1,(2*k-1)+1]       = (-1/np)*sum(cal_data_s$W)
    B[(2*k-1)+1,(2*k-1)]     = B[2*k-1,(2*k-1)+1]
    B[(2*k-1)+1,(2*k-1)+1]   = (-1/np)*sum(cal_data_s$W^2)
  }
  
  ### 8.b. B12 entries: 
  
  B12 = matrix(0, ncol=(1+nz), nrow=(2*nstud))
  
  for(k in (nref+1):nstud){
    
    data_s_h    = subset(data, S==k & H==1) # Calibration controls (int)
    data_s_g    = subset(data, S==k & G==1) # Cases paired with calibration controls (int)
    data_s_e    = subset(data, S==k & H==0 & G==0) # remaining subjects not in cal study (full cal)
    
    W1_s_g      = data_s_g$W
    Xd_s_h      = data_s_h$xd_int  
    BXd_s_h     = data_s_h$BXd
    Wd_s_e      = data_s_e$wd[seq(1, nrow(data_s_e), 2)]
    BXd_s_e     = data_s_e$BXd[seq(1, nrow(data_s_e), 2)]
    
    # Odd entries (a and betax), internalized only
    frac11                = 1/(1+exp(BXd_s_h))
    frac22                = -beta_hat*Xd_s_h*exp(BXd_s_h)/((1+exp(BXd_s_h))^2)
    B12[(2*k-1), 1]       = (1/np)*(sum(frac11)+sum(frac22))
    
    # a and betaz, which exists for internalized only
    
    for(i in 1:nz){ 
      col_num  = which(names(data) == covard[i])
      ZZ       = data_s_h[,col_num]
      frac     = -beta_hat*ZZ*exp(BXd_s_h)/(1+exp(BXd_s_h))^2
      
      B12[(2*k-1),(1+i)]  = (1/np)*sum(frac)
    }
    
    # Even entries (b and betax)
    
    ## For those in calibration subset
    frac1            = W1_s_g/(1+exp(BXd_s_h))
    frac2            = (-beta_hat*Xd_s_h*W1_s_g*exp(BXd_s_h))/((1+exp(BXd_s_h))^2)
    
    ## For those outside calibration subset
    frac3            = Wd_s_e/(1+exp(BXd_s_e))
    frac4            = exp(BXd_s_e)*(-beta_hat*b_hat[k]*Wd_s_e^2)/((1+exp(BXd_s_e))^2)
    
    ## Sum all components
    B12[(2*k), 1]   = (1/np)*(sum(frac1)+sum(frac2)+sum(frac3)+sum(frac4))
    
    
    
    ### Betaz and b : full calibration and internalized pairs
    
    for(i in 1:nz){ 
      col_num  = which(names(data_s_h) == covard[i])
      ZZ       = data_s_h[,col_num]
      
      # internalized
      frac1     = -beta_hat*ZZ*W1_s_g*exp(BXd_s_h)/(1+exp(BXd_s_h))^2
      # full calibration
      frac2     = -beta_hat*ZZ*Wd_s_e*exp(BXd_s_e)/(1+exp(BXd_s_e))^2
      
      B12[(2*k),(1+i)]  = (1/np)*(sum(frac1)+sum(frac2))
    }
  }
  
  #### Replace B12 and its transpose in the B matrix
  B[1:(2*nstud),(2*nstud+1):dim_sand] = B12
  B[(2*nstud+1):dim_sand,1:(2*nstud)] = t(B12)
  
  ### B22 entry : Completed component-wise #######
  
  B22    = matrix(NA, ncol=(nz+1), nrow=(nz+1))
  Xd     = data$xd_int[seq(1, nrow(data), 2)] # vector of all Xd 
  BXd    = data$BXd[seq(1, nrow(data), 2)] # vector of all BXd
  
  ## Betax^2 entry
  frac      = -(Xd^2)*exp(BXd)/(1+exp(BXd))^2 
  B22[1,1]  = (1/np)*sum(frac)
  
  ### Computations in B22 involving betax and betaz, betaz squared
  
  for(i in 1:nz){ 
    col_num  = which(names(data) == covard[i])
    ZZ       = data[,col_num]
    ZZ       = ZZ[seq(1, length(ZZ), 2)]
    fracXZ   = -Xd*ZZ*exp(BXd)/(1+exp(BXd))^2
    fracZZ   = -(ZZ^2)*exp(BXd)/(1+exp(BXd))^2
    
    B22[(1),(1+i)]      = (1/np)*sum(fracXZ)
    B22[(1+i),(1)]      = (1/np)*sum(fracXZ)
    B22[(1+i),(1+i)]    = (1/np)*sum(fracZZ) 
  }
  
  ### Computations involving cross terms of Z
  
  ## Interact the covariate psi terms (Z1 and Z2...) for all pairs
  
  pairs = combn(seq(1,nz,1),m=2)
  for(j in 1:choose(nz,2)){ #for each pairwise combination
    pair = pairs[,j]
    
    # Obtain appropriate matrix column vector based on pairwise selection
    col_num1 = which(names(data) == covard[pair[1]])
    ZZ1 = data[,col_num1]
    ZZ1 = ZZ1[seq(1, length(ZZ1), 2)] # get unique Zd obs
    
    col_num2 = which(names(data) == covar[pair[2]])
    ZZ2 = data[,col_num2]
    ZZ2 = ZZ2[seq(1, length(ZZ2), 2)] # get unique Zd obs
    
    #Fill matrix elements
    fracZZ = -ZZ1*ZZ2*exp(BXd)/(1+exp(BXd))^2
    B22[(1+pair[2]),(1+pair[1])] = (1/np)*sum(fracZZ)
    B22[(1+pair[1]),(1+pair[2])] = (1/np)*sum(fracZZ)
  }
  
  ## Now place B22 in the main B matrix
  B[(2*nstud+1):dim_sand, (2*nstud+1):dim_sand] = B22
  
  ## Eliminate extra entries from B associated with noncalibration studies
  B = B[(2*nref+1):dim_sand,(2*nref+1):dim_sand]
  
  ### Compute sandwich variance estimator and place variance estimate in output
  V                = solve(B)%*%A%*%t(solve(B))
  output_int[3]    = (1/np)*V[(2*nc+1),(2*nc+1)]
  
  #### 99. Return appropriate output
  return(output_int)
  
}

int_df_no_cov = function(data, X, S, H, G, W, Y, strata, nstud, nref){
  
  #### 1. Create matrix to store output
  output_int = matrix(NA, ncol=3, nrow=1)
  colnames(output_int) = c("Point estimate of betax","Estimated RR", "Estimated variance of betahatx")
  
  #### 2. Rename variables names given to function to ensure consistency
  data$X = data[[X]]
  data$S = data[[S]]
  data$H = data[[H]]
  data$G = data[[G]]
  data$W = data[[W]]
  data$Y = data[[Y]]
  data$strata = data[[strata]]
  
  #### 3. Compute other useful quantities
  np = dim(data)[1]/2 # total number of strata
  n1 = c(rep(NA,nstud)) # number of matched pairs in each study
  nc = nstud - nref # number of studies needing calibration
  for(k in 1:nstud){
    n1[k] = sum(data$S==k)/2
  }
  
  #### 4. Complete calibration studies and add appropriate ahat, bhat to data frame
  a_hat = c(rep(0,nref),rep(NA,nc))
  b_hat = c(rep(1,nref),rep(NA,nc))
  ncal  = c(rep(NA,nstud))
  
  for(k in (nref+1):nstud){
    cal_data_s  = subset(data, S==k & H==1)
    ncal[k]     = nrow(cal_data_s)
    fit         = lm(X~W, data=cal_data_s)
    a_hat[k]    = fit$coefficients[1]
    b_hat[k]    = fit$coefficients[2]
  }
  
  data$a_hat = a_hat[data$S] ## adding a_hat and b_hat to the dataframe
  data$b_hat = b_hat[data$S]
  
  #### 5. Create xhat_fc and xhat_int variable- use H==2 to indicate when using X ref lab
  data$xhat_fc  = ifelse(data$H==2, data$X, data$a_hat + data$b_hat*data$W)
  data$xhat_int = ifelse(data$H==0, data$xhat_fc, data$X) 
  #data$BXd      = rep(beta_hat*xd_int, each=2)
  
  #### 6. Obtain point estimate from standard logistic regression
  formula        = as.formula(paste("Y~xhat_int+strata(strata)",sep="+"))
  int_fit        = clogit(formula, data=data)
  beta_hat       = int_fit$coefficients[1]
  output_int[1]   = beta_hat
  output_int[2]   = exp(beta_hat)
  
  #### 6.b. Compute Xd for each observation
  cases        = subset(data, Y==1)
  controls     = subset(data, Y==0)
  xd_int       = cases$xhat_int - controls$xhat_int
  data$xd_int  = rep(xd_int, each=2)
  wd           = cases$W - controls$W
  data$wd      = rep(wd,each=2)
  data$BXd     = rep(beta_hat*xd_int, each=2)
  
  #### 7. Compute variance; prepare matrices 
  dim_sand  = 2*nstud+1
  A         = matrix(0, ncol = dim_sand, nrow = dim_sand)
  B         = matrix(0, ncol = dim_sand, nrow = dim_sand)
  
  ###################### A MATRIX #######################################
  
  #### 7.a.i: The upper block diagonals of A
  for(k in (nref+1):nstud){
    cal_data_s               = subset(data, (H==1 & S==k)) # Specific cal data with only controls
    A[2*k-1, 2*k-1]          = (1/np)*sum((cal_data_s$X-a_hat[k]-b_hat[k]*cal_data_s$W)^2)
    A[2*k-1,(2*k-1)+1]       = (1/np)*sum(cal_data_s$W*(cal_data_s$X-a_hat[k]-b_hat[k]*cal_data_s$W)^2)
    A[(2*k-1)+1,(2*k-1)]     = A[2*k-1, (2*k-1)+1]
    A[(2*k-1)+1,(2*k-1)+1]   = (1/np)*sum((cal_data_s$W^2)*(cal_data_s$X-a_hat[k]-b_hat[k]*cal_data_s$W)^2)
  }
  
  #### 7.a.ii. Upper right entries of A (A12)
  A12    = matrix(NA, ncol=1, nrow=(2*nstud))
  
  for(k in 1:nstud){
    data_s_h     = subset(data, S==k & H==1) # Calibration controls
    Xd_s_h       = data_s_h$xd_int# vector of all Xd in study
    psi1         = data_s_h$X-a_hat[k]-b_hat[k]*data_s_h$W # estimating equ 1
    psi2         = (data_s_h$X-a_hat[k]-b_hat[k]*data_s_h$W)*data_s_h$W # estimating equ 2
    psiR         = Xd_s_h/(exp(beta_hat*Xd_s_h)+1) # last estimating equation
    A12[(2*k-1)] = (1/np)*sum(psi1*psiR)
    A12[2*k]     = (1/np)*sum(psi2*psiR)
  }
  
  A[1:(dim_sand-1), dim_sand] = A12 # Fill appropriate piece of A matrix
  A21                         = t(A12) # Reflect entries for A21
  A[dim_sand,1:(dim_sand-1)]  = A21
  
  ## Replace appropriate elements of A matrix
  A[1:(2*nstud), (2*nstud+1):dim_sand] = A12 
  A21                                  = t(A12)
  A[(2*nstud+1):dim_sand,1:(2*nstud)]  = A21
  
  
  ## Now compute A22 entry
  Xd       = data$xd_int[seq(1, nrow(data), 2)] # vector of all Xd overall
  A[dim_sand,dim_sand] = (1/np)*sum( (Xd/(exp(beta_hat*Xd)+1))^2 )
  
  
  ## Eliminate extra entries from A associated with calibration studies
  A = A[(2*nref+1):dim_sand,(2*nref+1):dim_sand]
  
  #################### B Matrix ############################  
  
  ### 8.a.i. B11 entries: derivatives associated with calibration studies
  for(k in (nref+1):nstud){
    cal_data_s               = subset(data, S==k & H==1) # Specific cal data with only controls
    B[2*k-1,2*k-1]           = -1
    B[2*k-1,(2*k-1)+1]       = (-1/np)*sum(cal_data_s$W)
    B[(2*k-1)+1,(2*k-1)]     = B[2*k-1,(2*k-1)+1]
    B[(2*k-1)+1,(2*k-1)+1]   = (-1/np)*sum(cal_data_s$W^2)
  }
  
  ### 8.b. B21 and B12 entries: place directly in B matrix
  
  B12 = matrix(NA, nrow=(2*nstud),ncol=1)
  
  for(k in (nref+1):nstud){
    data_s_h    = subset(data, S==k & H==1) # Calibration controls (int)
    data_s_g    = subset(data, S==k & G==1) # Cases paired with calibration controls (int)
    if(ncal[k]!=n1[k]){data_s_e = subset(data, S==k & H==0 & G==0)} # remaining subjects not in cal study (full cal)
    
    W1_s_g      = data_s_g$W
    Xd_s_h      = data_s_h$xd_int  
    BXd_s_h     = data_s_h$BXd
    if(ncal[k]!=n1[k]){ # for those studies having a FC component
      Wd_s_e      = data_s_e$wd[seq(1, nrow(data_s_e), 2)]
      BXd_s_e     = data_s_e$BXd[seq(1, nrow(data_s_e), 2)]
    }
    
    # Odd entries (a and betax), internalized only
    frac11                = 1/(1+exp(BXd_s_h))
    frac22                = -beta_hat*Xd_s_h*exp(BXd_s_h)/((1+exp(BXd_s_h))^2)
    B12[(2*k-1), 1]         = (1/np)*(sum(frac11)+sum(frac22))
    
    # Even entries (b and betax)
    
    ## For those in calibration subset
    frac1            = W1_s_g/(1+exp(BXd_s_h))
    frac2            = (-beta_hat*Xd_s_h*W1_s_g*exp(BXd_s_h))/((1+exp(BXd_s_h))^2)
    
    ## For those outside calibration subset
    frac3            = Wd_s_e/(1+exp(BXd_s_e))
    frac4            = exp(BXd_s_e)*(-beta_hat*b_hat[k]*Wd_s_e^2)/((1+exp(BXd_s_e))^2)
    
    ## Sum all components
    B12[(2*k), 1]   = (1/np)*(sum(frac1)+sum(frac2)+sum(frac3)+sum(frac4))
    if(ncal[k]==n1[k]){  B12[(2*k), 1] = (1/np)*(sum(frac1)+sum(frac2))}
    
  }
  
  ### Replace B12 and B21 in main B matrix
  B[1:(2*nstud), dim_sand] = B12
  B[dim_sand, 1:(2*nstud)] = t(B12)
  
  ### Single corner entry for betax^2 term
  
  B[dim_sand,dim_sand] = (1/np)*sum( -(Xd^2)*exp(beta_hat*Xd) / (exp(beta_hat*Xd)+1)^2)
  
  ## Eliminate extra entries from B associated with calibration studies
  B = B[(2*nref+1):dim_sand,(2*nref+1):dim_sand]
  
  ### Compute sandwich variance estimator and place variance estimate in output
  V                = solve(B)%*%A%*%t(solve(B))
  output_int[3]    = (1/np)*V[(2*nc+1),(2*nc+1)]
  
  #### 99. Return appropriate output
  return(output_int)
  
}

fc_df = function(data, X, S, H, W, Y, strata, nstud, nref, covar=NA){
  
  #### 1. Create matrix to store output
  output_fc = matrix(NA, ncol=3, nrow=1)
  colnames(output_fc) = c("Point estimate of betax","Estimated RR", "Estimated variance of betahatx")
  
  #### 2. Rename variables names given to function to ensure consistency
  data$X = data[[X]]
  data$S = data[[S]]
  data$H = data[[H]]
  data$W = data[[W]]
  data$Y = data[[Y]]
  data$strata = data[[strata]]
  
  #### 3. Compute other useful quantities
  nz = length(covar)
  np = dim(data)[1]/2 # total number of strata
  n1 = c(rep(NA,nstud)) # number of matched pairs in each study
  nc = nstud - nref # number of studies needing calibration
  for(k in 1:nstud){
    n1[k] = sum(data$S==k)/2
  }
  
  #### 4. Complete calibration studies and add appropriate ahat, bhat to data frame
  a_hat = c(rep(0,nref),rep(NA,nc))
  b_hat = c(rep(1,nref),rep(NA,nc))
  
  for(k in (nref+1):nstud){
    cal_data_s  = subset(data, S==k & H==1)
    fit         = lm(X~W, data=cal_data_s)
    a_hat[k]    = fit$coefficients[1]
    b_hat[k]    = fit$coefficients[2]
  }
  
  data$a_hat = a_hat[data$S] ## adding a_hat and b_hat to the dataframe
  data$b_hat = b_hat[data$S]
  
  #### 5. Create xhat_fc variable- use H==2 to indicate when using X ref lab
  data$xhat_fc = ifelse(data$H==2, data$X, data$a_hat + data$b_hat*data$W)
  
  #### 6. Obtain point estimate from standard logistic regression
  formula        = as.formula(paste("Y~xhat_fc+strata(strata)",paste(covar,collapse ="+"),sep="+"))
  fc_fit         = clogit(formula, data=data)
  beta_hat       = fc_fit$coefficients[1] # betahat_x value
  betazhat       = fc_fit$coefficients[-1] # extract betazhat values
  output_fc[1]   = beta_hat
  output_fc[2]   = exp(beta_hat)
  
  #### 6.b. Compute Xd for each observation, Zd, and BX (linear combo)
  cases        = subset(data, Y==1)
  controls     = subset(data, Y==0)
  xd_fc        = cases$xhat_fc - controls$xhat_fc
  data$xd_fc   = rep(xd_fc, each=2) # put in data frame for later reference in variance computation
  wd           = cases$W - controls$W
  data$wd      = rep(wd,each=2)
  
  ## Z differences
  col_nums        = which(names(data) %in% covar)
  Zvec_cases      = as.matrix(cases[,col_nums])
  Zvec_controls   = as.matrix(controls[,col_nums])
  Zd              = Zvec_cases - Zvec_controls
  BZ              = Zd %*% betazhat
  data$BXd        = rep(beta_hat*xd_fc + BZ, each=2)
  
  ### Zd columns; update naming, add to data frame
  covdiffs        = Zd[rep(1:nrow(Zd), each = 2), ]
  covard          = paste(covar,"d",sep = "")
  colnames(covdiffs) = covard
  data            = cbind(data,covdiffs)
  
  #### 7. Compute variance; prepare matrices 
  dim_sand  = 2*nstud+1+nz
  A         = matrix(0, ncol = dim_sand, nrow = dim_sand)
  B         = matrix(0, ncol = dim_sand, nrow = dim_sand)
  
  ###################### A MATRIX #######################################
  
  #### 7.a.i: The upper block diagonals of A
  for(k in (nref+1):nstud){
    cal_data_s               = subset(data, (H==1 & S==k)) # Specific cal data with only controls
    A[2*k-1, 2*k-1]          = (1/np)*sum((cal_data_s$X-a_hat[k]-b_hat[k]*cal_data_s$W)^2)
    A[2*k-1,(2*k-1)+1]       = (1/np)*sum(cal_data_s$W*(cal_data_s$X-a_hat[k]-b_hat[k]*cal_data_s$W)^2)
    A[(2*k-1)+1,(2*k-1)]     = A[2*k-1, (2*k-1)+1]
    A[(2*k-1)+1,(2*k-1)+1]   = (1/np)*sum((cal_data_s$W^2)*(cal_data_s$X-a_hat[k]-b_hat[k]*cal_data_s$W)^2)
  }
  
  #### 7.a.ii. Upper right entries of A (A12)
  A12    = matrix(NA, ncol=(nz+1), nrow=(2*nstud))
  
  for(k in 1:nstud){
    
    ### Betax and a,b ###
    data_s_h     = subset(data, S==k & H==1) # Calibration controls
    Xd_s_h       = data_s_h$xd_fc
    BXd_s_h      = data_s_h$BXd
    psi1         = data_s_h$X-a_hat[k]-b_hat[k]*data_s_h$W # estimating equ 1
    psi2         = (data_s_h$X-a_hat[k]-b_hat[k]*data_s_h$W)*data_s_h$W # estimating equ 2
    psiR         = Xd_s_h/(exp(BXd_s_h)+1) # last estimating equation
    A12[(2*k-1),1]   = (1/np)*sum(psi1*psiR)
    A12[(2*k),1]     = (1/np)*sum(psi2*psiR)
    
    ### Betaz and a,b ###
    ## Grab the Z covariate vector of interest; do computations with that particular Z
    for(i in 1:nz){ 
      col_num  = which(names(data_s_h) == covard[i])
      ZZ       = data_s_h[,col_num]
      psiZ     = ZZ/(1+exp(BXd_s_h))
      
      A12[(2*k-1),(1+i)]  = (1/np)*sum(psi1*psiZ)
      A12[(2*k),(1+i)]    = (1/np)*sum(psi2*psiZ) 
    }
    
  }
  
  A[1:(2*nstud), (2*nstud+1):dim_sand] = A12 # Fill appropriate piece of A matrix
  A21                                  = t(A12) # Reflect entries for A21
  A[(2*nstud+1):dim_sand,1:(2*nstud)]  = A21
  
  
  ## Now compute A22 entries
  A22    = matrix(NA, ncol=(nz+1), nrow=(nz+1))
  Xd     = data$xd_fc[seq(1, nrow(data), 2)] # vector of all Xd in study
  BXd    = data$BXd[seq(1, nrow(data), 2)] # vector of all BXd in study
  
  ## Betax^2 entry
  psiX      = Xd/(exp(BXd)+1) 
  A22[1,1]  = (1/np)*sum(psiX^2)
  
  ### Computations in A22 involving betax and betaz, betaz squared
  
  for(i in 1:nz){ 
    col_num  = which(names(data) == covard[i])
    ZZ       = data[,col_num]
    ZZ       = ZZ[seq(1, length(ZZ), 2)]
    psiZ     = ZZ/(1+exp(BXd))
    
    A22[(1),(1+i)]  = (1/np)*sum(psiX*psiZ)
    A22[(1+i),(1)]  = (1/np)*sum(psiX*psiZ)
    A22[(1+i),(1+i)]    = (1/np)*sum(psiZ*psiZ) 
  }
  
  ### Computations involving cross terms of Z
  
  ## Interact the covariate psi terms (Z1 and Z2...) for all pairs
  
  pairs = combn(seq(1,nz,1),m=2)
  for(j in 1:choose(nz,2)){ #for each pairwise combination
    pair = pairs[,j]
    
    # Obtain appropriate matrix column vector based on pairwise selection
    col_num1 = which(names(data) == covard[pair[1]])
    ZZ1 = data[,col_num1]
    ZZ1 = ZZ1[seq(1, length(ZZ1), 2)] # get unique Zd obs
    psi_betaz1    = ZZ1/(1+exp(BXd))
    
    col_num2 = which(names(data) == covar[pair[2]])
    ZZ2 = data[,col_num2]
    ZZ2 = ZZ2[seq(1, length(ZZ2), 2)] # get unique Zd obs
    psi_betaz2    = ZZ2/(1+exp(BXd))
    
    #Fill matrix elements
    A22[(1+pair[2]),(1+pair[1])] = (1/np)*sum(psi_betaz1*psi_betaz2)
    A22[(1+pair[1]),(1+pair[2])] = (1/np)*sum(psi_betaz1*psi_betaz2)
  }
  
  ## Now place A22 in the main A matrix
  A[(2*nstud+1):dim_sand, (2*nstud+1):dim_sand] = A22
  
  ## Eliminate extra entries from A not associated with calibration studies
  A = A[(2*nref+1):dim_sand,(2*nref+1):dim_sand]
  
  #################### B Matrix ############################  
  
  ### 8.a.i. B11 entries: derivatives associated with calibration studies
  for(k in (nref+1):nstud){
    cal_data_s               = subset(data, S==k & H==1) # Specific cal data with only controls
    B[2*k-1,2*k-1]           = -1
    B[2*k-1,(2*k-1)+1]       = (-1/np)*sum(cal_data_s$W)
    B[(2*k-1)+1,(2*k-1)]     = B[2*k-1,(2*k-1)+1]
    B[(2*k-1)+1,(2*k-1)+1]   = (-1/np)*sum(cal_data_s$W^2)
  }
  
  ### 8.b. B12 entries: 
  
  B12 = matrix(0, ncol=(1+nz), nrow=(2*nstud))
  
  for(k in (nref+1):nstud){
    
    ### Betax and b
    data_s           = subset(data, S==k) # Study specific values
    Wd_s             = data_s$wd[seq(1, nrow(data_s), 2)]
    Xd_s             = data_s$xd_fc[seq(1, nrow(data_s), 2)]
    BXd_s            = data_s$BXd[seq(1, nrow(data_s), 2)]
    frac1            = Wd_s/(1+exp(BXd_s))
    frac2            = exp(BXd_s)*(-beta_hat*b_hat[k]*Wd_s^2)/((1+exp(BXd_s))^2)
    
    B12[(2*k),1]  = (1/np)*(sum(frac1)+sum(frac2))
    
    ### Betaz and b
    for(i in 1:nz){ 
      col_num  = which(names(data) == covard[i])
      ZZ       = data_s[,col_num]
      ZZ       = ZZ[seq(1, length(ZZ), 2)]
      frac     = -beta_hat*ZZ*Wd_s*exp(BXd_s)/(1+exp(BXd_s))^2
      
      B12[(2*k),(1+i)]  = (1/np)*sum(frac)
    }
  }
  
  #### Replace B12 and its transpose in the B matrix
  B[1:(2*nstud),(2*nstud+1):dim_sand] = B12
  B[(2*nstud+1):dim_sand,1:(2*nstud)] = t(B12)
  
  ### B22 entry : Completed component-wise #######
  
  B22    = matrix(NA, ncol=(nz+1), nrow=(nz+1))
  Xd     = data$xd_fc[seq(1, nrow(data), 2)] # vector of all Xd 
  BXd    = data$BXd[seq(1, nrow(data), 2)] # vector of all BXd
  
  ## Betax^2 entry
  frac      = -(Xd^2)*exp(BXd)/(1+exp(BXd))^2 
  B22[1,1]  = (1/np)*sum(frac)
  
  ### Computations in B22 involving betax and betaz, betaz squared
  
  for(i in 1:nz){ 
    col_num  = which(names(data) == covard[i])
    ZZ       = data[,col_num]
    ZZ       = ZZ[seq(1, length(ZZ), 2)]
    fracXZ   = -Xd*ZZ*exp(BXd)/(1+exp(BXd))^2
    fracZZ   = -(ZZ^2)*exp(BXd)/(1+exp(BXd))^2
    
    B22[(1),(1+i)]      = (1/np)*sum(fracXZ)
    B22[(1+i),(1)]      = (1/np)*sum(fracXZ)
    B22[(1+i),(1+i)]    = (1/np)*sum(fracZZ) 
  }
  
  ### Computations involving cross terms of Z
  
  ## Interact the covariate psi terms (Z1 and Z2...) for all pairs
  
  pairs = combn(seq(1,nz,1),m=2)
  for(j in 1:choose(nz,2)){ #for each pairwise combination
    pair = pairs[,j]
    
    # Obtain appropriate matrix column vector based on pairwise selection
    col_num1 = which(names(data) == covard[pair[1]])
    ZZ1 = data[,col_num1]
    ZZ1 = ZZ1[seq(1, length(ZZ1), 2)] # get unique Zd obs
    
    col_num2 = which(names(data) == covar[pair[2]])
    ZZ2 = data[,col_num2]
    ZZ2 = ZZ2[seq(1, length(ZZ2), 2)] # get unique Zd obs
    
    #Fill matrix elements
    fracZZ = -ZZ1*ZZ2*exp(BXd)/(1+exp(BXd))^2
    B22[(1+pair[2]),(1+pair[1])] = (1/np)*sum(fracZZ)
    B22[(1+pair[1]),(1+pair[2])] = (1/np)*sum(fracZZ)
  }
  
  ## Now place B22 in the main B matrix
  B[(2*nstud+1):dim_sand, (2*nstud+1):dim_sand] = B22
  
  ## Eliminate extra entries from B associated with noncalibration studies
  B = B[(2*nref+1):dim_sand,(2*nref+1):dim_sand]
  
  ### Compute sandwich variance estimator and place variance estimate in output
  V                = solve(B)%*%A%*%t(solve(B))
  output_fc[3]    = (1/np)*V[(2*nc+1),(2*nc+1)]
  
  #### 99. Return appropriate output
  return(output_fc)
  
}

fc_df_no_cov = function(data, X, S, H, W, Y, strata, nstud, nref){
  
  #### 1. Create matrix to store output
  output_fc = matrix(NA, ncol=3, nrow=1)
  colnames(output_fc) = c("Point estimate of betax","Estimated RR", "Estimated variance of betahatx")
  
  #### 2. Rename variables names given to function to ensure consistency
  data$X = data[[X]]
  data$S = data[[S]]
  data$H = data[[H]]
  data$W = data[[W]]
  data$Y = data[[Y]]
  data$strata = data[[strata]]
  
  #### 3. Compute other useful quantities
  np = dim(data)[1]/2 # total number of strata
  n1 = c(rep(NA,nstud)) # number of matched pairs in each study
  nc = nstud - nref # number of studies needing calibration
  for(k in 1:nstud){
    n1[k] = sum(data$S==k)/2
  }
  
  #### 4. Complete calibration studies and add appropriate ahat, bhat to data frame
  a_hat = c(rep(0,nref),rep(NA,nc))
  b_hat = c(rep(1,nref),rep(NA,nc))
  
  for(k in (nref+1):nstud){
    cal_data_s  = subset(data, S==k & H==1)
    fit         = lm(X~W, data=cal_data_s)
    a_hat[k]    = fit$coefficients[1]
    b_hat[k]    = fit$coefficients[2]
  }
  
  data$a_hat = a_hat[data$S] ## adding a_hat and b_hat to the dataframe
  data$b_hat = b_hat[data$S]
  
  #### 5. Create xhat_fc variable- use H==2 to indicate when using X ref lab
  data$xhat_fc = ifelse(data$H==2, data$X, data$a_hat + data$b_hat*data$W)
  
  #### 6. Obtain point estimate from standard logistic regression
  formula        = as.formula(paste("Y~xhat_fc+strata(strata)",sep="+"))
  fc_fit         = clogit(formula, data=data)
  beta_hat       = fc_fit$coefficients[1]
  output_fc[1]   = beta_hat
  output_fc[2]   = exp(beta_hat)
  
  #### 6.b. Compute Xd for each observation
  cases        = subset(data, Y==1)
  controls     = subset(data, Y==0)
  xd_fc        = cases$xhat_fc - controls$xhat_fc
  data$xd_fc   = rep(xd_fc, each=2)
  wd           = cases$W - controls$W
  data$wd      = rep(wd,each=2)
  
  #### 7. Compute variance; prepare matrices 
  dim_sand  = 2*nstud+1
  A         = matrix(0, ncol = dim_sand, nrow = dim_sand)
  B         = matrix(0, ncol = dim_sand, nrow = dim_sand)
  
  ###################### A MATRIX #######################################
  
  #### 7.a.i: The upper block diagonals of A
  for(k in (nref+1):nstud){
    cal_data_s               = subset(data, (H==1 & S==k)) # Specific cal data with only controls
    A[2*k-1, 2*k-1]          = (1/np)*sum((cal_data_s$X-a_hat[k]-b_hat[k]*cal_data_s$W)^2)
    A[2*k-1,(2*k-1)+1]       = (1/np)*sum(cal_data_s$W*(cal_data_s$X-a_hat[k]-b_hat[k]*cal_data_s$W)^2)
    A[(2*k-1)+1,(2*k-1)]     = A[2*k-1, (2*k-1)+1]
    A[(2*k-1)+1,(2*k-1)+1]   = (1/np)*sum((cal_data_s$W^2)*(cal_data_s$X-a_hat[k]-b_hat[k]*cal_data_s$W)^2)
  }
  
  #### 7.a.ii. Upper right entries of A (A12)
  A12    = matrix(NA, ncol=1, nrow=(2*nstud))
  
  for(k in 1:nstud){
    data_s_h     = subset(data, S==k & H==1) # Calibration controls
    Xd_s_h       = data_s_h$xd_fc# vector of all Xd in study
    psi1         = data_s_h$X-a_hat[k]-b_hat[k]*data_s_h$W # estimating equ 1
    psi2         = (data_s_h$X-a_hat[k]-b_hat[k]*data_s_h$W)*data_s_h$W # estimating equ 2
    psiR         = Xd_s_h/(exp(beta_hat*Xd_s_h)+1) # last estimating equation
    A12[(2*k-1)] = (1/np)*sum(psi1*psiR)
    A12[2*k]     = (1/np)*sum(psi2*psiR)
  }
  
  A[1:(dim_sand-1), dim_sand] = A12 # Fill appropriate piece of A matrix
  A21                         = t(A12) # Reflect entries for A21
  A[dim_sand,1:(dim_sand-1)]  = A21
  
  ## Replace appropriate elements of A matrix
  A[1:(2*nstud), (2*nstud+1):dim_sand] = A12 
  A21                                  = t(A12)
  A[(2*nstud+1):dim_sand,1:(2*nstud)]  = A21
  
  
  ## Now compute A22 entry
  Xd       = data$xd_fc[seq(1, nrow(data), 2)] # vector of all Xd overall
  A[dim_sand,dim_sand] = (1/np)*sum( (Xd/(exp(beta_hat*Xd)+1))^2 )
  
  
  ## Eliminate extra entries from A associated with calibration studies
  A = A[(2*nref+1):dim_sand,(2*nref+1):dim_sand]
  
  #################### B Matrix ############################  
  
  ### 8.a.i. B11 entries: derivatives associated with calibration studies
  for(k in (nref+1):nstud){
    cal_data_s               = subset(data, S==k & H==1) # Specific cal data with only controls
    B[2*k-1,2*k-1]           = -1
    B[2*k-1,(2*k-1)+1]       = (-1/np)*sum(cal_data_s$W)
    B[(2*k-1)+1,(2*k-1)]     = B[2*k-1,(2*k-1)+1]
    B[(2*k-1)+1,(2*k-1)+1]   = (-1/np)*sum(cal_data_s$W^2)
  }
  
  ### 8.b. B21 and B12 entries: place directly in B matrix
  
  for(k in (nref+1):nstud){
    data_s           = subset(data, S==k) # Study specific values
    Wd_s             = data_s$wd[seq(1, nrow(data_s), 2)]
    Xd_s             = data_s$xd_fc[seq(1, nrow(data_s), 2)]
    frac1            = Wd_s/(1+exp(beta_hat*Xd_s))
    frac2            = exp(beta_hat*Xd_s)*(-beta_hat*b_hat[k]*Wd_s^2)/((1+exp(beta_hat*Xd_s))^2)
    
    
    B[2*k, 2*nstud+1]  = (1/np)*(sum(frac1)+sum(frac2))
    B[2*nstud+1, 2*k]  = (1/np)*(sum(frac1)+sum(frac2))
    
  }
  
  ### Single corner entry for betax
  
  B[dim_sand,dim_sand] = (1/np)*sum( -(Xd^2)*exp(beta_hat*Xd) / (exp(beta_hat*Xd)+1)^2)
  
  ## Eliminate extra entries from B associated with calibration studies
  B = B[(2*nref+1):dim_sand,(2*nref+1):dim_sand]
  
  ### Compute sandwich variance estimator and place variance estimate in output
  V                = solve(B)%*%A%*%t(solve(B))
  output_fc[3]    = (1/np)*V[(2*nc+1),(2*nc+1)]
  
  #### 99. Return appropriate output
  return(output_fc)
  
}

twostage_df = function(data, X, S, H, W, Y, G, strata, covar=NA, nref, nstud){
  
  #### 1.  Create output matrix: returns betahax and its estimated variance
  ####     wts: betahatw, var(betahatw), inversevar, weights, bhat, varbhat, betahatxi, varbetahatxi (not in that order)
  ####     ests:betahatx, var(betahatx), inversevar, weights
  output_ts           = matrix(NA, ncol=3, nrow=1)
  colnames(output_ts) = c("Point estimate of betax","Estimated RR", "Estimated variance of betahatx")
  wts                 = matrix(NA, nrow=(nstud-nref), ncol=8) # for local labs
  ests                = matrix(NA, nrow=nref, ncol=4) # for reference labs
  
  #### 2. Rename variables names given to function to ensure consistency
  data$X = data[[X]]
  data$S = data[[S]]
  data$H = data[[H]]
  data$W = data[[W]]
  data$Y = data[[Y]]
  data$strata = data[[strata]]
  
  #### 3. Get point and variance estimate of b_hat for each study in cal subset
  for(k in (nref+1):nstud){
    cal_data_s         = subset(data, S==k & H==1)
    fit                = lm(X~W, data=cal_data_s)
    wts[(k-nref),1]    = fit$coefficients[2]
    wts[(k-nref),2]    = vcov(fit)[2,2]
  }
  
  #### 4. Compute betahatw, var(betahatw) for each study in noncal subset
  for(k in (nref+1):nstud){
    #data_s        = subset(data, S==k & H==0 & G==0) #study specific, uncalibrated pairs
    data_s        = subset(data, S==k) #study specific pairs (not disregarding cal obs so we have enough data)
    if(is.na(covar[1])){ # if no covariates were specified
      fitW          = clogit(Y~W+strata(strata), data=data_s)
    }else{ # formula fit when covariates were specified
      formula       = as.formula(paste("Y~W+strata(strata)",paste(covar,collapse ="+"),sep="+"))
      fitW          = clogit(formula, data=data_s)
    }
    wts[(k-nref),3]      = fitW$coefficients[1] #betahatw
    wts[(k-nref),4]      = vcov(fitW)[1,1] #var betahatw
  }
  
  #### 5. Regression calibration to get corrected pt/var estimates
  wts[,5] = wts[,3]/wts[,1] # betahat xi
  wts[,6] = wts[,4]/(wts[,1]^2) + (wts[,2]*wts[,3]^2)/(wts[,1]^4) #var betahat xi
  wts[,7] = 1/wts[,6] # var(betahat) xi ^ (-1)
  
  #### 6. Compute inverse variance weights
  total = sum(wts[,7])
  wts[,8] = wts[,7]/total
  
  #### 7. Compute final estimate of betahatx and var betahatxi from LOCAL studies
  betahatx_loc    = wts[,8]%*%wts[,5]
  varbetahatx_loc = (wts[,8]^2)%*%wts[,6]
  
  #### 8. Compute betahatx, var(betahatx) for studies using reference lab, if any
  if(nref>0){
    for(k in 1:nref){
      data_s        = subset(data, S==k) #study specific
      if(is.na(covar[1])){ # if no covariates were specified
        fitX          = clogit(Y~X+strata(strata),data=data_s)
      }else{ # formula fit when covariates were specified
        formula       = as.formula(paste("Y~X+strata(strata)",paste(covar,collapse ="+"),sep="+"))
        fitX          = clogit(formula, data=data_s)
      }
      ests[k,1]      = fitX$coefficients[1] #betahatx
      ests[k,2]      = vcov(fitX)[1,1] #var betahatx
      ests[k,3]      = 1/ests[k,2] #inverse var weight
    }
    
    
    #### 9. Inverse var weights to get one estimate from studies using reference lab
    total    = sum(ests[ ,3])
    ests[,4] = ests[,3]/total
    
    #### 10. Compute final estimate of betahatx and var betahatxi from reference lab studies
    betahatx_ref    = ests[,4]%*%ests[,1]
    varbetahatx_ref = (ests[,4]^2)%*%ests[,2]
    
    #### 11. Combine ref and loc estimates to provide final estimates of point, RR, var
    total = 1/varbetahatx_ref + 1/varbetahatx_loc
    w_ref = (1/varbetahatx_ref)/total
    w_loc = (1/varbetahatx_loc)/total
    output_ts[1] = betahatx_loc*w_loc + betahatx_ref*w_ref
    output_ts[2] = exp(betahatx_loc*w_loc + betahatx_ref*w_ref)
    output_ts[3] = varbetahatx_ref*w_ref^2+varbetahatx_loc*w_loc^2
    
  }else{ # estimates if only local studies were used
    output_ts[1] = betahatx_loc
    output_ts[2] = exp(betahatx_loc)
    output_ts[3] = varbetahatx_loc
  } 
  
  
  #####
  return(output_ts)
  
}

main_x = function(data, X, S, H, W, Y, G, strata, nref, nstud, covar=NA){
  
  if(is.na(covar[1])){ # if there are no covariates
    
    output = matrix(NA, nrow=3, ncol=3)
    colnames(output) = c("betax estimate","RR estimate","Var(betax) estimate")
    rownames(output) = c("Internalized", "Full calibration", "Two stage")
    
    ## Computations
    output[1,] = int_df_no_cov(data=data, X=X, S=S, H=H, W=W, Y=Y, G=G, strata=strata, nref=nref, nstud=nstud)
    output[2,] = fc_df_no_cov(data=data, X=X, S=S, H=H, W=W, Y=Y, strata=strata, nref=nref, nstud=nstud)
    output[3,] = twostage_df(data=data, X=X, S=S, H=H, W=W, Y=Y, strata=strata, nref=nref, nstud=nstud)
    
    ## Rounding
    output[,1:2] = round(output[,1:2],8)
    output[,3]   = round(output[,3],8)
    
    ## Return output
    return(output)
    
  }else{ # Covariates given
    
    output = matrix(NA, nrow=6, ncol=3)
    colnames(output) = c("betax estimate","RR estimate","Var(betax) estimate")
    rownames(output) = c("Internalized, unadjusted", "Internalized, adjusted","Full calibration, unadjusted",
                         "Full calibration,adjusted","Two stage, unadjusted", "Two stage, adjusted")
    ## Computations
    output[1,] = int_df_no_cov(data=data, X=X, S=S, H=H, W=W, Y=Y, G=G, strata=strata, nref=nref, nstud=nstud)
    output[2,] = int_df(data=data, X=X, S=S, H=H, W=W, covar=covar, Y=Y, G=G, strata=strata, nref=nref, nstud=nstud)
    output[3,] = fc_df_no_cov(data=data, X=X, S=S, H=H, W=W, Y=Y, strata=strata, nref=nref, nstud=nstud)
    output[4,] = fc_df(data=data, X=X, S=S, H=H, W=W, strata=strata, covar=covar, Y=Y, nref=nref, nstud=nstud)
    output[5,] = twostage_df(data=data, X=X, S=S, H=H, W=W, strata=strata, Y=Y, nref=nref, nstud=nstud)
    output[6,] = twostage_df(data=data, X=X, S=S, H=H, W=W, covar=covar, strata=strata, Y=Y, nref=nref, nstud=nstud)
    
    ## Rounding
    output[,1:2] = round(output[,1:2],8)
    output[,3]   = round(output[,3],8)
    
    ## Return output
    return(output)
    
  }
}

#######################################################
#######################################################

#######################################################
#######################################################
###   III. Running functions on actual data

## Example syntax: main_x(data=mydata, X="ref",S="S",H="H",W="local",Y="Y",nref=2, nstud=5,
##                       covar=c("age","bmi"))

########### The function arguments of main_x are as follows:
## mydata: a dataframe of data that you have formatted in advance
## nstud:  is the total number of studies contributing to the analysis
## nref:   is the total number of studies that used the reference lab initially (if none, then nref=0)
## H:      is an indicator variable saying whether that observation was part of the calibration subset (0=no, 1=yes, 2=NA because
           the associated study used reference lab for all measurements and thus does not calibrate)
## S:      is the study number (numbering begins at 1 and increases incrementally by 1. Studies using the reference lab initially
           must be numbered before those that did use the reference lab)
## X:      is the reference lab measurement (NA if not available)
## W:      is the local laboratory measurement (NA if all measurements taken at reference lab initially)
## Y:      Binary outcome data (1/0 for yes/no)
## covar:  List of covariate names, syntax as shown above


#######################################################
#######################################################


