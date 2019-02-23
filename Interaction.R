#######################################################
#######################################################
### R Code for 1-1 matching analyses, Interaction term
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
###      b. int_df_xv
###      c. int_df_no_cov_xv
###      d. fc_df_xv
###      e. fc_df_no_cov_xv
###      f. twostage_df_xv
###      g. main_xv

expit = function(x){return(exp(x)/(exp(x)+1))}


in_df_xv = function(data, X, S, H, W, Y, V, G, strata, covar=NA, nstud, nref){
  
  #### 1. Create matrix to store output
  output_in = matrix(NA, ncol=9, nrow=1)
  colnames(output_in) = c("betax","RRx", "VARx","betav","RRv", "VARv","betaxv","RRxv", "VARxv")
  
  #### 2. Rename variables names given to function to ensure consistency
  data$X = data[[X]]
  data$S = data[[S]]
  data$H = data[[H]]
  data$G = data[[G]]
  data$W = data[[W]]
  data$Y = data[[Y]]
  data$V = data[[V]]
  data$strata = data[[strata]]
  
  #### 2b. Sort data frame, in case it wasn't already
  data = data[with(data, order(strata, -Y)),]
  
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
  data$xhat_fc    = ifelse(data$H==2, data$X, data$a_hat + data$b_hat*data$W)
  data$xhat_in    = ifelse(data$H==0, data$xhat_fc, data$X) 
  data$xhat_in_v  = ifelse(data$H==0, data$xhat_fc*data$V, data$X*data$V)
  
  #### 6. Obtain point estimate from standard logistic regression
  formula            = as.formula(paste("Y~xhat_in+V+xhat_in_v+strata(strata)",paste(covar,collapse ="+"),sep="+"))
  in_fit             = clogit(formula, data=data)
  beta_hatx          = in_fit$coefficients[1]
  beta_hatv          = in_fit$coefficients[2]
  beta_hatxv         = in_fit$coefficients[3]
  betazhat           = in_fit$coefficients[-c(1,2,3)] # extract betazhat values
  output_in[c(1,2)]  = c(beta_hatx,exp(beta_hatx))
  output_in[c(4,5)]  = c(beta_hatv,exp(beta_hatv))
  output_in[c(7,8)]  = c(beta_hatxv,exp(beta_hatxv))
  
  #### 6.b. Compute necessary differences for sandwich computations
  
  ## Take differences for each observation and store in matrix
  
  data$WV      = data$W*data$V
  cases        = subset(data, Y==1)
  controls     = subset(data, Y==0)
  xd_in        = cases$xhat_in - controls$xhat_in
  vd           = cases$V - controls$V
  xvd_in       = cases$xhat_in_v - controls$xhat_in_v
  wd           = cases$W - controls$W
  wvd          = cases$WV - controls$WV
  
  data$xd_in   = rep(xd_in, each=2)
  data$vd      = rep(vd, each=2)
  data$xvd_in  = rep(xvd_in, each=2)
  data$wd      = rep(wd, each=2)
  data$wvd     = rep(wvd, each=2)
  
  ## Z differences and linear combination
  col_nums        = which(names(data) %in% covar)
  Zvec_cases      = as.matrix(cases[,col_nums])
  Zvec_controls   = as.matrix(controls[,col_nums])
  Zd              = Zvec_cases - Zvec_controls
  BZ              = Zd %*% betazhat
  data$BZ         = rep(BZ,each=2)
  data$BX         = beta_hatx*data$xd_in + beta_hatv*data$vd + beta_hatxv*data$xvd_in + data$BZ
  
  ### Zd columns; update naming, add to data frame
  covdiffs        = Zd[rep(1:nrow(Zd), each = 2), ]
  covard          = paste(covar,"d",sep = "")
  colnames(covdiffs) = covard
  data            = cbind(data,covdiffs)
  
  #### 7. Compute variance; prepare matrices 
  dim_sand  = 2*nstud+3+nz
  A         = matrix(0, ncol = dim_sand, nrow = dim_sand)
  B         = matrix(0, ncol = dim_sand, nrow = dim_sand)
  
  ###################### A MATRIX #######################################
  
  #### 4a: The upper block diagonals of A
  for(k in 1:nstud){
    cal_data_s               = subset(data, (H==1 & S==k)) # Specific cal data with only controls
    A[2*k-1, 2*k-1]          = (1/np)*sum((cal_data_s$X-a_hat[k]-b_hat[k]*cal_data_s$W)^2)
    A[2*k-1,(2*k-1)+1]       = (1/np)*sum(cal_data_s$W*(cal_data_s$X-a_hat[k]-b_hat[k]*cal_data_s$W)^2)
    A[(2*k-1)+1,(2*k-1)]     = A[2*k-1, (2*k-1)+1]
    A[(2*k-1)+1,(2*k-1)+1]   = (1/np)*sum((cal_data_s$W^2)*(cal_data_s$X-a_hat[k]-b_hat[k]*cal_data_s$W)^2)
  }
  
  #### 4b: Upper right entries of A (A12)
  A12    = matrix(NA, ncol=(3+nz), nrow=(2*nstud))
  
  for(k in 1:nstud){
    data_s_h     = subset(data, S==k & H==1) # Calibration controls
    Xd_s_h       = data_s_h$xd_in
    Vd_s_h       = data_s_h$vd
    XVd_s_h      = data_s_h$xvd_in
    BX_s_h       = data_s_h$BX
    psi1         = data_s_h$X-a_hat[k]-b_hat[k]*data_s_h$W # estimating equ 1
    psi2         = (data_s_h$X-a_hat[k]-b_hat[k]*data_s_h$W)*data_s_h$W # estimating equ 2
    
    psi_x         = Xd_s_h/(exp(BX_s_h)+1) # betax estimating equation
    psi_v         = Vd_s_h/(exp(BX_s_h)+1) # betav estimating equation
    psi_xv        = XVd_s_h/(exp(BX_s_h)+1) # betaxv estimating equation
    
    ## Betax and a,b
    A12[(2*k-1),1]   = (1/np)*sum(psi1*psi_x)
    A12[(2*k),1]     = (1/np)*sum(psi2*psi_x)
    
    ## Betav and a,b
    A12[(2*k-1),2]   = (1/np)*sum(psi1*psi_v)
    A12[(2*k),2]     = (1/np)*sum(psi2*psi_v)
    
    ## Betaxv and a,b
    A12[(2*k-1),3]   = (1/np)*sum(psi1*psi_xv)
    A12[(2*k),3]     = (1/np)*sum(psi2*psi_xv)
    
    ## Grab the Z covariate vector of interest; do computations with that particular Z
    ## Betaz and a,b
    for(i in 1:nz){ 
      col_num  = which(names(data_s_h) == covard[i])
      ZZ       = data_s_h[,col_num]
      psiZ     = ZZ/(1+exp(BX_s_h))
      
      A12[(2*k-1),(3+i)]  = (1/np)*sum(psi1*psiZ)
      A12[(2*k),(3+i)]    = (1/np)*sum(psi2*psiZ) 
    }
    
  }
  
  A[1:(2*nstud), (2*nstud+1):dim_sand] = A12 # Fill appropriate piece of A matrix
  A21                                  = t(A12) # Reflect entries for A21
  A[(2*nstud+1):dim_sand,1:(2*nstud)]  = A21
  
  #### 4c: Preparations for computations in A22
  A22      = matrix(NA,nrow=(3+nz),ncol=(3+nz))
  Xd       = data$xd_in[seq(1, nrow(data), 2)] 
  Vd       = data$vd[seq(1, nrow(data), 2)] 
  XVd      = data$xvd_in[seq(1, nrow(data), 2)] 
  BX       = data$BX[seq(1, nrow(data), 2)] 
  
  ## Squared terms of X,V,XV on diagonal of A22
  psiX     = Xd/(exp(BX)+1)
  psiV     = Vd/(exp(BX)+1)
  psiXV    = XVd/(exp(BX)+1)
  A22[1,1] = (1/np)*sum( psiX^2 )
  A22[2,2] = (1/np)*sum( psiV^2 )
  A22[3,3] = (1/np)*sum( psiXV^2 )
  
  ## Betax and betav of A22
  A22[1,2] = (1/np)*sum( psiX*psiV )
  A22[2,1] = A22[1,2]
  
  ## Betav and betaxv of A22
  A22[2,3] = (1/np)*sum( psiV*psiXV )
  A22[3,2] = A22[2,3]
  
  ## Betax and betaxv of A22
  A22[1,3] = (1/np)*sum( psiX*psiXV )
  A22[3,1] = A22[1,3]
  
  ######## Now consider A22 entries involving an estimating equation with betaz
  ## Select specific psiZ and compute products
  for(i in 1:nz){ 
    col_num  = which(names(data) == covard[i])
    ZZ       = data[,col_num]
    ZZ       = ZZ[seq(1, length(ZZ), 2)]
    psiZ     = ZZ/(1+exp(BX))
    
    A22[1,(3+i)]  = (1/np)*sum(psiX*psiZ)
    A22[(3+i),1]  = (1/np)*sum(psiX*psiZ)
    A22[2,(3+i)]  = (1/np)*sum(psiV*psiZ)
    A22[(3+i),2]  = (1/np)*sum(psiV*psiZ)
    A22[3,(3+i)]  = (1/np)*sum(psiXV*psiZ)
    A22[(3+i),3]  = (1/np)*sum(psiXV*psiZ)
    A22[(3+i),(3+i)]    = (1/np)*sum(psiZ*psiZ) 
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
    psi_betaz1    = ZZ1/(1+exp(BX))
    
    col_num2 = which(names(data) == covar[pair[2]])
    ZZ2 = data[,col_num2]
    ZZ2 = ZZ2[seq(1, length(ZZ2), 2)] # get unique Zd obs
    psi_betaz2    = ZZ2/(1+exp(BX))
    
    #Fill matrix elements
    A22[(3+pair[2]),(3+pair[1])] = (1/np)*sum(psi_betaz1*psi_betaz2)
    A22[(3+pair[1]),(3+pair[2])] = (1/np)*sum(psi_betaz1*psi_betaz2)
  }
  
  
  ## Replace A22 in main matrix of A
  A[(2*nstud+1):dim_sand,(2*nstud+1):dim_sand ] = A22
  
  
  ## Eliminate extra entries from A associated with calibration studies
  A = A[(2*nref+1):dim_sand,(2*nref+1):dim_sand]
  
  ###################### B MATRIX #######################################
  
  #### 5a: Block diagonal entries of B associated with the calibration studies.
  
  for(k in 1:nstud){
    cal_data_s               = subset(data, S==k & H==1) # Specific cal data with only controls
    B[2*k-1,2*k-1]           = -1
    B[2*k-1,(2*k-1)+1]       = (-1/np)*sum(cal_data_s$W)
    B[(2*k-1)+1,(2*k-1)]     = B[2*k-1,(2*k-1)+1]
    B[(2*k-1)+1,(2*k-1)+1]   = (-1/np)*sum(cal_data_s$W^2)
  }
  
  ### 5b: Off diagonal even entries of B 
  B12 = matrix(NA, ncol=(3+nz), nrow=(2*nstud))
  
  for(k in 1:nstud){
    data_s      = subset(data, S==k & H==0 & G==0) #FC-esque component
    Xd_s        = data_s$xd_in[seq(1, nrow(data_s), 2)]
    Wd_s        = data_s$wd[seq(1, nrow(data_s), 2)] 
    BX_s        = data_s$BX[seq(1, nrow(data_s), 2)] 
    Vd_s        = data_s$vd[seq(1, nrow(data_s), 2)]
    WVd_s       = data_s$wvd[seq(1, nrow(data_s), 2)] 
    XVd_s       = data_s$xvd_in[seq(1, nrow(data_s), 2)] 
    
    data_s_g    = subset(data, S==k & G==1) #cal study, cases from internalized component
    Xd_s_g      = data_s_g$xd_in
    Wd_s_g      = data_s_g$wd
    BX_s_g      = data_s_g$BX
    Vd_s_g      = data_s_g$vd
    WVd_s_g     = data_s_g$wvd
    XVd_s_g     = data_s_g$xvd_in
    
    ## a and betax
    frac1             = (-b_hat[k]*Wd_s*Vd_s*beta_hatxv*exp(BX_s))/ ((1+exp(BX_s))^2) #FC
    frac2             = 1/(1+exp(BX_s_g)) #Int
    frac3             = -Xd_s_g*(beta_hatx+beta_hatxv*data_s_g$V)*exp(BX_s_g)/( (1+exp(BX_s_g))^2 ) #Int
    B12[2*k-1, 1]     = (1/np)*(sum(frac1)+sum(frac2)+sum(frac3))
    
    ## b and betax
    frac1             = Wd_s/(1+exp(BX_s)) #FC
    frac2             = -b_hat[k]*Wd_s*(beta_hatx*Wd_s+beta_hatxv*WVd_s)*exp(BX_s)/ ((1+exp(BX_s))^2) #FC
    frac3             = data_s_g$W/(1+exp(BX_s_g)) #Int
    frac4             = -data_s_g$W*Xd_s_g*(beta_hatx+beta_hatxv*data_s_g$V)*exp(BX_s_g)/( (1+exp(BX_s_g))^2 ) #Int
    B12[2*k, 1]       = (1/np)*(sum(frac1)+sum(frac2)+sum(frac3)+sum(frac4))
    
    ## a and betav
    frac1             = (-(Vd_s^2)*beta_hatxv*exp(BX_s))/ ((1+exp(BX_s))^2) #fc
    frac2             = (-Vd_s_g*(beta_hatx+beta_hatxv*data_s_g$V)*exp(BX_s_g))/ ( (1+exp(BX_s_g))^2 ) #int
    B12[2*k-1, 2]     = (1/np)*(sum(frac1)+sum(frac2))
    
    ## b and betav
    frac1             = (-Vd_s*(beta_hatx*Wd_s+beta_hatxv*WVd_s)*exp(BX_s))/ ((1+exp(BX_s))^2) #fc
    frac2             = -data_s_g$W*Vd_s_g*(beta_hatx+beta_hatxv*data_s_g$V)*exp(BX_s_g) /( (1+exp(BX_s_g))^2 ) #int
    B12[2*k, 2]       = (1/np)*(sum(frac1)+sum(frac2))
    
    ## a and betaxv
    frac1            = Vd_s/(1+exp(BX_s)) #fc
    frac2            = -Vd_s*beta_hatxv*XVd_s*exp(BX_s) / ((1+exp(BX_s))^2) #fc
    frac3            = data_s_g$V/(1+exp(BX_s_g)) #int
    frac4            = -(beta_hatx+beta_hatxv*data_s_g$V)*XVd_s_g*exp(BX_s_g)/( (1+exp(BX_s_g))^2 ) #int
    B12[2*k-1, 3]    = (1/np)*(sum(frac1)+sum(frac2)+sum(frac3)+sum(frac4))
    
    
    ## b and betaxv
    frac1            = WVd_s/(1+exp(BX_s)) #fc
    frac2            = -XVd_s*(beta_hatx*Wd_s+beta_hatxv*WVd_s)*exp(BX_s) / ((1+exp(BX_s))^2) #fc
    frac3            = data_s_g$WV/(1+exp(BX_s_g)) #int
    frac4            = -data_s_g$W*(beta_hatx+beta_hatxv*data_s_g$V)*XVd_s_g*exp(BX_s_g)/( (1+exp(BX_s_g))^2 ) #int
    B12[2*k, 3]      = (1/np)*(sum(frac1)+sum(frac2)+sum(frac3)+sum(frac4))
    
    ### Betaz and a,b
    for(i in 1:nz){ 
      col_num    = which(names(data) == covard[i])
      ZZ_s       = data_s[,col_num] #FC
      ZZ_s       = ZZ_s[seq(1, length(ZZ_s), 2)]
      ZZ_s_g     = data_s_g[,col_num] #IN
      
      ## Betaz and a
      frac1     = -Vd_s*ZZ_s*beta_hatxv*exp(BX_s)/(1+exp(BX_s))^2 #FC
      frac2     = -ZZ_s_g*(beta_hatx + beta_hatxv*data_s_g$V)*exp(BX_s_g)/(1+exp(BX_s_g))^2  #IN
      B12[(2*k-1),(3+i)]  = (1/np)*(sum(frac1)+sum(frac2))
      
      ## Betaz and b
      frac1     = -ZZ_s*(beta_hatx*Wd_s+beta_hatxv*WVd_s)*exp(BX_s)/(1+exp(BX_s))^2 #FC
      frac2     = -ZZ_s_g*data_s_g$W*(beta_hatx+beta_hatxv*data_s_g$V)*exp(BX_s_g)/(1+exp(BX_s_g))^2 #IN
      B12[(2*k),(3+i)]  = (1/np)*(sum(frac1)+sum(frac2))
    }
    
    
    
  }
  
  
  ### 5bi. Place B21 entries in main B matrix
  B[1:(2*nstud),(2*nstud+1):dim_sand] = B12
  B[(2*nstud+1):dim_sand,1:(2*nstud)] = t(B12)
  
  
  ### 5c: Lowermost right entry of B
  B22  = matrix(NA, nrow=(3+nz), ncol=(3+nz))
  
  ## betax_sq, betav_sq, betaxv_sq
  B22[1,1]  = (1/np)*sum(-(Xd^2)*exp(BX)/(exp(BX)+1)^2)
  B22[2,2]  = (1/np)*sum(-(Vd^2)*exp(BX)/(exp(BX)+1)^2)
  B22[3,3]  = (1/np)*sum(-(XVd^2)*exp(BX)/(exp(BX)+1)^2)
  
  ## Cross terms
  B22[1,2]  =  (1/np)*sum(-(Xd*Vd)*exp(BX)/(exp(BX)+1)^2)
  B22[2,1]  =  (1/np)*sum(-(Xd*Vd)*exp(BX)/(exp(BX)+1)^2)
  B22[1,3]  =  (1/np)*sum(-(Xd*XVd)*exp(BX)/(exp(BX)+1)^2)
  B22[3,1]  =  (1/np)*sum(-(Xd*XVd)*exp(BX)/(exp(BX)+1)^2)
  B22[3,2]  =  (1/np)*sum(-(XVd*Vd)*exp(BX)/(exp(BX)+1)^2)
  B22[2,3]  =  (1/np)*sum(-(XVd*Vd)*exp(BX)/(exp(BX)+1)^2)
  
  
  ### Computations in B22 involving betax/v/xv and betaz, betaz squared
  for(i in 1:nz){ 
    col_num  = which(names(data) == covard[i])
    ZZ       = data[,col_num]
    ZZ       = ZZ[seq(1, length(ZZ), 2)]
    fracXZ   = -Xd*ZZ*exp(BX)/(1+exp(BX))^2
    fracVZ   = -Vd*ZZ*exp(BX)/(1+exp(BX))^2
    fracXVZ  = -XVd*ZZ*exp(BX)/(1+exp(BX))^2
    fracZZ   = -(ZZ^2)*exp(BX)/(1+exp(BX))^2
    ## Cross terms
    B22[1,(3+i)]      = (1/np)*sum(fracXZ)
    B22[(3+i),1]      = (1/np)*sum(fracXZ)
    B22[2,(3+i)]      = (1/np)*sum(fracVZ)
    B22[(3+i),2]      = (1/np)*sum(fracVZ)
    B22[3,(3+i)]      = (1/np)*sum(fracXVZ)
    B22[(3+i),3]      = (1/np)*sum(fracXVZ)
    ## Squared terms
    B22[(3+i),(3+i)]    = (1/np)*sum(fracZZ) 
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
    fracZZ = -ZZ1*ZZ2*exp(BX)/(1+exp(BX))^2
    B22[(3+pair[2]),(3+pair[1])] = (1/np)*sum(fracZZ)
    B22[(3+pair[1]),(3+pair[2])] = (1/np)*sum(fracZZ)
  }
  
  
  ## Replace B22 in main B matrix
  B[(2*nstud+1):dim_sand, (2*nstud+1):dim_sand] = B22
  
  ## Eliminate entries associated with the reference lab studies
  B = B[(2*nref+1):dim_sand,(2*nref+1):dim_sand]
  
  ##################### Final variance computations ###################
  
  #### 6a. Create V matrix and extract the variance element
  V                = solve(B)%*%A%*%t(solve(B))
  output_in[3]     = (1/np)*V[(dim_sand-2),(dim_sand-2)] #betax
  output_in[6]     = (1/np)*V[(dim_sand-1),(dim_sand-1)] #betav
  output_in[9]     = (1/np)*V[dim_sand,dim_sand] #betaxv
  
  #### 99. Return appropriate output
  return(output_in)
  
}


in_df_no_cov_xv = function(data, X, S, H, W, Y, V, G, strata, nstud, nref){
  
  #### 1. Create matrix to store output
  output_in = matrix(NA, ncol=9, nrow=1)
  colnames(output_in) = c("betax","RRx", "VARx","betav","RRv", "VARv","betaxv","RRxv", "VARxv")
  
  #### 2. Rename variables names given to function to ensure consistency
  data$X = data[[X]]
  data$S = data[[S]]
  data$H = data[[H]]
  data$G = data[[G]]
  data$W = data[[W]]
  data$Y = data[[Y]]
  data$V = data[[V]]
  data$strata = data[[strata]]
  
  #### 2b. Sort data frame, in case it wasn't already
  data = data[with(data, order(strata, -Y)),]
  
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
  data$xhat_fc    = ifelse(data$H==2, data$X, data$a_hat + data$b_hat*data$W)
  data$xhat_in    = ifelse(data$H==0, data$xhat_fc, data$X) 
  data$xhat_in_v  = ifelse(data$H==0, data$xhat_fc*data$V, data$X*data$V)
  
  #### 6. Obtain point estimate from standard logistic regression
  formula            = as.formula(paste("Y~xhat_in+V+xhat_in_v+strata(strata)",sep="+"))
  in_fit             = clogit(formula, data=data)
  beta_hatx          = in_fit$coefficients[1]
  beta_hatv          = in_fit$coefficients[2]
  beta_hatxv         = in_fit$coefficients[3]
  output_in[c(1,2)]  = c(beta_hatx,exp(beta_hatx))
  output_in[c(4,5)]  = c(beta_hatv,exp(beta_hatv))
  output_in[c(7,8)]  = c(beta_hatxv,exp(beta_hatxv))
  
  #### 6.b. Compute necessary differences for sandwich computations
  
  ## Take differences for each observation and store in matrix
  
  data$wv      = data$W*data$V
  cases        = subset(data, Y==1)
  controls     = subset(data, Y==0)
  xd_in        = cases$xhat_in - controls$xhat_in
  vd           = cases$V - controls$V
  xvd_in       = cases$xhat_in_v - controls$xhat_in_v
  wd           = cases$W - controls$W
  wvd          = cases$wv - controls$wv
  
  data$xd_in   = rep(xd_in, each=2)
  data$vd      = rep(vd, each=2)
  data$xvd_in  = rep(xvd_in, each=2)
  data$wd      = rep(wd, each=2)
  data$wvd     = rep(wvd, each=2)
  data$BX      = beta_hatx*data$xd_in + beta_hatv*data$vd + beta_hatxv*data$xvd_in
  
  #### 7. Compute variance; prepare matrices 
  dim_sand  = 2*nstud+3
  A         = matrix(0, ncol = dim_sand, nrow = dim_sand)
  B         = matrix(0, ncol = dim_sand, nrow = dim_sand)
  
  ###################### A MATRIX #######################################
  
  #### 4a: The upper block diagonals of A
  for(k in 1:nstud){
    cal_data_s               = subset(data, (H==1 & S==k)) # Specific cal data with only controls
    A[2*k-1, 2*k-1]          = (1/np)*sum((cal_data_s$X-a_hat[k]-b_hat[k]*cal_data_s$W)^2)
    A[2*k-1,(2*k-1)+1]       = (1/np)*sum(cal_data_s$W*(cal_data_s$X-a_hat[k]-b_hat[k]*cal_data_s$W)^2)
    A[(2*k-1)+1,(2*k-1)]     = A[2*k-1, (2*k-1)+1]
    A[(2*k-1)+1,(2*k-1)+1]   = (1/np)*sum((cal_data_s$W^2)*(cal_data_s$X-a_hat[k]-b_hat[k]*cal_data_s$W)^2)
  }
  
  #### 4b: Upper right entries of A (A12)
  A12    = matrix(NA, ncol=3, nrow=(2*nstud))
  
  for(k in 1:nstud){
    data_s_h     = subset(data, S==k & H==1) # Calibration controls
    Xd_s_h       = data_s_h$xd_in
    Vd_s_h       = data_s_h$vd
    XVd_s_h      = data_s_h$xvd_in
    BX_s_h       = data_s_h$BX
    psi1         = data_s_h$X-a_hat[k]-b_hat[k]*data_s_h$W # estimating equ 1
    psi2         = (data_s_h$X-a_hat[k]-b_hat[k]*data_s_h$W)*data_s_h$W # estimating equ 2
    
    psi_x         = Xd_s_h/(exp(BX_s_h)+1) # betax estimating equation
    psi_v         = Vd_s_h/(exp(BX_s_h)+1) # betav estimating equation
    psi_xv        = XVd_s_h/(exp(BX_s_h)+1) # betaxv estimating equation
    
    ## Betax and a,b
    A12[(2*k-1),1]   = (1/np)*sum(psi1*psi_x)
    A12[(2*k),1]     = (1/np)*sum(psi2*psi_x)
    
    ## Betav and a,b
    A12[(2*k-1),2]   = (1/np)*sum(psi1*psi_v)
    A12[(2*k),2]     = (1/np)*sum(psi2*psi_v)
    
    ## Betaxv and a,b
    A12[(2*k-1),3]   = (1/np)*sum(psi1*psi_xv)
    A12[(2*k),3]     = (1/np)*sum(psi2*psi_xv)
  }
  
  A[1:(2*nstud), (dim_sand-2):dim_sand] = A12 # Fill appropriate piece of A matrix
  A21                                   = t(A12) # Reflect entries for A21
  A[(dim_sand-2):dim_sand,1:(2*nstud)]  = A21
  
  #### 4c: A22 entry, is 3x3
  A22      = matrix(NA,nrow=3,ncol=3)
  Xd       = data$xd_in[seq(1, nrow(data), 2)] 
  Vd       = data$vd[seq(1, nrow(data), 2)] 
  XVd      = data$xvd_in[seq(1, nrow(data), 2)] 
  BX       = data$BX[seq(1, nrow(data), 2)] 
  
  ## Squared terms on diagonal of A22
  A22[1,1] = (1/np)*sum( (Xd/(exp(BX)+1))^2 )
  A22[2,2] = (1/np)*sum( (Vd/(exp(BX)+1))^2 )
  A22[3,3] = (1/np)*sum( (XVd/(exp(BX)+1))^2 )
  
  ## Betax and betav of A22
  A22[1,2] = (1/np)*sum( (Xd*Vd/(exp(BX)+1)^2) )
  A22[2,1] = A22[1,2]
  
  ## Betav and betaxv of A22
  A22[2,3] = (1/np)*sum( (Vd*XVd/(exp(BX)+1)^2) )
  A22[3,2] = A22[2,3]
  
  ## Betax and betaxv of A22
  A22[1,3] = (1/np)*sum( (Xd*XVd/(exp(BX)+1)^2) )
  A22[3,1] = A22[1,3]
  
  ## Replace A22 in main matrix of A
  A[(2*nstud+1):dim_sand,(2*nstud+1):dim_sand ] = A22
  
  
  ## Eliminate extra entries from A associated with calibration studies
  A = A[(2*nref+1):dim_sand,(2*nref+1):dim_sand]
  
  ###################### B MATRIX #######################################
  
  #### 5a: Block diagonal entries of B associated with the calibration studies.
  
  for(k in 1:nstud){
    cal_data_s               = subset(data, S==k & H==1) # Specific cal data with only controls
    B[2*k-1,2*k-1]           = -1
    B[2*k-1,(2*k-1)+1]       = (-1/np)*sum(cal_data_s$W)
    B[(2*k-1)+1,(2*k-1)]     = B[2*k-1,(2*k-1)+1]
    B[(2*k-1)+1,(2*k-1)+1]   = (-1/np)*sum(cal_data_s$W^2)
  }
  
  ### 5b: Off diagonal even entries of B 
  B12 = matrix(NA, ncol=3, nrow=(2*nstud))
  
  for(k in 1:nstud){
    data_s      = subset(data, S==k & H==0 & G==0) #FC-esque component
    Xd_s        = data_s$xd_in[seq(1, nrow(data_s), 2)]
    Wd_s        = data_s$wd[seq(1, nrow(data_s), 2)] 
    BX_s        = data_s$BX[seq(1, nrow(data_s), 2)] 
    Vd_s        = data_s$vd[seq(1, nrow(data_s), 2)]
    WVd_s       = data_s$wvd[seq(1, nrow(data_s), 2)] 
    XVd_s       = data_s$xvd_in[seq(1, nrow(data_s), 2)] 
    
    data_s_g    = subset(data, S==k & G==1) #cal study, cases from internalized component
    Xd_s_g      = data_s_g$xd_in
    Wd_s_g      = data_s_g$wd
    BX_s_g      = data_s_g$BX
    Vd_s_g      = data_s_g$vd
    WVd_s_g     = data_s_g$wvd
    XVd_s_g     = data_s_g$xvd_in
    
    ## a and betax
    frac1             = (-b_hat[k]*Wd_s*Vd_s*beta_hatxv*exp(BX_s))/ ((1+exp(BX_s))^2) #FC
    frac2             = 1/(1+exp(BX_s_g)) #Int
    frac3             = -Xd_s_g*(beta_hatx+beta_hatxv*data_s_g$V)*exp(BX_s_g)/( (1+exp(BX_s_g))^2 ) #Int
    B12[2*k-1, 1]     = (1/np)*(sum(frac1)+sum(frac2)+sum(frac3))
    
    ## b and betax
    frac1             = Wd_s/(1+exp(BX_s)) #FC
    frac2             = -b_hat[k]*Wd_s*(beta_hatx*Wd_s+beta_hatxv*WVd_s)*exp(BX_s)/ ((1+exp(BX_s))^2) #FC
    frac3             = data_s_g$W/(1+exp(BX_s_g)) #Int
    frac4             = -data_s_g$W*Xd_s_g*(beta_hatx+beta_hatxv*data_s_g$V)*exp(BX_s_g)/( (1+exp(BX_s_g))^2 ) #Int
    B12[2*k, 1]       = (1/np)*(sum(frac1)+sum(frac2)+sum(frac3)+sum(frac4))
    
    ## a and betav
    frac1             = (-(Vd_s^2)*beta_hatxv*exp(BX_s))/ ((1+exp(BX_s))^2) #fc
    frac2             = (-Vd_s_g*(beta_hatx+beta_hatxv*data_s_g$V)*exp(BX_s_g))/ ( (1+exp(BX_s_g))^2 ) #int
    B12[2*k-1, 2]     = (1/np)*(sum(frac1)+sum(frac2))
    
    ## b and betav
    frac1             = (-Vd_s*(beta_hatx*Wd_s+beta_hatxv*WVd_s)*exp(BX_s))/ ((1+exp(BX_s))^2) #fc
    frac2             = -data_s_g$W*Vd_s_g*(beta_hatx+beta_hatxv*data_s_g$V)*exp(BX_s_g) /( (1+exp(BX_s_g))^2 ) #int
    B12[2*k, 2]       = (1/np)*(sum(frac1)+sum(frac2))
    
    ## a and betaxv
    frac1            = Vd_s/(1+exp(BX_s)) #fc
    frac2            = -Vd_s*beta_hatxv*XVd_s*exp(BX_s) / ((1+exp(BX_s))^2) #fc
    frac3            = data_s_g$V/(1+exp(BX_s_g)) #int
    frac4            = -(beta_hatx+beta_hatxv*data_s_g$V)*XVd_s_g*exp(BX_s_g)/( (1+exp(BX_s_g))^2 ) #int
    B12[2*k-1, 3]    = (1/np)*(sum(frac1)+sum(frac2)+sum(frac3)+sum(frac4))
    
    
    ## b and betaxv
    frac1            = WVd_s/(1+exp(BX_s)) #fc
    frac2            = -XVd_s*(beta_hatx*Wd_s+beta_hatxv*WVd_s)*exp(BX_s) / ((1+exp(BX_s))^2) #fc
    frac3            = data_s_g$WV/(1+exp(BX_s_g)) #int
    frac4            = -data_s_g$W*(beta_hatx+beta_hatxv*data_s_g$V)*XVd_s_g*exp(BX_s_g)/( (1+exp(BX_s_g))^2 ) #int
    B12[2*k, 3]      = (1/np)*(sum(frac1)+sum(frac2)+sum(frac3)+sum(frac4))
    
  }
  
  
  ### 5bi. Place B21 entries in main B matrix
  B[1:(2*nstud), (dim_sand-2):dim_sand] = B12
  B21                                   = t(B12) # Reflect entries for B21
  B[(dim_sand-2):dim_sand,1:(2*nstud)]  = B21
  
  
  ### 5c: Lowermost right entry of B
  B22  = matrix(NA, nrow=3, ncol=3)
  
  ## betax_sq, betav_sq, betaxv_sq
  B22[1,1]  = (1/np)*sum(-(Xd^2)*exp(BX)/(exp(BX)+1)^2)
  B22[2,2]  = (1/np)*sum(-(Vd^2)*exp(BX)/(exp(BX)+1)^2)
  B22[3,3]  = (1/np)*sum(-(XVd^2)*exp(BX)/(exp(BX)+1)^2)
  
  ## Cross terms
  B22[1,2]  =  (1/np)*sum(-(Xd*Vd)*exp(BX)/(exp(BX)+1)^2)
  B22[2,1]  =  (1/np)*sum(-(Xd*Vd)*exp(BX)/(exp(BX)+1)^2)
  B22[1,3]  =  (1/np)*sum(-(Xd*XVd)*exp(BX)/(exp(BX)+1)^2)
  B22[3,1]  =  (1/np)*sum(-(Xd*XVd)*exp(BX)/(exp(BX)+1)^2)
  B22[3,2]  =  (1/np)*sum(-(XVd*Vd)*exp(BX)/(exp(BX)+1)^2)
  B22[2,3]  =  (1/np)*sum(-(XVd*Vd)*exp(BX)/(exp(BX)+1)^2)
  
  ## Replace B22 in main B matrix
  B[(2*nstud+1):dim_sand, (2*nstud+1):dim_sand] = B22
  
  ## Eliminate entries associated with the reference lab studies
  B = B[(2*nref+1):dim_sand,(2*nref+1):dim_sand]
  
  ##################### Final variance computations ###################
  
  #### 6a. Create V matrix and extract the variance element
  V                = solve(B)%*%A%*%t(solve(B))
  output_in[3]     = (1/np)*V[(dim_sand-2),(dim_sand-2)] #betax
  output_in[6]     = (1/np)*V[(dim_sand-1),(dim_sand-1)] #betav
  output_in[9]     = (1/np)*V[dim_sand,dim_sand] #betaxv
  
  #### 99. Return appropriate output
  return(output_in)
  
}


fc_df_xv = function(data, X, S, H, W, Y, V, strata, nstud, nref, covar=NA){
  
  #### 1. Create matrix to store output
  output_fc = matrix(NA, ncol=9, nrow=1)
  colnames(output_fc) = c("betax","RRx", "VARx","betav","RRv", "VARv","betaxv","RRxv", "VARxv")
  
  #### 2. Rename variables names given to function to ensure consistency
  data$X = data[[X]]
  data$S = data[[S]]
  data$H = data[[H]]
  data$W = data[[W]]
  data$Y = data[[Y]]
  data$V = data[[V]]
  data$strata = data[[strata]]
  
  #### 2b. Sort data frame, in case it wasn't already
  data = data[with(data, order(strata, -Y)),]
  
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
  data$xhat_fc   = ifelse(data$H==2, data$X, data$a_hat + data$b_hat*data$W)
  data$xhat_fc_v = ifelse(data$H==2, data$X*data$V, data$xhat_fc*data$V)
  
  #### 6. Obtain point estimate from standard logistic regression
  formula           = as.formula(paste("Y~xhat_fc+V+xhat_fc_v+strata(strata)",paste(covar,collapse ="+"),sep="+"))
  fc_fit            = clogit(formula, data=data)
  beta_hatx         = fc_fit$coefficients[1]
  beta_hatv         = fc_fit$coefficients[2]
  beta_hatxv        = fc_fit$coefficients[3]
  betazhat          = fc_fit$coefficients[-c(1,2,3)] # extract betazhat values
  output_fc[c(1,2)] = c(beta_hatx,exp(beta_hatx))
  output_fc[c(4,5)] = c(beta_hatv,exp(beta_hatv))
  output_fc[c(7,8)] = c(beta_hatxv,exp(beta_hatxv))
  
  #### 6.b. Compute necessary differences for sandwich computations
  
  ## Take differences for each observation and store in matrix
  
  data$WV      = data$W*data$V
  cases        = subset(data, Y==1)
  controls     = subset(data, Y==0)
  xd_fc        = cases$xhat_fc - controls$xhat_fc
  vd           = cases$V - controls$V
  xvd_fc       = cases$xhat_fc_v - controls$xhat_fc_v
  wd           = cases$W - controls$W
  wvd          = cases$WV - controls$WV
  
  data$xd_fc   = rep(xd_fc, each=2)
  data$vd      = rep(vd, each=2)
  data$xvd_fc  = rep(xvd_fc, each=2)
  data$wd      = rep(wd, each=2)
  data$wvd     = rep(wvd, each=2)
  
  ## Z differences and linear combination
  col_nums        = which(names(data) %in% covar)
  Zvec_cases      = as.matrix(cases[,col_nums])
  Zvec_controls   = as.matrix(controls[,col_nums])
  Zd              = Zvec_cases - Zvec_controls
  BZ              = Zd %*% betazhat
  data$BZ         = rep(BZ,each=2)
  data$BX         = beta_hatx*data$xd_fc + beta_hatv*data$vd + beta_hatxv*data$xvd_fc + data$BZ
  
  ### Zd columns; update naming, add to data frame
  covdiffs        = Zd[rep(1:nrow(Zd), each = 2), ]
  covard          = paste(covar,"d",sep = "")
  colnames(covdiffs) = covard
  data            = cbind(data,covdiffs)
  
  #### 7. Compute variance; prepare matrices 
  dim_sand  = 2*nstud+3+nz
  A         = matrix(0, ncol = dim_sand, nrow = dim_sand)
  B         = matrix(0, ncol = dim_sand, nrow = dim_sand)
  
  ###################### A MATRIX #######################################
  
  #### 4a: The upper block diagonals of A
  for(k in 1:nstud){
    cal_data_s               = subset(data, (H==1 & S==k)) # Specific cal data with only controls
    A[2*k-1, 2*k-1]          = (1/np)*sum((cal_data_s$X-a_hat[k]-b_hat[k]*cal_data_s$W)^2)
    A[2*k-1,(2*k-1)+1]       = (1/np)*sum(cal_data_s$W*(cal_data_s$X-a_hat[k]-b_hat[k]*cal_data_s$W)^2)
    A[(2*k-1)+1,(2*k-1)]     = A[2*k-1, (2*k-1)+1]
    A[(2*k-1)+1,(2*k-1)+1]   = (1/np)*sum((cal_data_s$W^2)*(cal_data_s$X-a_hat[k]-b_hat[k]*cal_data_s$W)^2)
  }
  
  #### 4b: Upper right entries of A (A12)
  A12    = matrix(NA, ncol=(3+nz), nrow=(2*nstud))
  
  for(k in 1:nstud){
    data_s_h     = subset(data, S==k & H==1) # Calibration controls
    Xd_s_h       = data_s_h$xd_fc
    Vd_s_h       = data_s_h$vd
    XVd_s_h      = data_s_h$xvd_fc
    BX_s_h       = data_s_h$BX
    psi1         = data_s_h$X-a_hat[k]-b_hat[k]*data_s_h$W # estimating equ 1
    psi2         = (data_s_h$X-a_hat[k]-b_hat[k]*data_s_h$W)*data_s_h$W # estimating equ 2
    
    psi_x         = Xd_s_h/(exp(BX_s_h)+1) # betax estimating equation
    psi_v         = Vd_s_h/(exp(BX_s_h)+1) # betav estimating equation
    psi_xv        = XVd_s_h/(exp(BX_s_h)+1) # betaxv estimating equation
    
    ## Betax and a,b
    A12[(2*k-1),1]   = (1/np)*sum(psi1*psi_x)
    A12[(2*k),  1]   = (1/np)*sum(psi2*psi_x)
    
    ## Betav and a,b
    A12[(2*k-1),2]   = (1/np)*sum(psi1*psi_v)
    A12[(2*k),  2]   = (1/np)*sum(psi2*psi_v)
    
    ## Betax and a,b
    A12[(2*k-1),3]   = (1/np)*sum(psi1*psi_xv)
    A12[(2*k),  3]   = (1/np)*sum(psi2*psi_xv)
    
    ## Grab the Z covariate vector of interest; do computations with that particular Z
    ## Betaz and a,b
    for(i in 1:nz){ 
      col_num  = which(names(data_s_h) == covard[i])
      ZZ       = data_s_h[,col_num]
      psiZ     = ZZ/(1+exp(BX_s_h))
      
      A12[(2*k-1),(3+i)]  = (1/np)*sum(psi1*psiZ)
      A12[(2*k),(3+i)]    = (1/np)*sum(psi2*psiZ) 
    }
    
  }
  
  A[1:(2*nstud), (2*nstud+1):dim_sand] = A12 # Fill appropriate piece of A matrix
  A21                                  = t(A12) # Reflect entries for A21
  A[(2*nstud+1):dim_sand,1:(2*nstud)]  = A21
  
  #### 4c: Preparations for computations in A22
  Xd       = data$xd_fc[seq(1, nrow(data), 2)] 
  Vd       = data$vd[seq(1, nrow(data), 2)] 
  XVd      = data$xvd_fc[seq(1, nrow(data), 2)] 
  BX       = data$BX[seq(1, nrow(data), 2)] 
  
  ### A22; need to take cross values of all beta est equ
  A22 = matrix(NA,nrow=(3+nz),ncol=(3+nz))
  
  ## Squared terms of X,V,XV on diagonal of A22
  psiX     = Xd/(exp(BX)+1)
  psiV     = Vd/(exp(BX)+1)
  psiXV    = XVd/(exp(BX)+1)
  A22[1,1] = (1/np)*sum( psiX^2 )
  A22[2,2] = (1/np)*sum( psiV^2 )
  A22[3,3] = (1/np)*sum( psiXV^2 )

  ## Betax and betav of A22
  A22[1,2] = (1/np)*sum( psiX*psiV )
  A22[2,1] = A22[1,2]
  
  ## Betav and betaxv of A22
  A22[2,3] = (1/np)*sum( psiV*psiXV )
  A22[3,2] = A22[2,3]
  
  ## Betax and betaxv of A22
  A22[1,3] = (1/np)*sum( psiX*psiXV )
  A22[3,1] = A22[1,3]
  
  ######## Now consider A22 entries involving an estimating equation with betaz
  ## Select specific psiZ and compute products
  for(i in 1:nz){ 
    col_num  = which(names(data) == covard[i])
    ZZ       = data[,col_num]
    ZZ       = ZZ[seq(1, length(ZZ), 2)]
    psiZ     = ZZ/(1+exp(BX))
    
    A22[1,(3+i)]  = (1/np)*sum(psiX*psiZ)
    A22[(3+i),1]  = (1/np)*sum(psiX*psiZ)
    A22[2,(3+i)]  = (1/np)*sum(psiV*psiZ)
    A22[(3+i),2]  = (1/np)*sum(psiV*psiZ)
    A22[3,(3+i)]  = (1/np)*sum(psiXV*psiZ)
    A22[(3+i),3]  = (1/np)*sum(psiXV*psiZ)
    A22[(3+i),(3+i)]    = (1/np)*sum(psiZ*psiZ) 
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
    psi_betaz1    = ZZ1/(1+exp(BX))
    
    col_num2 = which(names(data) == covar[pair[2]])
    ZZ2 = data[,col_num2]
    ZZ2 = ZZ2[seq(1, length(ZZ2), 2)] # get unique Zd obs
    psi_betaz2    = ZZ2/(1+exp(BX))
    
    #Fill matrix elements
    A22[(3+pair[2]),(3+pair[1])] = (1/np)*sum(psi_betaz1*psi_betaz2)
    A22[(3+pair[1]),(3+pair[2])] = (1/np)*sum(psi_betaz1*psi_betaz2)
  }
  
  ## Replace A22 in main matrix of A
  A[(2*nstud+1):dim_sand,(2*nstud+1):dim_sand ] = A22
  
  ## Eliminate extra entries from A associated with calibration studies
  A = A[(2*nref+1):dim_sand,(2*nref+1):dim_sand]
  
  ###################### B MATRIX #######################################
  
  #### 5a: Block diagonal entries of B associated with the calibration studies.
  for(k in 1:nstud){
    cal_data_s               = subset(data, S==k & H==1) # Specific cal data with only controls
    B[2*k-1,2*k-1]           = -1
    B[2*k-1,(2*k-1)+1]       = (-1/np)*sum(cal_data_s$W)
    B[(2*k-1)+1,(2*k-1)]     = B[2*k-1,(2*k-1)+1]
    B[(2*k-1)+1,(2*k-1)+1]   = (-1/np)*sum(cal_data_s$W^2)
  }
  
  
  ### 5b: Off diagonal even entries of B 
  B12 = matrix(NA, ncol=(3+nz), nrow=(2*nstud))
  
  for(k in 1:nstud){
    data_s     = subset(data, S==k) # Study-specific
    Wd_s       = data_s$wd[seq(1, nrow(data_s), 2)] 
    Vd_s       = data_s$vd[seq(1, nrow(data_s), 2)] 
    WVd_s      = data_s$wvd[seq(1, nrow(data_s), 2)]
    XVd_s      = data_s$xvd_fc[seq(1, nrow(data_s), 2)] 
    BX_s       = data_s$BX[seq(1, nrow(data_s), 2)] 
    
    ## a and betax
    frac             = (-b_hat[k]*Wd_s*Vd_s*beta_hatxv*exp(BX_s))/ ((1+exp(BX_s))^2)
    B12[2*k-1, 1]    = (1/np)*sum(frac)
    
    ## b and betax
    frac1            = Wd_s/(1+exp(BX_s))
    frac2            = -b_hat[k]*Wd_s*(beta_hatx*Wd_s+beta_hatxv*WVd_s)*exp(BX_s)/ ((1+exp(BX_s))^2)
    B12[2*k, 1]      = (1/np)*(sum(frac1)+sum(frac2))
    
    ## a and betav
    frac             = (-(Vd_s^2)*beta_hatxv*exp(BX_s))/ ((1+exp(BX_s))^2)
    B12[2*k-1, 2]    = (1/np)*sum(frac)
    
    ## b and betav
    frac             = (-Vd_s*(beta_hatx*Wd_s + beta_hatxv*WVd_s)*exp(BX_s))/ ((1+exp(BX_s))^2)
    B12[2*k, 2]      = (1/np)*sum(frac)
    
    ## a and betaxv
    frac1            = Vd_s/(1+exp(BX_s))
    frac2            = -Vd_s*beta_hatxv*XVd_s*exp(BX_s) / ((1+exp(BX_s))^2)
    B12[2*k-1, 3]    = (1/np)*(sum(frac1)+sum(frac2))
    
    ## b and betaxv
    frac1            = WVd_s/(1+exp(BX_s))
    frac2            = -XVd_s*(beta_hatx*Wd_s+beta_hatxv*WVd_s)*exp(BX_s) / ((1+exp(BX_s))^2)
    B12[2*k, 3]      = (1/np)*(sum(frac1)+sum(frac2))
    
    ### Betaz and a,b
    for(i in 1:nz){ 
      col_num  = which(names(data) == covard[i])
      ZZ       = data_s[,col_num]
      ZZ       = ZZ[seq(1, length(ZZ), 2)]
      
      ## Betaz and a
      frac     = -Vd_s*ZZ*beta_hatxv*exp(BX_s)/(1+exp(BX_s))^2
      B12[(2*k-1),(3+i)]  = (1/np)*sum(frac)
      
      ## Betaz and b
      frac     = -ZZ*(beta_hatx*Wd_s+beta_hatxv*WVd_s)*exp(BX_s)/(1+exp(BX_s))^2
      B12[(2*k),(3+i)]  = (1/np)*sum(frac)
    }
    
  }
  
  
  ### 5bi. Place B21 entries in main B matrix
  B[1:(2*nstud),(2*nstud+1):dim_sand] = B12
  B[(2*nstud+1):dim_sand,1:(2*nstud)] = t(B12)
  
  
  ### 5c: Lowermost right entry of B
  B22        = matrix(NA, nrow=(3+nz), ncol=(3+nz))
  
  ## betax_sq, betav_sq, betaxv_sq
  B22[1,1]  = (1/np)*sum(-(Xd^2)*exp(BX)/(exp(BX)+1)^2)
  B22[2,2]  = (1/np)*sum(-(Vd^2)*exp(BX)/(exp(BX)+1)^2)
  B22[3,3]  = (1/np)*sum(-(XVd^2)*exp(BX)/(exp(BX)+1)^2)
  
  ## Cross terms
  B22[1,2]  =  (1/np)*sum(-(Xd*Vd)*exp(BX)/(exp(BX)+1)^2)
  B22[2,1]  =  (1/np)*sum(-(Xd*Vd)*exp(BX)/(exp(BX)+1)^2)
  B22[1,3]  =  (1/np)*sum(-(Xd*XVd)*exp(BX)/(exp(BX)+1)^2)
  B22[3,1]  =  (1/np)*sum(-(Xd*XVd)*exp(BX)/(exp(BX)+1)^2)
  B22[3,2]  =  (1/np)*sum(-(XVd*Vd)*exp(BX)/(exp(BX)+1)^2)
  B22[2,3]  =  (1/np)*sum(-(XVd*Vd)*exp(BX)/(exp(BX)+1)^2)
  
  ### Computations in B22 involving betax/v/xv and betaz, betaz squared
  for(i in 1:nz){ 
    col_num  = which(names(data) == covard[i])
    ZZ       = data[,col_num]
    ZZ       = ZZ[seq(1, length(ZZ), 2)]
    fracXZ   = -Xd*ZZ*exp(BX)/(1+exp(BX))^2
    fracVZ   = -Vd*ZZ*exp(BX)/(1+exp(BX))^2
    fracXVZ  = -XVd*ZZ*exp(BX)/(1+exp(BX))^2
    fracZZ   = -(ZZ^2)*exp(BX)/(1+exp(BX))^2
    ## Cross terms
    B22[1,(3+i)]      = (1/np)*sum(fracXZ)
    B22[(3+i),1]      = (1/np)*sum(fracXZ)
    B22[2,(3+i)]      = (1/np)*sum(fracVZ)
    B22[(3+i),2]      = (1/np)*sum(fracVZ)
    B22[3,(3+i)]      = (1/np)*sum(fracXVZ)
    B22[(3+i),3]      = (1/np)*sum(fracXVZ)
    ## Squared terms
    B22[(3+i),(3+i)]    = (1/np)*sum(fracZZ) 
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
    fracZZ = -ZZ1*ZZ2*exp(BX)/(1+exp(BX))^2
    B22[(3+pair[2]),(3+pair[1])] = (1/np)*sum(fracZZ)
    B22[(3+pair[1]),(3+pair[2])] = (1/np)*sum(fracZZ)
  }
  
  ## Replace B22 in main B matrix
  B[(2*nstud+1):dim_sand, (2*nstud+1):dim_sand] = B22
  
  ## Eliminate entries associated with the reference lab studies
  B = B[(2*nref+1):dim_sand,(2*nref+1):dim_sand]
  
  ##################### Final variance computations ###################
  
  #### 6a. Create V matrix and extract the variance element
  V                = solve(B)%*%A%*%t(solve(B))
  output_fc[3]     = (1/np)*V[(2*nstud+1),(2*nstud+1)] #betax
  output_fc[6]     = (1/np)*V[(2*nstud+2),(2*nstud+2)] #betav
  output_fc[9]     = (1/np)*V[(2*nstud+3),(2*nstud+3)] #betaxv
  
  #### 99. Return appropriate output
  return(output_fc)
  
}


fc_df_no_cov_xv = function(data, X, S, H, W, Y, V, strata, nstud, nref){
  
  #### 1. Create matrix to store output
  output_fc = matrix(NA, ncol=9, nrow=1)
  colnames(output_fc) = c("betax","RRx", "VARx","betav","RRv", "VARv","betaxv","RRxv", "VARxv")
  
  #### 2. Rename variables names given to function to ensure consistency
  data$X = data[[X]]
  data$S = data[[S]]
  data$H = data[[H]]
  data$W = data[[W]]
  data$Y = data[[Y]]
  data$V = data[[V]]
  data$strata = data[[strata]]
  
  #### 2b. Sort data frame, in case it wasn't already
  data = data[with(data, order(strata, -Y)),]
  
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
  data$xhat_fc   = ifelse(data$H==2, data$X, data$a_hat + data$b_hat*data$W)
  data$xhat_fc_v = ifelse(data$H==2, data$X*data$V, data$xhat_fc*data$V)
  
  #### 6. Obtain point estimate from standard logistic regression
  formula           = as.formula(paste("Y~xhat_fc+V+xhat_fc_v+strata(strata)",sep="+"))
  fc_fit            = clogit(formula, data=data)
  beta_hatx         = fc_fit$coefficients[1]
  beta_hatv         = fc_fit$coefficients[2]
  beta_hatxv        = fc_fit$coefficients[3]
  output_fc[c(1,2)] = c(beta_hatx,exp(beta_hatx))
  output_fc[c(4,5)] = c(beta_hatv,exp(beta_hatv))
  output_fc[c(7,8)] = c(beta_hatxv,exp(beta_hatxv))
  
  #### 6.b. Compute necessary differences for sandwich computations
  
  ## Take differences for each observation and store in matrix

  data$WV      = data$W*data$V
  cases        = subset(data, Y==1)
  controls     = subset(data, Y==0)
  xd_fc        = cases$xhat_fc - controls$xhat_fc
  vd           = cases$V - controls$V
  xvd_fc       = cases$xhat_fc_v - controls$xhat_fc_v
  wd           = cases$W - controls$W
  wvd          = cases$WV - controls$WV
  
  data$xd_fc   = rep(xd_fc, each=2)
  data$vd      = rep(vd, each=2)
  data$xvd_fc  = rep(xvd_fc, each=2)
  data$wd      = rep(wd, each=2)
  data$wvd     = rep(wvd, each=2)
  data$BX      = beta_hatx*data$xd_fc + beta_hatv*data$vd + beta_hatxv*data$xvd_fc
  
  #### 7. Compute variance; prepare matrices 
  dim_sand  = 2*nstud+3
  A         = matrix(0, ncol = dim_sand, nrow = dim_sand)
  B         = matrix(0, ncol = dim_sand, nrow = dim_sand)
  
  ###################### A MATRIX #######################################
  
  #### 4a: The upper block diagonals of A
  for(k in 1:nstud){
    cal_data_s               = subset(data, (H==1 & S==k)) # Specific cal data with only controls
    A[2*k-1, 2*k-1]          = (1/np)*sum((cal_data_s$X-a_hat[k]-b_hat[k]*cal_data_s$W)^2)
    A[2*k-1,(2*k-1)+1]       = (1/np)*sum(cal_data_s$W*(cal_data_s$X-a_hat[k]-b_hat[k]*cal_data_s$W)^2)
    A[(2*k-1)+1,(2*k-1)]     = A[2*k-1, (2*k-1)+1]
    A[(2*k-1)+1,(2*k-1)+1]   = (1/np)*sum((cal_data_s$W^2)*(cal_data_s$X-a_hat[k]-b_hat[k]*cal_data_s$W)^2)
  }
  
  #### 4b: Upper right entries of A (A12)
  A12    = matrix(NA, ncol=3, nrow=(2*nstud))
  
  for(k in 1:nstud){
    data_s_h     = subset(data, S==k & H==1) # Calibration controls
    Xd_s_h       = data_s_h$xd_fc
    Vd_s_h       = data_s_h$vd
    XVd_s_h      = data_s_h$xvd_fc
    BX_s_h       = data_s_h$BX
    psi1         = data_s_h$X-a_hat[k]-b_hat[k]*data_s_h$W # estimating equ 1
    psi2         = (data_s_h$X-a_hat[k]-b_hat[k]*data_s_h$W)*data_s_h$W # estimating equ 2
    
    psi_x         = Xd_s_h/(exp(BX_s_h)+1) # betax estimating equation
    psi_v         = Vd_s_h/(exp(BX_s_h)+1) # betav estimating equation
    psi_xv        = XVd_s_h/(exp(BX_s_h)+1) # betaxv estimating equation
    
    ## Betax and a,b
    A12[(2*k-1),1]   = (1/np)*sum(psi1*psi_x)
    A12[(2*k),  1]   = (1/np)*sum(psi2*psi_x)
    
    ## Betav and a,b
    A12[(2*k-1),2]   = (1/np)*sum(psi1*psi_v)
    A12[(2*k),  2]   = (1/np)*sum(psi2*psi_v)
    
    ## Betax and a,b
    A12[(2*k-1),3]   = (1/np)*sum(psi1*psi_xv)
    A12[(2*k),  3]   = (1/np)*sum(psi2*psi_xv)
  }
  
  A[1:(2*nstud), (dim_sand-2):dim_sand] = A12 # Fill appropriate piece of A matrix
  A21                                   = t(A12) # Reflect entries for A21
  A[(dim_sand-2):dim_sand,1:(2*nstud)]  = A21
  
  #### 4c: A22 entry, is 3x3
  Xd       = data$xd_fc[seq(1, nrow(data), 2)] 
  Vd       = data$vd[seq(1, nrow(data), 2)] 
  XVd      = data$xvd_fc[seq(1, nrow(data), 2)] 
  BX       = data$BX[seq(1, nrow(data), 2)] 
  
  ### A22 is 3x3; need to take cross values of all beta est equ
  A22 = matrix(NA,nrow=3,ncol=3)
  
  ## Squared terms on diagonal of A22
  A22[1,1] = (1/np)*sum( (Xd/(exp(BX)+1))^2 )
  A22[2,2] = (1/np)*sum( (Vd/(exp(BX)+1))^2 )
  A22[3,3] = (1/np)*sum( (XVd/(exp(BX)+1))^2 )
  
  ## Betax and betav of A22
  A22[1,2] = (1/np)*sum( Xd*Vd/(exp(BX)+1)^2 )
  A22[2,1] = A22[1,2]
  
  ## Betav and betaxv of A22
  A22[2,3] = (1/np)*sum( Vd*XVd/(exp(BX)+1)^2 )
  A22[3,2] = A22[2,3]
  
  ## Betax and betaxv of A22
  A22[1,3] = (1/np)*sum( Xd*XVd/(exp(BX)+1)^2 )
  A22[3,1] = A22[1,3]
  
  ## Replace A22 in main matrix of A
  A[(2*nstud+1):dim_sand,(2*nstud+1):dim_sand ] = A22
  
  ## Eliminate extra entries from A associated with calibration studies
  A = A[(2*nref+1):dim_sand,(2*nref+1):dim_sand]
  
  ###################### B MATRIX #######################################
  
  #### 5a: Block diagonal entries of B associated with the calibration studies.
  for(k in 1:nstud){
    cal_data_s               = subset(data, S==k & H==1) # Specific cal data with only controls
    B[2*k-1,2*k-1]           = -1
    B[2*k-1,(2*k-1)+1]       = (-1/np)*sum(cal_data_s$W)
    B[(2*k-1)+1,(2*k-1)]     = B[2*k-1,(2*k-1)+1]
    B[(2*k-1)+1,(2*k-1)+1]   = (-1/np)*sum(cal_data_s$W^2)
  }
  
  
  ### 5b: Off diagonal even entries of B 
  B12 = matrix(NA, ncol=3, nrow=(2*nstud))
  
  for(k in 1:nstud){
    data_s     = subset(data, S==k) # Study-specific
    Wd_s       = data_s$wd[seq(1, nrow(data_s), 2)] 
    Vd_s       = data_s$vd[seq(1, nrow(data_s), 2)] 
    WVd_s      = data_s$wvd[seq(1, nrow(data_s), 2)]
    XVd_s      = data_s$xvd_fc[seq(1, nrow(data_s), 2)] 
    BX_s       = data_s$BX[seq(1, nrow(data_s), 2)] 
    
    ## a and betax
    frac             = (-b_hat[k]*Wd_s*Vd_s*beta_hatxv*exp(BX_s))/ ((1+exp(BX_s))^2)
    B12[2*k-1, 1]    = (1/np)*sum(frac)
    
    ## b and betax
    frac1            = Wd_s/(1+exp(BX_s))
    frac2            = -b_hat[k]*Wd_s*(beta_hatx*Wd_s+beta_hatxv*WVd_s)*exp(BX_s)/ ((1+exp(BX_s))^2)
    B12[2*k, 1]      = (1/np)*(sum(frac1)+sum(frac2))
    
    ## a and betav
    frac             = (-(Vd_s^2)*beta_hatxv*exp(BX_s))/ ((1+exp(BX_s))^2)
    B12[2*k-1, 2]    = (1/np)*sum(frac)
    
    ## b and betav
    frac             = (-Vd_s*(beta_hatx*Wd_s + beta_hatxv*WVd_s)*exp(BX_s))/ ((1+exp(BX_s))^2)
    B12[2*k, 2]      = (1/np)*sum(frac)
    
    ## a and betaxv
    frac1            = Vd_s/(1+exp(BX_s))
    frac2            = -Vd_s*beta_hatxv*XVd_s*exp(BX_s) / ((1+exp(BX_s))^2)
    B12[2*k-1, 3]    = (1/np)*(sum(frac1)+sum(frac2))
    
    ## b and betaxv
    frac1            = WVd_s/(1+exp(BX_s))
    frac2            = -XVd_s*(beta_hatx*Wd_s+beta_hatxv*WVd_s)*exp(BX_s) / ((1+exp(BX_s))^2)
    B12[2*k, 3]      = (1/np)*(sum(frac1)+sum(frac2))
    
  }
  
  
  ### 5bi. Place B21 entries in main B matrix
  B[1:(2*nstud), (dim_sand-2):dim_sand] = B12
  B21                                   = t(B12) # Reflect entries for B21
  B[(dim_sand-2):dim_sand,1:(2*nstud)]  = B21
  
  
  ### 5c: Lowermost right entry of B
  B22                  = matrix(NA, nrow=3, ncol=3)
  
  ## betax_sq, betav_sq, betaxv_sq
  B22[1,1]  = (1/np)*sum(-(Xd^2)*exp(BX)/(exp(BX)+1)^2)
  B22[2,2]  = (1/np)*sum(-(Vd^2)*exp(BX)/(exp(BX)+1)^2)
  B22[3,3]  = (1/np)*sum(-(XVd^2)*exp(BX)/(exp(BX)+1)^2)
  
  ## Cross terms
  B22[1,2]  =  (1/np)*sum(-(Xd*Vd)*exp(BX)/(exp(BX)+1)^2)
  B22[2,1]  =  (1/np)*sum(-(Xd*Vd)*exp(BX)/(exp(BX)+1)^2)
  B22[1,3]  =  (1/np)*sum(-(Xd*XVd)*exp(BX)/(exp(BX)+1)^2)
  B22[3,1]  =  (1/np)*sum(-(Xd*XVd)*exp(BX)/(exp(BX)+1)^2)
  B22[3,2]  =  (1/np)*sum(-(XVd*Vd)*exp(BX)/(exp(BX)+1)^2)
  B22[2,3]  =  (1/np)*sum(-(XVd*Vd)*exp(BX)/(exp(BX)+1)^2)
  
  ## Replace B22 in main B matrix
  B[(2*nstud+1):dim_sand, (2*nstud+1):dim_sand] = B22
  
  ## Eliminate entries associated with the reference lab studies
  B = B[(2*nref+1):dim_sand,(2*nref+1):dim_sand]
  
  ##################### Final variance computations ###################
  
  #### 6a. Create V matrix and extract the variance element
  V                = solve(B)%*%A%*%t(solve(B))
  output_fc[3]     = (1/np)*V[(dim_sand-2),(dim_sand-2)] #betax
  output_fc[6]     = (1/np)*V[(dim_sand-1),(dim_sand-1)] #betav
  output_fc[9]     = (1/np)*V[dim_sand,dim_sand] #betaxv
  
  #### 99. Return appropriate output
  return(output_fc)
  
}


twostage_df_xv = function(data, X, S, H, W, Y, V, strata, covar=NA, nref, nstud){
  
  #### 0. Rename variables names given to function to ensure consistency
  data$X = data[[X]]
  data$S = data[[S]]
  data$H = data[[H]]
  data$W = data[[W]]
  data$Y = data[[Y]]
  data$V = data[[V]]
  data$strata = data[[strata]]
  data$WV = data$W*data$V
  
  #### 1. Create storage matrix for point, variance estimate, and coverage rate
  
  output_ts   = matrix(NA, ncol=9, nrow=1) # betax pb var cr, betav pb var cr, betaxv pb var cr
  colnames(output_ts) = c("betax","RRx", "VARx","betav","RRv", "VARv","betaxv","RRxv", "VARxv")
  wts_x       = matrix(NA, nrow=nstud, ncol=8)
  wts_v       = matrix(NA, nrow=nstud, ncol=14)
  wts_xv      = matrix(NA, nrow=nstud, ncol=8) 
  
  #### 2. Calibration parameter estimates
  for(k in (nref+1):nstud){
    cal_data_s    = subset(data, S==k & H==1)
    fit           = lm(X~W, data=cal_data_s)
    
    ## Betax
    wts_x[k,1]    = fit$coefficients[2]
    wts_x[k,2]    = vcov(fit)[2,2]
    
    ## Betaxv
    wts_xv[k,1]   = fit$coefficients[2]
    wts_xv[k,2]   = vcov(fit)[2,2]
    
    ## Betav
    wts_v[k,1]    = fit$coefficients[1]
    wts_v[k,2]    = vcov(fit)[1,1]
    wts_v[k,3]    = fit$coefficients[2]
    wts_v[k,4]    = vcov(fit)[2,2]
    wts_v[k,5]    = vcov(fit)[1,2]
  }
  
  #### 3. Compute betahatw, var(betahatw) for each study in noncal subset
  
  for(k in (nref+1):nstud){
    data_s        = subset(data, S==k & H==0 & G==0) #study specific, uncalibrated pairs
    if(is.na(covar)[1]){ # no covariates
      formula = as.formula(paste("Y~W+V+WV+strata(strata)",sep="+"))
    }else{ # covariates named
      formula = as.formula(paste("Y~W+V+WV+strata(strata)",paste(covar,collapse ="+"),sep="+"))
    }
    fc_fit        = clogit(formula,data=data_s)
    
    wts_x[k,3]    = fc_fit$coefficients[1] #betahatw
    wts_x[k,4]    = vcov(fc_fit)[1,1] # var betahatw
    
    wts_xv[k,3]   = fc_fit$coefficients[3] #betahatwv
    wts_xv[k,4]   = vcov(fc_fit)[3,3] # var betahatwv
    
    wts_v[k,6]    = fc_fit$coefficients[2] #betahatv
    wts_v[k,7]    = vcov(fc_fit)[2,2] # var betahatv
    wts_v[k,13]   = vcov(fc_fit)[2,3] # cov betav, betawv
  }
  
  #### 4. Additional calculations for betav: E(a/b) and Var(a/b)
  wts_v[,8] = wts_v[,1]/wts_v[,3]
  wts_v[,9] = (wts_v[,8])^2 *( (wts_v[,2]/wts_v[,1]^2) + (wts_v[,4]/wts_v[,3]^2) - 2*wts_v[,5]/(wts_v[,1]*wts_v[,3]))
  
  
  #### 5. Regression calibration to get corrected pt/var estimates
  
  ### Beta hat x calculations
  wts_x[,5]  = wts_x[,3]/wts_x[,1] # betahat xk
  wts_x[,6]  = wts_x[,4]/(wts_x[,1]^2) + (wts_x[,2]*wts_x[,3]^2)/(wts_x[,1]^4) #var betahat xk
  wts_x[,7]  = 1/wts_x[,6] # var(betahat) xk ^ (-1)
  
  ### Beta hat xv calculations
  wts_xv[,5] = wts_xv[,3]/wts_xv[,1] # betahat xvk
  wts_xv[,6] = wts_xv[,4]/(wts_xv[,1]^2) + (wts_xv[,2]*wts_xv[,3]^2)/(wts_xv[,1]^4) #var betahat xvk
  wts_xv[,7] = 1/wts_xv[,6] # var(betahat) xvk ^ (-1)
  
  ### Beta hat v calculations
  wts_v[,10] = wts_v[,6] - wts_xv[,3]*wts_v[,8] # betahat vk
  wts_v[,11] = (wts_v[,7] + wts_v[,9]*wts_xv[,4] + wts_v[,9]*wts_xv[,3]^2 + wts_xv[,4]*wts_v[,8]^2
                - 2*wts_v[,8]*wts_v[,13]) #var betahat vk
  wts_v[,12] = 1/wts_v[,11] # var(betahat) vk ^ (-1)
  
  #### 6. Compute inverse variance weights
  total_x = sum(wts_x[,7])
  wts_x[,8] = wts_x[,7]/total_x
  
  total_xv = sum(wts_xv[,7])
  wts_xv[,8] = wts_xv[,7]/total_xv
  
  total_v = sum(wts_v[,12])
  wts_v[,14] = wts_v[,12]/total_v
  
  #### 7. Compute final estimate of betahatx and var betahatxi
  
  betahatx      = wts_x[,8]%*%wts_x[,5]
  varbetahatx   = (wts_x[,8]^2)%*%wts_x[,6]
  
  betahatxv     = wts_xv[,8]%*%wts_xv[,5]
  varbetahatxv  = (wts_xv[,8]^2)%*%wts_xv[,6]
  
  betahatv     = wts_v[,14]%*%wts_v[,10]
  varbetahatv  = (wts_v[,14]^2)%*%wts_v[,11]
  
  #### 8. Point, variance coverage rate
  output_ts[1] = betahatx
  output_ts[2] = exp(betahatx)
  output_ts[3] = varbetahatx
  
  output_ts[4] = betahatv
  output_ts[5] = exp(betahatv)
  output_ts[6] = varbetahatv
  
  output_ts[7] = betahatxv
  output_ts[8] = exp(betahatxv)
  output_ts[9] = varbetahatxv
  
  if(nref!=0){ # If we have reference labs, need to adjust for reference lab observations.
    data_ref     = subset(data, H==2)
    data_ref$XV  = data_ref$X*data_ref$V
    if(is.na(covar)[1]){ # no covariates
      formula = as.formula(paste("Y~X+V+XV+strata(strata)",sep="+"))
    }else{ # covariates named
      formula = as.formula(paste("Y~X+V+XV+strata(strata)",paste(covar,collapse ="+"),sep="+"))
    }
    
    ## Model with X available and estimated variances
    true_fit   = clogit(formula,data=data_ref)
    betaxref   = true_fit$coefficients[1]
    varxref    = vcov(true_fit)[1,1]
    betavref   = true_fit$coefficients[2]
    varvref    = vcov(true_fit)[2,2]
    betaxvref  = true_fit$coefficients[3]
    varxvref   = vcov(true_fit)[3,3]
    
    ## Compute weights to combine estimates
    w_ref_x    = (1/varxref)/(1/varxref + 1/varbetahatx)
    w_loc_x    = 1-w_ref_x
    w_ref_v    = (1/varvref)/(1/varvref + 1/varbetahatv)
    w_loc_v    = 1-w_ref_x
    w_ref_xv   = (1/varxvref)/(1/varxvref + 1/varbetahatxv)
    w_loc_xv   = 1-w_ref_xv
    
    ## Compute updated point and variance estimates; place in output
    output_ts[1] = w_ref_x*betaxref + w_loc_x*betahatx
    output_ts[2] = var(output_ts[1])
    output_ts[3] = (w_ref_x^2)*varxref + (w_loc_x^2)*varbetahatx
    output_ts[4] = w_ref_v*betavref + w_loc_v*betahatv
    output_ts[5] = var(output_ts[4])
    output_ts[6] = (w_ref_v^2)*varvref + (w_loc_v^2)*varbetahatv
    output_ts[7] = w_ref_v*betavref + w_loc_v*betahatv
    output_ts[8] = var(output_ts[7])
    output_ts[9] = (w_ref_xv^2)*varxvref + (w_loc_xv^2)*varbetahatxv
    
  }
  
  #### 9. Return point estimate in output_ts
  return(output_ts) 
}


main_xv = function(data, X, S, H, W, Y, G, V, strata, nref, nstud, covar=NA){  
  if(is.na(covar[1])){ # if there are no covariates
    
    output = matrix(NA, nrow=3, ncol=9)
    colnames(output) = c("betax","ORx","Var(betax)",
                         "betav","ORv","Var(betav)",
                         "betaxv","ORxv","Var(betaxv)")
    rownames(output) = c("Internalized", "Full calibration", "Two stage")
    
    ## Computations
    output[1,] = in_df_no_cov_xv(data=data, X=X, S=S, H=H, W=W, Y=Y, G=G, V=V, strata=strata, nref=nref, nstud=nstud)
    output[2,] = fc_df_no_cov_xv(data=data, X=X, S=S, H=H, W=W, Y=Y, V=V, strata=strata, nref=nref, nstud=nstud)
    output[3,] = twostage_df_xv(data=data, X=X, S=S, H=H, W=W, Y=Y, V=V, strata=strata, nref=nref, nstud=nstud)
    
    ## Rounding
    #output[,1:2] = round(output[,1:2],8)
    #output[,3]   = round(output[,3],8)
    output = round(output,5)
    
    ## Return output
    return(output)
    
  }else{ # Covariates given
    
    output = matrix(NA, nrow=6, ncol=9)
    colnames(output) = c("betax","ORx","Var(betax)",
                         "betav","ORv","Var(betav)",
                         "betaxv","ORxv","Var(betaxv)")
    rownames(output) = c("IN unadjusted", "IN adjusted","FC unadjusted",
                         "FC adjusted","TS unadjusted", "TS adjusted")
    ## Computations
    output[1,] = in_df_no_cov_xv(data=data, X=X, S=S, H=H, W=W, Y=Y, G=G, V=V, strata=strata, nref=nref, nstud=nstud)
    output[2,] = in_df_xv(data=data, X=X, S=S, H=H, W=W, covar=covar, Y=Y, G=G, V=V, strata=strata, nref=nref, nstud=nstud)
    output[3,] = fc_df_no_cov_xv(data=data, X=X, S=S, H=H, W=W, Y=Y, V=V, strata=strata, nref=nref, nstud=nstud)
    output[4,] = fc_df_xv(data=data, X=X, S=S, H=H, W=W, V=V, strata=strata, covar=covar, Y=Y, nref=nref, nstud=nstud)
    output[5,] = twostage_df_xv(data=data, X=X, S=S, H=H, W=W, V=V, strata=strata, Y=Y, nref=nref, nstud=nstud)
    output[6,] = twostage_df_xv(data=data, X=X, S=S, H=H, W=W, V=V, covar=covar, strata=strata, Y=Y, nref=nref, nstud=nstud)
    
    ## Rounding
    #output[,1:2] = round(output[,1:2],8)
    #output[,3]   = round(output[,3],8)
    output = round(output,5)
    
    ## Return output
    return(output)
    }
  }




#######################################################
#######################################################

#######################################################
#######################################################
###   III. Running functions on actual data

## Example syntax: main_xv(data=mydata, X="ref",S="S",V="V",G="G",H="H",W="local",Y="Y",nref=2, nstud=5,
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
## V:	   is the covariate with which you plan to create an interaction term
## G:      is an indicator for whether the observation is a case paired with a control in the calibration subset (1=yes, 0=no)
## covar:  List of covariate names, syntax as shown above

########## Returns estimates for three methods: internalized, full calibration, two stage; adjusted and unadjusted
## 1st column: betax estimate (parameter associated with continuous biomarker)
## 2nd column: RR estimate associated with betax
## 3rd column: Variance estimate of estimated betax
## Adjusted: adjusts for covariates included in analysis
## Unadjusted: does not adjust for covariates


#######################################################
#######################################################


