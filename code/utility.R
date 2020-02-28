
discretize_system<-function(par,dims,dt){
  # this function creates discrete system equation from parameter. 
  # input, parameter (Ci, Ce, Rie, Rea, simga_i, sigma_e in this order)
  # dims = dimension of x, u, y
  # dt time interval (in seconds) (1/4 hour= 15 minutes)
  
  nx=dims[1]
  nu=dims[2]
  ny=dims[3]
  # create continuous system matrix A, B, C,
  AA=matrix(0,nrow=nx,ncol=nx)
  BB=matrix(0,nrow=nx,ncol=nu)
  CC=matrix(0,nrow=ny,ncol=nx)
  # process noise matrix Q
  QQ=matrix(0,nrow=nx,ncol=nx)
  # parameters
  Ci=par[1];Ce=par[2]
  Rie=par[3];Rea=par[4];
  if(!is.na(par[5])){
    sigma_i=par[5]
    sigma_e=par[6]
    skip_Q=FALSE
  }else{
    Qd=0
    skip_Q=TRUE
  }
  # assign values into continuous scale
  AA[1,1:2]=c(-1/(Rie*Ci),1/(Rie*Ci)  ) 
  AA[2,1:2]=c(1/(Rie*Ce),-1/(Rie*Ce)-1/(Rea*Ce))
  BB[1,2]=1/Ci
  BB[2,1]=1/(Rea*Ce)
  CC[1,1]=1
  if(!skip_Q){
    QQ[1,1]=sigma_i^2
    QQ[2,2]=sigma_e^2
  }
  # Discretize A, B
  a_cols <- ncol(AA) #column of A matrix
  b_cols <- ncol(BB) #column of B matrix
  mat1 <- rbind( (cbind(AA, BB) * dt ), matrix(0, nrow = b_cols, ncol = a_cols + b_cols) )
  smat <- expm::expm(mat1) #or Matrix::expm(mat1)
  # Discretized matrix Ad, Bd
  Ad <- as.matrix(smat[1:a_cols, 1:a_cols, drop = FALSE])
  Bd <- as.matrix(smat[1:a_cols, (a_cols + 1):(a_cols + b_cols), drop = FALSE])
  if(!skip_Q){
    # Discretize Q
    QQ_temp=rbind(cbind(-AA,QQ),cbind(AA*0,t(AA)))*dt
    QQ_temp2=expm::expm(QQ_temp)
    # Discertized matrix Qd
    Qd=Ad%*%QQ_temp2[1:nx,((nx+1):(nx+nx))]
  }
  return(list(Ad=Ad,Bd=Bd,Cd=CC,Qd=Qd))
}


nstep<-function(Ad,Bd,Cd,x0,u_mat,...){
  # deterministic nstep ahead prediction
  t_K=ncol(u_mat) #total time t_K 
  nx=ncol(Ad)
  nu=ncol(Bd)
  ny=nrow(Cd)
  
  x_mat=matrix(0,nrow=nx,ncol=t_K)
  y_mat=matrix(0,nrow=ny,ncol=t_K)
  x_mat[,1]=x0 #initial states
  
  for (t_k in 2:t_K){
    x_mat[,t_k]=Ad%*%x_mat[,t_k-1]+Bd%*%u_mat[,t_k-1] #f_d(x,u,theta)
    y_mat[,t_k-1]=Cd%*%x_mat[,t_k-1] #g_d(x,u,theta)
  }
  y_mat[,t_K]=Cd%*%x_mat[,t_K] #g_d(x,u,theta)
  return(list(x_pred=x_mat,y_pred=y_mat))
}


create_data_set<-function(dsys,u_mat,x0,dims,dt,seed_num=1234,stochastic=FALSE){
  # create dataset
  nx=dims[1] 
  nu=dims[2]
  ny=dims[3]
  t_K=dim(u_mat)[2]
  
  Ad=dsys$Ad
  Bd=dsys$Bd
  Cd=dsys$Cd
  Qd=dsys$Qd
  
  sdy=0.25 #fixed value
  
  x_mat=matrix(0,nrow=nx,ncol=t_K)
  y_mat=matrix(0,nrow=ny,ncol=t_K)
  x_mat[,1]=x0
  
  if(stochastic){
    # stochastic data generation with noise parameter
    for (t_k in 2:t_K){
      set.seed(seed_num+t_k+200)
      x_mat[,t_k]=Ad%*%x_mat[,t_k-1]+Bd%*%u_mat[,t_k-1]+matrix(mvtnorm::rmvnorm(n=1,sigma=Qd),ncol=1)
      y_mat[,t_k-1]=Cd%*%x_mat[,t_k-1]+rnorm(ny,mean=0,sd=sdy)
      
    }
    set.seed(seed_num+t_K)
    y_mat[,t_K]=Cd%*%x_mat[,t_K]+rnorm(ny,mean=0,sd=sdy)
    print('stochastic')
    return(list(x_mat=x_mat,y_mat=y_mat))
    
  }else{
    for (t_k in 2:t_K){
      x_mat[,t_k]=Ad%*%x_mat[,t_k-1]+Bd%*%u_mat[,t_k-1]
      y_mat[,t_k-1]=Cd%*%x_mat[,t_k-1]
    }
    y_mat[,t_K]=Cd%*%x_mat[,t_K]
    print('non-stochastic')
    return(list(x_mat=x_mat,y_mat=y_mat))
  }
  
}


get_rmse<-function(par,dt,u_mat,y_mat,dims,output="cost",normalize=T,L_val=NA){
  # root mean squared error of n-step prediction
  if(normalize){
    if(is.na(L_val[1])){
      stop("L_val is NA")
    }
    par=(par)*L_val
  }
  par_temp=c(par[-c(1:2)])
  dsys=discretize_system(par=par_temp,dt=dt,dims=dims)

  x0=matrix(par[1:2],ncol=1)
  res=nstep(A=dsys$Ad,B=dsys$Bd,C=dsys$Cd,x0=x0,u_mat=u_mat)
  ## library(kfsang) Don't use this library. This package is written by me to have cpp function of nstep function for the speed. Pacakage is not ready-to-publish format.
  #res=kfsang::nstep_cpp(A=dsys$Ad,B=dsys$Bd,C=dsys$Cd,D=matrix(0,ncol=2),x0=x0,u_mat=u_mat)
  if(output=="cost"){
    return(sum(((as.vector(res$y_pred))-as.vector(y_mat))^2))
  }else{
    return(res)
  }
}

update_system_onestep<-function(par,dims,dt,x,u){
  # discrete system and one step ahead prediction
  par[5]=NA # to skip Qd calculation
  dsys=discretize_system(par,dims,dt)
  x=dsys$Ad%*%matrix(x,ncol=1)+dsys$Bd%*%matrix(u,ncol=1)
  y=dsys$Cd%*%matrix(x,ncol=1)#+DD%*%matrix(u,ncol=1)
  return(list(x=x,y=y))
}

discretize_Q<-function(par,dims,dt){
  # this function creates discrete system equation from parameter. 
  # input, parameter (Ci, Ce, Rie, Rea, simga_i, sigma_e in this order)
  # dims = dimension of x, u, y
  # dt time interval (in hour) (1/12 hour= 5 minutes)
  
  nx=dims[1]
  nu=dims[2]
  ny=dims[3]
  # create continuous system matrix A, B, C,
  AA=matrix(0,nrow=nx,ncol=nx)
  BB=matrix(0,nrow=nx,ncol=nu)
  CC=matrix(0,nrow=ny,ncol=nx)
  # process noise matrix Q
  QQ=matrix(0,nrow=nx,ncol=nx)
  
  # parameters
  Ci=par[1];Ce=par[2]
  Rie=par[3];Rea=par[4];
  sigma_i=par[5]
  sigma_e=par[6]
  
  # assign values into continuous scale
  AA[1,1:2]=c(-1/(Rie*Ci),1/(Rie*Ci)  ) # Ti, Te, Th 
  AA[2,1:2]=c(1/(Rie*Ce),-1/(Rie*Ce)-1/(Rea*Ce))
  BB[1,2]=1/Ci
  BB[2,1]=1/(Rea*Ce)

  CC[1,1]=1
  QQ[1,1]=sigma_i^2
  QQ[2,2]=sigma_e^2
  
  QQ_temp=matrix(0,nrow=nx+nx,ncol=nx+nx)
  QQ_temp[1:nx,1:nx]=-1*AA*dt
  QQ_temp[1:nx,(nx+1):(2*nx)]=QQ*dt
  
  QQ_temp[(nx+1):(2*nx),(nx+1):(2*nx)]=t(AA)*dt
  
  QQ_temp2=expm::expm(QQ_temp)
   
  Qd=QQ_temp2[(nx+1):(2*nx),(nx+1):(2*nx)]%*%QQ_temp2[1:nx,(nx+1):(2*nx)]
  Qd=(Qd+t(Qd))/2 #make it symmetric
  if(mean(abs(Qd))>100000|is.nan(mean(Qd))){
    # Qd sometimes gets unrealistic values during the matrix exponential.
    # The particles that gives this unrealistic value need to be dropped out when the filter moves.
    # However, it results in numerical error during the calculation, so we force those values to smaller numbers.
    # Through this approach, we can drop out unrealistic particles while prevent filter from stopping due to numerical error.
    Qd=diag(rep(20^2,nx))
  }
  return(list(Qd=Qd))
}



suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(lubridate))
suppressPackageStartupMessages(library(forcats))
suppressPackageStartupMessages(library(rlang))
suppressPackageStartupMessages(library(lhs))
suppressPackageStartupMessages(library(kfsang))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(mvtnorm))


r_theta_recover<-function(r_theta,th=1e-11,idx=5:6){
  # to reduce numerical error, negative sigma values are forced to positive numbers
  # if we use small number close to zero, it also can give numerical error.
  r_theta[idx,]=abs(r_theta[idx,])
  # replace negative parameters to very small number
  r_theta[r_theta<=0]<-th
  return(r_theta)
}

lw_ss_filter<-function(inputs,NP,pii,x,theta,delta=0.90,seed=13234,...){
  # Liu-West filter
  set.seed(seed)
  # input data
  y_mat=matrix(inputs$y_mat,nrow=1)
  dt=inputs$dt
  u_mat=(inputs$u_mat)
  par_names=inputs$par

  L_val=inputs$L_val
  xs_idx=1 # sensor measurement index in x verctor.
  sdy=0.25 # sigma_y (measurement noise)
  th_value=1e-9 #threshold that used in parameter unnormalization.
  # Detect NA in u
  u_na=t(u_mat)
  u_na_vec=apply(u_na,MARGIN=c(1),anyNA)
  # Detect NA in y
  y_na_vec=is.na(y_mat[1,])
  
  # dimensions
  nx=dim(x)[1];nu=dim(u_mat)[1]
  ny=dim(y_mat)[1] 
  dims=c(nx,nu,ny) #dimension of x, u, y
  t_K = dim(y_mat)[2] #total time
  
  num_par=length(par_names)
  num_state=nx

  sdti_idx=which(str_detect(par_names,"sdti")) # noise sigma_ti  
  sdte_idx=which(str_detect(par_names,"sdte")) # noise sigma_te

  sd_idx=c(sdti_idx,sdte_idx)
  xs = array(0,dim=c(num_state,NP,t_K)) # store filtered states
  thetas=array(0,dim=c(num_par,NP,t_K))  # store normalized parameters
  r_thetas=array(0,dim=c(num_par,NP,t_K)) #store unnormalized parameters
  
  thetas[,,1]=theta #normalized parameter
  xs[,,1]=x #initial states
  x_new=x*0 # initialize x_new's array space
  
  # Liu-Wst filter tuning parameter a and h
  a=(3*delta-1)/(2*delta)
  h=sqrt(1-a^2)
  
  r_theta=(theta+0.5)*L_val #unnormalize parameter (vectorized calculation)
  r_theta=r_theta_recover(r_theta,th=th_value,idx=sd_idx) 
  r_thetas[,,1]=r_theta
  
  pii=pii #initial weight of particles
  chol_fail=0 #count faulier in cholesky decomposition
  ll=0 #log-likelihood
  ll_prev=0 #previous log-likelihood to detect wrong filter update
  ll_vec=rep(0,t_K)  #store log-likelihood
  ll_run=TRUE # stop calculation when log-likelihood is unrealistic
  
  # parallel processing
  no_cores <- parallel::detectCores() - 1
    # Initiate cluster
  cl <- parallel::makeCluster(no_cores)
  envir=environment() 
  
  for (t_k in 2:t_K){
    #seed=1345
    set.seed(seed*t_k+1)
    if(u_na_vec[t_k-1]==FALSE&y_na_vec[t_k]==FALSE){
      # refer the algorithm shown in Table 2 in the paper.
      theta_prev=theta #store theta (not to update when the parameter is not involved in the process)
      # Line 5 
      thetabar = Rfast::rowmeans(theta)  
      # Line 6 
      cov_theta=(cov.wt(Rfast::transpose(theta),wt=pii)$cov)*(NP-1)/NP #weighted cov (not sample variance)
      cov_theta=(cov_theta+t(cov_theta))/2 #make symmetric
      Vtheta=h^2*cov_theta #h^2 is multiplied (see line 12)
      # prevent particle degeneration. see 5. Discussion in the paper
      diag(Vtheta)=if_else(diag(Vtheta)<1e-8,1e-8,0)+diag(Vtheta) 
      
      # Line 7
      # Parallel apply see below link
      # https://stackoverflow.com/questions/19467133/performance-of-clusterapply-deteriorates-when-called-inside-a-function
      clusterExport(cl, list("update_system_onestep","discretize_system","r_theta","dims","dt","x","u_mat","t_k"),envir=envir) #line8
      f00=function(ix){
        update_system_onestep(par =as.vector(r_theta[,ix]),dims = dims,dt = dt,x =matrix(x[,ix],ncol=1),u=matrix(u_mat[,t_k-1] ,ncol=1)  )$x
      }
      environment(f00)<-.GlobalEnv
      mu_x_pred=(parSapply(cl,1:NP,f00 )) 
      # Line 8
      mu_theta_pred=a*theta+(1-a)*thetabar
      # print intermediate results
      if (t_k%%5==0){
        cat("mu_theta_t_k|t_k-1")
        cat("\n")
        cat(Rfast::rowmeans((mu_theta_pred+0.5)*L_val))
        cat("\n")
        cat("x")
        cat("\n")
        cat((Rfast::rowmeans(x)))
        cat("\n")
        cat((diag(Vtheta)))
        cat("\n")
      }
      r_mu_theta_pred=(mu_theta_pred+0.5)*L_val #unnormalize parameter (vectorized calculation)
      r_mu_theta_pred=r_theta_recover(r_mu_theta_pred,th=th_value,idx=sd_idx) # remove negative values
      # line 9
      w = as.vector(dnorm(y_mat[1,t_k],mu_x_pred[xs_idx,],(sdy),log=T))+log(pii) 
      w=exp(w-max(w)) #prevent numerical overflow
      pii=w/sum(w) 
      
      #resampling
      m = sample(1:NP,size=NP,replace=TRUE,prob=pii) 
      m=sort(m)
      
      # 
      x_m = x[,m]    # used for line 13
      mu_theta_pred_m = mu_theta_pred[,m]  #used for line 12
      
      mu_x_pred_m=mu_x_pred[,m]  # used for line 15

      # line 12. trycatch in case of cholesky decomposition failure
      theta_new=tryCatch(
        (mu_theta_pred_m[1:num_par,]+Rfast::transpose(Rfast::rmvnorm(n=NP, mu=rep(0,num_par), sigma=Vtheta))),
        error=function(e){
          chol_fail=chol_fail+1 #count cholesky decomposition fails
          print("chol fail")
          mu_theta_pred_m[1:num_par,]+Rfast::transpose(mvtnorm::rmvnorm(NP,mean=rep(0,num_par),sigma = Vtheta))
        })
      
      r_theta_new=(theta_new+0.5)*L_val #unnormalize parameter (vectorized calculation)
      r_theta_new=r_theta_recover(r_theta_new,th=th_value,idx=sd_idx) # remove negative values
      
      # line 13
      clusterExport(cl, list("update_system_onestep","discretize_system","r_theta_new","dims","dt","x_m","u_mat","t_k"),envir=envir) #line8
      f01=function(ix){
        update_system_onestep(par =as.vector(r_theta_new[,ix]),dims = dims,dt = dt,x =matrix(x_m[,ix],ncol=1),u=matrix(u_mat[,t_k-1] ,ncol=1)  )$x
      }
      environment(f01)<-.GlobalEnv
      mu_x_new=(parSapply(cl,1:NP,f01 )) 
      
      # line 14
      set.seed(seed*t_k+2)
      clusterExport(cl, list("rmvnorm","mu_x_new","discretize_Q","r_theta_new","dt","dims"),envir=envir) #line8
      f02=function(ix){
        #update_system_lw_onestep_cpp(par =as.vector(r_theta_new[,ix]),dims = dims,dt = dt,x =matrix(x_m[,ix],ncol=1),u=matrix(u_mat[,t-1] ,ncol=1)  )$x
        rmvnorm(1,mu_x_new[,ix],sigma=discretize_Q(par=r_theta_new[,ix],dt=dt,dims=dims)$Qd)
      }
      environment(f02)<-.GlobalEnv
      x_new=(parSapply(cl,1:NP,f02 )) 
      
      # line 15
      w = as.vector(dnorm(y_mat[1,t_k],x_new[xs_idx,],(sdy),log=TRUE))-as.vector(dnorm(y_mat[1,t_k],mu_x_pred_m[xs_idx,],(sdy),log=TRUE)) 
      
      # remove extreme value
      # to remove extreme outlier particles 
      w_norm<-scale(w)
      w_outlier=which(abs(w_norm)>=3.090232) # qnorm(0.999)
      w[w_outlier]=mean(w[-w_outlier])
  
      ct=max(w)
      wstar=exp(w-ct)
      ll=ll_vec[t_k-1]+(-log(NP)+ct+log(sum(wstar))) #calculate log-likelihood
      
      # final weight of each particle
      pii=wstar/sum(wstar)
      
      # resampling 
      kk = sample(1:NP,size=NP,replace=TRUE,prob=pii) # o in algorithm
      kk=sort(kk)
      pii=rep(1/NP,NP) #since we do resampling, make equal weight again
      x= (x_new)[,kk]
      # store values 
      xs[,,t_k] =x #x_new
      theta_new=theta_new[,kk]
      
      r_theta_new=(theta_new+0.5)*L_val #unnormalize parameter (vectorized calculation)
      r_theta_new=r_theta_recover(r_theta_new,th=th_value,idx=sd_idx) # remove negative values
      
      theta= (theta_new)
      thetas[,,t_k] =theta
      r_theta=r_theta_new
      r_thetas[,,t_k] =r_theta_new
      cat(paste0("time is ",t_k," ll is ",ll,"\n"))
    }else{
      # skip training when there is NA data
      ll=ll_vec[t_k-1]
      ll_run=FALSE
      xs[,,t_k] =x
      thetas[,,t_k] =theta
      r_thetas[,,t_k] =r_theta
      cat(paste0("time is ",t_k," ll is ",ll,"\n"))
    }
    
    ll_vec[t_k]=ll
  }
  # stop parallel cluster
  parallel::stopCluster(cl)
  # return the results
  return(list(thetas=thetas,r_thetas=r_thetas,xs=xs,NP=NP,pii=pii,delta=delta,ll=ll_vec,L_val=L_val)) # lw_mean=lw_mean,lw_sd=lw_sd
  
}





