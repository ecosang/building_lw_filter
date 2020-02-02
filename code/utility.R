
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
  #res=kfsang::nstep_cpp(A=dsys$Ad,B=dsys$Bd,C=dsys$Cd,D=matrix(0,ncol=2),x0=x0,u_mat=u_mat)
  if(output=="cost"){
    return(sum((((as.vector(res$y_pred))-as.vector(y_mat))/c(1))^2))
  }else{
    return(res)
  }
}











update_system_onestep<-function(par,dims,dt,x,u){
  # discrete system and one step ahead prediction
  
  par[5]=NA # to skip Q calculation
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
  #QQ_temp2
  
  Qd=QQ_temp2[(nx+1):(2*nx),(nx+1):(2*nx)]%*%QQ_temp2[1:nx,(nx+1):(2*nx)]
  Qd=(Qd+t(Qd))/2
  if(mean(abs(Qd))>100000|is.nan(mean(Qd))){
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
  
  r_theta[idx,]=abs(r_theta[idx,])
  r_theta[r_theta<=0]<-th
  return(r_theta)
}

lw_ss_filter<-function(inputs,NP,pii,x,theta,delta=0.90,seed=13234,...){
  
  set.seed(seed)
  
  # input data
  y_mat=matrix(inputs$y_mat,nrow=1)
  dt=inputs$dt
  u_mat=(inputs$u_mat)
  par_names=inputs$par
  # Line xx,  
  # u_update_temp=t(u_mat)[,-1]
  # u_update=u_update_temp%>%as_tibble()%>%mutate_all(list(~if_else(.==0,F,T)))
  
  L_val=inputs$L_val
  
  xs_idx=1 #sensor measure
  sdy=0.25
  # Detect NA in u
  u_na=t(u_mat)
  u_na_vec=apply(u_na,MARGIN=c(1),anyNA)
  # Detect NA in y
  y_na_vec=is.na(y_mat[1,])
  
  # dimensions
  nx=dim(x)[1];nu=dim(u_mat)[1]
  ny=dim(y_mat)[1]
  dims=c(nx,nu,ny)
  t_K = dim(y_mat)[2]
  
  
  num_par=length(par_names)
  num_state=nx
  
  sdte_idx=which(str_detect(par_names,"sdte"))
  sdti_idx=which(str_detect(par_names,"sdti"))
  sd_idx=c(sdte_idx,sdti_idx)
  xs = array(0,dim=c(num_state,NP,t_K))
  thetas=array(0,dim=c(num_par,NP,t_K))
  r_thetas=array(0,dim=c(num_par,NP,t_K))
  
  thetas[,,1]=theta
  xs[,,1]=x
  
  r_theta=(theta+0.5)*L_val
  r_theta=r_theta_recover(r_theta,th=1e-9,idx=sd_idx)
  
  r_thetas[,,1]=r_theta
  
  
  pii=pii
  chol_fail=0
  ll=0
  ll_prev=0
  ll_vec=rep(0,t_K)
  
  no_cores <- parallel::detectCores() - 1
  
  # Initiate cluster
  cl <- parallel::makeCluster(no_cores)
  envir=environment()
  ll_run=TRUE
  
  x_new=x*0 # space declaration
  
  a=(3*delta-1)/(2*delta)
  h=sqrt(1-a^2)
  
  #t_k=2
  for (t_k in 2:t_K){
    #seed=1345
    set.seed(seed*t_k+1)
    if(u_na_vec[t_k-1]==FALSE&y_na_vec[t_k]==FALSE){
      
      theta_prev=theta #store theta (not to update when the parameter is not involved in the process)
      cov_theta=(cov.wt(Rfast::transpose(theta),wt=pii)$cov)*(NP-1)/NP #weighted cov (not sample variance)
      cov_theta=(cov_theta+t(cov_theta))/2 #make a symmetry
      
      thetabar = Rfast::rowmeans(theta) 
      Vtheta=h^2*cov_theta
      diag(Vtheta)=if_else(diag(Vtheta)<1e-7,1e-7,0)+diag(Vtheta)
      
      #update_system_lw_onestep_cpp_sx(par =as.vector(r_theta[t-1,]),dims = dims,dt = dt,x =matrix(r_x[t-1,],ncol=1),u=matrix(u_mat[t-1,] ,ncol=1)  )
      # https://stackoverflow.com/questions/19467133/performance-of-clusterapply-deteriorates-when-called-inside-a-function
      
      clusterExport(cl, list("update_system_onestep","discretize_system","r_theta","dims","dt","x","u_mat","t_k"),envir=envir) #line8
      f00=function(ix){
        update_system_onestep(par =as.vector(r_theta[,ix]),dims = dims,dt = dt,x =matrix(x[,ix],ncol=1),u=matrix(u_mat[,t_k-1] ,ncol=1)  )$x
      }
  
      #update_system_onestep(par =as.vector(r_theta[,ix]),dims = dims,dt = dt,x =matrix(x[,ix],ncol=1),u=matrix(u_mat[,t_k-1] ,ncol=1)  )$x
      
      
      environment(f00)<-.GlobalEnv
      mu_x_pred=(parSapply(cl,1:NP,f00 )) # line 12 of nimble
      
      mu_theta_pred=a*theta+(1-a)*thetabar
      
      # dim(theta)
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
      
      r_mu_theta_pred=(mu_theta_pred+0.5)*L_val
      r_mu_theta_pred=r_theta_recover(r_mu_theta_pred,th=1e-9,idx=sd_idx)
      # /sqrt(3) if y is mean value.
      
      w = as.vector(dnorm(y_mat[1,t_k],mu_x_pred[xs_idx,],(sdy),log=T))+log(pii) # w+1=wp(y+1|g(x),m(t)) #line 10, 11 in NIMBLE  # r_mu_theta_pred[,20]
      w=exp(w-max(w))
      pii=w/sum(w)
      
      #resampling
      m = sample(1:NP,size=NP,replace=TRUE,prob=pii) 
      
      m=sort(m)
      
      #line 12 first part. x_tilde_t-1, theta_tilde_t-1
      x_m = x[,m]         # resample   # line 12.  x_tilde_{t-1} from x_{t-1}  
      mu_theta_pred_m = mu_theta_pred[,m] # this is theta_tilde_{t|t-1} used in line 15 and line 13 (type) whichi is m(theta_tilde_t))
      
      mu_x_pred_m=mu_x_pred[,m]
      theta_m = theta[,m]
      
      r_theta_m=(theta_m+0.5)*L_val
      r_theta_m=r_theta_recover(r_theta_m,th=1e-9,idx=sd_idx)
      #r_theta_m[r_theta_m<=0]<-1e-7
      
      theta_new=tryCatch(
        (mu_theta_pred_m[1:num_par,]+Rfast::transpose(Rfast::rmvnorm(n=NP, mu=rep(0,num_par), sigma=Vtheta))),
        #(mu_theta_pred_m[1:num_dpar,]+Rfast::transpose(mvnfast::rmvn(n=NP, mu=rep(0,num_dpar), sigma=chol(dVtheta*10^8)/(10^4),isChol=TRUE))),
        #(mu_theta_pred_m[1:num_dpar,]+(Rfast::transpose((chol(dVtheta*1000000)/1000))%*%matrix(rnorm(num_dpar*NP),nrow=num_dpar))), #mu+A x   chol(cov) gives AT
        error=function(e){
          chol_fail=chol_fail+1 #count cholesky decomposition fails
          print("chol fail")
          mu_theta_pred_m[1:num_par,]+Rfast::transpose(mvtnorm::rmvnorm(NP,mean=rep(0,num_par),sigma = Vtheta))
          #(mu_theta_pred_m[1:num_dpar,]+armarmv(n=NP,mu=rep(0,num_dpar), sigma=dVtheta))
        })
      
      
      r_theta_new=(theta_new+0.5)*L_val
      r_theta_new=r_theta_recover(r_theta_new,th=1e-10,idx=sd_idx)
      
      clusterExport(cl, list("update_system_onestep","discretize_system","r_theta_new","dims","dt","x_m","u_mat","t_k"),envir=envir) #line8
      f01=function(ix){
        update_system_onestep(par =as.vector(r_theta_new[,ix]),dims = dims,dt = dt,x =matrix(x_m[,ix],ncol=1),u=matrix(u_mat[,t_k-1] ,ncol=1)  )$x
      }
      environment(f01)<-.GlobalEnv
      
      mu_x_new=(parSapply(cl,1:NP,f01 )) # line 12 of nimble
      set.seed(seed*t_k+2)
      
      clusterExport(cl, list("rmvnorm","mu_x_new","discretize_Q","r_theta_new","dt","dims"),envir=envir) #line8
      f02=function(ix){
        #update_system_lw_onestep_cpp(par =as.vector(r_theta_new[,ix]),dims = dims,dt = dt,x =matrix(x_m[,ix],ncol=1),u=matrix(u_mat[,t-1] ,ncol=1)  )$x
        rmvnorm(1,mu_x_new[,ix],sigma=discretize_Q(par=r_theta_new[,ix],dt=dt,dims=dims)$Qd)
      }
      environment(f02)<-.GlobalEnv
      #mu_xy_new=t(parLapply(cl,1:NP,f01 )) # line 12 of nimble
      x_new=(parSapply(cl,1:NP,f02 )) # line 12 of nimble
      
      w = as.vector(dnorm(y_mat[1,t_k],x_new[xs_idx,],(sdy ),log=TRUE))-as.vector(dnorm(y_mat[1,t_k],mu_x_pred_m[xs_idx,],(sdy ),log=TRUE)) 
      
      # remove extreme value
      w_norm<-scale(w)
      w_outlier=which(abs(w_norm)>=3.090232) # qnorm(0.999)
      w[w_outlier]=mean(w[-w_outlier])
      
      # r_theta_new[,20]         # r_mu_theta_pred_k[,20]
      ct=max(w)
      wstar=exp(w-ct)
      ll=ll_vec[t_k-1]+(-log(NP)+ct+log(sum(wstar)))
      
      
      pii=wstar/sum(wstar)
      
      kk = sample(1:NP,size=NP,replace=TRUE,prob=pii)
      kk=sort(kk)
      pii=rep(1/NP,NP)
      x= (x_new)[,kk]
      xs[,,t_k] =x #x_new
      theta_new=theta_new[,kk]
      
      r_theta_new=(theta_new+0.5)*L_val
      r_theta_new=r_theta_recover(r_theta_new,th=1e-10,idx=sd_idx)
      
      theta= (theta_new)
      thetas[,,t_k] =theta
      r_theta=r_theta_new
      r_thetas[,,t_k] =r_theta_new
      
      cat(paste0("time is ",t_k," ll is ",ll,"\n"))
    }else{
      ll=ll_vec[t_k-1]
      ll_run=FALSE
      xs[,,t_k] =x
      thetas[,,t_k] =theta
      r_thetas[,,t_k] =r_theta
      
      cat(paste0("time is ",t_k," ll is ",ll,"\n"))
      
    }
    
    
    ll_vec[t_k]=ll
    
  }
  
  parallel::stopCluster(cl)
  
  return(list(thetas=thetas,r_thetas=r_thetas,xs=xs,NP=NP,pii=pii,delta=delta,ll=ll_vec,L_val=L_val)) # lw_mean=lw_mean,lw_sd=lw_sd
  
}






