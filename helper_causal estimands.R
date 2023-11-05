# helper functions for paper 1

# ===================== mathematical derivation version ===================== #
# "calculate.effect_singlet" and "calculate.effect_allt" returns causal effects without CI
#      only recommended for quick computation for confirmation purpose
# "calculate.effect_allt_withCI" returns causal effects with CI and is recommanded to be used
#      tx is a binary vector/scalar, specified forward as (...,x_{t-2},x_{t-1},x_t)
# change name from "calculate.counterfactual_outcome_singlet"
# Updated on 11/03/2023
calculate.effect_singlet=function(t,tx,y_coeffi_table,c_coeffi_table,printFlag=T){
  if(printFlag){
    cat(blue(" =================================================================================================== \n"))
    cat(blue("This function calculates contemporaneous effect and q-lag effects directly from mathematical expressions of coefficients.\n"))
    cat(blue("  Caution: it only applied for recent treatment up to length 5 (or 4 time points ago)\n"))
    cat(red ("  The formulat for outcome regression has to be: y_t = intercept + y_{t-1} + x_t + x_{t-1} + c_t\n"))
    cat(red ("               for covariate regression has to be: c_t = intercept + c_{t-1} + x_{t-1} + c_{t-1}\n"))  
    cat(blue("  Inputs include: observational time t,\n"))
    cat(blue("                  a vector/scalar of recent tx=(...,x_{t-2},x_{t-1},x_t),\n"))
    cat(blue("                  a table of outcome regression coefficients of all time,\n"))
    cat(blue("                  a table of covariate regression coefficients of all time.\n"))
    cat(blue("  Output is: contemporaneous or q-lag effect at time t.\n"))
    cat(blue(" =================================================================================================== \n"))
  }
  if(!all(is.numeric(t) & t %% 1 == 0 & t>0 & t<=nrow(y_coeffi_table) & t<=nrow(c_coeffi_table))){stop("Chosen t a propriate positive integer between [1,T].")}
  if(!all(is.data.frame(y_coeffi_table) & is.data.frame(c_coeffi_table))){stop("y_coeffi_table and c_coeffi_table should be data.frame.")}
  if(!all(c("(Intercept)","y_1","x","x_1","c") %in% colnames(y_coeffi_table))){stop("Column names for y_coeffi_table should be: (Intercept), y_1, x, x_1, c.")}
  if(!all(c("(Intercept)","c_1","x_1","y_1")  %in% colnames(c_coeffi_table))){stop("Column names for c_coeffi_table should be: (Intercept), c_1, x_1, y_1.")}
  
  # extract coefficients to be more easy to access
  rho_y=y_coeffi_table$y_1
  beta_1=y_coeffi_table$x
  beta_2=y_coeffi_table$x_1
  beta_c=y_coeffi_table$c
  
  rho_c=c_coeffi_table$c_1
  mu_x=c_coeffi_table$x_1
  mu_y=c_coeffi_table$y_1
  
  if(length(tx)==1){
    if(printFlag){cat("Compute contemporanous effect.")}
    return(beta_1[t]*tx)
  }else if(length(tx)==2){
    if(t>1){
      if(printFlag){cat("Compute 1-lag effect (treatment of length 2).")}
      result=beta_1[t]*tx[2]+(beta_2[t]+beta_c[t]*mu_x[t]+rho_y[t]*beta_1[t-1]+beta_c[t]*mu_y[t]*beta_1[t-1])*tx[1]
      return(result)
    }else{return(NA)}
  }else if(length(tx)==3){
    if(t>2){
      if(printFlag){cat("Compute 2-lag effect (treatment of length 3).")}
      term_a=beta_c[t]*rho_c[t]+rho_y[t]*beta_c[t-1]+beta_c[t]*mu_y[t]*beta_c[t-1]
      term_b=(rho_y[t]+beta_c[t]*mu_y[t])*rho_y[t-1]+term_a*mu_y[t-1]
      result=beta_1[t]*tx[3]+
        (beta_2[t]+beta_c[t]*mu_x[t]+(rho_y[t]+beta_c[t]*mu_y[t])*beta_1[t-1])*tx[2]+
        ((rho_y[t]+beta_c[t]*mu_y[t])*beta_2[t-1]+term_a*mu_x[t-1]+term_b*beta_1[t-2])*tx[1]
      return(result) 
    }else{return(NA)}
  }else if(length(tx)==4){
    if(t>3){
      if(printFlag){cat("Compute 3-lag effect (treatment of length 4).")}
      term_a=beta_c[t]*rho_c[t]+rho_y[t]*beta_c[t-1]+beta_c[t]*mu_y[t]*beta_c[t-1]
      term_b=(rho_y[t]+beta_c[t]*mu_y[t])*rho_y[t-1]+term_a*mu_y[t-1]
      result=beta_1[t]*tx[4]+
        (beta_2[t]+beta_c[t]*mu_x[t]+(rho_y[t]+beta_c[t]*mu_y[t])*beta_1[t-1])*tx[3]+
        ((rho_y[t]+beta_c[t]*mu_y[t])*beta_2[t-1]+term_a*mu_x[t-1]+term_b*beta_1[t-2])*tx[2]+
        (term_b*beta_2[t-2]+term_a*rho_c[t-1]*mu_x[t-2]+term_b*beta_c[t-2]*mu_x[t-2]+(term_a*rho_c[t-1]+term_b*beta_c[t-2])*mu_y[t-2]*beta_1[t-3]+term_b*rho_y[t-2]*beta_1[t-3])*tx[1]
      return(result) 
    }else{return(NA)}
  }else if(length(tx)==5){
    if(t>4){
      if(printFlag){cat("Compute 4-lag effect (treatment of length 5).")}
      term_a=beta_c[t]*rho_c[t]+rho_y[t]*beta_c[t-1]+beta_c[t]*mu_y[t]*beta_c[t-1]
      term_b=(rho_y[t]+beta_c[t]*mu_y[t])*rho_y[t-1]+term_a*mu_y[t-1]
      term_c=term_a*rho_c[t-1]+term_b*beta_c[t-2]
      term_d=term_c*mu_y[t-2]+term_b*rho_y[t-2]
      term_e=term_c*rho_c[t-2]+term_d*beta_c[t-2]
      result=beta_1[t]*tx[5]+
        (beta_2[t]+beta_c[t]*mu_x[t]+(rho_y[t]+beta_c[t]*mu_y[t])*beta_1[t-1])*tx[4]+
        ((rho_y[t]+beta_c[t]*mu_y[t])*beta_2[t-1]+term_a*mu_x[t-1]+term_b*beta_1[t-2])*tx[3]+
        (term_b*beta_2[t-2]+term_c*mu_x[t-2]+term_d*beta_1[t-3])*tx[2]+
        (term_d*beta_2[t-3]+term_e*mu_x[t-3]+(term_e*mu_y[t-4]+term_d*rho_y[t-3])*beta_1[t-4])*tx[1]
      return(result)
    }else{return(NA)}
  }else{stop("This function temporarily takes tx of length from 1 to 5.")}
}
# change name from "calculate.effect_allt"
# Updated on 11/03/2023
calculate.effect_allt=function(tx,y_coeffi_table,c_coeffi_table,printFlag=T){
  if(printFlag){
    cat(blue(" =================================================================================================== \n"))
    cat(blue("This function calculates contemporaneous effect and q-lag effects directly from mathematical expressions of coefficients.\n"))
    cat(blue("  Caution: it only applied for recent treatment up to length 5 (or 4 time points ago)\n"))
    cat(red ("  The formulat for outcome regression has to be: y_t = intercept + y_{t-1} + x_t + x_{t-1} + c_t\n"))
    cat(red ("               for covariate regression has to be: c_t = intercept + c_{t-1} + x_{t-1} + c_{t-1}\n"))  
    cat(blue("  Inputs include: a vector/scalar of recent tx=(...,x_{t-2},x_{t-1},x_t),\n"))
    cat(blue("                  a table of outcome regression coefficients of all time,\n"))
    cat(blue("                  a table of covariate regression coefficients of all time.\n"))
    cat(blue("  Output is: estimate of a q-lag effect for all time t.\n"))
    cat(blue(" =================================================================================================== \n"))
  }
  if(!all(is.data.frame(y_coeffi_table) & is.data.frame(c_coeffi_table))){stop("y_coeffi_table and c_coeffi_table should be data.frame.")}
  if(!all(c("(Intercept)","y_1","x","x_1","c") %in% colnames(y_coeffi_table))){stop("Colum names for y_coeffi_table should be: (Intercept), y_1, x, x_1, c.")}
  if(!all(c("(Intercept)","c_1","x_1","y_1")  %in% colnames(c_coeffi_table))){stop("Colum names for c_coeffi_table should be: (Intercept), c_1, x_1, y_1.")}

  result=rep(NA,nrow(y_coeffi_table))
  for(i in 1:nrow(y_coeffi_table)){
    result[i]=calculate.effect_singlet(t=i,tx,y_coeffi_table,c_coeffi_table,printFlag=F)
  }
  return(result)
}
# Created on 11/03/2023
# this function is created to replace previous "calculate.upto5_lag_effect_allt_withCI",
#                                              "calculate.upto5_total_effect_allt_withCI",
#                                              "calculate.general_effect0101_allt_withCI",
#                                              "calculate.general_effect1001_allt_withCI".
#   so that it calculate contemporanous effect,
#                        lagged effect up to 4, tx=1, tx=c(1,0), tx=c(1,0,0), tx=c(1,0,0,0), tx=c(1,0,0,0,0);
#                        total effect up to 4, tx=1, tx=c(1,1), tx=c(1,1,1), tx=c(1,1,1,1), tx=c(1,1,1,1,1);
#                        general effect of 0101, tx=0, tx=c(0,1), tx=c(0,1,0), tx=c(0,1,0,1);
#                        general effect of 1001, tx=1, tx=c(1,0), tx=c(1,0,0), tx=c(1,0,0,1);
#   separately and more flexibly
calculate.effect_allt_withCI=function(tx,y_coeffi_table,y_coeffi_var_table,c_coeffi_table,c_coeffi_var_table,n_sim=5,printFlag=T){
  if(printFlag){
    cat(blue(" =================================================================================================== \n"))
    cat(blue("This function calculates contemporaneous effect and q-lag effects as well as simulated distribution"))
    cat(blue("  directly from mathematical expressions of coefficients.\n"))
    cat(blue("  Caution: it only applied for recent treatment up to length 5 (or 4 time points ago)\n"))
    cat(red ("  The formulat for outcome regression has to be: y_t = intercept + y_{t-1} + x_t + x_{t-1} + c_t\n"))
    cat(red ("               for covariate regression has to be: c_t = intercept + c_{t-1} + x_{t-1} + c_{t-1}\n"))  
    cat(blue("  Inputs include: [tx] a vector/scalar of recent tx, specified forward as (...,x_{t-2},x_{t-1},x_t)\n"))
    cat(blue("                  [y_coeffi_table] a data.frame of outcome regression coefficients of all times\n"))
    cat(blue("                  [y_coeffi_var_table] a list of outcome regression coefficients Variance Matrix of all times,\n"))
    cat(blue("                  [c_coeffi_table] a data.frame of covariate regression coefficients of all times.\n"))
    cat(blue("                  [c_coeffi_var_table] a list of covariate regression coefficients Variance Matrix of all times,\n"))
    cat(blue("                  [n_sim] a number specifying the numbef of random draws.\n"))
    cat(blue("  Output is: estimate of a q-lag effect and its distribution for all time t.\n"))
    cat(blue("             rows represent all time points and columns represent n_sim.\n"))
    cat(blue(" =================================================================================================== \n"))
  }
  if(!all(is.data.frame(y_coeffi_table) & is.data.frame(c_coeffi_table))){stop("y_coeffi_table and c_coeffi_table should be data.frame.")}
  if(!all(c("(Intercept)","y_1","x","x_1","c") %in% colnames(y_coeffi_table))){stop("Colum names for y_coeffi_table should be: (Intercept), y_1, x, x_1, c.")}
  if(!all(c("(Intercept)","c_1","x_1","y_1")  %in% colnames(c_coeffi_table))){stop("Colum names for c_coeffi_table should be: (Intercept), c_1, x_1, y_1.")}
  
  simulated_coeffi_y=list()
  simulated_coeffi_c=list()
  for(k in 1:nrow(y_coeffi_table)){
    simulated_coeffi_y[[k]]=as.data.frame(mvrnorm(n=n_sim,mu=as.numeric(y_coeffi_table[k,]),Sigma=y_coeffi_var_table[[k]]))
    simulated_coeffi_c[[k]]=as.data.frame(mvrnorm(n=n_sim,mu=as.numeric(c_coeffi_table[k,]),Sigma=c_coeffi_var_table[[k]]))
  }
  
  result_nsim=matrix(NA,ncol=n_sim,nrow=nrow(y_coeffi_table))
  for(i in 1:n_sim){
    y_coeffi_table_temp=list.rbind(lapply(simulated_coeffi_y,function(x){x[i,]}))
    colnames(y_coeffi_table_temp)=colnames(y_coeffi_table)
    c_coeffi_table_temp=list.rbind(lapply(simulated_coeffi_c,function(x){x[i,]}))
    colnames(c_coeffi_table_temp)=colnames(c_coeffi_table)
    result_temp=calculate.effect_allt(tx,y_coeffi_table_temp,c_coeffi_table_temp,printFlag=F)
    result_nsim[,i]=result_temp
  }
  return(result_nsim)
}
# Here listed the no longer used "calculate.upto5_lag_effect_allt_withCI",
#                                "calculate.upto5_total_effect_allt_withCI",
#                                "calculate.general_effect0101_allt_withCI",
#                                "calculate.general_effect1001_allt_withCI".
{
  # calculate.upto5_lag_effect_allt_withCI=function(y_coeffi_table,y_coeffi_var_table,c_coeffi_table,c_coeffi_var_table,n_sim=5,printFlag=T){
  #   if(printFlag){
  #     cat(blue(" =================================================================================================== \n"))
  #     cat(blue("This function calculates q-lag effects for q=0,1,2,3,4 directly from mathematical expressions of coefficients.\n"))
  #     cat(blue("  Caution: it only applied for tx length of 1 to 5.\n"))
  #     cat(blue("                  [y_coeffi_table] a data.frame of outcome regression coefficients of all times\n"))
  #     cat(blue("                  [y_coeffi_var_table] a list of outcome regression coefficients Variance Matrix of all times,\n"))
  #     cat(blue("                  [c_coeffi_table] a data.frame of covariate regression coefficients of all times.\n"))
  #     cat(blue("                  [c_coeffi_var_table] a list of covariate regression coefficients Variance Matrix of all times,\n"))
  #     cat(blue("                  [n_sim] a number specifying the numbef of random draws.\n"))
  #     cat(blue("  Output is: estimate of a distribution of (n_sim copies of) q-lag effects at time t for q=0,1,2,3,4.\n"))
  #     cat(blue(" =================================================================================================== \n"))
  #   }
  #   if(!all(is.data.frame(y_coeffi_table) & is.data.frame(c_coeffi_table))){stop("y_coeffi_table and c_coeffi_table should be data.frame.")}
  #   if(!all(c("(Intercept)","y_1","x","x_1","c") %in% colnames(y_coeffi_table))){stop("Colum names for y_coeffi_table should be: (Intercept), y_1, x, x_1, c.")}
  #   if(!all(c("(Intercept)","c_1","x_1","y_1")  %in% colnames(c_coeffi_table))){stop("Colum names for c_coeffi_table should be: (Intercept), c_1, x_1, y_1.")}
  #   
  #   simulated_coeffi_y=list()
  #   simulated_coeffi_c=list()
  #   for(k in 1:nrow(y_coeffi_table)){
  #     simulated_coeffi_y[[k]]=as.data.frame(mvrnorm(n=n_sim,mu=as.numeric(y_coeffi_table[k,]),Sigma=y_coeffi_var_table[[k]]))
  #     simulated_coeffi_c[[k]]=as.data.frame(mvrnorm(n=n_sim,mu=as.numeric(c_coeffi_table[k,]),Sigma=c_coeffi_var_table[[k]]))
  #   }
  #   
  #   result_nsim=list()
  #   for(i in 1:n_sim){
  #     y_coeffi_table_temp=list.rbind(lapply(simulated_coeffi_y,function(x){x[i,]}))
  #     colnames(y_coeffi_table_temp)=colnames(y_coeffi_table)
  #     c_coeffi_table_temp=list.rbind(lapply(simulated_coeffi_c,function(x){x[i,]}))
  #     colnames(c_coeffi_table_temp)=colnames(c_coeffi_table)
  #     
  #     result_temp=matrix(NA,ncol=5,nrow(y_coeffi_table))
  #     result_temp[,1]=calculate.effect_allt(tx=1,y_coeffi_table_temp,c_coeffi_table_temp,printFlag=F)
  #     result_temp[,2]=calculate.effect_allt(tx=c(1,0),y_coeffi_table_temp,c_coeffi_table_temp,printFlag=F)
  #     result_temp[,3]=calculate.effect_allt(tx=c(1,0,0),y_coeffi_table_temp,c_coeffi_table_temp,printFlag=F)
  #     result_temp[,4]=calculate.effect_allt(tx=c(1,0,0,0),y_coeffi_table_temp,c_coeffi_table_temp,printFlag=F)
  #     result_temp[,5]=calculate.effect_allt(tx=c(1,0,0,0,0),y_coeffi_table_temp,c_coeffi_table_temp,printFlag=F)
  #     colnames(result_temp)=c("contemporaneous","1-lag","2-lag","3-lag","4-lag")
  #     
  #     result_nsim[[i]]=result_temp
  #   }
  #   return(result_nsim)
  # }
  # calculate.upto5_total_effect_allt_withCI=function(y_coeffi_table,y_coeffi_var_table,c_coeffi_table,c_coeffi_var_table,n_sim=5,printFlag=T){
  #   if(printFlag){
  #     cat(blue(" =================================================================================================== \n"))
  #     cat(blue("This function calculates q-lag effects for q=0,1,2,3,4 directly from mathematical expressions of coefficients.\n"))
  #     cat(blue("  Caution: it only applied for tx length of 1 to 5.\n"))
  #     cat(blue("                  [y_coeffi_table] a data.frame of outcome regression coefficients of all times\n"))
  #     cat(blue("                  [y_coeffi_var_table] a list of outcome regression coefficients Variance Matrix of all times,\n"))
  #     cat(blue("                  [c_coeffi_table] a data.frame of covariate regression coefficients of all times.\n"))
  #     cat(blue("                  [c_coeffi_var_table] a list of covariate regression coefficients Variance Matrix of all times,\n"))
  #     cat(blue("                  [n_sim] a number specifying the numbef of random draws.\n"))
  #     cat(blue("  Output is: estimate of a distribution of (n_sim copies of) q-lag effects at time t for q=0,1,2,3,4.\n"))
  #     cat(blue(" =================================================================================================== \n"))
  #   }
  #   if(!all(is.data.frame(y_coeffi_table) & is.data.frame(c_coeffi_table))){stop("y_coeffi_table and c_coeffi_table should be data.frame.")}
  #   if(!all(c("(Intercept)","y_1","x","x_1","c") %in% colnames(y_coeffi_table))){stop("Colum names for y_coeffi_table should be: (Intercept), y_1, x, x_1, c.")}
  #   if(!all(c("(Intercept)","c_1","x_1","y_1")  %in% colnames(c_coeffi_table))){stop("Colum names for c_coeffi_table should be: (Intercept), c_1, x_1, y_1.")}
  #   
  #   simulated_coeffi_y=list()
  #   simulated_coeffi_c=list()
  #   for(k in 1:nrow(y_coeffi_table)){
  #     simulated_coeffi_y[[k]]=as.data.frame(mvrnorm(n=n_sim,mu=as.numeric(y_coeffi_table[k,]),Sigma=y_coeffi_var_table[[k]]))
  #     simulated_coeffi_c[[k]]=as.data.frame(mvrnorm(n=n_sim,mu=as.numeric(c_coeffi_table[k,]),Sigma=c_coeffi_var_table[[k]]))
  #   }
  #   
  #   result_nsim=list()
  #   for(i in 1:n_sim){
  #     y_coeffi_table_temp=list.rbind(lapply(simulated_coeffi_y,function(x){x[i,]}))
  #     colnames(y_coeffi_table_temp)=colnames(y_coeffi_table)
  #     c_coeffi_table_temp=list.rbind(lapply(simulated_coeffi_c,function(x){x[i,]}))
  #     colnames(c_coeffi_table_temp)=colnames(c_coeffi_table)
  #     
  #     result_temp=matrix(NA,ncol=5,nrow(y_coeffi_table))
  #     result_temp[,1]=calculate.effect_allt(tx=1,y_coeffi_table_temp,c_coeffi_table_temp,printFlag=F)
  #     result_temp[,2]=calculate.effect_allt(tx=rep(1,2),y_coeffi_table_temp,c_coeffi_table_temp,printFlag=F)
  #     result_temp[,3]=calculate.effect_allt(tx=rep(1,3),y_coeffi_table_temp,c_coeffi_table_temp,printFlag=F)
  #     result_temp[,4]=calculate.effect_allt(tx=rep(1,4),y_coeffi_table_temp,c_coeffi_table_temp,printFlag=F)
  #     result_temp[,5]=calculate.effect_allt(tx=rep(1,5),y_coeffi_table_temp,c_coeffi_table_temp,printFlag=F)
  #     colnames(result_temp)=c("contemporaneous","1-step","2-step","3-step","4-step")
  #     
  #     result_nsim[[i]]=result_temp
  #   }
  #   return(result_nsim)
  # }
  # calculate.general_effect0101_allt_withCI=function(y_coeffi_table,y_coeffi_var_table,c_coeffi_table,c_coeffi_var_table,n_sim=5,printFlag=T){
  #   if(printFlag){
  #     cat(blue(" =================================================================================================== \n"))
  #     cat(blue("This function calculates q-lag effects for q=0,1,2,3,4 directly from mathematical expressions of coefficients.\n"))
  #     cat(blue("  Caution: it only applied for tx length of 1 to 5.\n"))
  #     cat(blue("                  [y_coeffi_table] a data.frame of outcome regression coefficients of all times\n"))
  #     cat(blue("                  [y_coeffi_var_table] a list of outcome regression coefficients Variance Matrix of all times,\n"))
  #     cat(blue("                  [c_coeffi_table] a data.frame of covariate regression coefficients of all times.\n"))
  #     cat(blue("                  [c_coeffi_var_table] a list of covariate regression coefficients Variance Matrix of all times,\n"))
  #     cat(blue("                  [n_sim] a number specifying the numbef of random draws.\n"))
  #     cat(blue("  Output is: estimate of a distribution of (n_sim copies of) q-lag effects at time t for q=0,1,2,3,4.\n"))
  #     cat(blue(" =================================================================================================== \n"))
  #   }
  #   if(!all(is.data.frame(y_coeffi_table) & is.data.frame(c_coeffi_table))){stop("y_coeffi_table and c_coeffi_table should be data.frame.")}
  #   if(!all(c("(Intercept)","y_1","x","x_1","c") %in% colnames(y_coeffi_table))){stop("Colum names for y_coeffi_table should be: (Intercept), y_1, x, x_1, c.")}
  #   if(!all(c("(Intercept)","c_1","x_1","y_1")  %in% colnames(c_coeffi_table))){stop("Colum names for c_coeffi_table should be: (Intercept), c_1, x_1, y_1.")}
  #   
  #   simulated_coeffi_y=list()
  #   simulated_coeffi_c=list()
  #   for(k in 1:nrow(y_coeffi_table)){
  #     simulated_coeffi_y[[k]]=as.data.frame(mvrnorm(n=n_sim,mu=as.numeric(y_coeffi_table[k,]),Sigma=y_coeffi_var_table[[k]]))
  #     simulated_coeffi_c[[k]]=as.data.frame(mvrnorm(n=n_sim,mu=as.numeric(c_coeffi_table[k,]),Sigma=c_coeffi_var_table[[k]]))
  #   }
  #   
  #   result_nsim=list()
  #   for(i in 1:n_sim){
  #     y_coeffi_table_temp=list.rbind(lapply(simulated_coeffi_y,function(x){x[i,]}))
  #     colnames(y_coeffi_table_temp)=colnames(y_coeffi_table)
  #     c_coeffi_table_temp=list.rbind(lapply(simulated_coeffi_c,function(x){x[i,]}))
  #     colnames(c_coeffi_table_temp)=colnames(c_coeffi_table)
  #     # result_temp=calculate.effect_allt(tx=c(0,1,0,1),y_coeffi_table_temp,c_coeffi_table_temp,printFlag=F)
  #     # result_nsim[[i]]=result_temp
  #     
  #     result_temp=matrix(NA,ncol=4,nrow(y_coeffi_table))
  #     result_temp[,1]=calculate.effect_allt(tx=0,y_coeffi_table_temp,c_coeffi_table_temp,printFlag=F)
  #     result_temp[,2]=calculate.effect_allt(tx=c(0,1),y_coeffi_table_temp,c_coeffi_table_temp,printFlag=F)
  #     result_temp[,3]=calculate.effect_allt(tx=c(0,1,0),y_coeffi_table_temp,c_coeffi_table_temp,printFlag=F)
  #     result_temp[,4]=calculate.effect_allt(tx=c(0,1,0,1),y_coeffi_table_temp,c_coeffi_table_temp,printFlag=F)
  #     colnames(result_temp)=c("1-step","2-step","3-step","4-step")
  #     result_nsim[[i]]=result_temp
  #   }
  #   return(result_nsim)
  # }
  # calculate.general_effect1001_allt_withCI=function(y_coeffi_table,y_coeffi_var_table,c_coeffi_table,c_coeffi_var_table,n_sim=5,printFlag=T){
  #   if(printFlag){
  #     cat(blue(" =================================================================================================== \n"))
  #     cat(blue("This function calculates q-lag effects for q=0,1,2,3,4 directly from mathematical expressions of coefficients.\n"))
  #     cat(blue("  Caution: it only applied for tx length of 1 to 5.\n"))
  #     cat(blue("                  [y_coeffi_table] a data.frame of outcome regression coefficients of all times\n"))
  #     cat(blue("                  [y_coeffi_var_table] a list of outcome regression coefficients Variance Matrix of all times,\n"))
  #     cat(blue("                  [c_coeffi_table] a data.frame of covariate regression coefficients of all times.\n"))
  #     cat(blue("                  [c_coeffi_var_table] a list of covariate regression coefficients Variance Matrix of all times,\n"))
  #     cat(blue("                  [n_sim] a number specifying the numbef of random draws.\n"))
  #     cat(blue("  Output is: estimate of a distribution of (n_sim copies of) q-lag effects at time t for q=0,1,2,3,4.\n"))
  #     cat(blue(" =================================================================================================== \n"))
  #   }
  #   if(!all(is.data.frame(y_coeffi_table) & is.data.frame(c_coeffi_table))){stop("y_coeffi_table and c_coeffi_table should be data.frame.")}
  #   if(!all(c("(Intercept)","y_1","x","x_1","c") %in% colnames(y_coeffi_table))){stop("Colum names for y_coeffi_table should be: (Intercept), y_1, x, x_1, c.")}
  #   if(!all(c("(Intercept)","c_1","x_1","y_1")  %in% colnames(c_coeffi_table))){stop("Colum names for c_coeffi_table should be: (Intercept), c_1, x_1, y_1.")}
  #   
  #   simulated_coeffi_y=list()
  #   simulated_coeffi_c=list()
  #   for(k in 1:nrow(y_coeffi_table)){
  #     simulated_coeffi_y[[k]]=as.data.frame(mvrnorm(n=n_sim,mu=as.numeric(y_coeffi_table[k,]),Sigma=y_coeffi_var_table[[k]]))
  #     simulated_coeffi_c[[k]]=as.data.frame(mvrnorm(n=n_sim,mu=as.numeric(c_coeffi_table[k,]),Sigma=c_coeffi_var_table[[k]]))
  #   }
  #   
  #   result_nsim=list()
  #   for(i in 1:n_sim){
  #     y_coeffi_table_temp=list.rbind(lapply(simulated_coeffi_y,function(x){x[i,]}))
  #     colnames(y_coeffi_table_temp)=colnames(y_coeffi_table)
  #     c_coeffi_table_temp=list.rbind(lapply(simulated_coeffi_c,function(x){x[i,]}))
  #     colnames(c_coeffi_table_temp)=colnames(c_coeffi_table)
  #     # result_temp=calculate.effect_allt(tx=c(1,0,0,1),y_coeffi_table_temp,c_coeffi_table_temp,printFlag=F)
  #     # result_nsim[[i]]=result_temp
  #     result_temp=matrix(NA,ncol=4,nrow(y_coeffi_table))
  #     result_temp[,1]=calculate.effect_allt(tx=1,y_coeffi_table_temp,c_coeffi_table_temp,printFlag=F)
  #     result_temp[,2]=calculate.effect_allt(tx=c(1,0),y_coeffi_table_temp,c_coeffi_table_temp,printFlag=F)
  #     result_temp[,3]=calculate.effect_allt(tx=c(1,0,0),y_coeffi_table_temp,c_coeffi_table_temp,printFlag=F)
  #     result_temp[,4]=calculate.effect_allt(tx=c(1,0,0,1),y_coeffi_table_temp,c_coeffi_table_temp,printFlag=F)
  #     colnames(result_temp)=c("1-step","2-step","3-step","4-step")
  #     result_nsim[[i]]=result_temp
  #     
  #   }
  #   return(result_nsim)
  # }
}
# change name from "calculate.qlag_controlled_direct_effect_allt_withCI"
calculate.controlled_direct_effect_allt_withCI=function(y_coeffi_table,y_coeffi_var_table,n_sim=5,printFlag=T){
  if(printFlag){
    cat(blue(" =================================================================================================== \n"))
    cat(blue("This function calculates 1-lag controlled direct effect as well as simulated distribution.\n"))
    cat(blue("  Caution: it only applied for 1-lag controlled direct effect, as all other q-lag controlled direct effects for q>1 is 0"))
    cat(red ("  The formulat for outcome regression has to be: y_t = intercept + y_{t-1} + x_t + x_{t-1} + c_t\n"))
    cat(red ("               for covariate regression has to be: c_t = intercept + c_{t-1} + x_{t-1} + c_{t-1}\n"))  
    cat(blue("  Inputs include: [tx] a vector/scalar of recent tx, specified forward as (...,x_{t-2},x_{t-1},x_t)\n"))
    cat(blue("                  [y_coeffi_table] a data.frame of outcome regression coefficients of all times\n"))
    cat(blue("                  [y_coeffi_var_table] a list of outcome regression coefficients Variance Matrix of all times,\n"))
    cat(blue("                  [n_sim] a number specifying the numbef of random draws.\n"))
    cat(blue("  Output is: estimate of a distribution of (n_sim copies of) 1-lag controlled direct effect at all times.\n"))
    cat(blue("             rows represent all time points and columns represent n_sim.\n"))
    cat(blue(" =================================================================================================== \n"))
  }
  if(!is.data.frame(y_coeffi_table)){stop("y_coeffi_table should be data.frame.")}
  if(!all(c("(Intercept)","y_1","x","x_1","c") %in% colnames(y_coeffi_table))){stop("Colum names for y_coeffi_table should be: (Intercept), y_1, x, x_1, c.")}

  simulated_coeffi_y=list()
  for(k in 1:nrow(y_coeffi_table)){
    simulated_coeffi_y[[k]]=as.data.frame(mvrnorm(n=n_sim,mu=as.numeric(y_coeffi_table[k,]),Sigma=y_coeffi_var_table[[k]]))
  }
  
  result_nsim=matrix(NA,ncol=n_sim,nrow=nrow(y_coeffi_table))
  for(i in 1:n_sim){
    y_coeffi_table_temp=list.rbind(lapply(simulated_coeffi_y,function(x){x[i,]}))
    colnames(y_coeffi_table_temp)=colnames(y_coeffi_table)
    result_nsim[,i]=y_coeffi_table_temp[,"x_1"]
  }
  return(result_nsim)
}


# ========================== simulation version ============================ #
# tx is a binary vector/scalar, specified forward as (...,x_{t-2},x_{t-1},x_t)
# change name from "simulate.counterfactual_outcome_singlet"
# Updated on 11/03/2023
simulate.counterfactual_singlet=function(t,tx,y_coeffi_table,c_coeffi_table,raw_data,printFlag=T){
  if(printFlag){
    cat(blue(" =================================================================================================== \n"))
    cat(blue("This function simulates DGP to obtain potential outcomes at time t.\n"))
    cat(red ("  The formulat for outcome regression has to be: y_t = intercept + y_{t-1} + x_t + x_{t-1} + c_t\n"))
    cat(red ("               for covariate regression has to be: c_t = intercept + c_{t-1} + x_{t-1} + c_{t-1}\n"))  
    cat(blue("  Inputs include:  [t] observational time t,\n"))
    cat(blue("                   [tx] a vector/scalar of recent tx=(...,x_{t-2},x_{t-1},x_t),\n"))
    cat(blue("                   [y_coeffi_table] a data.frame of outcome regression coefficients of all times,\n"))
    cat(blue("                   [c_coeffi_table] a data.frame of covariate regression coefficients of all times.\n"))
    cat(blue("                   [raw_data] a data.frame of raw data of (y,x,c).\n"))
    cat(blue("  Output is: entire sequence of simulated t from starting t to current t.\n"))
    cat(blue(" =================================================================================================== \n"))
  }
  if(!all(is.numeric(t) & t %% 1 == 0 & t>0 & t<=nrow(y_coeffi_table) & t<=nrow(c_coeffi_table))){stop("Chosen t a propriate positive integer between [1,T].")}
  if(!all(is.data.frame(y_coeffi_table) & is.data.frame(c_coeffi_table))){stop("y_coeffi_table and c_coeffi_table should be data.frame.")}
  if(!all(c("(Intercept)","y_1","x","x_1","c") %in% colnames(y_coeffi_table))){stop("Column names for y_coeffi_table should contain: (Intercept), y_1, x, x_1, c.")}
  if(!all(c("(Intercept)","c_1","x_1","y_1") %in% colnames(c_coeffi_table))){stop("Column names for c_coeffi_table should contain: (Intercept), c_1, x_1, y_1.")}
  if(!is.data.frame(raw_data)){stop("raw_data should be data.frame.")}
  if(!all(c("Date","y","x","c") %in% colnames(raw_data))){stop("Colume names for raw_data should be: Date, y, x, c.")}
  
  # --------------- generate data --------------- #
  periods = length(tx)
  start_t=t-length(tx)+1
  if(start_t<1){stop("The time has been backward beyond 0.")}
  
  # parameters for c
  parameters_for_c=list(baseline=c_coeffi_table$`(Intercept)`[start_t:t],
                        c=c_coeffi_table$c_1[start_t:t],
                        x=c_coeffi_table$x_1[start_t:t],
                        y=c_coeffi_table$y_1[start_t:t],
                        sd=extend(0,periods)) # set 0 variation along the DGP
  model_for_c=c("baseline","c","x","y")
  # it is difficult to generate the first c, let's give it its true value
  initial.c=ifelse(!is.na(raw_data$c[start_t]),raw_data$c[start_t],mean(raw_data$c,na.rm = T)) 
  
  # parameters for x
  #  no initial value of x, no simulation of x
  #  x is entirely given
  if(all(raw_data$x %in% c(0,1)) & all(tx %in% c(0,1))){type_for_x="binary"}else{type_for_x="continuous"}
  
  # parameters for y
  #  no initial value of y
  #  y is entirely simulated from start_t to t
  lag_x_on_y=2;lag_y_on_y=1;lag_c_on_y=1
  model_for_y=c("baseline","c","x","y")
  baseline_y=y_coeffi_table$`(Intercept)`[start_t:t]
  # correction on the generation of first y:
  #     should be y_{st} = intercept + y_{st-1} + x_{st} + x_{st-1} + c_{st}
  #     however, function "sim_y_x_c_univariate" only handles the part of intercept + x_{st} + c_{st} in the generation of 1st obs
  #              therefore, an alternative solution is to add y_{st-1} + x_{st-1} onto the first intercept as compensation
  baseline_y[1]=baseline_y[1]+
    y_coeffi_table$y_1[start_t]*ifelse(start_t-1>0 & !is.na(raw_data$y[start_t-1]),raw_data$y[start_t-1],mean(raw_data$y,na.rm=T)) + # y_{st -1} is used here
    y_coeffi_table$x_1[start_t]*ifelse(start_t-1>0 & !is.na(raw_data$x[start_t-1]),raw_data$x[start_t-1],mean(raw_data$x,na.rm=T))   # x_{st -1} is used here
  parameters_for_y=list(baseline=as.matrix(baseline_y),
                        y=as.matrix(y_coeffi_table$y_1[start_t:t]),
                        x=as.matrix(y_coeffi_table[start_t:t,c("x","x_1")]),
                        c=as.matrix(y_coeffi_table$c[start_t:t]),
                        sd=extend(0,periods))
  colnames(parameters_for_y$x)=paste("t",0:(ncol(parameters_for_y$x)-1),sep="-")
  colnames(parameters_for_y$y)=paste("t",1:ncol(parameters_for_y$y),sep="-")
  colnames(parameters_for_y$c)=paste("t",0:(ncol(parameters_for_y$c)-1),sep="-")
  
  simulated_result=sim_y_x_c_univariate(given.c=F,parameters_for_c=parameters_for_c,model_for_c=model_for_c,initial.c = initial.c,
                                        given.x=T,input.x=tx,type_for_x=type_for_x,
                                        given.y=F,parameters_for_y=parameters_for_y,model_for_y=model_for_y,
                                        lag_y_on_y=lag_y_on_y,lag_x_on_y=lag_x_on_y,lag_c_on_y=lag_c_on_y,
                                        n=sum(periods)) 
  return(simulated_result$y)
}
# change name from "simulate.counterfactual_outcome_singlet_withCI"
# Updated on 11/03/2023
simulate.counterfactual_singlet_withCI=function(t,tx,y_coeffi_table,y_coeffi_var_table,c_coeffi_table,c_coeffi_var_table,
                                                raw_data,n_sim=5,printFlag=T){
  if(printFlag){
    cat(blue(" =================================================================================================== \n"))
    cat(blue("This function simulates DGP to obtain potential outcomes  at time t as well as its distribution .\n"))
    cat(red ("  The formulat for outcome regression has to be: y_t = intercept + y_{t-1} + x_t + x_{t-1} + c_t\n"))
    cat(red ("               for covariate regression has to be: c_t = intercept + c_{t-1} + x_{t-1} + c_{t-1}\n"))  
    cat(blue("  Inputs include: [t] a observational time t,\n"))
    cat(blue("                  [tx] a vector/scalar of recent tx, specified forward as (...,x_{t-2},x_{t-1},x_t)\n"))
    cat(blue("                  [y_coeffi_table] a data.frame of outcome regression coefficients of all times\n"))
    cat(blue("                  [y_coeffi_var_table] a list of outcome regression coefficients Variance Matrix of all times,\n"))
    cat(blue("                  [c_coeffi_table] a data.frame of covariate regression coefficients of all times.\n"))
    cat(blue("                  [c_coeffi_var_table] a list of covariate regression coefficients Variance Matrix of all times,\n"))
    cat(blue("                  [raw_data] a data.frame of raw_data of (y,x,c),\n"))
    cat(blue("                  [n_sim] a number specifying the numbef of random draws.\n"))
    cat(blue("  Output is: entire sequence of simulated outcome and its distribution from starting t to current t.\n"))
    cat(blue("             each row represents each simulated outcome from starting t to current t.\n"))
    cat(blue(" =================================================================================================== \n"))
  }
  if(!all(is.numeric(t) & t %% 1 == 0 & t>=0 & t<=nrow(y_coeffi_table) & t<=nrow(c_coeffi_table))){stop("Chosen t a propriate positive integer between [1,T].")}
  if(!all(is.data.frame(y_coeffi_table) & is.data.frame(c_coeffi_table))){stop("y_coeffi_table and c_coeffi_table should be data.frame.")}
  if(!all(c("(Intercept)","y_1","x","x_1","c") %in% colnames(y_coeffi_table))){stop("Colume names for y_coeffi_table should contain: (Intercept), y_1, x, x_1, c.")}
  if(!all(c("(Intercept)","c_1","x_1","y_1") %in% colnames(c_coeffi_table))){stop("Colume names for c_coeffi_table should contain: (Intercept), c_1, x_1, y_1.")}
  if(!is.data.frame(raw_data)){stop("raw_data should be data.frame.")}
  if(!all(c("Date","y","x","c") %in% colnames(raw_data))){stop("Colume names for raw_data should be: Date, y, x, c.")}
  start_t=t-length(tx)+1
  if(start_t<1){stop("The time has been backward beyond 0.")}

  # sample n_sim copies of coefficient tables for each time from start_t to t
  simulated_coeffi_y=list()
  simulated_coeffi_c=list()
  for(k in start_t:t){
    simulated_coeffi_y[[k]]=as.data.frame(mvrnorm(n=n_sim,mu=as.numeric(y_coeffi_table[k,]),Sigma=y_coeffi_var_table[[k]]))
    simulated_coeffi_c[[k]]=as.data.frame(mvrnorm(n=n_sim,mu=as.numeric(c_coeffi_table[k,]),Sigma=c_coeffi_var_table[[k]]))
  }
  
  result=matrix(NA,nrow=n_sim,ncol=length(tx))
  for(i in 1:n_sim){
    y_coeffi_table_temp=as.data.frame(matrix(NA,nrow=nrow(y_coeffi_table),ncol=ncol(y_coeffi_table))) # empty data frame meeting the requirement of "simulate.counterfactual_singlet"
    y_coeffi_table_temp[start_t:t,]=as.matrix(rbindlist(lapply(simulated_coeffi_y[start_t:t],function(x){x[i,]}))) # update the relevant proportion
    colnames(y_coeffi_table_temp)=colnames(y_coeffi_table) # give column names as required by "simulate.counterfactual_singlet"
    c_coeffi_table_temp=as.data.frame(matrix(NA,nrow=nrow(c_coeffi_table),ncol=ncol(c_coeffi_table))) # empty data frame meeting the requirement of "simulate.counterfactual_singlet"
    c_coeffi_table_temp[start_t:t,]=as.matrix(rbindlist(lapply(simulated_coeffi_c[start_t:t],function(x){x[i,]})))  # update the relevant proportion
    colnames(c_coeffi_table_temp)=colnames(c_coeffi_table) # give column names as required by "simulate.counterfactual_singlet"
    result[i,]=simulate.counterfactual_singlet(t,tx,y_coeffi_table = y_coeffi_table_temp,c_coeffi_table = c_coeffi_table_temp,
                                               raw_data,printFlag=F)
  }
  return(result)
}
# plot the results from function "simulate.counterfactual_singlet_withCI"
plot_simulatedCI=function(simulated_dist,probs=c(0.05,0.95),printFlag=T,ylim=c(0,2)){
  if(printFlag){
    cat(blue(" =================================================================================================== \n"))
    cat(blue("This function returns the alpha% lower and upper bounds and mean level of simulated conterfactuals .\n"))
    cat(blue("  Inputs include: [simulated_dist] result of function simulate.counterfactual_singlet_withCI,\n"))
    cat(blue("                                   each colunmn represents a t, each row represents a copy.\n"))
    cat(blue("                  [probs] lower and upper bound quantiles.\n"))
    cat(blue(" =================================================================================================== \n"))
  }
  result=as.data.frame(matrix(NA,ncol=3,nrow=ncol(simulated_dist)))
  colnames(result)=c("lower","upper","mean")
  for(j in 1:ncol(simulated_dist)){
    result[j,]=c(quantile(simulated_dist[,j],probs,na.rm = T),mean(simulated_dist[,j]))
  }
  if(printFlag==T){
    plot(0:(nrow(result)-1),result$mean,type="l",lty=3,ylim=ylim)
    points(0:(nrow(result)-1),result$mean)
    points(0:(nrow(result)-1),result$upper,col="blue",lty=3,type="l")
    points(0:(nrow(result)-1),result$lower,col="blue",lty=3,type="l")
  }
  return(result)
}


test.positivity=function(tx,data,length){
  if(!is.data.frame(data)){stop("data should be a data.frame.")}
  if(!tx %in% colnames(data)){stop("selected tx are not in the data.")}
  if(!(length>0 & length%%1 ==0)){stop("length should be a positive integer.")}
  
  param=list(list(lagged_param=list(variables=tx,param=length)))
  data_temp=add_variables_procedures(data,param)
  
  details=list()
  percentage=rep(NA,length)
  for(k in 1:length){
    if(k==1){selected=tx}else if(k>1){selected=c(tx,paste(tx,"_",1:(k-1),sep=""))}
    n_groups_of_0=data_temp %>% group_by(.dots=lapply(selected,as.symbol)) %>% summarize(count = n()) %>% as.data.frame()
    n_groups_of_0=n_groups_of_0[rowSums(is.na(n_groups_of_0))==0,]
    colnames(n_groups_of_0)=gsub(tx,"tx",colnames(n_groups_of_0))
    details[[k]]=n_groups_of_0
    percentage[k]=length(n_groups_of_0$count>0)/(2^k)
  }
  return(list(percentage=percentage,details=details))
}

