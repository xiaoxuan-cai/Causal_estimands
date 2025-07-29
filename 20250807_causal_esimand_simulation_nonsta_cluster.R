# File: 20250708_causal_estimand_simulation_non.R
# Date: 2025.07.08
source("/home/cai.1083/paper_causal_estimand_simu/Summary 1_packages.R")
source("/home/cai.1083/paper_causal_estimand_simu/Summary 2_helper functions.R")
source("/home/cai.1083/paper_causal_estimand_simu/helper_causal estimands.R")
source("/home/cai.1083/paper_causal_estimand_simu/helper_SSM_version3.R")
source("/home/cai.1083/paper_causal_estimand_simu/helper_other.R")

for(seed in 1:n_seed){
  cat(seed,"\n")
  # ------------------------------#
  #       Simulation Data         #
  # ------------------------------#
  # non-stationary case: 
  #         Y_t = random walk (var=0.01) + 0.5 Y_{t-1} + -2.5 or -5 X_t - 2 X_{t-1} - C_t + noise (var=0.1)
  #         X_t = 0 + 0.1 C_{t} + 0.3 X_{t-1} + 0.05 Y_{t-1} + noise (var=1)
  #         C_t = 10 + 0.3 C_{t-1} - 0.2 Y_{t-1} + 1 or 0.5 X_{t-1} + noise (var=1)
  # simulate 1100 time points and only keep 1000
  set.seed(seed)
  periods = c(500,300,300)
  # parameters for c
  parameters_for_c=list(baseline=extend(10,periods),
                        c=extend(0.3,periods),
                        x=extend(c(1,0.5,1),periods),
                        y=extend(-0.2,periods),
                        sd=extend(sqrt(1),periods))
  model_for_c=c("baseline","c","x","y")
  # parameters for x
  parameters_for_x=list(baseline=extend(0,periods),
                        c=extend(0.1,periods),
                        x=extend(0.3,periods),
                        y=extend(0.05,periods),
                        sd=extend(sqrt(1),periods))
  model_for_x=c("baseline","c","x","y")
  # parameters for y
  lag_x_on_y=2;lag_c_on_y=1;lag_y_on_y=1
  treatment_effects=matrix(c(-2.5,-2,-5,-2,-2.5,-2),ncol=lag_x_on_y,byrow=T)
  model_for_y=c("baseline","c","x","y")
  set.seed(seed+1, kind = "Mersenne-Twister", normal.kind = "Inversion")
  baseline=sim_randomwalk(baseline=40,sd=0.1,periods=sum(periods))
  parameters_for_y=list(baseline=baseline,
                        y=as.matrix(extend(0.5,periods)),
                        x=extend_into_matrix(treatment_effects,periods),
                        c=as.matrix(extend(-1,periods)),
                        sd=extend(sqrt(0.1),periods))
  colnames(parameters_for_y$x)=paste("t",0:(ncol(parameters_for_y$x)-1),sep="-")
  colnames(parameters_for_y$y)=paste("t",1:ncol(parameters_for_y$y),sep="-")
  colnames(parameters_for_y$c)=paste("t",0:(ncol(parameters_for_y$c)-1),sep="-")
  simulation=sim_y_x_c_univariate(given.c=F,input.c=NULL,parameters_for_c=parameters_for_c,model_for_c=model_for_c,initial.c=NULL,
                                  given.x=F,input.x=NULL,parameters_for_x=parameters_for_x,model_for_x=model_for_x,initial.x=NULL,type_for_x="continuous",
                                  given.y=F,input.y=NULL,parameters_for_y=parameters_for_y,model_for_y=model_for_y,initial.y=NULL,
                                  lag_y_on_y=lag_y_on_y,lag_x_on_y=lag_x_on_y,lag_c_on_y=lag_c_on_y,interaction_pairs=NULL,
                                  n=sum(periods),printFlag=F)
  data=data.frame(Date=Sys.Date()-days(length(simulation$y):1),y=simulation$y,x=simulation$x,c=simulation$c,simulation$noise_y)
  param=list(list(lagged_param=list(variables=c("y","x","c"),param=rep(1,3))))
  data=add_variables_procedures(data,param)
  data=data[-c(1:100),]
  
  # ------------------------------#
  #         Fit SSM model         #
  # ------------------------------#
  # fit state space model for y
  data_ss_y=data[,c("Date","y","y_1","x","x_1","c")]
  formula_y="y~y_1+x+x_1+c"
  ss_param_y=list(inits=c(log(0.1),log(0.01)),m0=c(40,0.5,-2.5,-2,-1),C0=diag(rep(10^3),5),
                  AR1_coeffi=NULL,rw_coeffi="intercept",v_cp_param=NULL,
                  w_cp_param=list(list(variable="x",segments=3,changepoints=c(400,700),fixed_cpts=F)))
  SSM_y_result=run.SSM(data_ss=data_ss_y,formula=formula_y,ss_param_temp=ss_param_y,max_iteration=100,
                       cpt_learning_param=list(cpt_method="mean",burnin=1/10,mergeband=20,convergence_cri=10),
                       cpt_merge_option = "unanimous",
                       dlm_option="smooth",
                       dlm_cpt_learning_option="smooth",bandwidth=20, cpt_V=1, printFlag=F)
  # SSM_y_result$out_filter$mod$W # true = 0.01
  # SSM_y_result$out_filter$mod$V # true = 0.1
  outF_y=SSM_y_result$out_filter
  outS_y=SSM_y_result$out_smooth
  
  # fit state space model for c
  data_ss_c=data[,c("c","c_1","x_1","y_1")]
  formula_c="c~c_1+x_1+y_1"
  ss_param_c=list(inits=c(log(1)),m0=c(10,0.3,-0.2,1),C0=diag(rep(10^3),4),
                  AR1_coeffi=NULL,rw_coeffi=NULL,v_cp_param=NULL,
                  w_cp_param=list(list(variable="x_1",segments=3,changepoints=c(400,700),fixed_cpts=F)))
  SSM_c_result=run.SSM(data_ss=data_ss_c,formula=formula_c,ss_param_temp=ss_param_c,max_iteration=100,
                       cpt_learning_param=list(cpt_method="mean",burnin=1/10,mergeband=20,convergence_cri=10),
                       cpt_merge_option = "unanimous",
                       dlm_option="smooth",
                       dlm_cpt_learning_option="smooth",bandwidth=20, cpt_V=5, printFlag=F)
  # SSM_c_result$out_filter$mod$V # true = 1
  outF_c=SSM_c_result$out_filter
  outS_c=SSM_c_result$out_smooth
  
  y_coeffi_table=as.data.frame(outS_y$s[-1,]);colnames(y_coeffi_table)=c("(Intercept)","y_1","x","x_1","c")
  # y_coeffi_var_table=lapply(1:(nrow(outS_y$s)-1),function(i){dlmSvd2var(outS_y$U.S[[i+1]],outS_y$D.S[i+1,])})
  c_coeffi_table=as.data.frame(outS_c$s[-1,]);colnames(c_coeffi_table)=c("(Intercept)","c_1","x_1","y_1")
  # c_coeffi_var_table=lapply(1:(nrow(outS_c$s)-1),function(i){dlmSvd2var(outS_c$U.S[[i+1]],outS_c$D.S[i+1,])})
  
  # ------------------------------#
  #      Calculate Estimand       #
  # ------------------------------#
  # Simulation result:
  # 1.contemporaneous effect, 1 lag controlled direct effect, 1-lag effect, 1 step total effect of all time points
  # 2.impulse impact plot and step response plot at t=600 for 10 days
  
  # the following causal estimands are calculated analytically without CI (faster than simulated method)
  # 1.contemporaneous effect, 1 lag controlled direct effect, 1-lag effect, 1 step total effect of all time points
  time=1:(nrow(outF_y$m)-1)
  # contemporaneous effect: Y(1) - Y(0)
  contemporaneous[,seed]=calculate.causaleffect(t=time,tx=1,y_coeffi_table,c_coeffi_table,printFlag=F)-
    calculate.causaleffect(t=time,tx=0,y_coeffi_table,c_coeffi_table,printFlag=F)
  # 1-lag effect: Y(1,0) - Y(0,0)
  lag1[,seed]=calculate.causaleffect(t=time,tx=c(1,0),y_coeffi_table,c_coeffi_table,printFlag=F)-
    calculate.causaleffect(t=time,tx=c(0,0),y_coeffi_table,c_coeffi_table,printFlag=F)
  # 2-lag effect: Y(1,0,0) - Y(0,0,0)
  lag2[,seed]=calculate.causaleffect(t=time,tx=c(1,0,0),y_coeffi_table,c_coeffi_table,printFlag=F)-calculate.causaleffect(t=time,tx=c(0,0,0),y_coeffi_table,c_coeffi_table,printFlag=F)
  # 3-lag effect: Y(1,0,0,0) - Y(0,0,0,0): -2.085623
  lag3[,seed]=calculate.causaleffect(t=time,tx=c(1,0,0,0),y_coeffi_table,c_coeffi_table,printFlag=F)-calculate.causaleffect(t=time,tx=c(0,0,0,0),y_coeffi_table,c_coeffi_table,printFlag=F)
  # 4-lag effect: Y(1,0,0,0,0) - Y(0,0,0,0,0): -1.615715
  lag4[,seed]=calculate.causaleffect(t=time,tx=c(1,0,0,0,0),y_coeffi_table,c_coeffi_table,printFlag=F)-calculate.causaleffect(t=time,tx=c(0,0,0,0,0),y_coeffi_table,c_coeffi_table,printFlag=F)
  # 1 lag controlled direct effect
  lag1_direct[,seed]=calculate.controlled_direct_effect(t=time,y_coeffi_table,printFlag=F)
  # 1 step total effect: Y(1,1) - Y(0,0)
  total1[,seed]=calculate.causaleffect(t=time,tx=c(1,1),y_coeffi_table,c_coeffi_table,printFlag=F)-
    calculate.causaleffect(t=time,tx=c(0,0),y_coeffi_table,c_coeffi_table,printFlag=F)
  
  # 2.impulse impact plot and step response plot at t=500 for 10 days
  one_time=500
  # # Method 1: calculate time-invariant causal estimands (withoutCI) from analytical expressions
  # # Y(1) - Y(0)
  # estimand1=calculate.causaleffect(t=one_time,tx=1,y_coeffi_table,c_coeffi_table)-calculate.causaleffect(t=one_time,tx=0,y_coeffi_table,c_coeffi_table);estimand1
  # # Y(1,0) - Y(0,0)
  # estimand2=calculate.causaleffect(t=one_time,tx=c(1,0),y_coeffi_table,c_coeffi_table)-calculate.causaleffect(t=one_time,tx=c(0,0),y_coeffi_table,c_coeffi_table);estimand2
  # # Y(1,0,0) - Y(0,0,0)
  # estimand3=calculate.causaleffect(t=one_time,tx=c(1,0,0),y_coeffi_table,c_coeffi_table)-calculate.causaleffect(t=one_time,tx=c(0,0,0),y_coeffi_table,c_coeffi_table);estimand3
  # # Y(1,0,0,0) - Y(0,0,0,0)
  # estimand4=calculate.causaleffect(t=one_time,tx=c(1,0,0,0),y_coeffi_table,c_coeffi_table)-calculate.causaleffect(t=one_time,tx=c(0,0,0,0),y_coeffi_table,c_coeffi_table);estimand4
  # # Y(1,0,0,0,0) - Y(0,0,0,0,0)
  # estimand5=calculate.causaleffect(t=one_time,tx=c(1,0,0,0,0),y_coeffi_table,c_coeffi_table)-calculate.causaleffect(t=one_time,tx=c(0,0,0,0,0),y_coeffi_table,c_coeffi_table);estimand5
  # plot(0:4,c(estimand1,estimand2,estimand3,estimand4,estimand5),type="l",ylim=c(-5.1,2),xlab="# lags",ylab="causal estimands")
  # abline(h=0,col="red")
  # Method 2: from simulated version
  estimand_1to7[,seed]=simulate.counterfactual_path_singlet(t=one_time,tx=c(1,0,0,0,0,0,0),y_coeffi_table,c_coeffi_table,raw_data=data,printFlag=F)-
    simulate.counterfactual_path_singlet(t=one_time,tx=c(0,0,0,0,0,0,0),y_coeffi_table,c_coeffi_table,raw_data=data,printFlag=F)
}
# save(contemporaneous,lag1,lag2,lag3,lag4,lag1_direct,total1,estimand_1to7,
#      file="simu_causal_nonsta.Rdata")
save(contemporaneous,lag1,lag2,lag3,lag4,lag1_direct,total1,estimand_1to7,
     file="/home/cai.1083/paper_causal_estimand_simu/simu_causal_nonsta.Rdata")

