# File: 20210801_causal_estimand_simulation.R
# Date: 2021-08-01
# Update: 2025-07-08
# Note: This script was created for the paper revision.
#       Previously, we included only stationary-case simulations to test the algorithm.
#       These results were not presented in the paper or appendix.
#       In response to reviewer feedback, we added a new simulation example for the non-stationary case, reflecting the structure of the real data application.
# Results: 1. Various causal estimands for appendix, simulation section
#          2. Impulse impact plot and Step response plot for for appendix, simulation section

# Structure:
# Part I (Stationary Case):
#     This section serves as a comprehensive demonstration of the full set of implemented functions.
#     It includes both analytical and simulation-based methods for estimating causal effects with and without CI.
#     Estimation without CI does not require setting the seed or specifying coefficient variance.
#     These results are NOT included in the paper or appendix and are only for illustrative purposes.
# Part II (Non-Stationary Case):
#     This section is used in the appendix and is designed to reflect the time-varying nature of the real data application.
#     It makes use of only a subset of the full functionality:
#     For estimating causal effect trajectories across all time points, we use analytical methods *without* CI, due to
#     - computation efficiency (for up to 4 lags) as it is faster than simulated version and needs fewer input
#     - CI is not calculated as CI coverage is not reported in the result
#     Simulations in this section were run on a computing cluster, and the outputs were used to generate summary plots.
#     Features of the non-stationary case:
#     - Time-varying coefficients
#     - Feedback between treatment, outcome, and covariates
#     - Continuous treatment, outcome, and covariate variables

# Workflow:
# 1. Run simulations on the cluster
# 2. Calculate true causal effects
# 3. Generate summary plots

library(MASS) 
library(crayon) # for cat colors, OSU
library(lubridate) # for dates, OSU
library(mice) # for multiple imputation, OSU
library(forecast) # for most the time series functions: auto.arima, OSU
library(changepoint) # for finding change points, OSU
library(dlm) # for statespace model, OSU
library(data.table) # for simulate.counterfactual_path_singlet to combine items of list
source("/Users/xiaoxuancai/Documents/GitHub/mHealth_data_processing/Summary 2_helper functions.R")
source("/Users/xiaoxuancai/Documents/GitHub/Causal_estimands/helper_causal estimands.R")
source("/Users/xiaoxuancai/Documents/GitHub/SSMimpute2/helper_SSM_version3.R")
source("/Users/xiaoxuancai/Documents/GitHub/SSMimpute2/helper_other.R")

#####################################################################################
###########                     Stationary case                           ###########
###########      (you may skip as it's not used for paper/appendix)       ###########
#####################################################################################
# Y_t= 35 + 0.5 Y_{t-1} - 2 X_t - X_{t-1} - C_t + noise (var=0.1)
# X_t= 0  + 0.1 C_{t} + 0.3 X_{t-1} + 0.05 Y_{t-1} + noise (var=0.1)
# C_t= 7  + 0.2 C_{t-1} - 0.2 Y_{t-1} + X_{t-1} + noise (var=0.1)
set.seed(1)
periods=1000
# parameters for c
parameters_for_c=list(baseline=extend(7,periods),
                      c=extend(0.2,periods),
                      y=extend(-0.2,periods),
                      x=extend(1,periods),
                      sd=extend(sqrt(0.1),periods))
model_for_c=c("baseline","c","x","y")
# parameters for x
parameters_for_x=list(baseline=extend(0,periods),
                      c=extend(0.1,periods),
                      x=extend(0.3,periods),
                      y=extend(0.05,periods),
                      sd=extend(sqrt(0.1),periods))
model_for_x=c("baseline","c","x","y")
# parameters for y
lag_x_on_y=2;lag_c_on_y=1;lag_y_on_y=1
treatment_effects=matrix(c(-2,-1),ncol=lag_x_on_y,byrow=T)
model_for_y=c("baseline","c","x","y")
parameters_for_y=list(baseline=extend(35,periods),
                      y=as.matrix(extend(0.5,periods)),
                      x=extend_into_matrix(treatment_effects,periods),
                      c=as.matrix(extend(-1,periods)),
                      sd=extend(sqrt(0.1),periods))
colnames(parameters_for_y$x)=paste("t",0:(ncol(parameters_for_y$x)-1),sep="-")
colnames(parameters_for_y$y)=paste("t",1:ncol(parameters_for_y$y),sep="-")
colnames(parameters_for_y$c)=paste("t",0:(ncol(parameters_for_y$c)-1),sep="-")
# generate data
simulation=sim_y_x_c_univariate(given.c=F,input.c=NULL,parameters_for_c=parameters_for_c,model_for_c=model_for_c,initial.c=NULL,
                                given.x=F,input.x=NULL,parameters_for_x=parameters_for_x,model_for_x=model_for_x,initial.x=NULL,type_for_x="continuous",
                                given.y=F,input.y=NULL,parameters_for_y=parameters_for_y,model_for_y=model_for_y,initial.y=NULL,
                                lag_y_on_y=lag_y_on_y,lag_x_on_y=lag_x_on_y,lag_c_on_y=lag_c_on_y,interaction_pairs=NULL,
                                n=sum(periods),printFlag=T)
data=data.frame(Date=Sys.Date()-days(length(simulation$y):1),y=simulation$y,x=simulation$x,c=simulation$c,simulation$noise_y)
param=list(list(lagged_param=list(variables=c("y","x","c"),param=rep(1,3))))
data=add_variables_procedures(data,param)

result_y=summary(lm(y~y_1+x+x_1+c,data=data[-c(1),]));result_y
# Coefficients:
#               Estimate Std. Error t value Pr(>|t|)    
#   (Intercept) 34.120015   0.409160   83.39   <2e-16 ***
#   y_1          0.521032   0.009281   56.14   <2e-16 ***
#   x           -1.970228   0.031350  -62.84   <2e-16 ***
#   x_1         -1.070711   0.041688  -25.68   <2e-16 ***
#   c           -0.959736   0.031287  -30.68   <2e-16 ***
result_c=summary(lm(c~c_1+x_1+y_1,data=data[-c(1),]));result_c
# Coefficients:
#               Estimate Std. Error t value Pr(>|t|)    
# (Intercept)    7.32109    0.77945   9.393  < 2e-16 ***
#   c_1          0.17448    0.03815   4.573 5.41e-06 ***
#   x_1          0.97918    0.03555  27.541  < 2e-16 ***
#   y_1         -0.20456    0.01399 -14.622  < 2e-16 ***

y_coeffi_table=as.data.frame(matrix(rep(result_y$coefficients[,1],sum(periods)),nrow=sum(periods),byrow=T))
colnames(y_coeffi_table)=names(result_y$coefficients[,1])
y_coeffi_var_table=lapply(1:sum(periods),function(i){vcov(result_y)})
c_coeffi_table=as.data.frame(matrix(rep(result_c$coefficients[,1],sum(periods)),nrow=sum(periods),byrow=T))
colnames(c_coeffi_table)=names(result_c$coefficients[,1])
c_coeffi_var_table=lapply(1:sum(periods),function(i){vcov(result_c)})

# Calculate causal estimands for t = 10
# (Note: the choice of time point is arbitrary, as the coefficients are time-invariant)

# Method 1: Analytical calculation of causal estimands (with/without CI)
# ---------------------------------------------------------------
# - Confidence intervals (CI) can be optionally computed.
# - The unified function `calculate.causaleffect` replaces earlier versions:
#     * `calculate.effect_singlet`
#     * `calculate.effect_allt`
#     * `calculate.effect_allt_withCI`
# - This new function yields identical results but offers greater flexibility:
#     * Allows estimation at multiple time point(s)
#     * Supports both with and without CI in a single function
#     * Can only calculate contemporaneous effect and lag effect up to 4 lags
time=10
# Y(1) - Y(0): -1.970228
estimand1=calculate.causaleffect(t=time,tx=1,y_coeffi_table,c_coeffi_table)-calculate.causaleffect(t=time,tx=0,y_coeffi_table,c_coeffi_table);estimand1
# calculate estimand1 with CI
# estimand1=calculate.causaleffect(t=time,tx=1,y_coeffi_table,c_coeffi_table,CI=T,n_sim=10,y_coeffi_var_table,c_coeffi_var_table,seed=1)-
#           calculate.causaleffect(t=time,tx=0,y_coeffi_table,c_coeffi_table,CI=T,n_sim=10,y_coeffi_var_table,c_coeffi_var_table,seed=1);estimand1
# Y(1,0) - Y(0,0): -3.423824
estimand2=calculate.causaleffect(t=time,tx=c(1,0),y_coeffi_table,c_coeffi_table)-calculate.causaleffect(t=time,tx=c(0,0),y_coeffi_table,c_coeffi_table);estimand2
# Y(1,0,0) - Y(0,0,0): -2.687571
estimand3=calculate.causaleffect(t=time,tx=c(1,0,0),y_coeffi_table,c_coeffi_table)-calculate.causaleffect(t=time,tx=c(0,0,0),y_coeffi_table,c_coeffi_table);estimand3
# Y(1,0,0,0) - Y(0,0,0,0): -2.085623
estimand4=calculate.causaleffect(t=time,tx=c(1,0,0,0),y_coeffi_table,c_coeffi_table)-calculate.causaleffect(t=time,tx=c(0,0,0,0),y_coeffi_table,c_coeffi_table);estimand4
# Y(1,0,0,0,0) - Y(0,0,0,0,0): -1.615715
estimand5=calculate.causaleffect(t=time,tx=c(1,0,0,0,0),y_coeffi_table,c_coeffi_table)-calculate.causaleffect(t=time,tx=c(0,0,0,0,0),y_coeffi_table,c_coeffi_table);estimand5
plot(0:4,c(estimand1,estimand2,estimand3,estimand4,estimand5),type="l",ylim=c(-4,0),xlab="# lags",ylab="causal estimands")
abline(h=0,col="red")

# Method 2: Simulation-based calculation of causal estimands (with/without CI)
# ---------------------------------------------------------------------
# - Simulates the full path of lagged causal effects ending at a time point (e.g., t = 10)
# - Unlike the analytical method (which simulates 1-lag effects across multiple time points),
#   this approach simulates *all* lag effects (e.g., contemporaneous effect and 1 to 4 lag effects) ending at a *single* time point
one_time=10
estimand_1to5=simulate.counterfactual_path_singlet(t=one_time,tx=c(1,0,0,0,0),y_coeffi_table,c_coeffi_table,raw_data=data[,c(1:4)])-
  simulate.counterfactual_path_singlet(t=one_time,tx=c(0,0,0,0,0),y_coeffi_table,c_coeffi_table,raw_data=data[,c(1:4)])
estimand_1to5 # the same as c(estimand1,estimand2,estimand3,estimand4,estimand5)
estimand_1to5_withCI=simulate.counterfactual_path_singlet(t=one_time,tx=c(1,0,0,0,0),y_coeffi_table,c_coeffi_table,raw_data=data[,c(1:4)],CI=T,n_sim=100,y_coeffi_var_table,c_coeffi_var_table,seed=1)-
  simulate.counterfactual_path_singlet(t=one_time,tx=c(0,0,0,0,0),y_coeffi_table,c_coeffi_table,raw_data=data[,c(1:4)],CI=T,n_sim=100,y_coeffi_var_table,c_coeffi_var_table,seed=1)
boxplot(estimand_1to5_withCI)
points(estimand_1to5,col="blue",pch=16,type="l",lty=3)

#####################################################################################
###########                   Non-stationary case                         ###########
#####################################################################################
# Cluster Execution Instructions:
# 1. Save R commands (between "start: code on cluster" and "end: code on cluster”)
#    as a new script: "20250807_causal_estimand_simulation_nonsta_cluster.R"
# 2. Copy all R source files referenced into the same folder on cluster.
# 3. Create a job script "causal_estimand_simu_stat.sh" with the following content:
# --------------------------------------------------------------
# #!/usr/bin/env bash
# #SBATCH --partition=stat,batch
# #SBATCH --time=01:00:00
# #SBATCH --nodes=1
# #SBATCH --ntasks-per-node=8
# #SBATCH --mem=1G
# #SBATCH --mail-type=FAIL
# #SBATCH --mail-user=cai.1083@osu.edu
#
# # Navigate to the submission directory so Rscript can find the files
# cd $SLURM_SUBMIT_DIR
#
# # Load necessary modules
# ml gnu/9.1.0
# ml R/4.4.0
#
# # Execute the R script
# Rscript "20250807_causal_estimand_simulation_nonsta_cluster.R"
# --------------------------------------------------------------
# 4. Submit the job using:
# sbatch causal_estimand_simu_stat.sh


# --------------------- start: code on cluster --------------------- #
library(MASS) 
library(crayon) # for cat colors, OSU
library(lubridate) # for dates, OSU
library(mice) # for multiple imputation, OSU
library(forecast) # for most the time series functions: auto.arima, OSU
library(changepoint) # for finding change points, OSU
library(dlm) # for statespace model, OSU
source("/home/cai.1083/paper_causal_estimand_simu/Summary 2_helper functions.R")
source("/home/cai.1083/paper_causal_estimand_simu/helper_causal estimands.R")
source("/home/cai.1083/paper_causal_estimand_simu/helper_SSM_version3.R")
source("/home/cai.1083/paper_causal_estimand_simu/helper_other.R")

n_seed=500
contemporaneous=lag1=lag2=lag3=lag4=lag1_direct=total1=matrix(NA,nrow=1000,ncol=n_seed)
estimand_qlag=estimand_qstep_total=matrix(NA,nrow=20,ncol=n_seed)
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
  estimand_qlag[,seed]=simulate.counterfactual_path_singlet(t=one_time,tx=c(1,rep(0,19)),y_coeffi_table,c_coeffi_table,raw_data=data,printFlag=F)-
    simulate.counterfactual_path_singlet(t=one_time,tx=rep(0,20),y_coeffi_table,c_coeffi_table,raw_data=data,printFlag=F)
  estimand_qstep_total[,seed]=simulate.counterfactual_path_singlet(t=one_time,tx=rep(1,20),y_coeffi_table,c_coeffi_table,raw_data=data,printFlag=F)-
    simulate.counterfactual_path_singlet(t=one_time,tx=rep(0,20),y_coeffi_table,c_coeffi_table,raw_data=data,printFlag=F)
}
save(contemporaneous,lag1,lag2,lag3,lag4,lag1_direct,total1,estimand_qlag,estimand_qstep_total,
     file="/home/cai.1083/paper_causal_estimand_simu/simu_causal_nonsta.Rdata")
# --------------------- end: code on cluster --------------------- #

# calculate true values
y_coeffi_table_new=as.data.frame(list.cbind(parameters_for_y)[-c(1:100),-6]);colnames(y_coeffi_table_new)=c("(Intercept)","y_1","x","x_1","c")
c_coeffi_table_new=as.data.frame(list.cbind(parameters_for_c)[-c(1:100),-5]);colnames(c_coeffi_table_new)=c("(Intercept)","c_1","x_1","y_1")
lag1_ture=calculate.causaleffect(t=1:1000,tx=c(1,0),y_coeffi_table_new,c_coeffi_table_new,printFlag=F)-calculate.causaleffect(t=1:1000,tx=c(0,0),y_coeffi_table_new,c_coeffi_table_new,printFlag=F)
lag2_ture=calculate.causaleffect(t=1:1000,tx=c(1,0,0),y_coeffi_table_new,c_coeffi_table_new,printFlag=F)-calculate.causaleffect(t=1:1000,tx=c(0,0,0),y_coeffi_table_new,c_coeffi_table_new,printFlag=F)
lag3_ture=calculate.causaleffect(t=1:1000,tx=c(1,0,0,0),y_coeffi_table_new,c_coeffi_table_new,printFlag=F)-calculate.causaleffect(t=1:1000,tx=c(0,0,0,0),y_coeffi_table_new,c_coeffi_table_new,printFlag=F)
lag4_ture=calculate.causaleffect(t=1:1000,tx=c(1,0,0,0,0),y_coeffi_table_new,c_coeffi_table_new,printFlag=F)-calculate.causaleffect(t=1:1000,tx=c(0,0,0,0,0),y_coeffi_table_new,c_coeffi_table_new,printFlag=F)
estimand_qlag_ture=simulate.counterfactual_path_singlet(t=500,tx=c(1,0,0,0,0,0,0),y_coeffi_table_new,c_coeffi_table_new,raw_data=data,printFlag=F)-
  simulate.counterfactual_path_singlet(t=500,tx=c(0,0,0,0,0,0,0),y_coeffi_table_new,c_coeffi_table_new,raw_data=data,printFlag=F)
estimand_qstep_total_ture=simulate.counterfactual_path_singlet(t=500,tx=rep(1,7),y_coeffi_table_new,c_coeffi_table_new,raw_data=data,printFlag=F)-
  simulate.counterfactual_path_singlet(t=500,tx=rep(0,7),y_coeffi_table_new,c_coeffi_table_new,raw_data=data,printFlag=F)

# plots
library(rlist)
load("/Users/xiaoxuancai/Documents/GitHub/Causal_estimands/simu_causal_nonsta.Rdata")
address_appendix = "/Users/xiaoxuancai/Dropbox/MHealthPsychSummerProject2020/Xiaoxuan_Cai/[Paper 1] Causal estimands for time series data/Graphs_appendix/"

contemporaneous_summary=data.frame(mean=rowMeans(contemporaneous),t(apply(contemporaneous,1,quantile,probs=c(0.025,0.975))))
colnames(contemporaneous_summary)=c("mean","2.5%","97.5%")
pdf(file=paste(address_appendix,"contemp_simu.pdf",sep=""),width=11,height=8)
par(mar = c(5, 6, 1, 1))
plot(1:1000,contemporaneous_summary$mean,ylim=c(-7,-1),type="l",
     bty="n",cex.lab=3,cex.axis=3,
     ylab="contemporaneous",xlab="time")
polygon(c(1:1000,rev(1:1000)),c(contemporaneous_summary$`2.5%`,rev(contemporaneous_summary$`97.5%`)),col="grey90",border="grey")
points(1:1000,contemporaneous_summary$mean,type="l")
points(1:400,rep(-2.5,400),col="red",type="l",lty=2,lwd=2)
points(401:700,rep(-5,300),col="red",type="l",lty=2,lwd=2)
points(701:1000,rep(-2.5,300),col="red",type="l",lty=2,lwd=2)
legend("bottomright",legend=c("estimate", "95% CI","true"),
       pch = c(NA,15,NA), lty=c(1,1,2), col=c("black","grey","red"),
       bty = "n", # remove the bounder of the legend
       lwd = c(1,NA,2), pt.bg = c(NA,"grey90",NA),cex=2.5)
dev.off()

lag1_summary=data.frame(mean=rowMeans(lag1[-1,]),t(apply(lag1[-1,],1,quantile,probs=c(0.025,0.975))))
colnames(lag1_summary)=c("mean","2.5%","97.5%")
pdf(file=paste(address_appendix,"lag1_simu.pdf",sep=""),width=11,height=8)
par(mar = c(5, 6, 1, 1))
plot(2:1000,lag1_summary$mean,ylim=c(-7,-1),type="l",
     bty="n",cex.lab=3,cex.axis=3,
     ylab="1-lag effect",xlab="time")
polygon(c(2:1000,rev(2:1000)),c(lag1_summary$`2.5%`,rev(lag1_summary$`97.5%`)),col="grey90",border="grey")
points(2:1000,lag1_summary$mean,type="l")
points(1:1000,lag1_ture,col="red",type="l",lty=2,lwd=2)
legend("topright",legend=c("estimate", "95% CI","true"),
       pch = c(NA,15,NA), lty=c(1,1,2), col=c("black","grey","red"),
       bty = "n", # remove the bounder of the legend
       lwd = c(1,NA,2), pt.bg = c(NA,"grey90",NA),cex=2.5)
dev.off()

lag2_summary=data.frame(mean=rowMeans(lag2[-c(1:2),]),t(apply(lag2[-c(1:2),],1,quantile,probs=c(0.025,0.975))))
colnames(lag2_summary)=c("mean","2.5%","97.5%")
pdf(file=paste(address_appendix,"lag2_simu.pdf",sep=""),width=11,height=8)
par(mar = c(5, 6, 1, 1))
plot(3:1000,lag2_summary$mean,ylim=c(-7,-1),type="l",
     bty="n",cex.lab=3,cex.axis=3,
     ylab="2-lag effect",xlab="time")
polygon(c(3:1000,rev(3:1000)),c(lag2_summary$`2.5%`,rev(lag2_summary$`97.5%`)),col="grey90",border="grey")
points(3:1000,lag2_summary$mean,type="l")
points(1:1000,lag2_ture,col="red",type="l",lty=2,lwd=2)
legend("bottomright",legend=c("estimate", "95% CI","true"),
       pch = c(NA,15,NA), lty=c(1,1,2), col=c("black","grey","red"),
       bty = "n", # remove the bounder of the legend
       lwd = c(1,NA,2), pt.bg = c(NA,"grey90",NA),cex=2.5)
dev.off()

lag3_summary=data.frame(mean=rowMeans(lag3[-c(1:3),]),t(apply(lag3[-c(1:3),],1,quantile,probs=c(0.025,0.975))))
colnames(lag3_summary)=c("mean","2.5%","97.5%")
pdf(file=paste(address_appendix,"lag3_simu.pdf",sep=""),width=11,height=8)
par(mar = c(5, 6, 1, 1))
plot(4:1000,lag3_summary$mean,ylim=c(-7,-1),type="l",
     bty="n",cex.lab=3,cex.axis=3,
     ylab="3-lag effect",xlab="time")
polygon(c(4:1000,rev(4:1000)),c(lag3_summary$`2.5%`,rev(lag3_summary$`97.5%`)),col="grey90",border="grey")
points(4:1000,lag3_summary$mean,type="l")
points(1:1000,lag3_ture,col="red",type="l",lty=2,lwd=2)
legend("bottomright",legend=c("estimate", "95% CI","true"),
       pch = c(NA,15,NA), lty=c(1,1,2), col=c("black","grey","red"),
       bty = "n", # remove the bounder of the legend
       lwd = c(1,NA,2), pt.bg = c(NA,"grey90",NA),cex=2.5)
dev.off()

lag4_summary=data.frame(mean=rowMeans(lag4[-c(1:4),]),t(apply(lag4[-c(1:4),],1,quantile,probs=c(0.025,0.975))))
colnames(lag4_summary)=c("mean","2.5%","97.5%")
pdf(file=paste(address_appendix,"lag4_simu.pdf",sep=""),width=11,height=8)
par(mar = c(5, 6, 1, 1))
plot(5:1000,lag4_summary$mean,ylim=c(-7,-1),type="l",
     bty="n",cex.lab=3,cex.axis=3,
     ylab="4-lag effect",xlab="time")
polygon(c(5:1000,rev(5:1000)),c(lag4_summary$`2.5%`,rev(lag4_summary$`97.5%`)),col="grey90",border="grey")
points(5:1000,lag4_summary$mean,type="l")
points(1:1000,lag4_ture,col="red",type="l",lty=2,lwd=2)
legend("bottomright",legend=c("estimate", "95% CI","true"),
       pch = c(NA,15,NA), lty=c(1,1,2), col=c("black","grey","red"),
       bty = "n", # remove the bounder of the legend
       lwd = c(1,NA,2), pt.bg = c(NA,"grey90",NA),cex=2.5)
dev.off()

estimand_qlag_summary=data.frame(mean=rowMeans(estimand_qlag),t(apply(estimand_qlag,1,quantile,probs=c(0.025,0.975))))
colnames(estimand_qlag_summary)=c("mean","2.5%","97.5%")
pdf(file=paste(address_appendix,"impulse_impact_simu.pdf",sep=""),width=11,height=8)
par(mar = c(5, 6, 1, 1))
nlag=7
plot(1:nlag,estimand_qlag_summary$mean[1:nlag],type="l",
     bty="n",cex.lab=3,cex.axis=3,
     ylab="q-lag effect at t=500",xlab="# lags")
polygon(c(1:nlag,rev(1:nlag)),c(estimand_qlag_summary$`2.5%`[1:nlag],rev(estimand_qlag_summary$`97.5%`[1:nlag])),col="grey90",border="grey")
points(1:nlag,estimand_qlag_summary$mean[1:nlag],type="l")
points(1:nlag,estimand_qlag_ture,col="red",type="l",lty=2,lwd=2)
legend("bottomright",legend=c("estimate", "95% CI","true"),
       pch = c(NA,15,NA), lty=c(1,1,2), col=c("black","grey","red"),
       bty = "n", # remove the bounder of the legend
       lwd = c(1,NA,2), pt.bg = c(NA,"grey90",NA),cex=2.5)
dev.off()

estimand_qstep_total_summary=data.frame(mean=rowMeans(estimand_qstep_total),t(apply(estimand_qstep_total,1,quantile,probs=c(0.025,0.975))))
colnames(estimand_qstep_total_summary)=c("mean","2.5%","97.5%")
pdf(file=paste(address_appendix,"step_response_simu.pdf",sep=""),width=11,height=8)
par(mar = c(5, 6, 1, 1))
nstep=7
plot(1:nstep,estimand_qstep_total_summary$mean[1:nstep],type="l",
     bty="n",cex.lab=3,cex.axis=3,
     ylab="q-step total effect at t=500",xlab="# lags")
polygon(c(1:nstep,rev(1:nstep)),c(estimand_qstep_total_summary$`2.5%`[1:nstep],rev(estimand_qstep_total_summary$`97.5%`[1:nstep])),col="grey90",border="grey")
points(1:nstep,estimand_qstep_total_summary$mean[1:nstep],type="l")
points(1:nstep,estimand_qstep_total_ture,col="red",type="l",lty=2,lwd=2)
legend("topright",legend=c("estimate", "95% CI","true"),
       pch = c(NA,15,NA), lty=c(1,1,2), col=c("black","grey","red"),
       bty = "n", # remove the bounder of the legend
       lwd = c(1,NA,2), pt.bg = c(NA,"grey90",NA),cex=2.5)
dev.off()

