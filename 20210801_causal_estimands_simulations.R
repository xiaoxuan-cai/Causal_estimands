# File: outcome regression using statespace model
# Date: 2021.08.01
source("/Users/xiaoxuancai/Documents/GitHub/mHealth_data_processing/Summary 1_packages.R")
source("/Users/xiaoxuancai/Documents/GitHub/mHealth_data_processing/Summary 2_helper functions.R")
source("/Users/xiaoxuancai/Documents/GitHub/Causal_estimands/helper_causal estimands.R")

#################################################################################################################
##        stable coefficients + treatment/outcome-covariate feedback + continuous tx/outcome/covariate         ##
#################################################################################################################
# case 2: Y_t= 35 + 0.5 Y_{t-1} - 2 X_t - X_{t-1} - C_t + noise (var=0.1)
#         X_t= 0  + 0.1 C_{t} + 0.3 X_{t-1} + 0.05 Y_{t-1} + noise (var=0.1)
#         C_t= 7  + 0.2 C_{t-1} - 0.2 Y_{t-1} + X_{t-1} + noise (var=0.1)
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
print(head(data))
plot(simulation$c);plot(simulation$x);plot(simulation$y);plot(simulation$noise_y)
param=list(list(lagged_param=list(variables=c("y","x","c"),param=rep(7,3))))
data=add_variables_procedures(data,param)

## ------------------------------------------- ##
##  time-invariant linear regression analysis  ##
## ------------------------------------------- ##
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
head(y_coeffi_table);head(y_coeffi_var_table)
c_coeffi_table=as.data.frame(matrix(rep(result_c$coefficients[,1],sum(periods)),nrow=sum(periods),byrow=T))
colnames(c_coeffi_table)=names(result_c$coefficients[,1])
c_coeffi_var_table=lapply(1:sum(periods),function(i){vcov(result_c)})
head(c_coeffi_table);head(c_coeffi_var_table)

# Calculate time-invariant causal estimands for t=10
# (the choice of time point is irrelevant, as the coefficients are time-invariate)
# 
# Method 1: calculate causal estimands (withoutCI) using analytical expressions
# Note:
# - Confidence intervals (CI) can also be computed if needed.
# - The new unified function `calculate.causaleffect` replaces the older versions:
#     `calculate.effect_singlet`, `calculate.effect_allt`, and `calculate.effect_allt_withCI`
# - This new version provides identical results but with improved flexibility,
#   allowing estimation at arbitrary time point(s), with or without CI, in a single function.
time=10
# Y(1) - Y(0): -1.970228
estimand1=calculate.causaleffect(t=time,tx=1,y_coeffi_table,c_coeffi_table)-calculate.causaleffect(t=time,tx=0,y_coeffi_table,c_coeffi_table);estimand1
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

# Method 2: from simulated version
one_time=10
estimand_1to5=simulate.counterfactual_path_singlet(t=one_time,tx=c(1,0,0,0,0),y_coeffi_table,c_coeffi_table,raw_data=data[,c(1:4)])-
  simulate.counterfactual_path_singlet(t=one_time,tx=c(0,0,0,0,0),y_coeffi_table,c_coeffi_table,raw_data=data[,c(1:4)])
estimand_1to5 # the same as c(estimand1,estimand2,estimand3,estimand4,estimand5)
estimand_1to5_withCI=simulate.counterfactual_path_singlet(t=one_time,tx=c(1,0,0,0,0),y_coeffi_table,c_coeffi_table,raw_data=data[,c(1:4)],CI=T,n_sim=100,y_coeffi_var_table,c_coeffi_var_table,seed=1)-
  simulate.counterfactual_path_singlet(t=one_time,tx=c(0,0,0,0,0),y_coeffi_table,c_coeffi_table,raw_data=data[,c(1:4)],CI=T,n_sim=100,y_coeffi_var_table,c_coeffi_var_table,seed=1)
boxplot(estimand_1to5_withCI)
points(estimand_1to5,col="blue",pch=16,type="l",lty=3)


# plots effects with CI for contemporaneous effect and q-lag effects for q=1,2,3,4 for all time points
qlag_result_alltimepoints=lapply(1:nrow(data),function(ts){simulate.counterfactual_path_singlet(t=ts,tx=c(1,0,0,0,0),y_coeffi_table,c_coeffi_table,raw_data=data[,c(1:4)],CI=T,n_sim=10,y_coeffi_var_table,c_coeffi_var_table,seed=1,printFlag=F)-
    simulate.counterfactual_path_singlet(t=ts,tx=c(0,0,0,0,0),y_coeffi_table,c_coeffi_table,raw_data=data[,c(1:4)],CI=T,n_sim=10,y_coeffi_var_table,c_coeffi_var_table,seed=1,printFlag=F)})

# plots effects with CI for contemporaneous effect and q-step total effects for q=1,2,3,4 for all time points
{
  lag1_analytical_qstep=list.cbind(lapply(result_alltimepoints_analytical_qstep,function(x){x[,2]}))
  lag1_CIband_analytical_qstep=plot_simulatedCI(t(lag1_analytical_qstep),probs=c(0.05,0.95),printFlag=F)
  plot(1:1000,lag1_CIband_analytical_qstep$mean,type="l",ylab="1-step total effect",xlab="Date",bty="n",cex.axis=2,cex.lab=2,ylim=c(-15,1.5))
  polygon(c(1:1000,rev(1:1000)),c(lag1_CIband_analytical_qstep$upper,rev(lag1_CIband_analytical_qstep$lower)),col="grey90",border="grey")
  points(1:1000,lag1_CIband_analytical_qstep$mean,type="l")
  abline(h=0,lty=3,lwd=2)
  legend("bottomright",legend=c("estimate", "95% CI"),
         pch = c(NA,15), lty=c(1,NA), col=c("black","grey"),
         bty = "n", # remove the bounder of the legend
         lwd = c(1,NA), pt.bg = c(NA,"grey90"),cex=2)
  
  lag2_analytical_qstep=list.cbind(lapply(result_alltimepoints_analytical_qstep,function(x){x[,3]}))
  lag2_CIband_analytical_qstep=plot_simulatedCI(t(lag2_analytical_qstep),probs=c(0.05,0.95),printFlag=F)
  plot(1:1000,lag2_CIband_analytical_qstep$mean,type="l",ylab="2-step total effect",xlab="",bty="n",cex.axis=2,cex.lab=2,ylim=c(-15,1.5))
  polygon(c(1:1000,rev(1:1000)),c(lag2_CIband_analytical_qstep$upper,rev(lag2_CIband_analytical_qstep$lower)),col="grey90",border="grey")
  points(1:1000,lag2_CIband_analytical_qstep$mean,type="l")
  abline(h=0,lty=3,lwd=2)
  legend("bottomright",legend=c("estimate", "95% CI"),
         pch = c(NA,15), lty=c(1,NA), col=c("black","grey"),
         bty = "n", # remove the bounder of the legend
         lwd = c(1,NA), pt.bg = c(NA,"grey90"),cex=2)
  

  t=600;n_sim=100
  set.seed(1)
  qstep_effects_withCI_part1_ssm_600=simulate.counterfactual_outcome_singlet_withCI(t,tx=rep(1,21),y_coeffi_table,y_coeffi_var_table,c_coeffi_table,c_coeffi_var_table,raw_data,n_sim=n_sim,printFlag=T)
  set.seed(1)
  qstep_effects_withCI_part2_ssm_600=simulate.counterfactual_outcome_singlet_withCI(t,tx=rep(0,21),y_coeffi_table,y_coeffi_var_table,c_coeffi_table,c_coeffi_var_table,raw_data,n_sim=n_sim,printFlag=T)
  stepeffects_CIband_simulated_600=plot_simulatedCI(qstep_effects_withCI_part1_ssm_600-qstep_effects_withCI_part2_ssm_600,probs=c(0.05,0.95),printFlag=F)

  stepeffects_CIband_simulated_600$mean[16]/stepeffects_CIband_simulated_600$mean[20]
  
  plot(0:20,stepeffects_CIband_simulated_600$mean,type="l",ylab="q-step total effect",xlab="# steps",bty="n",cex.axis=2,cex.lab=2,ylim=c(-20,0.5))
  polygon(c(0:20,rev(0:20)),c(stepeffects_CIband_simulated_600$upper,rev(stepeffects_CIband_simulated_600$lower)),col="grey90",border="grey")
  points(0:20,stepeffects_CIband_simulated_600$mean,type="l")
  points(0:20,stepeffects_CIband_simulated_600$mean,pch=16)
  abline(h=0,lty=3,lwd=2)
  abline(v=c(7,11,15),col="red",lty=3,lwd=2)
  text(7,-1.5,"85%",cex=2,col="red")
  text(11,-1.5,"95%",cex=2,col="red")
  text(15,-1.5,"99%",cex=2,col="red")
  legend("bottomright",legend=c("estimate", "95% CI"),
         pch = c(NA,15), lty=c(1,NA), col=c("black","grey"),
         bty = "n", # remove the bounder of the legend
         lwd = c(1,NA), pt.bg = c(NA,"grey90"),cex=2)
  
  
}
t=600;n_sim=100
set.seed(1)
qstep1001001_ssm_600=simulate.counterfactual_outcome_singlet_withCI(t,tx=c(1,0,0,1,0,0,1,rep(0,14)),y_coeffi_table,y_coeffi_var_table,c_coeffi_table,c_coeffi_var_table,raw_data,n_sim=n_sim,printFlag=T)
set.seed(1)
qstep0101010_ssm_600=simulate.counterfactual_outcome_singlet_withCI(t,tx=c(0,1,0,1,0,1,0,rep(0,14)),y_coeffi_table,y_coeffi_var_table,c_coeffi_table,c_coeffi_var_table,raw_data,n_sim=n_sim,printFlag=T)
set.seed(1)
qstep1110000_ssm_600=simulate.counterfactual_outcome_singlet_withCI(t,tx=c(1,1,1,0,0,0,0,rep(0,14)),y_coeffi_table,y_coeffi_var_table,c_coeffi_table,c_coeffi_var_table,raw_data,n_sim=n_sim,printFlag=T)
set.seed(1)
qstep0000000_ssm_600=simulate.counterfactual_outcome_singlet_withCI(t,tx=c(0,0,0,0,0,0,0,rep(0,14)),y_coeffi_table,y_coeffi_var_table,c_coeffi_table,c_coeffi_var_table,raw_data,n_sim=n_sim,printFlag=T)
qstep1001001_CIband_simulated_600=plot_simulatedCI(qstep1001001_ssm_600-qstep0000000_ssm_600,probs=c(0.05,0.95),printFlag=F)
qstep0101010_CIband_simulated_600=plot_simulatedCI(qstep0101010_ssm_600-qstep0000000_ssm_600,probs=c(0.05,0.95),printFlag=F)
qstep1110000_CIband_simulated_600=plot_simulatedCI(qstep1110000_ssm_600-qstep0000000_ssm_600,probs=c(0.05,0.95),printFlag=F)
plot(0:20,qstep1110000_CIband_simulated_600$mean,type="l",ylab="general effect of 1-week intervention",xlab="# days",bty="n",cex.axis=2,cex.lab=1.5,ylim=c(-10,0.5))
points(0:20,qstep1110000_CIband_simulated_600$mean,pch=15)
points(0:20,qstep0101010_CIband_simulated_600$mean,type="l",pch=17,lty=2, col="brown")
points(0:20,qstep0101010_CIband_simulated_600$mean,pch=16,lty=3, col="brown")
points(0:20,qstep1001001_CIband_simulated_600$mean,type="l",pch=15,lty=4,col="blue")
points(0:20,qstep1001001_CIband_simulated_600$mean,pch=17,lty=2,col="blue")
abline(v=6,lty=2,lwd=2,col="red")
legend("bottomright",legend=c("tx=(1,1,1,0,0,0,0)","tx=(0,1,0,1,0,1,0)","tx=(1,0,0,1,0,0,1)"),
       pch = c(15,16,17), lty=c(1,2,4), col=c("black","brown","blue"),
       bty = "n", # remove the bounder of the legend
       lwd = c(1,NA), pt.bg = c(NA,"grey90"),cex=2)











## ------------------------------------------- ##
##  time-variant state space model analysis    ##
## ------------------------------------------- ##
# fit state space model for y
X=data[-1,c("y_1","x","x_1","c")];head(X)
build=function(u){dlmModReg(X,dV=u[1],m0=c(35,0.5,-2,-1,0.5),dW=rep(0,5))}
outMLE=dlmMLE(data$y[-1],0.1,build)
outF=dlmFilter(data$y[-1],build(outMLE$par))
outS=dlmSmooth(outF)
y_coeffi_table=as.data.frame(outS$s);colnames(y_coeffi_table)=c("(Intercept)","y_1","x","x_1","c")
# fit state space model for c
C=data[-1,c("c_1","x_1","y_1")];head(C)
build_c=function(u){dlmModReg(C,dV=u[1],m0=c(7,0.2,1,-0.2),dW=rep(0,4))}
outMLE_c=dlmMLE(data$c[-1],0.1,build_c)
outF_c=dlmFilter(data$c[-1],build_c(outMLE_c$par))
outS_c=dlmSmooth(outF_c)
c_coeffi_table=as.data.frame(outS_c$s);colnames(c_coeffi_table)=c("(Intercept)","c_1","x_1","y_1")
# causal estimands at time t=10
estimand_1to5_ssm=simulate.counterfactual_outcome_singlet(t=one_time,tx=c(1,0,0,0,0),y_coeffi_table,c_coeffi_table,raw_data=data)-
  simulate.counterfactual_outcome_singlet(t=one_time,tx=c(0,0,0,0,0),y_coeffi_table,c_coeffi_table,raw_data=data)
estimand_1to5_ssm;estimand_1to5
# causal estimands at time t=10 with CI
y_coeffi_var_table=lapply(1:nrow(outS$s),function(i){dlmSvd2var(outS$U.S[[i]],outS$D.S[i,])})
c_coeffi_var_table=lapply(1:nrow(outS_c$s),function(i){dlmSvd2var(outS_c$U.S[[i]],outS_c$D.S[i,])})
set.seed(1)
estimand_1to5_withCI_part1_ssm=simulate.counterfactual_outcome_singlet_withCI(t=one_time,tx=c(1,0,0,0,0),y_coeffi_table,y_coeffi_var_table,c_coeffi_table,c_coeffi_var_table,raw_data,n_sim=n_sim,printFlag=T)
set.seed(1)
estimand_1to5_withCI_part2_ssm=simulate.counterfactual_outcome_singlet_withCI(t=one_time,tx=c(0,0,0,0,0),y_coeffi_table,y_coeffi_var_table,c_coeffi_table,c_coeffi_var_table,raw_data,n_sim=n_sim,printFlag=T)
par(mfrow=c(1,1))
boxplot(estimand_1to5_withCI_part1-estimand_1to5_withCI_part2)
points(estimand_1to5,col="blue",pch=16,type="l",lty=3)
boxplot(estimand_1to5_withCI_part1_ssm-estimand_1to5_withCI_part2_ssm)
points(estimand_1to5_ssm,col="blue",pch=16,type="l",lty=3)

