# File: 20230810_5BT65_causalestimands.R
# Date: 2023.08.10
# Updated: 2023.11.05
# Goal: fit both covariate and outcome SSM regression model -> calculate various causal estimands
# Contents: consider two candidate exposures: 
#             <keycontacts_call_totaldegree> (binary) + <keycontacts_text_reciprocity_degree> (binary)
#             for each of them, fit proper SSM without any imputation (ignore case)
# Explanation: we only consider binary exposure in the paper, therefore, we transform exposure into binary
#              for the pair of <keycontacts_call_totaldegree> + <keycontacts_text_reciprocity_degree>, which both range 0~3
#                 we merge level 3 to level 2 since there are two few of them
# Plots: Figure of raw outcome/exposures/covariates
source("/Users/xiaoxuancai/Documents/GitHub/mHealth_data_processing/helper_extract_combine_data.R")
source("/Users/xiaoxuancai/Documents/GitHub/mHealth_data_processing/Summary 1_packages.R")
source("/Users/xiaoxuancai/Documents/GitHub/mHealth_data_processing/Summary 2_helper functions.R")
source("/Users/xiaoxuancai/Documents/GitHub/SSMimpute/helper_multipleimputation.R")
source("/Users/xiaoxuancai/Documents/GitHub/Causal_estimands/helper_causal estimands.R")
load("/Users/xiaoxuancai/Documents/analysis/subject_5BT65/Data/processed/combined_all_5BT65.RData")
library(xtable)
# -------------------------------------------------------------------- #
#             Organize and add variables needed for analysis           #
# -------------------------------------------------------------------- #
# change values for certain variables
TAM_phone_modified=data_5BT65$TAM_phone
TAM_phone_modified[which(TAM_phone_modified==1)]=0.9999
TAM_phone_modified[which(TAM_phone_modified==0)]=0.0001

table(data_5BT65$keycontacts_call_totaldegree) # merge level 2 to level 1 -> turn into a binary variable
table(data_5BT65$keycontacts_text_reciprocity_degree) # merge level 2&3 to level 1 -> turn into a binary variable
keycontacts_call_totaldegree_binary=data_5BT65$keycontacts_call_totaldegree;keycontacts_call_totaldegree_binary[keycontacts_call_totaldegree_binary==2]=1;table(keycontacts_call_totaldegree_binary)
keycontacts_text_reciprocity_degree_binary=data_5BT65$keycontacts_text_reciprocity_degree;keycontacts_text_reciprocity_degree_binary[keycontacts_text_reciprocity_degree_binary>1]=1;table(keycontacts_text_reciprocity_degree_binary)
data=data.frame(Date=data_5BT65$Date,negative_total=data_5BT65$negative_total,
                keycontacts_call_totaldegree_binary,
                keycontacts_text_reciprocity_degree_binary,
                logit_TAM_phone=log(TAM_phone_modified/(1-TAM_phone_modified)))

param=list(list(lagged_param=list(variables=colnames(data)[-1],param=rep(3,length(colnames(data)[-1])))))
data=add_variables_procedures(data,param);
data=data.frame(intercept=1,data);colnames(data)
# plots of outcome/exposure/covariates (Figure 3)
{
  par(mfrow=c(2,2))
  par(mar = c(4, 5, 0.1, 0.5))
  plot(1:708,data$negative_total,type="l",bty="n",xlab="Time (days)",ylab="Negative mood",cex.lab=1.8,cex.axis=1.4)
  text(10,20.2, "a)",cex=2)
  # abline(v=c(1:708)[is.na(data$negative_total)],col="grey")
  plot(1:708,data$keycontacts_call_totaldegree_binary,ylim=c(0,1.1),bty="n",type="l",xlab="Time (days)",ylab="Call connectivity",cex.lab=1.8,cex.axis=1.4)
  text(10,1.06, "b)",cex=2)
  plot(1:708,data$keycontacts_text_reciprocity_degree_binary,ylim=c(0,1.1),bty="n",type="h",xlab="Time (days)",ylab="Text connectivity",cex.lab=1.8,cex.axis=1.4)
  text(10,1.06, "c)",cex=2)
  plot(1:708,data_5BT65$TAM_phone,ylim=c(0,1.1),bty="n",type="l",xlab="Time (days)",ylab="Phone mobility",cex.lab=1.8,cex.axis=1.4)
  text(10,1.06, "d)",cex=2)
}

##################################################################
######          SSM ignore modeling of covariates            #####
##################################################################
# This paper cannot handle missing data for both covariates and outcome
#   thus, we use SSM ignore temporarily
# covariate regression: logit_TAM_phone ~ intercept + logit_TAM_phone_1  
#                                         + keycontacts_call_totaldegree_binary
#                                         + keycontacts_text_reciprocity_degree_binary
#                                         + negative_total
{
  ######################## Part I:  initial guess of ssm structure ########################
  n_variables=5
  ssm_c_init=SSModel(logit_TAM_phone ~ -1 + SSMregression(~ intercept + logit_TAM_phone_1+  
                                                            keycontacts_call_totaldegree_binary+
                                                            keycontacts_text_reciprocity_degree_binary+
                                                            negative_total,  
                                                            data=data,Q=diag(rep(NA,n_variables))),
                    data,H=NA) # baseline variance constant
  ssm_c_init_fit=fitSSM(ssm_c_init, inits =rep(0,n_variables+1), method = "BFGS") # n_variables + one NA in H
  ssm_c_init_out=KFS(ssm_c_init_fit$model)
  plot.KFS(ssm_c_init_out,range=30:708)
  ssm_c_init_out$model$H;ssm_c_init_out$model$Q 
  # check intercept(0.01162), logit_TAM_phone_1(0.00050), keycontacts_text_reciprocity_degree_binary (0.00051)
  
  ######################## Part II: statespace model without imputation ########################
  formula="logit_TAM_phone ~ intercept + logit_TAM_phone_1 + keycontacts_call_totaldegree_binary + keycontacts_text_reciprocity_degree_binary + negative_total"
  all_var=unlist(strsplit(gsub(" ","",formula),"\\+|\\~"))
  outcome_var=strsplit(gsub(" ","",formula),"\\~")[[1]][1]
  formula_var=all_var[!(all_var %in% c(outcome_var,"intercept"))]
  data_ss_ignore=data[,all_var]
  colnames(data_ss_ignore)=gsub(outcome_var,"y",colnames(data_ss_ignore))
  formula_var=gsub(outcome_var,"y",formula_var)
  data_ss_ignore=data_ss_ignore[-1,]
  data_ss_ignore$y[rowSums(is.na(data_ss_ignore))!=0]=NA
  for(l in formula_var){
    data_ss_ignore[is.na(data_ss_ignore[,l]),l]=mean(data_ss_ignore[,l],na.rm=T)
  }
  head(data_ss_ignore);dim(data_ss_ignore)
  
  ss_param_keycontacts=list(inits=c(log(0.05),log(5)),m0=ssm_c_init_out$att[708,],C0=diag(rep(10^3),5),
                            AR1_coeffi=NULL,rw_coeffi=c("intercept"),v_cp_param=NULL,
                            w_cp_param=list(list(variable="y_1",segments=3,changepoints=c(410,450),fixed_cpts=F),
                                            list(variable="keycontacts_text_reciprocity_degree_binary",segments=2,changepoints=c(415),fixed_cpts=F)))
  ssm_c=run.SSM_partial_converged_cpts(data_ss=data_ss_ignore,formula_var=formula_var,ss_param=ss_param_keycontacts,max_iteration=100,
                                       cpt_learning_param=list(cpt_method="mean",burnin=1/10,mergeband=10,convergence_cri=5),
                                       dlm_cpt_learning_option="filter",
                                       dlm_option="smooth",printFlag=T,bandwidth=5,cpt_V =1)
  ssm_c_smooth=dlmSmooth(ssm_c$out_filter)
  ssm_c$estimated_cpts
  plot.dlm(ssm_c$out_filter,benchmark = rep(0,5),option="filtered_state",range=30:708)
  ssm_c_all=cbind(round(ssm_c$result,digits=3),
                  paste("(",round(ssm_c$result$Estimate-1.965*ssm_c$result$Std.Error,digits=3),",",
                        round(ssm_c$result$Estimate+1.965*ssm_c$result$Std.Error,digits=3),")",sep=""),
                  paste("(",round(ssm_c$result$Estimate-1.645*ssm_c$result$Std.Error,digits=3),",",
                        round(ssm_c$result$Estimate+1.645*ssm_c$result$Std.Error,digits=3),")",sep=""))
  rownames(ssm_c_all)=rownames(ssm_c$result)
  colnames(ssm_c_all)=c("Estimate","Std.Error","95% CI", "90% CI")
  ssm_c_all
  # Estimate Std.Error          95% CI          90% CI
  # (Intercept)                                            1.297     0.600   (0.117,2.476)   (0.309,2.284)
  # y_1(period1)                                           0.066     0.098  (-0.125,0.258)  (-0.094,0.227)
  # y_1(period2)                                           0.178     0.092  (-0.003,0.358)   (0.027,0.329)
  # y_1(period3)                                          -0.109     0.047 (-0.201,-0.018) (-0.186,-0.033)
  # keycontacts_call_totaldegree_binary                    0.012     0.292  (-0.561,0.585)  (-0.468,0.492)
  # keycontacts_text_reciprocity_degree_binary(period1)   -0.059     0.285  (-0.619,0.502)   (-0.528,0.41)
  # keycontacts_text_reciprocity_degree_binary(period2)   -0.781     0.304 (-1.379,-0.183)  (-1.282,-0.28)
  # negative_total                                        -0.014     0.026  (-0.066,0.038)  (-0.058,0.029)
  xtable(ssm_c_all)
}
# ssm_c plots
{  
  # random walk intercept
  par(mar = c(5, 5, .5, .5))
  plot(1:708,ssm_c$out_filter$m[,1],type="l",ylab=expression(paste(mu["0,t"]," (intercept)")),xlab="Time (days)",bty="n",cex.lab=2,cex.axis=2,ylim=c(-5,5))
  sd=unlist(lapply(1:708,function(i){sqrt(diag(dlmSvd2var(ssm_c$out_filter$U.C[[i]],ssm_c$out_filter$D.C[i,]))[1])}))
  polygon(c(1:708,rev(1:708)),c(ssm_c$out_filter$m[,1]+1.65*sd,rev(ssm_c$out_filter$m[,1]-1.65*sd)),col="grey90",border="grey")
  points(1:708,ssm_c$out_filter$m[,1],type="l")
  abline(h=0,col="black",lty=3,lwd=2)
  legend("bottomright",legend=c("estimate", "95% CI"),
         pch = c(NA,15), lty=c(1,NA), col=c("black","grey"),
         bty = "n", # remove the bounder of the legend
         lwd = c(1,NA), pt.bg = c(NA,"grey90"),cex=2)
  
  # 3 period autocorrelation
  a=cpt.mean(ssm_c$out_filter$m[,2],penalty="Manual",Q=4,method="BinSeg")
  plot(1:708,ssm_c$out_filter$m[,2],type="l",ylab=expression(paste(rho["pa,t"]," (auto-correlation)")),xlab="Time (days)",bty="n",cex.lab=2,cex.axis=2,xaxt="n")
  sd=unlist(lapply(1:708,function(i){sqrt(diag(dlmSvd2var(ssm_c$out_filter$U.C[[i]],ssm_c$out_filter$D.C[i,]))[2])}))
  polygon(c(1:708,rev(1:708)),c(ssm_c$out_filter$m[,2]+1.65*sd,rev(ssm_c$out_filter$m[,2]-1.65*sd)),col="grey90",border="grey")
  points(1:708,ssm_c$out_filter$m[,2],type="l")
  lines(c(1:708)[c(1,cpts(a)[3])],rep(param.est(a)$mean[3],2),col="red")
  lines(c(1:708)[c(cpts(a)[3]+1,cpts(a)[4])],rep(param.est(a)$mean[4],2),col="red")
  lines(c(1:708)[c(cpts(a)[4]+1,708)],rep(param.est(a)$mean[5],2),col="red")
  abline(h=0,col="black",lty=3,lwd=2)
  legend("topright",legend=c("estimate", "95% CI"),
         pch = c(NA,15), lty=c(1,NA), col=c("black","grey"),
         bty = "n", # remove the bounder of the legend
         lwd = c(1,NA), pt.bg = c(NA,"grey90"),cex=2)
  
  # time-invariant coeffi for degree of calls with keycontacts
  plot(1:708,ssm_c$out_filter$m[,3],type="l",ylab=expression(paste(mu["1,t"]," (degree of calls)")),xlab="Time (days)",bty="n",cex.lab=2,cex.axis=2,ylim=c(-1,1))
  sd=unlist(lapply(1:708,function(i){sqrt(diag(dlmSvd2var(ssm_c$out_filter$U.C[[i]],ssm_c$out_filter$D.C[i,]))[3])}))
  polygon(c(1:708,rev(1:708)),c(ssm_c$out_filter$m[,3]+1.65*sd,rev(ssm_c$out_filter$m[,3]-1.65*sd)),col="grey90",border="grey")
  points(1:708,ssm_c$out_filter$m[,3],type="l")
  abline(h=0,col="black",lty=3,lwd=2)
  legend("bottomright",legend=c("estimate", "95% CI"),
         pch = c(NA,15), lty=c(1,NA), col=c("black","grey"),
         bty = "n", # remove the bounder of the legend
         lwd = c(1,NA), pt.bg = c(NA,"grey90"),cex=2)
  
  # periodic-stable coeff for degree of texts with keycontacts
  a=cpt.mean(ssm_c$out_filter$m[,4],penalty="Manual",Q=2,method="BinSeg")
  plot(1:708,ssm_c$out_filter$m[,4],type="l",ylab=expression(paste(mu["2,t"]," (degree of texts)")),xlab="Time (days)",bty="n",cex.lab=2,cex.axis=2, ylim=c(-5,1))
  sd=unlist(lapply(1:708,function(i){sqrt(diag(dlmSvd2var(ssm_c$out_filter$U.C[[i]],ssm_c$out_filter$D.C[i,]))[4])}))
  polygon(c(1:708,rev(1:708)),c(ssm_c$out_filter$m[,4]+1.65*sd,rev(ssm_c$out_filter$m[,4]-1.65*sd)),col="grey90",border="grey")
  points(1:708,ssm_c$out_filter$m[,4],type="l")
  lines(c(1:708)[c(1,cpts(a)[2])],rep(param.est(a)$mean[2],2),col="red")
  lines(c(1:708)[c(cpts(a)[2]+1,708)],rep(param.est(a)$mean[3],2),col="red")
  abline(h=0,col="black",lty=3,lwd=2)
  legend("bottomright",legend=c("estimate", "95% CI"),
         pch = c(NA,15), lty=c(1,NA), col=c("black","grey"),
         bty = "n", # remove the bounder of the legend
         lwd = c(1,NA), pt.bg = c(NA,"grey90"),cex=2)
  
  # time-invariant coeffi for negative mood
  plot(1:708,ssm_c$out_filter$m[,5],type="l",ylab=expression(paste(mu["3,t"]," (negative mood)")),xlab="Time (days)",bty="n",cex.lab=2,cex.axis=2,ylim=c(-0.7,0.5))
  sd=unlist(lapply(1:708,function(i){sqrt(diag(dlmSvd2var(ssm_c$out_filter$U.C[[i]],ssm_c$out_filter$D.C[i,]))[5])}))
  polygon(c(1:708,rev(1:708)),c(ssm_c$out_filter$m[,5]+1.65*sd,rev(ssm_c$out_filter$m[,5]-1.65*sd)),col="grey90",border="grey")
  points(1:708,ssm_c$out_filter$m[,5],type="l")
  abline(h=0,col="black",lty=3,lwd=2)
  legend("bottomright",legend=c("estimate", "95% CI"),
         pch = c(NA,15), lty=c(1,NA), col=c("black","grey"),
         bty = "n", # remove the bounder of the legend
         lwd = c(1,NA), pt.bg = c(NA,"grey90"),cex=2)
}
# outcome regression: negative_total~ intercept + negative_total_1 
#                                     + keycontacts_call_totaldegree_binary + keycontacts_call_totaldegree_binary_1
#                                     + keycontacts_text_reciprocity_degree_binary + keycontacts_text_reciprocity_degree_binary
#                                     + logit_TAM_phone_1
{
  ######################## Part I: initial guess of ssm structure ########################
  n_variables=7
  ssm_y_init=SSModel(negative_total ~ -1 + SSMregression(~ intercept + negative_total_1+  
                                                           keycontacts_call_totaldegree_binary+keycontacts_call_totaldegree_binary_1+
                                                           keycontacts_text_reciprocity_degree_binary+keycontacts_text_reciprocity_degree_binary_1+
                                                           logit_TAM_phone_1,  
                                                           data=data,Q=diag(rep(NA,n_variables))),
                        data,H=NA) # baseline variance constant
  ssm_y_init_fit=fitSSM(ssm_y_init, inits =rep(0,n_variables+1), method = "BFGS") # n_variables + one NA in H
  ssm_y_init_out=KFS(ssm_y_init_fit$model)
  plot.KFS(ssm_y_init_out,range=30:708)
  ssm_y_init_out$model$Q; ssm_y_init_out$model$H
  # check: intercept(0.16839), keycontacts_text_reciprocity_degree_binary (0.001768), keycontacts_text_reciprocity_degree_binary_1 (0.00064)
  
  ######################## Part II: statespace model without imputation ########################
  formula=" negative_total ~ intercept + negative_total_1 + keycontacts_call_totaldegree_binary + keycontacts_call_totaldegree_binary_1 + keycontacts_text_reciprocity_degree_binary + keycontacts_text_reciprocity_degree_binary_1+ logit_TAM_phone_1"
  all_var=unlist(strsplit(gsub(" ","",formula),"\\+|\\~"))
  outcome_var=strsplit(gsub(" ","",formula),"\\~")[[1]][1]
  formula_var=all_var[!(all_var %in% c(outcome_var,"intercept"))]
  data_ss_ignore_y=data[,all_var]
  colnames(data_ss_ignore_y)=gsub(outcome_var,"y",colnames(data_ss_ignore_y))
  formula_var=gsub(outcome_var,"y",formula_var)
  data_ss_ignore_y=data_ss_ignore_y[-1,]
  data_ss_ignore_y$y[rowSums(is.na(data_ss_ignore_y))!=0]=NA
  for(l in formula_var){
    data_ss_ignore_y[is.na(data_ss_ignore_y[,l]),l]=mean(data_ss_ignore_y[,l],na.rm=T)
  }
  
  ss_param_temp=list(inits=c(log(0.2),log(2.5)),m0=ssm_y_init$att[708,],C0=diag(rep(10^3),7),
                     AR1_coeffi=NULL,rw_coeffi=c("intercept"),v_cp_param=NULL,
                     w_cp_param=list(list(variable="keycontacts_text_reciprocity_degree_binary",segments=3,changepoints=c(560,640),fixed_cpts=F),
                                     list(variable="keycontacts_text_reciprocity_degree_binary_1",segments=2,changepoints=c(460),fixed_cpts=F)))
  ssm_y=run.SSM_partial_converged_cpts(data_ss=data_ss_ignore_y,formula_var=formula_var,ss_param=ss_param_temp,max_iteration=100,
                                       cpt_learning_param=list(cpt_method="mean",burnin=1/10,mergeband=20,convergence_cri=10),
                                       dlm_cpt_learning_option="filter",
                                       dlm_option="smooth",printFlag=T,bandwidth = 5, cpt_V = 1)
  ssm_y$estimated_cpts
  ssm_y_smooth=dlmSmooth(ssm_y$out_filter)
  plot.dlm(ssm_y$out_filter,benchmark = rep(0,7),option="filtered_state",range=30:708)
  ssm_y_all=cbind(round(ssm_y$result,digits=3),
                  paste("(",round(ssm_y$result$Estimate-1.965*ssm_y$result$Std.Error,digits=3),",",
                         round(ssm_y$result$Estimate+1.965*ssm_y$result$Std.Error,digits=3),")",sep=""),
                  paste("(",round(ssm_y$result$Estimate-1.645*ssm_y$result$Std.Error,digits=3),",",
                         round(ssm_y$result$Estimate+1.645*ssm_y$result$Std.Error,digits=3),")",sep=""))
  rownames(ssm_y_all)=rownames(ssm_y$result)
  colnames(ssm_y_all)=c("Estimate","Std.Error","95% CI", "90% CI")
  xtable(ssm_y_all)
  ssm_y_all
  #                                                       Estimate Std.Error          95% CI          90% CI
  # (Intercept)                                              5.004     0.880   (3.276,6.733)   (3.557,6.451)
  # y_1                                                      0.633     0.038   (0.559,0.707)   (0.571,0.695)
  # keycontacts_call_totaldegree_merged                     -0.218     0.254   (-0.717,0.28)  (-0.636,0.199)
  # keycontacts_call_totaldegree_merged_1                    0.088     0.257  (-0.418,0.593)  (-0.335,0.511)
  # keycontacts_text_reciprocity_degree_merged(period1)     -0.087     0.210  (-0.499,0.326)  (-0.432,0.258)
  # keycontacts_text_reciprocity_degree_merged(period2)     -1.154     0.608  (-2.349,0.041) (-2.154,-0.154)
  # keycontacts_text_reciprocity_degree_merged(period3)     -0.739     0.583  (-1.886,0.408)  (-1.699,0.221)
  # keycontacts_text_reciprocity_degree_merged_1(period1)   -0.165     0.238  (-0.632,0.302)  (-0.556,0.226)
  # keycontacts_text_reciprocity_degree_merged_1(period2)   -0.724     0.312 (-1.338,-0.111) (-1.238,-0.211)
  # logit_TAM_phone_1                                       -0.012     0.036  (-0.083,0.058)  (-0.072,0.047)
}
# ssm_y plots
{
  # random walk intercept
  par(mar = c(5, 5, .5, .5))
  plot(1:708,ssm_y$out_filter$m[,1],type="l",ylab=expression(paste(beta["0,t"]," (intercept)")),xlab="Time (days)" ,bty="n",cex.lab=2,cex.axis=2, ylim=c(-3,8))
  sd=unlist(lapply(1:708,function(i){sqrt(diag(dlmSvd2var(ssm_y$out_filter$U.C[[i]],ssm_y$out_filter$D.C[i,]))[1])}))
  polygon(c(1:708,rev(1:708)),c(ssm_y$out_filter$m[,1]+1.65*sd,rev(ssm_y$out_filter$m[,1]-1.65*sd)),col="grey90",border="grey")
  points(1:708,ssm_y$out_filter$m[,1],type="l")
  abline(h=0,col="black",lty=3,lwd=2)
  legend("bottomright",legend=c("estimate", "95% CI"),
         pch = c(NA,15), lty=c(1,NA), col=c("black","grey"),
         bty = "n", # remove the bounder of the legend
         lwd = c(1,NA), pt.bg = c(NA,"grey90"),cex=2)
  
  # time-invariant autocorrelation
  plot(1:708,ssm_y$out_filter$m[,2],type="l",ylab=expression(paste(rho[t]," (auto-correlation)")),xlab="Time (days)" ,bty="n",cex.lab=2,cex.axis=2, ylim=c(0.0,1))
  sd=unlist(lapply(1:708,function(i){sqrt(diag(dlmSvd2var(ssm_y$out_filter$U.C[[i]],ssm_y$out_filter$D.C[i,]))[2])}))
  polygon(c(1:708,rev(1:708)),c(ssm_y$out_filter$m[,2]+1.65*sd,rev(ssm_y$out_filter$m[,2]-1.65*sd)),col="grey90",border="grey")
  points(1:708,ssm_y$out_filter$m[,2],type="l")
  abline(h=0,col="black",lty=3,lwd=2)
  legend("bottomright",legend=c("estimate", "95% CI"),
         pch = c(NA,15), lty=c(1,NA), col=c("black","grey"),
         bty = "n", # remove the bounder of the legend
         lwd = c(1,NA), pt.bg = c(NA,"grey90"),cex=2)
  
  # time-invariant coeffi for degree of calls with keycontacts
  plot(1:708,ssm_y$out_filter$m[,3],type="l",ylab=expression(paste(beta["11,t"]," (degree of calls)")),xlab="Time (days)" ,bty="n",cex.lab=2,cex.axis=2, ylim=c(-1,1))
  sd=unlist(lapply(1:708,function(i){sqrt(diag(dlmSvd2var(ssm_y$out_filter$U.C[[i]],ssm_y$out_filter$D.C[i,]))[3])}))
  polygon(c(1:708,rev(1:708)),c(ssm_y$out_filter$m[,3]+1.65*sd,rev(ssm_y$out_filter$m[,3]-1.65*sd)),col="grey90",border="grey")
  points(1:708,ssm_y$out_filter$m[,3],type="l")
  abline(h=0,col="black",lty=3,lwd=2)
  legend("bottomright",legend=c("estimate", "95% CI"),
         pch = c(NA,15), lty=c(1,NA), col=c("black","grey"),
         bty = "n", # remove the bounder of the legend
         lwd = c(1,NA), pt.bg = c(NA,"grey90"),cex=2)
  
  # time-invariant coeffi for degree of calls yesterday with keycontacts
  plot(1:708,ssm_y$out_filter$m[,4],type="l",ylab=expression(paste(beta["12,t"]," (degree of calls previous day)")),xlab="Time (days)" ,bty="n",cex.lab=2,cex.axis=2, ylim=c(-1,1))
  sd=unlist(lapply(1:708,function(i){sqrt(diag(dlmSvd2var(ssm_y$out_filter$U.C[[i]],ssm_y$out_filter$D.C[i,]))[4])}))
  polygon(c(1:708,rev(1:708)),c(ssm_y$out_filter$m[,4]+1.65*sd,rev(ssm_y$out_filter$m[,4]-1.65*sd)),col="grey90",border="grey")
  points(1:708,ssm_y$out_filter$m[,4],type="l")
  abline(h=0,col="black",lty=3,lwd=2)
  legend("bottomright",legend=c("estimate", "95% CI"),
         pch = c(NA,15), lty=c(1,NA), col=c("black","grey"),
         bty = "n", # remove the bounder of the legend
         lwd = c(1,NA), pt.bg = c(NA,"grey90"),cex=2)
  
  # periodic-stable coeff for degree of texts  
  a=cpt.mean(ssm_y$out_filter$m[,5],penalty="Manual",Q=2,method="BinSeg")
  plot(1:708,ssm_y$out_filter$m[,5],type="l",ylab=expression(paste(beta["21,t"]," (degree of texts)")),xlab="Time (days)" ,bty="n",cex.lab=2,cex.axis=2, ylim=c(-4,1))
  sd=unlist(lapply(1:708,function(i){sqrt(diag(dlmSvd2var(ssm_y$out_filter$U.C[[i]],ssm_y$out_filter$D.C[i,]))[5])}))
  polygon(c(1:708,rev(1:708)),c(ssm_y$out_filter$m[,5]+1.65*sd,rev(ssm_y$out_filter$m[,5]-1.65*sd)),col="grey90",border="grey")
  points(1:708,ssm_y$out_filter$m[,5],type="l")
  lines(c(1:708)[c(1,cpts(a)[1])],rep(param.est(a)$mean[1],2),col="red")
  lines(c(1:708)[c(cpts(a)[1]+1,cpts(a)[2])],rep(param.est(a)$mean[2],2),col="red")
  lines(c(1:708)[c(cpts(a)[2]+1,708)],rep(param.est(a)$mean[3],2),col="red")
  abline(h=0,col="black",lty=3,lwd=2)
  legend("bottomright",legend=c("estimate", "95% CI"),
         pch = c(NA,15), lty=c(1,NA), col=c("black","grey"),
         bty = "n", # remove the bounder of the legend
         lwd = c(1,NA), pt.bg = c(NA,"grey90"),cex=2)
  
  # time-invariant coeffi for degree of texts yesterday
  a=cpt.mean(ssm_y$out_filter$m[,6],penalty="Manual",Q=1,method="BinSeg")
  plot(1:708,ssm_y$out_filter$m[,6],type="l",ylab=expression(paste(beta["22,t"]," (degree of outgoing texts previous day)")),xlab="Time (days)" ,bty="n",cex.lab=1.5,cex.axis=2, ylim=c(-3,2))
  sd=unlist(lapply(1:708,function(i){sqrt(diag(dlmSvd2var(ssm_y$out_filter$U.C[[i]],ssm_y$out_filter$D.C[i,]))[6])}))
  polygon(c(1:708,rev(1:708)),c(ssm_y$out_filter$m[,6]+1.65*sd,rev(ssm_y$out_filter$m[,6]-1.65*sd)),col="grey90",border="grey")
  points(1:708,ssm_y$out_filter$m[,6],type="l")
  lines(c(1:708)[c(1,cpts(a)[1])],rep(param.est(a)$mean[1],2),col="red")
  lines(c(1:708)[c(cpts(a)[1]+1,708)],rep(param.est(a)$mean[2],2),col="red")
  abline(h=0,col="black",lty=3,lwd=2)
  legend("bottomright",legend=c("estimate", "95% CI"),
         pch = c(NA,15), lty=c(1,NA), col=c("black","grey"),
         bty = "n", # remove the bounder of the legend
         lwd = c(1,NA), pt.bg = c(NA,"grey90"),cex=2)
  
  # time-invariant coeffi for physical activity
  plot(1:708,ssm_y$out_filter$m[,7],type="l",ylab=expression(paste(beta["PA,t"]," (daily physical activity)")),xlab="Time (days)" ,bty="n",cex.lab=2,cex.axis=2, ylim=c(-1,1))
  sd=unlist(lapply(1:708,function(i){sqrt(diag(dlmSvd2var(ssm_y$out_filter$U.C[[i]],ssm_y$out_filter$D.C[i,]))[7])}))
  polygon(c(1:708,rev(1:708)),c(ssm_y$out_filter$m[,7]+1.65*sd,rev(ssm_y$out_filter$m[,7]-1.65*sd)),col="grey90",border="grey")
  points(1:708,ssm_y$out_filter$m[,7],type="l")
  abline(h=0,col="black",lty=3,lwd=2)
  legend("bottomright",legend=c("estimate", "95% CI"),
         pch = c(NA,15), lty=c(1,NA), col=c("black","grey"),
         bty = "n", # remove the bounder of the legend
         lwd = c(1,NA), pt.bg = c(NA,"grey90"),cex=2)
}

############################################################################
######     simulate/compute causal estimands of calls with CI          #####
############################################################################
# jump directly to load the data without recalculation!!
y_coeffi_var_table=lapply(1:nrow(ssm_y$out_smooth$s),function(i){dlmSvd2var(ssm_y$out_smooth$U.S[[i]],ssm_y$out_smooth$D.S[i,])});head(y_coeffi_var_table)
c_coeffi_var_table=lapply(1:nrow(ssm_c$out_smooth$s),function(i){dlmSvd2var(ssm_c$out_smooth$U.S[[i]],ssm_c$out_smooth$D.S[i,])});head(c_coeffi_var_table)
raw_data=data[,c(2,3,4,6)];colnames(raw_data)=c("Date","y","x","c");head(raw_data)
# first, calculate for calls
y_coeffi_table=as.data.frame(ssm_y$out_smooth$s);colnames(y_coeffi_table)[c(1:4,7)]=c("(Intercept)","y_1","x","x_1","c");head(y_coeffi_table)
c_coeffi_table=as.data.frame(ssm_c$out_smooth$s);colnames(c_coeffi_table)[c(1,2,3,5)]=c("(Intercept)","c_1","x_1","y_1");head(c_coeffi_table)
# calculate call's contemporaneous effect directly with CI (save call_contemporaneous)
set.seed(1)
call_contemporaneous = calculate.effect_allt_withCI(tx=c(1),y_coeffi_table=y_coeffi_table,y_coeffi_var_table=y_coeffi_var_table,c_coeffi_table=c_coeffi_table,c_coeffi_var_table=c_coeffi_var_table,n_sim=1000,printFlag=T)
# calculate call's 1-lag effect directly with CI
set.seed(2)
call_1_lag = calculate.effect_allt_withCI(tx=c(1,0),y_coeffi_table=y_coeffi_table,y_coeffi_var_table=y_coeffi_var_table,c_coeffi_table=c_coeffi_table,c_coeffi_var_table=c_coeffi_var_table,n_sim=1000,printFlag=T)
# calculate call's 1-lag structural direct effect directly with CI
set.seed(3)
call_1_lag_structural_direct = calculate.controlled_direct_effect_allt_withCI(y_coeffi_table,y_coeffi_var_table,n_sim=1000,printFlag=T)
# calculate call's 2-step total effect directly with CI
set.seed(4)
call_2_step = calculate.effect_allt_withCI(tx=c(1,1),y_coeffi_table=y_coeffi_table,y_coeffi_var_table=y_coeffi_var_table,c_coeffi_table=c_coeffi_table,c_coeffi_var_table=c_coeffi_var_table,n_sim=1000,printFlag=T)
# calculate call's 3-step general effect directly with CI
set.seed(4)
call_3_step_general_101 = calculate.effect_allt_withCI(tx=c(1,0,1),y_coeffi_table=y_coeffi_table,y_coeffi_var_table=y_coeffi_var_table,c_coeffi_table=c_coeffi_table,c_coeffi_var_table=c_coeffi_var_table,n_sim=1000,printFlag=T)
# impulse impact graph at t=600
t=600;n_sim=1000
set.seed(1)
call_qlag_effects_part1=simulate.counterfactual_singlet_withCI(t,tx=c(1,rep(0,10)),y_coeffi_table,y_coeffi_var_table,c_coeffi_table,c_coeffi_var_table,raw_data,n_sim=n_sim,printFlag=T)
set.seed(1)
call_qlag_effects_part2=simulate.counterfactual_singlet_withCI(t,tx=c(0,rep(0,10)),y_coeffi_table,y_coeffi_var_table,c_coeffi_table,c_coeffi_var_table,raw_data,n_sim=n_sim,printFlag=T)
call_effects_vs_qlag=call_qlag_effects_part1-call_qlag_effects_part2

save(call_contemporaneous,
     call_1_lag,call_1_lag_structural_direct,
     call_2_step,call_3_step_general_101,
     call_effects_vs_qlag,
     file="/Users/xiaoxuancai/Dropbox (Personal)/MHealthPsychSummerProject2020/Xiaoxuan_Cai/[Paper 1] Causal estimands and Graphical representation for time series data/result_calls.Rdata")
load("/Users/xiaoxuancai/Dropbox (Personal)/MHealthPsychSummerProject2020/Xiaoxuan_Cai/[Paper 1] Causal estimands and Graphical representation for time series data/result_calls.Rdata")
# plots
{
  # call's contemporaneous effect
  par(mar = c(2.5, 5, .5, .5))
  call_contemporaneous_CIband=plot_simulatedCI(t(call_contemporaneous),probs=c(0.05,0.95),printFlag=F)
  plot(1:708,call_contemporaneous_CIband$mean,type="l",ylab="contemporaneous effect (calls)",xlab="",bty="n", cex.axis=2,cex.lab=2,ylim=c(-1,0.5))
  polygon(c(1:708,rev(1:708)),c(call_contemporaneous_CIband$upper,rev(call_contemporaneous_CIband$lower)),col="grey90",border="grey")
  points(1:708,call_contemporaneous_CIband$mean,type="l")
  abline(h=0,lty=3,lwd=2)
  legend("bottomright",legend=c("estimate", "95% CI"),
         pch = c(NA,15), lty=c(1,NA), col=c("black","grey"),
         bty = "n", # remove the bounder of the legend
         lwd = c(1,NA), pt.bg = c(NA,"grey90"),cex=2)
  
  # call's 1-lag effect 
  par(mar = c(2.5, 5, .5, .5))
  call_1_lag_CIband=plot_simulatedCI(t(call_1_lag),probs=c(0.05,0.95),printFlag=F)
  plot(1:708,call_1_lag_CIband$mean,type="l",ylab="1-lag effect (calls)",xlab="Date",bty="n", cex.axis=2,cex.lab=2,ylim=c(-1,0.5))
  polygon(c(1:708,rev(1:708)),c(call_1_lag_CIband$upper,rev(call_1_lag_CIband$lower)),col="grey90",border="grey")
  points(1:708,call_1_lag_CIband$mean,type="l")
  abline(h=0,lty=3,lwd=2)
  legend("bottomright",legend=c("estimate", "95% CI"),
         pch = c(NA,15), lty=c(1,NA), col=c("black","grey"),
         bty = "n", # remove the bounder of the legend
         lwd = c(1,NA), pt.bg = c(NA,"grey90"),cex=2)
  
  # unused code for plotting call's q-lag effects for q=2,3,4 for all time points
  # lag2_analytical_calls=list.cbind(lapply(result_alltimepoints_analytical_calls,function(x){x[,3]}))
  # lag2_CIband_analytical_calls=plot_simulatedCI(t(lag2_analytical_calls),probs=c(0.05,0.95),printFlag=F)
  # plot(1:708,lag2_CIband_analytical_calls$mean,type="l",ylab="2-lag effect (calls)",xlab="Date",bty="n", cex.axis=2,cex.lab=2,ylim=c(-1,0.5))
  # polygon(c(1:708,rev(1:708)),c(lag2_CIband_analytical_calls$upper,rev(lag2_CIband_analytical_calls$lower)),col="grey90",border="grey")
  # points(1:708,lag2_CIband_analytical_calls$mean,type="l")
  # abline(h=0,lty=3,lwd=2)
  # legend("bottomright",legend=c("estimate", "95% CI"),
  #        pch = c(NA,15), lty=c(1,NA), col=c("black","grey"),
  #        bty = "n", # remove the bounder of the legend
  #        lwd = c(1,NA), pt.bg = c(NA,"grey90"),cex=2)
  # 
  # lag3_analytical_calls=list.cbind(lapply(result_alltimepoints_analytical_calls,function(x){x[,4]}))
  # lag3_CIband_analytical_calls=plot_simulatedCI(t(lag3_analytical_calls),probs=c(0.05,0.95),printFlag=F)
  # plot(1:708,lag3_CIband_analytical_calls$mean,type="l",ylab="3-lag effect (calls)",xlab="Date",bty="n", cex.axis=2,cex.lab=2,ylim=c(-1,0.5))
  # polygon(c(1:708,rev(1:708)),c(lag3_CIband_analytical_calls$upper,rev(lag3_CIband_analytical_calls$lower)),col="grey90",border="grey")
  # points(1:708,lag3_CIband_analytical_calls$mean,type="l")
  # abline(h=0,lty=3,lwd=2)
  # legend("bottomright",legend=c("estimate", "95% CI"),
  #        pch = c(NA,15), lty=c(1,NA), col=c("black","grey"),
  #        bty = "n", # remove the bounder of the legend
  #        lwd = c(1,NA), pt.bg = c(NA,"grey90"),cex=2)
  # 
  # lag4_analytical_calls=list.cbind(lapply(result_alltimepoints_analytical_calls,function(x){x[,5]}))
  # lag4_CIband_analytical_calls=plot_simulatedCI(t(lag4_analytical_calls),probs=c(0.05,0.95),printFlag=F)
  # plot(1:708,lag4_CIband_analytical_calls$mean,type="l",ylab="4-lag effect (texts)",xlab="Date",bty="n", cex.axis=2,cex.lab=2,ylim=c(-1,0.5))
  # polygon(c(1:708,rev(1:708)),c(lag4_CIband_analytical_calls$upper,rev(lag4_CIband_analytical_calls$lower)),col="grey90",border="grey")
  # points(1:708,lag4_CIband_analytical$mean,type="l")
  # abline(h=0,lty=3,lwd=2)
  # legend("bottomright",legend=c("estimate", "95% CI"),
  #        pch = c(NA,15), lty=c(1,NA), col=c("black","grey"),
  #        bty = "n", # remove the bounder of the legend
  #        lwd = c(1,NA), pt.bg = c(NA,"grey90"),cex=2)
  
  # call's 1-lag structural direct effect
  par(mar = c(2.5, 5, .5, .5))
  call_1_lag_structural_direct_CIband = plot_simulatedCI(t(call_1_lag_structural_direct),probs=c(0.05,0.95),printFlag=F)
  plot(1:708,call_1_lag_structural_direct_CIband$mean,type="l",ylab="1-lag controlled direct effect (calls)",xlab="# lags",bty="n",cex.axis=1.5,cex.lab=1.5,ylim=c(-.7,0.7))
  polygon(c(1:708,rev(1:708)),c(call_1_lag_structural_direct_CIband$upper,rev(call_1_lag_structural_direct_CIband$lower)),col="grey90",border="grey")
  points(1:708,call_1_lag_structural_direct_CIband$mean,type="l")
  abline(h=0,lty=3,lwd=2)
  legend("bottomright",legend=c("estimate", "95% CI"),
         pch = c(NA,15), lty=c(1,NA), col=c("black","grey"),
         bty = "n", # remove the bounder of the legend
         lwd = c(1,NA), pt.bg = c(NA,"grey90"),cex=2)
  
  # call's 2-step total effect
  par(mar = c(2.5, 5, .5, .5))
  call_2_step_CIband = plot_simulatedCI(t(call_2_step),probs=c(0.05,0.95),printFlag=F)
  plot(1:708,call_2_step_CIband$mean,type="l",ylab="2-step total effect (calls)",xlab="# lags",bty="n",cex.axis=1.5,cex.lab=1.5,ylim=c(-2,1))
  polygon(c(1:708,rev(1:708)),c(call_2_step_CIband$upper,rev(call_2_step_CIband$lower)),col="grey90",border="grey")
  points(1:708,call_2_step_CIband$mean,type="l")
  abline(h=0,lty=3,lwd=2)
  legend("bottomright",legend=c("estimate", "95% CI"),
         pch = c(NA,15), lty=c(1,NA), col=c("black","grey"),
         bty = "n", # remove the bounder of the legend
         lwd = c(1,NA), pt.bg = c(NA,"grey90"),cex=2)
  
  # call's 3-step general effect
  par(mar = c(2.5, 5, .5, .5))
  call_3_step_general_101_CIband = plot_simulatedCI(t(call_3_step_general_101),probs=c(0.05,0.95),printFlag=F)
  plot(1:708,call_3_step_general_101_CIband$mean,type="l",ylab="3-step general effect (calls)",xlab="# lags",bty="n",cex.axis=1.5,cex.lab=1.5,ylim=c(-2,1))
  polygon(c(1:708,rev(1:708)),c(call_3_step_general_101_CIband$upper,rev(call_3_step_general_101_CIband$lower)),col="grey90",border="grey")
  points(1:708,call_3_step_general_101_CIband$mean,type="l")
  abline(h=0,lty=3,lwd=2)
  legend("bottomright",legend=c("estimate", "95% CI"),
         pch = c(NA,15), lty=c(1,NA), col=c("black","grey"),
         bty = "n", # remove the bounder of the legend
         lwd = c(1,NA), pt.bg = c(NA,"grey90"),cex=2)
  
  # impulse impact graph at t=600
  call_effects_vs_qlag_CIband=plot_simulatedCI(call_qlag_effects,probs=c(0.05,0.95),printFlag=F)
  plot(0:10,call_effects_vs_qlag_CIband$mean,type="l",ylab="q-lag effect at t=600 (calls)",xlab="# lags",bty="n",cex.axis=2,cex.lab=2,ylim=c(-1,0.5))
  polygon(c(0:10,rev(0:10)),c(call_effects_vs_qlag_CIband$upper,rev(call_effects_vs_qlag_CIband$lower)),col="grey90",border="grey")
  points(0:10,call_effects_vs_qlag_CIband$mean,type="l")
  points(0:10,call_effects_vs_qlag_CIband$mean,pch=16)
  abline(h=0,lty=3,lwd=2)
  legend("bottomright",legend=c("estimate", "95% CI"),
         pch = c(NA,15), lty=c(1,NA), col=c("black","grey"),
         bty = "n", # remove the bounder of the legend
         lwd = c(1,NA), pt.bg = c(NA,"grey90"),cex=2)
}


#########################################################################
######          simulate causal estimands of texts with CI          #####
#########################################################################
# jump directly to load the data without recalculation!!
# covariates has been extracted before in the calculation of calls
# y_coeffi_var_table=lapply(1:nrow(ssm_y$out_smooth$s),function(i){dlmSvd2var(ssm_y$out_smooth$U.S[[i]],ssm_y$out_smooth$D.S[i,])});head(y_coeffi_var_table)
# c_coeffi_var_table=lapply(1:nrow(ssm_c$out_smooth$s),function(i){dlmSvd2var(ssm_c$out_smooth$U.S[[i]],ssm_c$out_smooth$D.S[i,])});head(c_coeffi_var_table)
# second, calculate for texts
y_coeffi_table=as.data.frame(ssm_y$out_smooth$s);colnames(y_coeffi_table)[c(1:2,5:7)]=c("(Intercept)","y_1","x","x_1","c");head(y_coeffi_table)
c_coeffi_table=as.data.frame(ssm_c$out_smooth$s);colnames(c_coeffi_table)[c(1,2,4,5)]=c("(Intercept)","c_1","x_1","y_1");head(c_coeffi_table)
raw_data=data[,c(2,3,5,7)];colnames(raw_data)=c("Date","y","x","c");head(raw_data)

# calculate text's contemporaneous effect directly with CI (save call_contemporaneous)
set.seed(1)
text_contemporaneous = calculate.effect_allt_withCI(tx=c(1),y_coeffi_table=y_coeffi_table,y_coeffi_var_table=y_coeffi_var_table,c_coeffi_table=c_coeffi_table,c_coeffi_var_table=c_coeffi_var_table,n_sim=1000,printFlag=T)
# calculate text's 1-lag effect directly with CI
set.seed(2)
text_1_lag = calculate.effect_allt_withCI(tx=c(1,0),y_coeffi_table=y_coeffi_table,y_coeffi_var_table=y_coeffi_var_table,c_coeffi_table=c_coeffi_table,c_coeffi_var_table=c_coeffi_var_table,n_sim=1000,printFlag=T)
# calculate text's 1-lag structural direct effect directly with CI
set.seed(3)
text_1_lag_structural_direct = calculate.controlled_direct_effect_allt_withCI(y_coeffi_table,y_coeffi_var_table,n_sim=1000,printFlag=T)
# calculate text's 2-step total effect directly with CI
set.seed(4)
text_1_step = calculate.effect_allt_withCI(tx=c(1,1),y_coeffi_table=y_coeffi_table,y_coeffi_var_table=y_coeffi_var_table,c_coeffi_table=c_coeffi_table,c_coeffi_var_table=c_coeffi_var_table,n_sim=1000,printFlag=T)
# calculate text's 3-step general effect directly with CI
set.seed(5)
text_2_step_general_101 = calculate.effect_allt_withCI(tx=c(1,0,1),y_coeffi_table=y_coeffi_table,y_coeffi_var_table=y_coeffi_var_table,c_coeffi_table=c_coeffi_table,c_coeffi_var_table=c_coeffi_var_table,n_sim=1000,printFlag=T)
# calculate text's 4-step general effect directly with CI
set.seed(6)
text_3_step_general_0101 = calculate.effect_allt_withCI(tx=c(0,1,0,1),y_coeffi_table=y_coeffi_table,y_coeffi_var_table=y_coeffi_var_table,c_coeffi_table=c_coeffi_table,c_coeffi_var_table=c_coeffi_var_table,n_sim=1000,printFlag=T)

# impulse impact graph at t=600
t=600;n_sim=1000
set.seed(1)
text_qlag_effects_part1=simulate.counterfactual_singlet_withCI(t,tx=c(1,rep(0,10)),y_coeffi_table,y_coeffi_var_table,c_coeffi_table,c_coeffi_var_table,raw_data,n_sim=n_sim,printFlag=T)
set.seed(1)
text_qlag_effects_part2=simulate.counterfactual_singlet_withCI(t,tx=c(0,rep(0,10)),y_coeffi_table,y_coeffi_var_table,c_coeffi_table,c_coeffi_var_table,raw_data,n_sim=n_sim,printFlag=T)
text_effects_vs_qlag=text_qlag_effects_part1-text_qlag_effects_part2
set.seed(1)
text_qstep_effects_part1=simulate.counterfactual_singlet_withCI(t,tx=c(1,rep(1,10)),y_coeffi_table,y_coeffi_var_table,c_coeffi_table,c_coeffi_var_table,raw_data,n_sim=n_sim,printFlag=T)
set.seed(1)
text_qstep_effects_part2=simulate.counterfactual_singlet_withCI(t,tx=c(0,rep(0,10)),y_coeffi_table,y_coeffi_var_table,c_coeffi_table,c_coeffi_var_table,raw_data,n_sim=n_sim,printFlag=T)
text_effects_vs_qstep=text_qstep_effects_part1-text_qstep_effects_part2

save(text_contemporaneous,
     text_1_lag,text_1_lag_structural_direct,
     text_1_step,text_2_step_general_101,
     text_3_step_general_0101,
     text_effects_vs_qlag,
     text_effects_vs_qstep,
     file="/Users/xiaoxuancai/Dropbox (Personal)/MHealthPsychSummerProject2020/Xiaoxuan_Cai/[Paper 1] Causal estimands and Graphical representation for time series data/result_texts.Rdata")
load("/Users/xiaoxuancai/Dropbox (Personal)/MHealthPsychSummerProject2020/Xiaoxuan_Cai/[Paper 1] Causal estimands and Graphical representation for time series data/result_texts.Rdata")


# plots
{
  # text's contemporaneous effect
  # Chosen in maintext of the paper (Figure 4a)
  location="/Users/xiaoxuancai/Dropbox (Personal)/MHealthPsychSummerProject2020/Xiaoxuan_Cai/[Paper 1] Causal estimands and Graphical representation for time series data/"
  pdf(file = paste(location,"contemp_texts.pdf",sep=""),
      width = 10, # The width of the plot in inches
      height = 6) # The height of the plot in inches
  par(mar = c(2, 5, .1, .1))
  par(mfrow=c(1,1))
  text_contemporaneous_CIband=plot_simulatedCI(t(text_contemporaneous),probs=c(0.05,0.95),printFlag=F)
  plot(1:708,text_contemporaneous_CIband$mean,type="l",ylab="contemporaneous (texts)",xlab="",bty="n", cex.axis=2.5,cex.lab=2.5,ylim=c(-3,0.5))
  polygon(c(1:708,rev(1:708)),c(text_contemporaneous_CIband$upper,rev(text_contemporaneous_CIband$lower)),col="grey90",border="grey")
  points(1:708,text_contemporaneous_CIband$mean,type="l")
  abline(h=0,lty=3,lwd=2)
  text(10,0.5, "a)",cex=2.5)
  legend("bottomright",legend=c("estimate", "95% CI"),
         pch = c(NA,15), lty=c(1,NA), col=c("black","grey"),
         bty = "n", # remove the bounder of the legend
         lwd = c(1,NA), pt.bg = c(NA,"grey90"),cex=2.5)
  dev.off()
  
  # text's 1-lag effect  (Figure 4b)
  pdf(file = paste(location,"lag1_texts.pdf",sep=""),
      width = 10, # The width of the plot in inches
      height = 6) # The height of the plot in inches
  par(mar = c(2, 5, .1, .1))
  par(mfrow=c(1,1))
  text_1_lag_CIband=plot_simulatedCI(t(text_1_lag),probs=c(0.05,0.95),printFlag=F)
  plot(1:708,text_1_lag_CIband$mean,type="l",ylab="1-lag effect (texts)",xlab="Date",bty="n", cex.axis=2.5,cex.lab=2.5,ylim=c(-3,0.5))
  polygon(c(1:708,rev(1:708)),c(text_1_lag_CIband$upper,rev(text_1_lag_CIband$lower)),col="grey90",border="grey")
  points(1:708,text_1_lag_CIband$mean,type="l")
  abline(h=0,lty=3,lwd=2)
  text(10,0.5, "b)",cex=2.5)
  legend("bottomright",legend=c("estimate", "95% CI"),
         pch = c(NA,15), lty=c(1,NA), col=c("black","grey"),
         bty = "n", # remove the bounder of the legend
         lwd = c(1,NA), pt.bg = c(NA,"grey90"),cex=2.5)
  dev.off()
  
  # unused code for plotting text's q-lag effects for q=2,3,4 for all time points
  # lag2_analytical=list.cbind(lapply(result_alltimepoints_analytical,function(x){x[,3]}))
  # lag2_CIband_analytical=plot_simulatedCI(t(lag2_analytical),probs=c(0.05,0.95),printFlag=F)
  # plot(data$Date,lag2_CIband_analytical$mean,type="l",ylab="2-lag effect (texts)",xlab="Date",bty="n",xaxt="n",cex.axis=2,cex.lab=2,ylim=c(-3,0.3))
  # axis(side=1,at=xlabs,labels=as.character(xlabs),cex.axis=1.5)
  # polygon(c(data$Date,rev(data$Date)),c(lag2_CIband_analytical$upper,rev(lag2_CIband_analytical$lower)),col="grey90",border="grey")
  # points(data$Date,lag2_CIband_analytical$mean,type="l")
  # abline(h=0,lty=3,lwd=2)
  # legend("bottomright",legend=c("estimate", "95% CI"),
  #        pch = c(NA,15), lty=c(1,NA), col=c("black","grey"),
  #        bty = "n", # remove the bounder of the legend
  #        lwd = c(1,NA), pt.bg = c(NA,"grey90"),cex=2)
  # 
  # lag3_analytical=list.cbind(lapply(result_alltimepoints_analytical,function(x){x[,4]}))
  # lag3_CIband_analytical=plot_simulatedCI(t(lag3_analytical),probs=c(0.05,0.95),printFlag=F)
  # plot(data$Date,lag3_CIband_analytical$mean,type="l",ylab="3-lag effect (texts)",xlab="Date",bty="n",xaxt="n",cex.axis=2,cex.lab=2,ylim=c(-3,0.3))
  # axis(side=1,at=xlabs,labels=as.character(xlabs),cex.axis=1.5)
  # polygon(c(data$Date,rev(data$Date)),c(lag3_CIband_analytical$upper,rev(lag3_CIband_analytical$lower)),col="grey90",border="grey")
  # points(data$Date,lag3_CIband_analytical$mean,type="l")
  # abline(h=0,lty=3,lwd=2)
  # legend("bottomright",legend=c("estimate", "95% CI"),
  #        pch = c(NA,15), lty=c(1,NA), col=c("black","grey"),
  #        bty = "n", # remove the bounder of the legend
  #        lwd = c(1,NA), pt.bg = c(NA,"grey90"),cex=2)
  # 
  # lag4_analytical=list.cbind(lapply(result_alltimepoints_analytical,function(x){x[,5]}))
  # lag4_CIband_analytical=plot_simulatedCI(t(lag4_analytical),probs=c(0.05,0.95),printFlag=F)
  # plot(data$Date,lag4_CIband_analytical$mean,type="l",ylab="4-lag effect (texts)",xlab="Date",bty="n",xaxt="n",cex.axis=2,cex.lab=2,ylim=c(-3,0.3))
  # axis(side=1,at=xlabs,labels=as.character(xlabs),cex.axis=1.5)
  # polygon(c(data$Date,rev(data$Date)),c(lag4_CIband_analytical$upper,rev(lag4_CIband_analytical$lower)),col="grey90",border="grey")
  # points(data$Date,lag4_CIband_analytical$mean,type="l")
  # abline(h=0,lty=3,lwd=2)
  # legend("bottomright",legend=c("estimate", "95% CI"),
  #        pch = c(NA,15), lty=c(1,NA), col=c("black","grey"),
  #        bty = "n", # remove the bounder of the legend
  #        lwd = c(1,NA), pt.bg = c(NA,"grey90"),cex=2)
  
  
  # text's 1-lag structural direct effect  (Figure 4c)
  pdf(file = paste(location,"lag1_controlled_direct_texts.pdf",sep=""),
      width = 10, # The width of the plot in inches
      height = 6) # The height of the plot in inches
  par(mar = c(2, 5, .1, .1))
  par(mfrow=c(1,1))
  text_1_lag_structural_direct_CIband = plot_simulatedCI(t(text_1_lag_structural_direct),probs=c(0.05,0.95),printFlag=F)
  plot(1:708,text_1_lag_structural_direct_CIband$mean,type="l",ylab="1-lag controlled direct (texts)",xlab="# lags",bty="n",cex.axis=2.5,cex.lab=2.3,ylim=c(-3,0.5))
  polygon(c(1:708,rev(1:708)),c(text_1_lag_structural_direct_CIband$upper,rev(text_1_lag_structural_direct_CIband$lower)),col="grey90",border="grey")
  points(1:708,text_1_lag_structural_direct_CIband$mean,type="l")
  abline(h=0,lty=3,lwd=2)
  text(10,0.5, "c)",cex=2.5)
  legend("bottomright",legend=c("estimate", "95% CI"),
         pch = c(NA,15), lty=c(1,NA), col=c("black","grey"),
         bty = "n", # remove the bounder of the legend
         lwd = c(1,NA), pt.bg = c(NA,"grey90"),cex=2.5)
  dev.off()
  
  # text's 1-step total effect  (Figure 4d)
  pdf(file = paste(location,"step1_total_texts.pdf",sep=""),
      width = 10, # The width of the plot in inches
      height = 6) # The height of the plot in inches
  par(mar = c(2, 5, .1, .1))
  par(mfrow=c(1,1))
  text_1_step_CIband = plot_simulatedCI(t(text_1_step),probs=c(0.05,0.95),printFlag=F)
  plot(1:708,text_1_step_CIband$mean,type="l",ylab="1-step total (texts)",xlab="# lags",bty="n",cex.axis=2.5,cex.lab=2.5,ylim=c(-6,0.5))
  polygon(c(1:708,rev(1:708)),c(text_1_step_CIband$upper,rev(text_1_step_CIband$lower)),col="grey90",border="grey")
  points(1:708,text_1_step_CIband$mean,type="l")
  abline(h=0,lty=3,lwd=2)
  text(10,0.5, "d)",cex=2.5)
  legend("bottomright",legend=c("estimate", "95% CI"),
         pch = c(NA,15), lty=c(1,NA), col=c("black","grey"),
         bty = "n", # remove the bounder of the legend
         lwd = c(1,NA), pt.bg = c(NA,"grey90"),cex=2.5)
  dev.off()
  
  # text's 3-step general effect (Figure 6a)
  pdf(file = paste(location,"step3_general_0101_texts.pdf",sep=""),
      width = 10, # The width of the plot in inches
      height = 6) # The height of the plot in inches
  par(mar = c(5, 5, .1, .1))
  par(mfrow=c(1,1))
  text_3_step_general_0101_CIband = plot_simulatedCI(t(text_3_step_general_0101),probs=c(0.05,0.95),printFlag=F)
  plot(1:708,text_3_step_general_0101_CIband$mean,type="l",ylab="4-step general effect (texts)",xlab="days",bty="n",cex.axis=2.5,cex.lab=2.5,ylim=c(-4,0.5))
  polygon(c(1:708,rev(1:708)),c(text_3_step_general_0101_CIband$upper,rev(text_3_step_general_0101_CIband$lower)),col="grey90",border="grey")
  points(1:708,text_3_step_general_0101_CIband$mean,type="l")
  abline(h=0,lty=3,lwd=2)
  text(10,0.5, "a)",cex=2.5)
  legend("bottomright",legend=c("estimate", "95% CI"),
         pch = c(NA,15), lty=c(1,NA), col=c("black","grey"),
         bty = "n", # remove the bounder of the legend
         lwd = c(1,NA), pt.bg = c(NA,"grey90"),cex=2.5)
  dev.off()
  
  # impulse impact graph at t=600 (Figure 5a)
  pdf(file = paste(location,"impulse_impact_texts.pdf",sep=""),
      width = 10, # The width of the plot in inches
      height = 6) # The height of the plot in inches
  par(mar = c(5, 5, .1, .1))
  par(mfrow=c(1,1))
  text_effects_vs_qlag_CIband=plot_simulatedCI(text_effects_vs_qlag,probs=c(0.05,0.95),printFlag=F)
  plot(0:10,text_effects_vs_qlag_CIband$mean,type="l",ylab="q-lag at t=600 (texts)",xlab="# lags",bty="n",cex.axis=2.5,cex.lab=2.5,ylim=c(-2.5,0.6))
  polygon(c(0:10,rev(0:10)),c(text_effects_vs_qlag_CIband$upper,rev(text_effects_vs_qlag_CIband$lower)),col="grey90",border="grey")
  points(0:10,text_effects_vs_qlag_CIband$mean,type="l")
  points(0:10,text_effects_vs_qlag_CIband$mean,pch=16)
  abline(h=0,lty=3,lwd=2)
  text(0.2,0.6, "a)",cex=2.5)
  legend("bottomright",legend=c("estimate", "95% CI"),
         pch = c(NA,15), lty=c(1,NA), col=c("black","grey"),
         bty = "n", # remove the bounder of the legend
         lwd = c(1,NA), pt.bg = c(NA,"grey90"),cex=2)
  dev.off()
  
  # impulse response graph at t=600 (Figure 5b)
  pdf(file = paste(location,"impulse_response_texts.pdf",sep=""),
      width = 10, # The width of the plot in inches
      height = 6) # The height of the plot in inches
  par(mar = c(5, 5, .1, .1))
  par(mfrow=c(1,1))
  text_effects_vs_qstep_CIband=plot_simulatedCI(text_effects_vs_qstep,probs=c(0.05,0.95),printFlag=F)
  plot(0:10,text_effects_vs_qstep_CIband$mean,type="l",ylab="q-step total at t=600 (texts)",xlab="# lags",bty="n",cex.axis=2.5,cex.lab=2.5,ylim=c(-9.5,0.6))
  polygon(c(0:10,rev(0:10)),c(text_effects_vs_qstep_CIband$upper,rev(text_effects_vs_qstep_CIband$lower)),col="grey90",border="grey")
  points(0:10,text_effects_vs_qstep_CIband$mean,type="l")
  points(0:10,text_effects_vs_qstep_CIband$mean,pch=16)
  abline(h=0,lty=3,lwd=2)
  text(0.2,0.6, "b)",cex=2.5)
  legend("bottomright",legend=c("estimate", "95% CI"),
         pch = c(NA,15), lty=c(1,NA), col=c("black","grey"),
         bty = "n", # remove the bounder of the legend
         lwd = c(1,NA), pt.bg = c(NA,"grey90"),cex=2.5)
  dev.off()
}

########################################################################
######         personalized intervention effect (text)             #####
########################################################################
t=600;n_sim=1000
set.seed(1)
qstep1001001=simulate.counterfactual_singlet(t,tx=c(1,0,0,1,0,0,1,rep(0,10)),y_coeffi_table,c_coeffi_table,raw_data,printFlag=T)
set.seed(1)
qstep0101010=simulate.counterfactual_singlet(t,tx=c(0,1,0,1,0,1,0,rep(0,10)),y_coeffi_table,c_coeffi_table,raw_data,printFlag=T)
set.seed(1)
qstep1110000=simulate.counterfactual_singlet(t,tx=c(1,1,1,0,0,0,0,rep(0,10)),y_coeffi_table,c_coeffi_table,raw_data,printFlag=T)
set.seed(1)
qstep0000000=simulate.counterfactual_singlet(t,tx=c(0,0,0,0,0,0,0,rep(0,10)),y_coeffi_table,c_coeffi_table,raw_data,printFlag=T)

qstep1001001_CIband_simulated_600=plot_simulatedCI(qstep1001001-qstep0000000,probs=c(0.05,0.95),printFlag=F)
qstep0101010_CIband_simulated_600=plot_simulatedCI(qstep0101010-qstep0000000,probs=c(0.05,0.95),printFlag=F)
qstep1110000_CIband_simulated_600=plot_simulatedCI(qstep1110000-qstep0000000,probs=c(0.05,0.95),printFlag=F)

# Figure 6b
pdf(file = paste(location,"optimal_tx.pdf",sep=""),
    width = 10, # The width of the plot in inches
    height = 6) # The height of the plot in inches
par(mar = c(5, 5, 1.5, .1))
par(mfrow=c(1,1))
plot(0:16,qstep1110000-qstep0000000,type="l",ylab="7-step general effect(texts)",xlab="days",bty="n",cex.axis=2.5,cex.lab=2.5,ylim=c(-4,0.7))
polygon(x = c(-0.1,6,6,-0.1), 
        y = c(-4,-4,0.7,0.7),
        col = "lightgrey",border ="NA")
text(2.9,0.3,"1-week of intervention",cex=1.8)
points(0:16,qstep1110000-qstep0000000,type="l")
points(0:16,qstep1110000-qstep0000000,pch=15,cex=2)
points(0:16,qstep0101010-qstep0000000,type="l",pch=17,lty=2, col="brown",cex=2)
points(0:16,qstep0101010-qstep0000000,pch=16,lty=3, col="brown",cex=2)
points(0:16,qstep1001001-qstep0000000,type="l",pch=15,lty=4,col="blue",cex=2)
points(0:16,qstep1001001-qstep0000000,pch=17,lty=2,col="blue",cex=2)
abline(v=6,lty=2,lwd=2,col="red")
text(0,.7,"b)",cex=2.5)
legend("bottomright",legend=c("tx=(1,1,1,0,0,0,0)","tx=(0,1,0,1,0,1,0)","tx=(1,0,0,1,0,0,1)"),
       pch = c(15,16,17), lty=c(1,2,4), col=c("black","brown","blue"),
       bty = "n", # remove the bounder of the legend
       lwd = c(1,1,1), pt.bg = c(NA,"grey90"),cex=2.5)
dev.off()

##############################################################
######         test positivity assumption                #####
##############################################################
positivity_calls=test.positivity(tx="keycontacts_call_totaldegree_merged",data=data[,c(1:7)],length=10)
positivity_texts=test.positivity(tx="keycontacts_text_reciprocity_degree_merged",data=data[,c(1:7)],length=10)

nrow(positivity_texts$details[[7]][rowSums(positivity_texts$details[[7]][,-8])==3,])

par(mar = c(5, 5, 2.5, 0.5))
barplot(height=positivity_calls$percentage,names=as.character(1:10),ylab="positivity % (calls)",xlab="exposure duration (days)",type="h",bty="n",cex.axis=1.5,cex.lab=2,las=1)
barplot(height=positivity_texts$percentage,names=as.character(1:10),ylab="positivity % (texts)",xlab="exposure duration (days)",type="h",bty="n",cex.axis=1.5,cex.lab=2,las=1)





