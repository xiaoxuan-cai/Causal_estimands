# File: 20230810_5BT65_causalestimands.R
# Date: 2023.08.10
# Goal: fit both covariate and outcome SSM regression model -> calculate various causal estimands
# Contents: consider two candidate sets of exposures: 
#                     <call_outdegree> (binary) + <text_outdegree> (binary)
#                     <keycontacts_call_outdegree> (binary) + <keycontacts_text_reciprocity_degree> (binary)
#           for each of them, fit proper SSM without any imputation (ignore case)
# Explanation: we only consider binary exposure in the paper, therefore, we transform exposure into binary
#              for the pair of <keycontacts_call_outdegree> + <keycontacts_text_reciprocity_degree>, which both range 0~3
#                 we merge level 3 to level 2 since there are two few of them
source("/Users/xiaoxuancai/Documents/GitHub/mHealth_data_processing/helper_extract_combine_data.R")
source("/Users/xiaoxuancai/Documents/GitHub/mHealth_data_processing/Summary 1_packages.R")
source("/Users/xiaoxuancai/Documents/GitHub/mHealth_data_processing/Summary 2_helper functions.R")
source("/Users/xiaoxuancai/Documents/GitHub/SSMimpute/helper_multipleimputation.R")
source("/Users/xiaoxuancai/Dropbox (Personal)/MHealthPsychSummerProject2020/Xiaoxuan_Cai/[Paper 1] Causal estimands and Graphical representation for time series data/Rcode/helper_causal estimands.R")
load("/Users/xiaoxuancai/Documents/analysis/subject_5BT65/Data/processed/combined_all_5BT65.RData")

# -------------------------------------------------------------------- #
#             Organize and add variables needed for analysis           #
# -------------------------------------------------------------------- #
# change values for certain variables
TAM_phone_modified=data_5BT65$TAM_phone
TAM_phone_modified[which(TAM_phone_modified==1)]=0.9999
TAM_phone_modified[which(TAM_phone_modified==0)]=0.0001

table(data_5BT65$keycontacts_call_totaldegree) # merge level 2 to level 1 -> turn into a binary variable
table(data_5BT65$keycontacts_text_totaldegree) # merge level 2&3 to level 1 -> turn into a binary variable
table(data_5BT65$keycontacts_text_reciprocity_degree) # merge level 2&3 to level 1 -> turn into a binary variable
keycontacts_call_totaldegree_binary=data_5BT65$keycontacts_call_totaldegree;keycontacts_call_totaldegree_binary[keycontacts_call_totaldegree_binary==2]=1;table(keycontacts_call_totaldegree_binary)
keycontacts_text_totaldegree_binary=data_5BT65$keycontacts_text_totaldegree;keycontacts_text_totaldegree_binary[keycontacts_text_totaldegree_binary>1]=1;table(keycontacts_text_totaldegree_binary)
keycontacts_text_reciprocity_degree_binary=data_5BT65$keycontacts_text_reciprocity_degree;keycontacts_text_reciprocity_degree_binary[keycontacts_text_reciprocity_degree_binary>1]=1;table(keycontacts_text_reciprocity_degree_binary)

data=data.frame(Date=data_5BT65$Date,negative_total=data_5BT65$negative_total,
                keycontacts_call_totaldegree_binary,
                keycontacts_text_totaldegree_binary,
                logit_TAM_phone=log(TAM_phone_modified/(1-TAM_phone_modified)))
param=list(list(lagged_param=list(variables=colnames(data)[-1],param=rep(3,length(colnames(data)[-1])))))
data=add_variables_procedures(data,param);
data=data.frame(intercept=1,data);colnames(data)
# plots of outcome/exposure/covariates
{
  par(mar = c(5, 5, 0.5, 0.5))
  par(mar = c(5, 5, 0.5, 0.5))
  plot(1:708,data$negative_total,type="l",bty="n",xlab="Time (days)",ylab="Negative mood",cex.lab=1.8,cex.axis=1.4)
  # abline(v=c(1:708)[is.na(data$negative_total)],col="grey")
  plot(1:708,data$keycontacts_call_totaldegree_binary,bty="n",type="h",xlab="Time (days)",ylab="Degree of outgoing calls",cex.lab=1.8,cex.axis=1.4)
  plot(1:708,data$keycontacts_text_totaldegree_binary,bty="n",type="h",xlab="Time (days)",ylab="Degree of outgoing texts",cex.lab=1.8,cex.axis=1.4)
  plot(1:708,data$logit_TAM_phone,bty="n",type="l",xlab="Time (days)",ylab="log(Phone Mobility)",cex.lab=2,cex.axis=1.4)
}


##################################################################
######          SSM ignore modeling of covariates            #####
##################################################################
# This paper cannot handle missing data for both covariates and outcome
#   thus, we use SSM ignore temporarily
# covariate regression: logit_TAM_phone ~ intercept + logit_TAM_phone_1  
#                                         + keycontacts_call_totaldegree_binary
#                                         + keycontacts_text_totaldegree_binary + negative_total
{
  ######################## Part I:  initial guess of ssm structure ########################
  n_variables=5
  ssm_c_init=SSModel(logit_TAM_phone ~ -1 + SSMregression(~ intercept + logit_TAM_phone_1+  
                                                            keycontacts_call_totaldegree_binary+
                                                            keycontacts_text_totaldegree_binary+
                                                            negative_total,  
                                                            data=data,Q=diag(rep(NA,n_variables))),
                       data,H=NA) # baseline variance constant
  ssm_c_init_fit=fitSSM(ssm_c_init, inits =rep(0,n_variables+1), method = "BFGS") # n_variables + one NA in H
  ssm_c_init_out=KFS(ssm_c_init_fit$model)
  plot.KFS(ssm_c_init_out,range=30:708)
  ssm_c_init_out$model$H;ssm_c_init_out$model$Q 
  # check intercept(.01240), logit_TAM_phone_1(0.00043), keycontacts_text_outdegree_merged (0.00242)
  
  ######################## Part II: statespace model without imputation ########################
  formula="logit_TAM_phone ~ intercept + logit_TAM_phone_1 + keycontacts_call_totaldegree_binary + keycontacts_text_totaldegree_binary + negative_total"
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
                                            list(variable="keycontacts_text_totaldegree_binary",segments=2,changepoints=c(415),fixed_cpts=F)))
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
  # (Intercept)                                     1.614     0.638   (0.361,2.867)   (0.565,2.663)
  # y_1(period1)                                    0.052     0.097   (-0.14,0.243)  (-0.109,0.212)
  # y_1(period2)                                    0.175     0.091  (-0.004,0.353)   (0.025,0.324)
  # y_1(period3)                                   -0.112     0.046 (-0.203,-0.021) (-0.188,-0.036)
  # keycontacts_call_totaldegree_binary             0.003     0.290  (-0.566,0.572)  (-0.474,0.479)
  # keycontacts_text_totaldegree_binary(period1)    0.096     0.257   (-0.41,0.601)  (-0.327,0.518)
  # keycontacts_text_totaldegree_binary(period2)   -1.054     0.319 (-1.681,-0.428)  (-1.579,-0.53)
  # negative_total                                 -0.010     0.027  (-0.063,0.043)  (-0.054,0.034)
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
#                                     + keycontacts_text_totaldegree_binary + keycontacts_text_totaldegree_binary_1 + logit_TAM_phone_1
{
  ######################## Part I: initial guess of ssm structure ########################
  n_variables=7
  ssm_y_init=SSModel(negative_total ~ -1 + SSMregression(~ intercept + negative_total_1+  
                                                           keycontacts_call_totaldegree_binary+keycontacts_call_totaldegree_binary_1+
                                                           keycontacts_text_totaldegree_binary+keycontacts_text_totaldegree_binary_1+
                                                           logit_TAM_phone_1,  
                                                           data=data,Q=diag(rep(NA,n_variables))),
                        data,H=NA) # baseline variance constant
  ssm_y_init_fit=fitSSM(ssm_y_init, inits =rep(0,n_variables+1), method = "BFGS") # n_variables + one NA in H
  ssm_y_init_out=KFS(ssm_y_init_fit$model)
  plot.KFS(ssm_y_init_out,range=30:708)
  ssm_y_init_out$model$Q; ssm_y_init_out$model$H
  # check: intercept(0.1474), text_totaldegree(0.0007), text_totaldegree_1 (0.00007)
  
  ######################## Part II: statespace model without imputation ########################
  formula=" negative_total ~ intercept + negative_total_1 + keycontacts_call_totaldegree_binary + keycontacts_call_totaldegree_binary_1 + keycontacts_text_totaldegree_binary + keycontacts_text_totaldegree_binary_1+ logit_TAM_phone_1"
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
                     w_cp_param=list(list(variable="keycontacts_text_totaldegree_binary",segments=3,changepoints=c(560,640),fixed_cpts=F),
                                     list(variable="keycontacts_text_totaldegree_binary_1",segments=2,changepoints=c(460),fixed_cpts=F)))
  ssm_y=run.SSM_partial_converged_cpts(data_ss=data_ss_ignore_y,formula_var=formula_var,ss_param=ss_param_temp,max_iteration=100,
                                       cpt_learning_param=list(cpt_method="mean",burnin=1/10,mergeband=20,convergence_cri=10),
                                       dlm_cpt_learning_option="filter",
                                       dlm_option="smooth",printFlag=T,bandwidth = 5, cpt_V = 1)
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

#########################################################################
######          simulate causal estimands of calls with CI          #####
#########################################################################
y_coeffi_table=as.data.frame(ssm_y$out_smooth$s);colnames(y_coeffi_table)[c(1:4,7)]=c("(Intercept)","y_1","x","x_1","c");head(y_coeffi_table)
c_coeffi_table=as.data.frame(ssm_c$out_smooth$s);colnames(c_coeffi_table)[c(1,2,3,5)]=c("(Intercept)","c_1","x_1","y_1");head(c_coeffi_table)
y_coeffi_var_table=lapply(1:nrow(ssm_y$out_smooth$s),function(i){dlmSvd2var(ssm_y$out_smooth$U.S[[i]],ssm_y$out_smooth$D.S[i,])});head(y_coeffi_var_table)
c_coeffi_var_table=lapply(1:nrow(ssm_c$out_smooth$s),function(i){dlmSvd2var(ssm_c$out_smooth$U.S[[i]],ssm_c$out_smooth$D.S[i,])});head(c_coeffi_var_table)
raw_data=data[,c(2,3,4,7)];colnames(raw_data)=c("Date","y","x","c");head(raw_data)

# simulate CI for contemporaneous effect and q-lag effects for q=1,2,3,4 for all time points
result_alltimepoints_analytical_calls=calculate.upto5_lag_effect_allt_withCI(y_coeffi_table,y_coeffi_var_table,c_coeffi_table,c_coeffi_var_table,n_sim=1000,printFlag=T)
save(result_alltimepoints_analytical_calls,file="/Users/xiaoxuancai/Dropbox (Personal)/MHealthPsychSummerProject2020/Xiaoxuan_Cai/[Paper 1] Causal estimands and Graphical representation for time series data/Rcode/result_alltimepoints_calls.Rdata")

# plots effects with CI for contemporaneous effect and q-lag effects for q=1,2,3,4 for all time points
{
  par(mar = c(2.5, 5, .5, .5))
  
  contemporenous_analytical_calls=list.cbind(lapply(result_alltimepoints_analytical_calls,function(x){x[,1]}))
  CIband_analytical_calls=plot_simulatedCI(t(contemporenous_analytical_calls),probs=c(0.05,0.95),printFlag=F)
  plot(1:708,CIband_analytical_calls$mean,type="l",ylab="contemporaneous effect (calls)",xlab="",bty="n", cex.axis=2,cex.lab=2,ylim=c(-1,0.5))
  polygon(c(1:708,rev(1:708)),c(CIband_analytical_calls$upper,rev(CIband_analytical_calls$lower)),col="grey90",border="grey")
  points(1:708,CIband_analytical_calls$mean,type="l")
  abline(h=0,lty=3,lwd=2)
  legend("bottomright",legend=c("estimate", "95% CI"),
         pch = c(NA,15), lty=c(1,NA), col=c("black","grey"),
         bty = "n", # remove the bounder of the legend
         lwd = c(1,NA), pt.bg = c(NA,"grey90"),cex=2)
  
  lag1_analytical_calls=list.cbind(lapply(result_alltimepoints_analytical_calls,function(x){x[,2]}))
  lag1_CIband_analytical_calls=plot_simulatedCI(t(lag1_analytical_calls),probs=c(0.05,0.95),printFlag=F)
  plot(1:708,lag1_CIband_analytical_calls$mean,type="l",ylab="1-lag effect (calls)",xlab="Date",bty="n", cex.axis=2,cex.lab=2,ylim=c(-1,0.5))
  polygon(c(1:708,rev(1:708)),c(lag1_CIband_analytical_calls$upper,rev(lag1_CIband_analytical_calls$lower)),col="grey90",border="grey")
  points(1:708,lag1_CIband_analytical_calls$mean,type="l")
  abline(h=0,lty=3,lwd=2)
  legend("bottomright",legend=c("estimate", "95% CI"),
         pch = c(NA,15), lty=c(1,NA), col=c("black","grey"),
         bty = "n", # remove the bounder of the legend
         lwd = c(1,NA), pt.bg = c(NA,"grey90"),cex=2)
  
  lag2_analytical_calls=list.cbind(lapply(result_alltimepoints_analytical_calls,function(x){x[,3]}))
  lag2_CIband_analytical_calls=plot_simulatedCI(t(lag2_analytical_calls),probs=c(0.05,0.95),printFlag=F)
  plot(1:708,lag2_CIband_analytical_calls$mean,type="l",ylab="2-lag effect (calls)",xlab="Date",bty="n", cex.axis=2,cex.lab=2,ylim=c(-1,0.5))
  polygon(c(1:708,rev(1:708)),c(lag2_CIband_analytical_calls$upper,rev(lag2_CIband_analytical_calls$lower)),col="grey90",border="grey")
  points(1:708,lag2_CIband_analytical_calls$mean,type="l")
  abline(h=0,lty=3,lwd=2)
  legend("bottomright",legend=c("estimate", "95% CI"),
         pch = c(NA,15), lty=c(1,NA), col=c("black","grey"),
         bty = "n", # remove the bounder of the legend
         lwd = c(1,NA), pt.bg = c(NA,"grey90"),cex=2)
  
  lag3_analytical_calls=list.cbind(lapply(result_alltimepoints_analytical_calls,function(x){x[,4]}))
  lag3_CIband_analytical_calls=plot_simulatedCI(t(lag3_analytical_calls),probs=c(0.05,0.95),printFlag=F)
  plot(1:708,lag3_CIband_analytical_calls$mean,type="l",ylab="3-lag effect (calls)",xlab="Date",bty="n", cex.axis=2,cex.lab=2,ylim=c(-1,0.5))
  polygon(c(1:708,rev(1:708)),c(lag3_CIband_analytical_calls$upper,rev(lag3_CIband_analytical_calls$lower)),col="grey90",border="grey")
  points(1:708,lag3_CIband_analytical_calls$mean,type="l")
  abline(h=0,lty=3,lwd=2)
  legend("bottomright",legend=c("estimate", "95% CI"),
         pch = c(NA,15), lty=c(1,NA), col=c("black","grey"),
         bty = "n", # remove the bounder of the legend
         lwd = c(1,NA), pt.bg = c(NA,"grey90"),cex=2)
  
  lag4_analytical_calls=list.cbind(lapply(result_alltimepoints_analytical_calls,function(x){x[,5]}))
  lag4_CIband_analytical_calls=plot_simulatedCI(t(lag4_analytical_calls),probs=c(0.05,0.95),printFlag=F)
  plot(1:708,lag4_CIband_analytical_calls$mean,type="l",ylab="4-lag effect (texts)",xlab="Date",bty="n", cex.axis=2,cex.lab=2,ylim=c(-1,0.5))
  polygon(c(1:708,rev(1:708)),c(lag4_CIband_analytical_calls$upper,rev(lag4_CIband_analytical_calls$lower)),col="grey90",border="grey")
  points(1:708,lag4_CIband_analytical$mean,type="l")
  abline(h=0,lty=3,lwd=2)
  legend("bottomright",legend=c("estimate", "95% CI"),
         pch = c(NA,15), lty=c(1,NA), col=c("black","grey"),
         bty = "n", # remove the bounder of the legend
         lwd = c(1,NA), pt.bg = c(NA,"grey90"),cex=2)
  
}
# plots of effects vs q-lags
{
  t=600;n_sim=1000
  set.seed(1)
  qlag_effects_withCI_part1_ssm_600_calls=simulate.counterfactual_outcome_singlet_withCI(t,tx=c(1,rep(0,10)),y_coeffi_table,y_coeffi_var_table,c_coeffi_table,c_coeffi_var_table,raw_data,n_sim=n_sim,printFlag=T)
  set.seed(1)
  qlag_effects_withCI_part2_ssm_600_calls=simulate.counterfactual_outcome_singlet_withCI(t,tx=c(0,rep(0,10)),y_coeffi_table,y_coeffi_var_table,c_coeffi_table,c_coeffi_var_table,raw_data,n_sim=n_sim,printFlag=T)
  lageffects_CIband_simulated_600_calls=plot_simulatedCI(qlag_effects_withCI_part1_ssm_600_calls-qlag_effects_withCI_part2_ssm_600_calls,probs=c(0.05,0.95),printFlag=F)
  t=200
  set.seed(1)
  qlag_effects_withCI_part1_ssm_200_calls=simulate.counterfactual_outcome_singlet_withCI(t,tx=c(1,rep(0,10)),y_coeffi_table,y_coeffi_var_table,c_coeffi_table,c_coeffi_var_table,raw_data,n_sim=n_sim,printFlag=T)
  set.seed(1)
  qlag_effects_withCI_part2_ssm_200_calls=simulate.counterfactual_outcome_singlet_withCI(t,tx=c(0,rep(0,10)),y_coeffi_table,y_coeffi_var_table,c_coeffi_table,c_coeffi_var_table,raw_data,n_sim=n_sim,printFlag=T)
  lageffects_CIband_simulated_200_calls=plot_simulatedCI(qlag_effects_withCI_part1_ssm_200_calls-qlag_effects_withCI_part2_ssm_200_calls,probs=c(0.05,0.95),printFlag=F)
  
  plot(0:10,lageffects_CIband_simulated_200_calls$mean,type="l",ylab="q-lag effect at t=200 (calls)",xlab="# lags",bty="n",cex.axis=2,cex.lab=2,ylim=c(-1,0.5))
  polygon(c(0:10,rev(0:10)),c(lageffects_CIband_simulated_200_calls$upper,rev(lageffects_CIband_simulated_200_calls$lower)),col="grey90",border="grey")
  points(0:10,lageffects_CIband_simulated_200_calls$mean,type="l")
  points(0:10,lageffects_CIband_simulated_200_calls$mean,pch=16)
  abline(h=0,lty=3,lwd=2)
  legend("bottomright",legend=c("estimate", "95% CI"),
         pch = c(NA,15), lty=c(1,NA), col=c("black","grey"),
         bty = "n", # remove the bounder of the legend
         lwd = c(1,NA), pt.bg = c(NA,"grey90"),cex=2)
  
  plot(0:10,lageffects_CIband_simulated_600_calls$mean,type="l",ylab="q-lag effect at t=600 (calls)",xlab="# lags",bty="n",cex.axis=2,cex.lab=2,ylim=c(-1,0.5))
  polygon(c(0:10,rev(0:10)),c(lageffects_CIband_simulated_600_calls$upper,rev(lageffects_CIband_simulated_600_calls$lower)),col="grey90",border="grey")
  points(0:10,lageffects_CIband_simulated_600_calls$mean,type="l")
  points(0:10,lageffects_CIband_simulated_600_calls$mean,pch=16)
  legend("bottomright",legend=c("estimate", "95% CI"),
         pch = c(NA,15), lty=c(1,NA), col=c("black","grey"),
         bty = "n", # remove the bounder of the legend
         lwd = c(1,NA), pt.bg = c(NA,"grey90"),cex=2)
}
