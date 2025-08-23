# File: 20230810_5BT65_causalestimands.R
# Date: 2023.08.10
# Updated: 2025.08.23
# This file is created for the real data analysis of the paper "Causal estimands and identification of time-varying effects in non-stationary time series from N-of-1 mobile device data"
# Goal: fit time-varying SSM model for both outcome and covariate -> to calculate various causal estimands
# Note: 1. the methodology of the paper only applies to binary exposure and # obs on level 3 is too few for proper inference
#          therefore, we merge <keycontacts_call_totaldegree>(range 0~3) and <keycontacts_text_reciprocity_degree> (range 0~3)
#          to binary variables
#          we consider the following binary variable as exposure of interest:
#          <keycontacts_call_totaldegree> (binary) and <keycontacts_text_reciprocity_degree> (binary)
#       2. with no proper strategies for missing data, we use ignore case analysis for now with additional work to be developed to deal with missing data
#       3. we also update SSM fitting function to their new versions
#          old version "GitHub/SSMimpute/helper_multipleimputation.R" -> new version
#       4. keycontacts_text_totaldegree or keycontacts_text_reciprocity_degree?
#       5. take a logit transformation on TAM_phone, to make it normally distributed from -infty to +infty
# Results: 1. SSM modeling fitting results in main paper
#          2. SSM fitting plots in both main paper and appendix
#          3. Various causal estimands for both main paper and appendix
#          4. Impulse impact plot and Step response plot for main paper
#          5. test positivity for main paper
source("/Users/xiaoxuancai/Documents/GitHub/mHealth_data_processing/Summary 1_packages.R")
source("/Users/xiaoxuancai/Documents/GitHub/mHealth_data_processing/Summary 2_helper functions.R")
source("/Users/xiaoxuancai/Documents/GitHub/Causal_estimands/helper_causal estimands.R")
source("/Users/xiaoxuancai/Documents/GitHub/SSMimpute2/helper_SSM_version3.R") 
source("/Users/xiaoxuancai/Documents/GitHub/SSMimpute/helper_multipleimputation.R") # include merge point used for SSM and plot functions
load("/Users/xiaoxuancai/Documents/analysis/subject_5BT65/Data/processed/combined_all_5BT65.RData")
library(data.table) # for simulate.counterfactual_path_singlet to combine items of list
library(rlist)
address_main = "/Users/xiaoxuancai/Dropbox/MHealthPsychSummerProject2020/Xiaoxuan_Cai/[Paper 1] Causal estimands for time series data/Graphs_main/"
address_appendix = "/Users/xiaoxuancai/Dropbox/MHealthPsychSummerProject2020/Xiaoxuan_Cai/[Paper 1] Causal estimands for time series data/Graphs_appendix/"

# -------------------------------------------------------------------- #
#             Organize and add variables needed for analysis           #
# -------------------------------------------------------------------- #
# change values for certain variables
TAM_phone_modified=data_5BT65$TAM_phone
TAM_phone_modified[which(TAM_phone_modified==1)]=0.9999
TAM_phone_modified[which(TAM_phone_modified==0)]=0.0001
table(data_5BT65$keycontacts_call_totaldegree) # merge level 2 to level 1 -> turn into a binary variable
keycontacts_call_totaldegree_binary=data_5BT65$keycontacts_call_totaldegree;keycontacts_call_totaldegree_binary[keycontacts_call_totaldegree_binary==2]=1;table(keycontacts_call_totaldegree_binary)
table(data_5BT65$keycontacts_text_reciprocity_degree) # merge level 2&3 to level 1 -> turn into a binary variable
keycontacts_text_reciprocity_degree_binary=data_5BT65$keycontacts_text_reciprocity_degree;keycontacts_text_reciprocity_degree_binary[keycontacts_text_reciprocity_degree_binary>1]=1;table(keycontacts_text_reciprocity_degree_binary)
data=data.frame(Date=data_5BT65$Date,negative_total=data_5BT65$negative_total,
                keycontacts_call_totaldegree_binary,
                keycontacts_text_reciprocity_degree_binary,
                logit_TAM_phone=log(TAM_phone_modified/(1-TAM_phone_modified))) 
param=list(list(lagged_param=list(variables=colnames(data)[-1],param=rep(3,length(colnames(data)[-1])))))
data=add_variables_procedures(data,param);
data=data.frame(intercept=1,data);colnames(data)
# plots of outcome/exposure/covariates ("raw_950_600.png" used in Figure 3)
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

##############################################################################
######             ARIMA (used in Main paper Table 3)                  #######
##############################################################################
formula_y=" negative_total ~ negative_total_1 + keycontacts_call_totaldegree_binary + keycontacts_call_totaldegree_binary_1 + keycontacts_text_reciprocity_degree_binary + keycontacts_text_reciprocity_degree_binary_1+ logit_TAM_phone_1"
all_var_y=unlist(strsplit(gsub(" ","",formula_y),"\\+|\\~"))
data_ARIMAreg_y=data[,all_var_y]
result_ARIMAreg_y=auto.arima(y=data_ARIMAreg_y$negative_total,
                           xreg=as.matrix(data_ARIMAreg_y)[,all_var_y[-1]],max.d=0)
#auto.arima(y=data_ARIMAreg_y$negative_total,xreg=as.matrix(data_ARIMAreg_y)[,all_var_y[-c(1,2)]],max.d=0,ic="bic")
digits=4
result_ARIMAreg_estimate_y=cbind("Estimate"=round(result_ARIMAreg_y$coef,digits=digits),
                               "Std. Error"=round(sqrt(diag(result_ARIMAreg_y$var.coef)),digits=digits),
                               "99% CI"=paste("(",round(result_ARIMAreg_y$coef-qnorm(1-0.01/2)*sqrt(diag(result_ARIMAreg_y$var.coef)),digits=digits),",",
                                              round(result_ARIMAreg_y$coef+qnorm(1-0.01/2)*sqrt(diag(result_ARIMAreg_y$var.coef)),digits=digits),")",sep=""),
                               "95% CI"=paste("(",round(result_ARIMAreg_y$coef-qnorm(1-0.05/2)*sqrt(diag(result_ARIMAreg_y$var.coef)),digits=digits),",",
                                              round(result_ARIMAreg_y$coef+qnorm(1-0.05/2)*sqrt(diag(result_ARIMAreg_y$var.coef)),digits=digits),")",sep=""),
                               "90% CI"=paste("(",round(result_ARIMAreg_y$coef-qnorm(1-0.1/2)*sqrt(diag(result_ARIMAreg_y$var.coef)),digits=digits),",",
                                              round(result_ARIMAreg_y$coef+qnorm(1-0.1/2)*sqrt(diag(result_ARIMAreg_y$var.coef)),digits=digits),")",sep=""))
result_ARIMAreg_estimate_y
#                                              Estimate   Std. Error 99% CI           95% CI           90% CI        
# intercept                                    "0.751"**  "0.196"    "(0.246,1.256)"  "(0.367,1.135)"  "(0.428,1.073)"  
# negative_total_1                             "0.936"**  "0.017"    "(0.892,0.98)"   "(0.903,0.969)"  "(0.908,0.964)"  
# keycontacts_call_totaldegree_binary          "-0.444"+  "0.268"    "(-1.135,0.248)" "(-0.97,0.083)"  "(-0.885,-0.002)"
# keycontacts_call_totaldegree_binary_1        "-0.093"   "0.269"    "(-0.786,0.6)"   "(-0.62,0.434)"  "(-0.536,0.349)" 
# keycontacts_text_reciprocity_degree_binary   "-0.200"   "0.193"    "(-0.698,0.298)" "(-0.579,0.179)" "(-0.518,0.118)" 
# keycontacts_text_reciprocity_degree_binary_1 "-0.285"   "0.194"    "(-0.785,0.214)" "(-0.666,0.095)" "(-0.604,0.034)" 
# logit_TAM_phone_1                            "-0.018"   "0.033"    "(-0.103,0.067)" "(-0.083,0.046)" "(-0.073,0.036)"

###############################################################################################################
######          SSM ignore modeling of covariates  (used in Main paper Table 2&3 and Appendix)            #####
###############################################################################################################
# This paper cannot handle missing data for both covariates and outcome
#   thus, we use SSM ignore temporarily
# covariate regression: logit_TAM_phone ~ intercept + logit_TAM_phone_1  
#                                         + keycontacts_call_totaldegree_binary
#                                         + keycontacts_text_reciprocity_degree_binary
#                                         + negative_total
# ssm_c model (used in the Main paper Table 2)
{
  ######################## Part I:  initial guess of ssm structure ########################
  n_variables=5
  ssm_c_init=SSModel(logit_TAM_phone ~ -1 + SSMregression(~ intercept + logit_TAM_phone_1 +  
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
  formula="logit_TAM_phone ~ logit_TAM_phone_1 + keycontacts_call_totaldegree_binary + keycontacts_text_reciprocity_degree_binary + negative_total"
  all_var=unlist(strsplit(gsub(" ","",formula),"\\+|\\~"))
  outcome_var=strsplit(gsub(" ","",formula),"\\~")[[1]][1]
  formula_var=all_var[!(all_var %in% c(outcome_var,"intercept"))]
  data_ss_ignore=data[-1,all_var]
  data_ss_ignore[rowSums(is.na(data_ss_ignore))!=0,outcome_var]=NA
  for(l in formula_var){
    data_ss_ignore[is.na(data_ss_ignore[,l]),l]=mean(data_ss_ignore[,l],na.rm=T)
  }
  ss_param_keycontacts=list(inits=c(log(0.05),log(5)),m0=ssm_c_init_out$att[708,],C0=diag(rep(10^3),5),
                            AR1_coeffi=NULL,rw_coeffi=c("intercept"),v_cp_param=NULL,
                            w_cp_param=list(list(variable="logit_TAM_phone_1",segments=3,changepoints=c(410,450),fixed_cpts=F),
                                            list(variable="keycontacts_text_reciprocity_degree_binary",segments=2,changepoints=c(415),fixed_cpts=F)))
  ssm_c=run.SSM(data_ss=data_ss_ignore,formula=formula,ss_param_temp=ss_param_keycontacts,max_iteration=100,
                                   cpt_learning_param=list(cpt_method="mean",burnin=1/10,mergeband=10,convergence_cri=5),
                                   cpt_merge_option = "separate",dlm_cpt_learning_option="filter",
                                   bandwidth=5,cpt_V =1,printFlag=T)
  #                                                         Estimate  Std.Error
  # (Intercept)                                          1.29666707 0.60011083
  # logit_TAM_phone_1(period1)                           0.06622593 0.09750778
  # logit_TAM_phone_1(period2)                           0.17776074 0.09184185
  # logit_TAM_phone_1(period3)                          -0.10914657 0.04657227
  # keycontacts_call_totaldegree_binary                  0.01220309 0.29165595
  # keycontacts_text_reciprocity_degree_binary(period1) -0.05881142 0.28519091
  # keycontacts_text_reciprocity_degree_binary(period2) -0.78101412 0.30446712
  # negative_total                                      -0.01410104 0.02640353
  plot.dlm(ssm_c$out_filter,benchmark = rep(0,5),option="filtered_state",range=30:708)
  digits=4
  ssm_c_all=cbind(round(ssm_c$result,digits=digits),
                  paste("(",round(ssm_c$result$Estimate-qnorm(0.995)*ssm_c$result$Std.Error,digits=digits),",",
                        round(ssm_c$result$Estimate+qnorm(0.995)*ssm_c$result$Std.Error,digits=digits),")",sep=""),
                  paste("(",round(ssm_c$result$Estimate-qnorm(0.975)*ssm_c$result$Std.Error,digits=digits),",",
                        round(ssm_c$result$Estimate+qnorm(0.975)*ssm_c$result$Std.Error,digits=digits),")",sep=""),
                  paste("(",round(ssm_c$result$Estimate-qnorm(0.95)*ssm_c$result$Std.Error,digits=digits),",",
                        round(ssm_c$result$Estimate+qnorm(0.95)*ssm_c$result$Std.Error,digits=digits),")",sep=""))
  rownames(ssm_c_all)=rownames(ssm_c$result)
  colnames(ssm_c_all)=c("Estimate","Std.Error","99% CI", "95% CI", "90% CI")
  ssm_c_all
  #                                                       Estimate  Std.Error  99% CI             95% CI            90% CI
  # (Intercept)                                           1.2967    0.6001     (-0.2491,2.8425)   (0.1205,2.4729)   (0.3096,2.2838)
  # logit_TAM_phone_1(period1)                            0.0662    0.0975     (-0.1849,0.3174)   (-0.1249,0.2573)  (-0.0942,0.2266)
  # logit_TAM_phone_1(period2)                            0.1778+   0.0918     (-0.0588,0.4143)   (-0.0022,0.3578)  (0.0267,0.3288)
  # logit_TAM_phone_1(period3)                           -0.1091*   0.0466     (-0.2291,0.0108)   (-0.2004,-0.0179) (-0.1858,-0.0325)
  # keycontacts_call_totaldegree_binary                   0.0122    0.2917     (-0.7391,0.7635)   (-0.5594,0.5838)  (-0.4675,0.4919)
  # keycontacts_text_reciprocity_degree_binary(period1)  -0.0588    0.2852     (-0.7934,0.6758)   (-0.6178,0.5002)  (-0.5279,0.4103)
  # keycontacts_text_reciprocity_degree_binary(period2)  -0.7810*   0.3045     (-1.5653,0.0032)   (-1.3778,-0.1843) (-1.2818,-0.2802)
  # negative_total                                       -0.0141    0.0264     (-0.0821,0.0539)   (-0.0659,0.0376)  (-0.0575,0.0293)
}
# ssm_c plots (used in the Appendix)
{  
  # random walk intercept
  pdf(file=paste(address_appendix,"ssmc_intercept.pdf",sep=""),width = 12, height = 6)
  par(mar = c(5, 5, .5, .5))
  plot(1:708,ssm_c$out_filter$m[,1],type="l",ylab=expression(paste(mu["0,t"]," (intercept)")),xlab="Time (days)",bty="n",cex.lab=2,cex.axis=2,ylim=c(-5,5))
  sd=unlist(lapply(1:708,function(i){sqrt(diag(dlmSvd2var(ssm_c$out_filter$U.C[[i]],ssm_c$out_filter$D.C[i,]))[1])}))
  polygon(c(1:708,rev(1:708)),c(ssm_c$out_filter$m[,1]+qnorm(0.975)*sd,rev(ssm_c$out_filter$m[,1]-qnorm(0.975)*sd)),col="grey90",border="grey")
  points(1:708,ssm_c$out_filter$m[,1],type="l")
  abline(h=0,col="black",lty=3,lwd=2)
  legend("bottomright",legend=c("estimate", "95% CI"),
         pch = c(NA,15), lty=c(1,NA), col=c("black","grey"),
         bty = "n", # remove the bounder of the legend
         lwd = c(1,NA), pt.bg = c(NA,"grey90"),cex=2)
  dev.off()
  
  # 3 period autocorrelation
  pdf(file=paste(address_appendix,"ssmc_rho.pdf",sep=""),width = 12, height = 6)
  par(mar = c(5, 5, .5, .5))
  a=cpt.mean(ssm_c$out_filter$m[,2],penalty="Manual",Q=4,method="BinSeg")
  plot(1:708,ssm_c$out_filter$m[,2],type="l",ylab=expression(paste(rho["pm,t"]," (auto-correlation)")),xlab="Time (days)",bty="n",cex.lab=2,cex.axis=2)
  sd=unlist(lapply(1:708,function(i){sqrt(diag(dlmSvd2var(ssm_c$out_filter$U.C[[i]],ssm_c$out_filter$D.C[i,]))[2])}))
  polygon(c(1:708,rev(1:708)),c(ssm_c$out_filter$m[,2]+qnorm(0.975)*sd,rev(ssm_c$out_filter$m[,2]-qnorm(0.975)*sd)),col="grey90",border="grey")
  points(1:708,ssm_c$out_filter$m[,2],type="l")
  lines(c(1:708)[c(1,cpts(a)[3])],rep(param.est(a)$mean[3],2),col="red")
  lines(c(1:708)[c(cpts(a)[3]+1,cpts(a)[4])],rep(param.est(a)$mean[4],2),col="red")
  lines(c(1:708)[c(cpts(a)[4]+1,708)],rep(param.est(a)$mean[5],2),col="red")
  abline(h=0,col="black",lty=3,lwd=2)
  legend("topright",legend=c("estimate", "95% CI"),
         pch = c(NA,15), lty=c(1,NA), col=c("black","grey"),
         bty = "n", # remove the bounder of the legend
         lwd = c(1,NA), pt.bg = c(NA,"grey90"),cex=2)
  dev.off()
  
  # time-invariant coeffi for degree of calls with keycontacts
  pdf(file=paste(address_appendix,"ssmc_calls.pdf",sep=""),width = 12, height = 6)
  par(mar = c(5, 5, .5, .5))
  plot(1:708,ssm_c$out_filter$m[,3],type="l",ylab=expression(paste(mu["1,t"]," (degree of calls)")),xlab="Time (days)",bty="n",cex.lab=2,cex.axis=2,ylim=c(-1,1))
  sd=unlist(lapply(1:708,function(i){sqrt(diag(dlmSvd2var(ssm_c$out_filter$U.C[[i]],ssm_c$out_filter$D.C[i,]))[3])}))
  polygon(c(1:708,rev(1:708)),c(ssm_c$out_filter$m[,3]+qnorm(0.975)*sd,rev(ssm_c$out_filter$m[,3]-qnorm(0.975)*sd)),col="grey90",border="grey")
  points(1:708,ssm_c$out_filter$m[,3],type="l")
  abline(h=0,col="black",lty=3,lwd=2)
  legend("bottomright",legend=c("estimate", "95% CI"),
         pch = c(NA,15), lty=c(1,NA), col=c("black","grey"),
         bty = "n", # remove the bounder of the legend
         lwd = c(1,NA), pt.bg = c(NA,"grey90"),cex=2)
  dev.off()
  
  # periodic-stable coeff for degree of texts with keycontacts
  pdf(file=paste(address_appendix,"ssmc_texts.pdf",sep=""),width = 12, height = 6)
  par(mar = c(5, 5, .5, .5))
  a=cpt.mean(ssm_c$out_filter$m[,4],penalty="Manual",Q=2,method="BinSeg")
  plot(1:708,ssm_c$out_filter$m[,4],type="l",ylab=expression(paste(mu["2,t"]," (degree of texts)")),xlab="Time (days)",bty="n",cex.lab=2,cex.axis=2, ylim=c(-5,1))
  sd=unlist(lapply(1:708,function(i){sqrt(diag(dlmSvd2var(ssm_c$out_filter$U.C[[i]],ssm_c$out_filter$D.C[i,]))[4])}))
  polygon(c(1:708,rev(1:708)),c(ssm_c$out_filter$m[,4]+qnorm(0.975)*sd,rev(ssm_c$out_filter$m[,4]-qnorm(0.975)*sd)),col="grey90",border="grey")
  points(1:708,ssm_c$out_filter$m[,4],type="l")
  lines(c(1:708)[c(1,cpts(a)[1])],rep(param.est(a)$mean[1],2),col="red")
  lines(c(1:708)[c(cpts(a)[1]+1,708)],rep(param.est(a)$mean[3],2),col="red")
  abline(h=0,col="black",lty=3,lwd=2)
  legend("bottomright",legend=c("estimate", "95% CI"),
         pch = c(NA,15), lty=c(1,NA), col=c("black","grey"),
         bty = "n", # remove the bounder of the legend
         lwd = c(1,NA), pt.bg = c(NA,"grey90"),cex=2)
  dev.off()
  
  # time-invariant coeffi for negative mood
  pdf(file=paste(address_appendix,"ssmc_negative.pdf",sep=""),width = 12, height = 6)
  par(mar = c(5, 5, .5, .5))
  plot(1:708,ssm_c$out_filter$m[,5],type="l",ylab=expression(paste(mu["3,t"]," (negative mood)")),xlab="Time (days)",bty="n",cex.lab=2,cex.axis=2,ylim=c(-0.7,0.5))
  sd=unlist(lapply(1:708,function(i){sqrt(diag(dlmSvd2var(ssm_c$out_filter$U.C[[i]],ssm_c$out_filter$D.C[i,]))[5])}))
  polygon(c(1:708,rev(1:708)),c(ssm_c$out_filter$m[,5]+qnorm(0.975)*sd,rev(ssm_c$out_filter$m[,5]-qnorm(0.975)*sd)),col="grey90",border="grey")
  points(1:708,ssm_c$out_filter$m[,5],type="l")
  abline(h=0,col="black",lty=3,lwd=2)
  legend("bottomright",legend=c("estimate", "95% CI"),
         pch = c(NA,15), lty=c(1,NA), col=c("black","grey"),
         bty = "n", # remove the bounder of the legend
         lwd = c(1,NA), pt.bg = c(NA,"grey90"),cex=2)
  dev.off()
}
# outcome regression: negative_total~ intercept + negative_total_1 
#                                     + keycontacts_call_totaldegree_binary + keycontacts_call_totaldegree_binary_1
#                                     + keycontacts_text_reciprocity_degree_binary + keycontacts_text_reciprocity_degree_binary
#                                     + logit_TAM_phone_1
# ssm_y model (used in the Main paper Tables 2&3)
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
  formula=" negative_total ~ negative_total_1 + keycontacts_call_totaldegree_binary + keycontacts_call_totaldegree_binary_1 + keycontacts_text_reciprocity_degree_binary + keycontacts_text_reciprocity_degree_binary_1+ logit_TAM_phone_1"
  all_var=unlist(strsplit(gsub(" ","",formula),"\\+|\\~"))
  outcome_var=strsplit(gsub(" ","",formula),"\\~")[[1]][1]
  formula_var=all_var[!(all_var %in% c(outcome_var,"intercept"))]
  data_ss_ignore_y=data[-1,all_var]
  data_ss_ignore_y[rowSums(is.na(data_ss_ignore_y))!=0,outcome_var]=NA
  for(l in formula_var){
    data_ss_ignore_y[is.na(data_ss_ignore_y[,l]),l]=mean(data_ss_ignore_y[,l],na.rm=T)
  }
  ss_param_y=list(inits=c(log(0.2),log(2.5)),m0=ssm_y_init$att[708,],C0=diag(rep(10^3),7),
                     AR1_coeffi=NULL,rw_coeffi=c("intercept"),v_cp_param=NULL,
                     w_cp_param=list(list(variable="keycontacts_text_reciprocity_degree_binary",segments=3,changepoints=c(560,640),fixed_cpts=F),
                                     list(variable="keycontacts_text_reciprocity_degree_binary_1",segments=2,changepoints=c(460),fixed_cpts=F)))
  ssm_y=run.SSM(data_ss=data_ss_ignore_y,formula=formula,ss_param=ss_param_y,max_iteration=100,
                cpt_learning_param=list(cpt_method="mean",burnin=1/10,mergeband=20,convergence_cri=10),
                cpt_merge_option="separate",dlm_cpt_learning_option="filter",
                bandwidth=5,cpt_V =1,printFlag=T)
  plot.dlm(ssm_y$out_filter,benchmark = rep(0,7),option="filtered_state",range=30:708)
  digits=3
  ssm_y_all=cbind(round(ssm_y$result,digits=digits),
                  paste("(",round(ssm_y$result$Estimate-qnorm(1-0.01/2)*ssm_y$result$Std.Error,digits=digits),",",
                        round(ssm_y$result$Estimate+qnorm(1-0.01/2)*ssm_y$result$Std.Error,digits=digits),")",sep=""),
                  paste("(",round(ssm_y$result$Estimate-qnorm(1-0.05/2)*ssm_y$result$Std.Error,digits=digits),",",
                         round(ssm_y$result$Estimate+qnorm(1-0.05/2)*ssm_y$result$Std.Error,digits=digits),")",sep=""),
                  paste("(",round(ssm_y$result$Estimate-qnorm(1-0.1/2)*ssm_y$result$Std.Error,digits=digits),",",
                         round(ssm_y$result$Estimate+qnorm(1-0.1/2)*ssm_y$result$Std.Error,digits=digits),")",sep=""))
  rownames(ssm_y_all)=rownames(ssm_y$result)
  colnames(ssm_y_all)=c("Estimate","Std.Error","99% CI","95% CI", "90% CI")
  ssm_y_all
  #                                                          Estimate  Std.Error  99% CI          95% CI          90% CI
  # (Intercept)                                              5.004     0.880      (2.738,7.269)   (3.28,6.728)    (3.557,6.451)
  # negative_total_1                                         0.633**   0.038      (0.536,0.731)   (0.559,0.707)   (0.571,0.695)
  # keycontacts_call_totaldegree_binary                     -0.218     0.254      (-0.872,0.436)  (-0.715,0.279)  (-0.635,0.199)
  # keycontacts_call_totaldegree_binary_1                    0.088     0.257      (-0.575,0.751)  (-0.416,0.592)  (-0.335,0.511)
  # keycontacts_text_reciprocity_degree_binary(period1)     -0.086     0.210      (-0.627,0.454)  (-0.498,0.325)  (-0.431,0.259)
  # keycontacts_text_reciprocity_degree_binary(period2)     -1.154+    0.608      (-2.72,0.412)   (-2.346,0.038)  (-2.154,-0.154)
  # keycontacts_text_reciprocity_degree_binary(period3)     -0.739     0.583      (-2.242,0.764)  (-1.883,0.405)  (-1.699,0.221)
  # keycontacts_text_reciprocity_degree_binary_1(period1)   -0.165     0.238      (-0.777,0.447)  (-0.631,0.301)  (-0.556,0.226)
  # keycontacts_text_reciprocity_degree_binary_1(period2)   -0.724*    0.312      (-1.529,0.08)   (-1.336,-0.113) (-1.238,-0.211)
  # logit_TAM_phone_1                                       -0.012     0.036      (-0.105,0.08)   (-0.083,0.058)  (-0.072,0.047)
}
# ssm_y plots (used in the Appendix)
{
  # random walk intercept
  pdf(file=paste(address_appendix,"ssmy_intercept.pdf",sep=""),width = 12, height = 6)
  par(mar = c(5, 5, .5, .5))
  plot(1:708,ssm_y$out_filter$m[,1],type="l",ylab=expression(paste(beta["0,t"]," (intercept)")),xlab="Time (days)" ,bty="n",cex.lab=2,cex.axis=2, ylim=c(-3,8))
  sd=unlist(lapply(1:708,function(i){sqrt(diag(dlmSvd2var(ssm_y$out_filter$U.C[[i]],ssm_y$out_filter$D.C[i,]))[1])}))
  polygon(c(1:708,rev(1:708)),c(ssm_y$out_filter$m[,1]+qnorm(0.975)*sd,rev(ssm_y$out_filter$m[,1]-qnorm(0.975)*sd)),col="grey90",border="grey")
  points(1:708,ssm_y$out_filter$m[,1],type="l")
  abline(h=0,col="black",lty=3,lwd=2)
  legend("bottomright",legend=c("estimate", "95% CI"),
         pch = c(NA,15), lty=c(1,NA), col=c("black","grey"),
         bty = "n", # remove the bounder of the legend
         lwd = c(1,NA), pt.bg = c(NA,"grey90"),cex=2)
  dev.off()
  
  # time-invariant autocorrelation
  pdf(file=paste(address_appendix,"ssmy_rho.pdf",sep=""),width = 12, height = 6)
  par(mar = c(5, 5, .5, .5))
  plot(1:708,ssm_y$out_filter$m[,2],type="l",ylab=expression(paste(rho[t]," (auto-correlation)")),xlab="Time (days)" ,bty="n",cex.lab=2,cex.axis=2, ylim=c(0.0,1))
  sd=unlist(lapply(1:708,function(i){sqrt(diag(dlmSvd2var(ssm_y$out_filter$U.C[[i]],ssm_y$out_filter$D.C[i,]))[2])}))
  polygon(c(1:708,rev(1:708)),c(ssm_y$out_filter$m[,2]+qnorm(0.975)*sd,rev(ssm_y$out_filter$m[,2]-qnorm(0.975)*sd)),col="grey90",border="grey")
  points(1:708,ssm_y$out_filter$m[,2],type="l")
  abline(h=0,col="black",lty=3,lwd=2)
  legend("bottomright",legend=c("estimate", "95% CI"),
         pch = c(NA,15), lty=c(1,NA), col=c("black","grey"),
         bty = "n", # remove the bounder of the legend
         lwd = c(1,NA), pt.bg = c(NA,"grey90"),cex=2)
  dev.off()
  
  # time-invariant coeffi for degree of calls with keycontacts
  pdf(file=paste(address_appendix,"ssmy_calls.pdf",sep=""),width = 12, height = 6)
  par(mar = c(5, 5, .5, .5))
  plot(1:708,ssm_y$out_filter$m[,3],type="l",ylab=expression(paste(beta["11,t"]," (degree of calls)")),xlab="Time (days)" ,bty="n",cex.lab=2,cex.axis=2, ylim=c(-1,1))
  sd=unlist(lapply(1:708,function(i){sqrt(diag(dlmSvd2var(ssm_y$out_filter$U.C[[i]],ssm_y$out_filter$D.C[i,]))[3])}))
  polygon(c(1:708,rev(1:708)),c(ssm_y$out_filter$m[,3]+qnorm(0.975)*sd,rev(ssm_y$out_filter$m[,3]-qnorm(0.975)*sd)),col="grey90",border="grey")
  points(1:708,ssm_y$out_filter$m[,3],type="l")
  abline(h=0,col="black",lty=3,lwd=2)
  legend("bottomright",legend=c("estimate", "95% CI"),
         pch = c(NA,15), lty=c(1,NA), col=c("black","grey"),
         bty = "n", # remove the bounder of the legend
         lwd = c(1,NA), pt.bg = c(NA,"grey90"),cex=2)
  dev.off()
  
  # time-invariant coeffi for degree of calls yesterday with keycontacts
  pdf(file=paste(address_appendix,"ssmy_calls2.pdf",sep=""),width = 12, height = 6)
  par(mar = c(5, 5, .5, .5))
  plot(1:708,ssm_y$out_filter$m[,4],type="l",ylab=expression(paste(beta["12,t"]," (degree of calls previous day)")),xlab="Time (days)" ,bty="n",cex.lab=2,cex.axis=2, ylim=c(-1,1))
  sd=unlist(lapply(1:708,function(i){sqrt(diag(dlmSvd2var(ssm_y$out_filter$U.C[[i]],ssm_y$out_filter$D.C[i,]))[4])}))
  polygon(c(1:708,rev(1:708)),c(ssm_y$out_filter$m[,4]+qnorm(0.975)*sd,rev(ssm_y$out_filter$m[,4]-qnorm(0.975)*sd)),col="grey90",border="grey")
  points(1:708,ssm_y$out_filter$m[,4],type="l")
  abline(h=0,col="black",lty=3,lwd=2)
  legend("bottomright",legend=c("estimate", "95% CI"),
         pch = c(NA,15), lty=c(1,NA), col=c("black","grey"),
         bty = "n", # remove the bounder of the legend
         lwd = c(1,NA), pt.bg = c(NA,"grey90"),cex=2)
  dev.off()
  
  # periodic-stable coeff for degree of texts
  pdf(file=paste(address_appendix,"ssmy_texts.pdf",sep=""),width = 12, height = 6)
  par(mar = c(5, 5, .5, .5))
  a=cpt.mean(ssm_y$out_filter$m[,5],penalty="Manual",Q=2,method="BinSeg")
  plot(1:708,ssm_y$out_filter$m[,5],type="l",ylab=expression(paste(beta["21,t"]," (degree of texts)")),xlab="Time (days)" ,bty="n",cex.lab=2,cex.axis=2, ylim=c(-4,1))
  sd=unlist(lapply(1:708,function(i){sqrt(diag(dlmSvd2var(ssm_y$out_filter$U.C[[i]],ssm_y$out_filter$D.C[i,]))[5])}))
  polygon(c(1:708,rev(1:708)),c(ssm_y$out_filter$m[,5]+qnorm(0.975)*sd,rev(ssm_y$out_filter$m[,5]-qnorm(0.975)*sd)),col="grey90",border="grey")
  points(1:708,ssm_y$out_filter$m[,5],type="l")
  lines(c(1:708)[c(1,cpts(a)[1])],rep(param.est(a)$mean[1],2),col="red")
  lines(c(1:708)[c(cpts(a)[1]+1,cpts(a)[2])],rep(param.est(a)$mean[2],2),col="red")
  lines(c(1:708)[c(cpts(a)[2]+1,708)],rep(param.est(a)$mean[3],2),col="red")
  abline(h=0,col="black",lty=3,lwd=2)
  legend("bottomright",legend=c("estimate", "95% CI"),
         pch = c(NA,15), lty=c(1,NA), col=c("black","grey"),
         bty = "n", # remove the bounder of the legend
         lwd = c(1,NA), pt.bg = c(NA,"grey90"),cex=2)
  dev.off()
  
  # time-invariant coeffi for degree of texts yesterday
  pdf(file=paste(address_appendix,"ssmy_texts2.pdf",sep=""),width = 12, height = 6)
  par(mar = c(5, 5, .5, .5))
  a=cpt.mean(ssm_y$out_filter$m[,6],penalty="Manual",Q=1,method="BinSeg")
  plot(1:708,ssm_y$out_filter$m[,6],type="l",ylab=expression(paste(beta["22,t"]," (degree of outgoing texts previous day)")),xlab="Time (days)" ,bty="n",cex.lab=1.5,cex.axis=2, ylim=c(-3,2))
  sd=unlist(lapply(1:708,function(i){sqrt(diag(dlmSvd2var(ssm_y$out_filter$U.C[[i]],ssm_y$out_filter$D.C[i,]))[6])}))
  polygon(c(1:708,rev(1:708)),c(ssm_y$out_filter$m[,6]+qnorm(0.975)*sd,rev(ssm_y$out_filter$m[,6]-qnorm(0.975)*sd)),col="grey90",border="grey")
  points(1:708,ssm_y$out_filter$m[,6],type="l")
  lines(c(1:708)[c(1,cpts(a)[1])],rep(param.est(a)$mean[1],2),col="red")
  lines(c(1:708)[c(cpts(a)[1]+1,708)],rep(param.est(a)$mean[2],2),col="red")
  abline(h=0,col="black",lty=3,lwd=2)
  legend("bottomright",legend=c("estimate", "95% CI"),
         pch = c(NA,15), lty=c(1,NA), col=c("black","grey"),
         bty = "n", # remove the bounder of the legend
         lwd = c(1,NA), pt.bg = c(NA,"grey90"),cex=2)
  dev.off()
  
  # time-invariant coeffi for physical activity
  pdf(file=paste(address_appendix,"ssmy_pa.pdf",sep=""),width = 12, height = 6)
  par(mar = c(5, 5, .5, .5))
  plot(1:708,ssm_y$out_filter$m[,7],type="l",ylab=expression(paste(beta["PA,t"]," (daily physical activity)")),xlab="Time (days)" ,bty="n",cex.lab=2,cex.axis=2, ylim=c(-1,1))
  sd=unlist(lapply(1:708,function(i){sqrt(diag(dlmSvd2var(ssm_y$out_filter$U.C[[i]],ssm_y$out_filter$D.C[i,]))[7])}))
  polygon(c(1:708,rev(1:708)),c(ssm_y$out_filter$m[,7]+qnorm(0.975)*sd,rev(ssm_y$out_filter$m[,7]-qnorm(0.975)*sd)),col="grey90",border="grey")
  points(1:708,ssm_y$out_filter$m[,7],type="l")
  abline(h=0,col="black",lty=3,lwd=2)
  legend("bottomright",legend=c("estimate", "95% CI"),
         pch = c(NA,15), lty=c(1,NA), col=c("black","grey"),
         bty = "n", # remove the bounder of the legend
         lwd = c(1,NA), pt.bg = c(NA,"grey90"),cex=2)
  dev.off()
}

##############################################################################################
######     simulate/compute causal estimands of calls with CI (used in Appendix)         #####
##############################################################################################
y_coeffi_table_call=as.data.frame(ssm_y$out_smooth$s);colnames(y_coeffi_table_call)[c(1:4,7)]=c("(Intercept)","y_1","x","x_1","c");head(y_coeffi_table_call)
c_coeffi_table_call=as.data.frame(ssm_c$out_smooth$s);colnames(c_coeffi_table_call)[c(1,2,3,5)]=c("(Intercept)","c_1","x_1","y_1");head(c_coeffi_table_call)
raw_data_call=data[,c(2,3,4,6)];colnames(raw_data_call)=c("Date","y","x","c");head(raw_data_call)
y_coeffi_var_table_call=lapply(1:nrow(ssm_y$out_smooth$s),function(i){dlmSvd2var(ssm_y$out_smooth$U.S[[i]],ssm_y$out_smooth$D.S[i,])});head(y_coeffi_var_table_call)
c_coeffi_var_table_call=lapply(1:nrow(ssm_c$out_smooth$s),function(i){dlmSvd2var(ssm_c$out_smooth$U.S[[i]],ssm_c$out_smooth$D.S[i,])});head(c_coeffi_var_table_call)
time=1:nrow(ssm_y$out_smooth$s)
# calculate call's contemporaneous effect with CI (save call_contemporaneous)
call_contemporaneous = calculate.causaleffect(t=time,tx=1,y_coeffi_table=y_coeffi_table_call,c_coeffi_table=c_coeffi_table_call,
                                 CI=T,n_sim=1000,y_coeffi_var_table=y_coeffi_var_table_call,c_coeffi_var_table=c_coeffi_var_table_call,seed=1,printFlag=T)
# calculate call's 1-lag effect with CI (save call_1_lag)
call_1_lag = calculate.causaleffect(t=time,tx=c(1,0),y_coeffi_table=y_coeffi_table_call,c_coeffi_table=c_coeffi_table_call,
                                    CI=T,n_sim=1000,y_coeffi_var_table=y_coeffi_var_table_call,c_coeffi_var_table=c_coeffi_var_table_call,seed=2,printFlag=T)
# calculate call's 2-lag effect with CI (save call_2_lag)
call_2_lag = calculate.causaleffect(t=time,tx=c(1,0,0),y_coeffi_table=y_coeffi_table_call,c_coeffi_table=c_coeffi_table_call,
                                    CI=T,n_sim=1000,y_coeffi_var_table=y_coeffi_var_table_call,c_coeffi_var_table=c_coeffi_var_table_call,seed=2,printFlag=T)
# calculate call's 3-lag effect with CI (save call_3_lag)
call_3_lag = calculate.causaleffect(t=time,tx=c(1,0,0,0),y_coeffi_table=y_coeffi_table_call,c_coeffi_table=c_coeffi_table_call,
                                    CI=T,n_sim=1000,y_coeffi_var_table=y_coeffi_var_table_call,c_coeffi_var_table=c_coeffi_var_table_call,seed=2,printFlag=T)
# calculate call's 1-lag structural effect with CI (save call_1_lag_structural_direct)
call_1_lag_structural_direct = calculate.controlled_direct_effect(t=time,y_coeffi_table=y_coeffi_table_call,
                                         CI=T,n_sim=1000,y_coeffi_var_table=y_coeffi_var_table_call,seed=3,printFlag=T)
# calculate call's 2-step total effect with CI (save call_2_step)
call_1_step = calculate.causaleffect(t=time,tx=c(1,1),y_coeffi_table=y_coeffi_table_call,c_coeffi_table=c_coeffi_table_call,
                       CI=T,n_sim=1000,y_coeffi_var_table=y_coeffi_var_table_call,c_coeffi_var_table=c_coeffi_var_table_call,seed=4,printFlag=T)
# calculate call's 3-step general effect with CI (save call_3_step_general_101)
call_2_step_general_101 = calculate.causaleffect(t=time,tx=c(1,0,1),y_coeffi_table=y_coeffi_table_call,c_coeffi_table=c_coeffi_table_call,
                       CI=T,n_sim=1000,y_coeffi_var_table=y_coeffi_var_table_call,c_coeffi_var_table=c_coeffi_var_table_call,seed=4,printFlag=T)

# impulse impact graph at t=600 without CI (save call_effects_vs_qlag)
call_qlag_effects_part1 = simulate.counterfactual_path_singlet(t=600,tx=c(1,rep(0,10)),y_coeffi_table=y_coeffi_table_call,c_coeffi_table=c_coeffi_table_call,raw_data=raw_data_call,CI=F)
call_qlag_effects_part2 = simulate.counterfactual_path_singlet(t=600,tx=c(0,rep(0,10)),y_coeffi_table=y_coeffi_table_call,c_coeffi_table=c_coeffi_table_call,raw_data=raw_data_call,CI=F)
# Note: tried n_sim=1000 and n_sim=10000 before, results not accurate enough -> end up useing n_sim
# call_qlag_effects_part1 = simulate.counterfactual_path_singlet(t=600,tx=c(1,rep(0,10)),y_coeffi_table=y_coeffi_table_call,c_coeffi_table=c_coeffi_table_call,raw_data=raw_data_call,
#                                               CI=T,n_sim=100000,y_coeffi_var_table=y_coeffi_var_table_call,c_coeffi_var_table=c_coeffi_var_table_call,seed=1,printFlag=T)
# call_qlag_effects_part2 = simulate.counterfactual_path_singlet(t=600,tx=c(0,rep(0,10)),y_coeffi_table=y_coeffi_table_call,c_coeffi_table=c_coeffi_table_call,raw_data=raw_data_call,
#                                               CI=T,n_sim=100000,y_coeffi_var_table=y_coeffi_var_table_call,c_coeffi_var_table=c_coeffi_var_table_call,seed=1,printFlag=T)
call_effects_vs_qlag = call_qlag_effects_part1-call_qlag_effects_part2
save(call_contemporaneous,
     call_1_lag,call_2_lag,call_3_lag,call_1_lag_structural_direct,
     call_1_step,call_2_step_general_101,
     call_effects_vs_qlag,
     file="/Users/xiaoxuancai/Documents/GitHub/Causal_estimands/result_calls.Rdata")

# compare results from analytical solution with CI, Monte Carlo simulation without CI, and true value (not used, just for checking)
par(mfrow=c(1,1),mar=c(4, 4, 2, 0.5))
plot(0:10,call_effects_vs_qlag,main="Analytical vs Monte Carlo vs Ture",ylab="q-lag effect",xlab="#lags",cex=1.5,cex.lab=1.5,cex.axis= 1.5, bty="n",col="black",pch=5)
points(0:3,c(rowMeans(call_contemporaneous)[600],rowMeans(call_1_lag)[600],rowMeans(call_2_lag)[600],rowMeans(call_3_lag)[600]),col="blue",cex=1.5,pch=16)
points(0:3,c(calculate.causaleffect(t=600,tx=1,y_coeffi_table=y_coeffi_table_call,c_coeffi_table=c_coeffi_table_call,CI=F),
             calculate.causaleffect(t=600,tx=c(1,0),y_coeffi_table=y_coeffi_table_call,c_coeffi_table=c_coeffi_table_call,CI=F),
             calculate.causaleffect(t=600,tx=c(1,0,0),y_coeffi_table=y_coeffi_table_call,c_coeffi_table=c_coeffi_table_call,CI=F),
             calculate.causaleffect(t=600,tx=c(1,0,0,0),y_coeffi_table=y_coeffi_table_call,c_coeffi_table=c_coeffi_table_call,CI=F)),
       col="red",pch=15,cex=1.5)
legend("bottomright",legend=c("Monte Carlo without CI","Analytical with CI","Truth"),
       pch = c(5,2,15), lty=c(1,1,1), col=c("black","blue","red"),
       bty = "n", # remove the bounder of the legend
       lwd = c(NA,NA,NA),cex=1.5)

# plots (used in the Appendix)
{
  # jump directly to load the data without recalculation!!
  load("/Users/xiaoxuancai/Documents/GitHub/Causal_estimands/result_calls.Rdata")

  # call's contemporaneous effect (used in appendix)
  pdf(file=paste(address_appendix,"contemp_calls.pdf",sep=""),width = 9, height = 6)
  par(mar = c(2.5, 5, .5, .5))
  call_contemporaneous_CIband=plot_simulatedCI(call_contemporaneous,probs=c(0.05,0.95),printFlag=F)
  plot(1:708,call_contemporaneous_CIband$mean,type="l",ylab="contemporaneous effect (calls)",xlab="",bty="n", cex.axis=2,cex.lab=2,ylim=c(-1,0.5))
  polygon(c(1:708,rev(1:708)),c(call_contemporaneous_CIband$upper,rev(call_contemporaneous_CIband$lower)),col="grey90",border="grey")
  points(1:708,call_contemporaneous_CIband$mean,type="l")
  abline(h=0,lty=3,lwd=2)
  legend("bottomright",legend=c("estimate", "95% CI"),
         pch = c(NA,15), lty=c(1,NA), col=c("black","grey"),
         bty = "n", # remove the bounder of the legend
         lwd = c(1,NA), pt.bg = c(NA,"grey90"),cex=2)
  dev.off()
  
  # call's 1-lag effect (used in appendix)
  pdf(file=paste(address_appendix,"lag1_calls.pdf",sep=""),width = 9, height = 6)
  par(mar = c(2.5, 5, .5, .5))
  call_1_lag_CIband=plot_simulatedCI(call_1_lag,probs=c(0.05,0.95),printFlag=F)
  plot(1:708,call_1_lag_CIband$mean,type="l",ylab="1-lag effect (calls)",xlab="Date",bty="n", cex.axis=2,cex.lab=2,ylim=c(-1,0.5))
  polygon(c(1:708,rev(1:708)),c(call_1_lag_CIband$upper,rev(call_1_lag_CIband$lower)),col="grey90",border="grey")
  points(1:708,call_1_lag_CIband$mean,type="l")
  abline(h=0,lty=3,lwd=2)
  legend("bottomright",legend=c("estimate", "95% CI"),
         pch = c(NA,15), lty=c(1,NA), col=c("black","grey"),
         bty = "n", # remove the bounder of the legend
         lwd = c(1,NA), pt.bg = c(NA,"grey90"),cex=2)
  dev.off()
  
  # call's 2-lag effect (used in appendix)
  pdf(file=paste(address_appendix,"lag2_calls.pdf",sep=""),width = 9, height = 6)
  par(mar = c(2.5, 5, .5, .5))
  call_2_lag_CIband=plot_simulatedCI(call_2_lag,probs=c(0.05,0.95),printFlag=F)
  plot(1:708,call_2_lag_CIband$mean,type="l",ylab="2-lag effect (calls)",xlab="Date",bty="n", cex.axis=2,cex.lab=2,ylim=c(-1,0.5))
  polygon(c(1:708,rev(1:708)),c(call_2_lag_CIband$upper,rev(call_2_lag_CIband$lower)),col="grey90",border="grey")
  points(1:708,call_2_lag_CIband$mean,type="l")
  abline(h=0,lty=3,lwd=2)
  legend("bottomright",legend=c("estimate", "95% CI"),
         pch = c(NA,15), lty=c(1,NA), col=c("black","grey"),
         bty = "n", # remove the bounder of the legend
         lwd = c(1,NA), pt.bg = c(NA,"grey90"),cex=2)
  dev.off()
  
  # call's 3-lag effect (used in appendix)
  pdf(file=paste(address_appendix,"lag3_calls.pdf",sep=""),width = 9, height = 6)
  par(mar = c(2.5, 5, .5, .5))
  call_3_lag_CIband=plot_simulatedCI(call_3_lag,probs=c(0.05,0.95),printFlag=F)
  plot(1:708,call_3_lag_CIband$mean,type="l",ylab="3-lag effect (calls)",xlab="Date",bty="n", cex.axis=2,cex.lab=2,ylim=c(-1,0.5))
  polygon(c(1:708,rev(1:708)),c(call_3_lag_CIband$upper,rev(call_3_lag_CIband$lower)),col="grey90",border="grey")
  points(1:708,call_3_lag_CIband$mean,type="l")
  abline(h=0,lty=3,lwd=2)
  legend("bottomright",legend=c("estimate", "95% CI"),
         pch = c(NA,15), lty=c(1,NA), col=c("black","grey"),
         bty = "n", # remove the bounder of the legend
         lwd = c(1,NA), pt.bg = c(NA,"grey90"),cex=2)
  dev.off()
  
  # call's 1-lag structural direct effect
  pdf(file=paste(address_appendix,"lag1_controlled_direct_effect_calls.pdf",sep=""),width = 9, height = 6)
  par(mar = c(2.5, 5, .5, .5))
  call_1_lag_structural_direct_CIband = plot_simulatedCI(call_1_lag_structural_direct,probs=c(0.05,0.95),printFlag=F)
  plot(1:708,call_1_lag_structural_direct_CIband$mean,type="l",ylab="1-lag controlled direct effect (calls)",xlab="# lags",bty="n",cex.axis=2,cex.lab=2,ylim=c(-.7,0.7))
  polygon(c(1:708,rev(1:708)),c(call_1_lag_structural_direct_CIband$upper,rev(call_1_lag_structural_direct_CIband$lower)),col="grey90",border="grey")
  points(1:708,call_1_lag_structural_direct_CIband$mean,type="l")
  abline(h=0,lty=3,lwd=2)
  legend("bottomright",legend=c("estimate", "95% CI"),
         pch = c(NA,15), lty=c(1,NA), col=c("black","grey"),
         bty = "n", # remove the bounder of the legend
         lwd = c(1,NA), pt.bg = c(NA,"grey90"),cex=2)
  dev.off()
  
  # call's 2-step total effect
  pdf(file=paste(address_appendix,"step1_total_calls.pdf",sep=""),width = 9, height = 6)
  par(mar = c(2.5, 5, .5, .5))
  call_1_step_CIband = plot_simulatedCI(call_1_step,probs=c(0.05,0.95),printFlag=F)
  plot(1:708,call_1_step_CIband$mean,type="l",ylab="1-step total effect (calls)",xlab="# lags",bty="n",cex.axis=2,cex.lab=2,ylim=c(-2,1))
  polygon(c(1:708,rev(1:708)),c(call_1_step_CIband$upper,rev(call_1_step_CIband$lower)),col="grey90",border="grey")
  points(1:708,call_1_step_CIband$mean,type="l")
  abline(h=0,lty=3,lwd=2)
  legend("bottomright",legend=c("estimate", "95% CI"),
         pch = c(NA,15), lty=c(1,NA), col=c("black","grey"),
         bty = "n", # remove the bounder of the legend
         lwd = c(1,NA), pt.bg = c(NA,"grey90"),cex=2)
  dev.off()
  
  # call's 3-step general effect
  pdf(file=paste(address_appendix,"step2_general_101_calls.pdf",sep=""),width = 9, height = 6)
  par(mar = c(2.5, 5, .5, .5))
  call_2_step_general_101_CIband = plot_simulatedCI(call_2_step_general_101,probs=c(0.05,0.95),printFlag=F)
  plot(1:708,call_2_step_general_101_CIband$mean,type="l",ylab="2-step general effect (calls)",xlab="# lags",bty="n",cex.axis=2,cex.lab=2,ylim=c(-2,1))
  polygon(c(1:708,rev(1:708)),c(call_2_step_general_101_CIband$upper,rev(call_2_step_general_101_CIband$lower)),col="grey90",border="grey")
  points(1:708,call_2_step_general_101_CIband$mean,type="l")
  abline(h=0,lty=3,lwd=2)
  legend("bottomright",legend=c("estimate", "95% CI"),
         pch = c(NA,15), lty=c(1,NA), col=c("black","grey"),
         bty = "n", # remove the bounder of the legend
         lwd = c(1,NA), pt.bg = c(NA,"grey90"),cex=2)
  dev.off()
}

##############################################################################################################
######          simulate causal estimands of texts with CI  (used in the Main Paper and Appendix)        #####
##############################################################################################################
y_coeffi_table_text=as.data.frame(ssm_y$out_smooth$s);colnames(y_coeffi_table_text)[c(1:2,5:7)]=c("(Intercept)","y_1","x","x_1","c");head(y_coeffi_table_text)
c_coeffi_table_text=as.data.frame(ssm_c$out_smooth$s);colnames(c_coeffi_table_text)[c(1,2,4,5)]=c("(Intercept)","c_1","x_1","y_1");head(c_coeffi_table_text)
raw_data_text=data[,c(2,3,5,7)];colnames(raw_data_text)=c("Date","y","x","c");head(raw_data_text)
y_coeffi_var_table_text=y_coeffi_var_table_call
c_coeffi_var_table_text=c_coeffi_var_table_call
time=1:nrow(ssm_y$out_smooth$s)
# calculate text's contemporaneous effect with CI (save text_contemporaneous)
text_contemporaneous = calculate.causaleffect(t=time,tx=1,y_coeffi_table=y_coeffi_table_text,c_coeffi_table=c_coeffi_table_text,
                                              CI=T,n_sim=1000,y_coeffi_var_table=y_coeffi_var_table_text,c_coeffi_var_table=c_coeffi_var_table_text,seed=1,printFlag=T)
# calculate text's 1-lag effect with CI (save text_1_lag)
text_1_lag = calculate.causaleffect(t=time,tx=c(1,0),y_coeffi_table=y_coeffi_table_text,c_coeffi_table=c_coeffi_table_text,
                                    CI=T,n_sim=1000,y_coeffi_var_table=y_coeffi_var_table_text,c_coeffi_var_table=c_coeffi_var_table_text,seed=2,printFlag=T)
# calculate text's 2-lag effect with CI (save text_2_lag)
text_2_lag = calculate.causaleffect(t=time,tx=c(1,0,0),y_coeffi_table=y_coeffi_table_text,c_coeffi_table=c_coeffi_table_text,
                                    CI=T,n_sim=1000,y_coeffi_var_table=y_coeffi_var_table_text,c_coeffi_var_table=c_coeffi_var_table_text,seed=2,printFlag=T)
# calculate text's 3-lag effect with CI (save text_3_lag)
text_3_lag = calculate.causaleffect(t=time,tx=c(1,0,0,0),y_coeffi_table=y_coeffi_table_text,c_coeffi_table=c_coeffi_table_text,
                                    CI=T,n_sim=1000,y_coeffi_var_table=y_coeffi_var_table_text,c_coeffi_var_table=c_coeffi_var_table_text,seed=2,printFlag=T)
# calculate text's 1-lag structural effect with CI (save text_1_lag_structural_direct)
text_1_lag_structural_direct = calculate.controlled_direct_effect(t=time,y_coeffi_table=y_coeffi_table_text,
                                                                  CI=T,n_sim=1000,y_coeffi_var_table=y_coeffi_var_table_text,seed=3,printFlag=T)
# calculate text's 2-step total effect with CI (save text_1_step)
text_1_step = calculate.causaleffect(t=time,tx=c(1,1),y_coeffi_table=y_coeffi_table_text,c_coeffi_table=c_coeffi_table_text,
                                     CI=T,n_sim=1000,y_coeffi_var_table=y_coeffi_var_table_text,c_coeffi_var_table=c_coeffi_var_table_text,seed=4,printFlag=T)
# calculate text's 2-step general effect with CI
text_2_step_general_101 = calculate.causaleffect(t=time,tx=c(1,0,1),y_coeffi_table=y_coeffi_table_text,c_coeffi_table=c_coeffi_table_text,
                                              CI=T,n_sim=1000,y_coeffi_var_table=y_coeffi_var_table_text,c_coeffi_var_table=c_coeffi_var_table_text,seed=5,printFlag=T)
# calculate text's 3-step general effect with CI
text_3_step_general_0101 = calculate.causaleffect(t=time,tx=c(0,1,0,1),y_coeffi_table=y_coeffi_table_text,c_coeffi_table=c_coeffi_table_text,
                                                  CI=T,n_sim=1000,y_coeffi_var_table=y_coeffi_var_table_text,c_coeffi_var_table=c_coeffi_var_table_text,seed=6,printFlag=T)
# impulse impact graph at t=600 (save text_effects_vs_qlag)
text_qlag_effects_part1 = simulate.counterfactual_path_singlet(t=600,tx=c(1,rep(0,10)),y_coeffi_table=y_coeffi_table_text,c_coeffi_table=c_coeffi_table_text,raw_data=raw_data_text,
                                                              CI=T,n_sim=10000,y_coeffi_var_table=y_coeffi_var_table_text,c_coeffi_var_table=c_coeffi_var_table_text,seed=1,
                                                              printFlag=T)
text_qlag_effects_part2 = simulate.counterfactual_path_singlet(t=600,tx=c(0,rep(0,10)),y_coeffi_table=y_coeffi_table_text,c_coeffi_table=c_coeffi_table_text,raw_data=raw_data_text,
                                                               CI=T,n_sim=10000,y_coeffi_var_table=y_coeffi_var_table_text,c_coeffi_var_table=c_coeffi_var_table_text,seed=1,
                                                               printFlag=T)
text_effects_vs_qlag = text_qlag_effects_part1-text_qlag_effects_part2
# step response plot at t=600
text_qstep_effects_part1 = simulate.counterfactual_path_singlet(t=600,tx=c(1,rep(1,10)),y_coeffi_table=y_coeffi_table_text,c_coeffi_table=c_coeffi_table_text,raw_data=raw_data_text,
                                                               CI=T,n_sim=10000,y_coeffi_var_table=y_coeffi_var_table_text,c_coeffi_var_table=c_coeffi_var_table_text,seed=1,
                                                               printFlag=T)
text_qstep_effects_part2 = simulate.counterfactual_path_singlet(t=600,tx=c(0,rep(0,10)),y_coeffi_table=y_coeffi_table_text,c_coeffi_table=c_coeffi_table_text,raw_data=raw_data_text,
                                                                CI=T,n_sim=10000,y_coeffi_var_table=y_coeffi_var_table_text,c_coeffi_var_table=c_coeffi_var_table_text,seed=1,
                                                                printFlag=T)
text_effects_vs_qstep = text_qstep_effects_part1-text_qstep_effects_part2
save(text_contemporaneous,
     text_1_lag,text_2_lag,text_3_lag,text_1_lag_structural_direct,
     text_1_step,text_2_step_general_101,text_3_step_general_0101,
     text_effects_vs_qlag,
     text_effects_vs_qstep,
     file="/Users/xiaoxuancai/Documents/GitHub/Causal_estimands/result_texts.Rdata")

# compare results from analytical solution with CI, Monte Carlo simulation with CI, and true value (not used, just for checking)
# Check results of q-lag effects
par(mfrow=c(1,1),mar=c(4, 4, 2, 0.5))
plot(0:10,colMeans(text_effects_vs_qlag),main="Analytical vs Monte Carlo vs Ture",ylab="q-lag effect",xlab="#lags",cex=1.5,cex.lab=1.5,cex.axis= 1.5, bty="n",col="black",pch=5)
points(0:3,c(rowMeans(text_contemporaneous)[600],rowMeans(text_1_lag)[600],rowMeans(text_2_lag)[600],rowMeans(text_3_lag)[600]),col="blue",cex=1.5,pch=16)
points(0:3,c(calculate.causaleffect(t=600,tx=1,y_coeffi_table=y_coeffi_table_text,c_coeffi_table=c_coeffi_table_text,CI=F),
             calculate.causaleffect(t=600,tx=c(1,0),y_coeffi_table=y_coeffi_table_text,c_coeffi_table=c_coeffi_table_text,CI=F),
             calculate.causaleffect(t=600,tx=c(1,0,0),y_coeffi_table=y_coeffi_table_text,c_coeffi_table=c_coeffi_table_text,CI=F),
             calculate.causaleffect(t=600,tx=c(1,0,0,0),y_coeffi_table=y_coeffi_table_text,c_coeffi_table=c_coeffi_table_text,CI=F)),
       col="red",pch=15,cex=1.5)
legend("bottomright",legend=c("Monte Carlo with CI","Analytical with CI","Truth"),
       pch = c(5,16,15), lty=c(1,1,1), col=c("black","blue","red"),
       bty = "n", # remove the bounder of the legend
       lwd = c(NA,NA,NA),cex=1.5)
# Check results of q-step effects
plot(0:10,colMeans(text_effects_vs_qstep),main="Analytical vs Monte Carlo vs Ture",ylab="q-step effect",xlab="#lags",cex=1.5,cex.lab=1.5,cex.axis= 1.5, bty="n",col="black",pch=5)
points(0:1,c(rowMeans(text_contemporaneous)[600],rowMeans(text_1_step)[600]),col="blue",cex=1.5,pch=16)
points(0:3,c(calculate.causaleffect(t=600,tx=1,y_coeffi_table=y_coeffi_table_text,c_coeffi_table=c_coeffi_table_text,CI=F),
             calculate.causaleffect(t=600,tx=c(1,1),y_coeffi_table=y_coeffi_table_text,c_coeffi_table=c_coeffi_table_text,CI=F),
             calculate.causaleffect(t=600,tx=c(1,1,1),y_coeffi_table=y_coeffi_table_text,c_coeffi_table=c_coeffi_table_text,CI=F),
             calculate.causaleffect(t=600,tx=c(1,1,1,1),y_coeffi_table=y_coeffi_table_text,c_coeffi_table=c_coeffi_table_text,CI=F)),
       col="red",pch=15,cex=1.5)
legend("topright",legend=c("Monte Carlo with CI","Analytical with CI","Truth"),
       pch = c(5,16,15), lty=c(1,1,1), col=c("black","blue","red"),
       bty = "n", # remove the bounder of the legend
       lwd = c(NA,NA,NA),cex=1.5)

# plots (used in the Main Paper & Appendix)
{
  load("/Users/xiaoxuancai/Documents/GitHub/Causal_estimands/result_texts.Rdata")
  # text's contemporaneous effect
  # Chosen in main text of the paper (Figure 4a)
  pdf(file = paste(address_main,"contemp_texts.pdf",sep=""),
      width = 10, # The width of the plot in inches
      height = 6) # The height of the plot in inches
  par(mar = c(2, 5, .1, .1))
  par(mfrow=c(1,1))
  text_contemporaneous_CIband=plot_simulatedCI(text_contemporaneous,probs=c(0.05,0.95),printFlag=F)
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
  pdf(file = paste(address_main,"lag1_texts.pdf",sep=""),
      width = 10, # The width of the plot in inches
      height = 6) # The height of the plot in inches
  par(mar = c(2, 5, .1, .1))
  par(mfrow=c(1,1))
  text_1_lag_CIband=plot_simulatedCI(text_1_lag,probs=c(0.05,0.95),printFlag=F)
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
  
  # text's 2-lag effect  (Appendix)
  pdf(file = paste(address_appendix,"lag2_texts.pdf",sep=""),
      width = 10, # The width of the plot in inches
      height = 6) # The height of the plot in inches
  par(mar = c(2, 5, .1, .1))
  par(mfrow=c(1,1))
  text_2_lag_CIband=plot_simulatedCI(text_2_lag,probs=c(0.05,0.95),printFlag=F)
  plot(1:708,text_2_lag_CIband$mean,type="l",ylab="2-lag effect (texts)",xlab="Date",bty="n", cex.axis=2.5,cex.lab=2.5,ylim=c(-3,0.5))
  polygon(c(1:708,rev(1:708)),c(text_2_lag_CIband$upper,rev(text_2_lag_CIband$lower)),col="grey90",border="grey")
  points(1:708,text_2_lag_CIband$mean,type="l")
  abline(h=0,lty=3,lwd=2)
  legend("bottomright",legend=c("estimate", "95% CI"),
         pch = c(NA,15), lty=c(1,NA), col=c("black","grey"),
         bty = "n", # remove the bounder of the legend
         lwd = c(1,NA), pt.bg = c(NA,"grey90"),cex=2.5)
  dev.off()
  
  # text's 3-lag effect  (Appendix)
  pdf(file = paste(address_appendix,"lag3_texts.pdf",sep=""),
      width = 10, # The width of the plot in inches
      height = 6) # The height of the plot in inches
  par(mar = c(2, 5, .1, .1))
  par(mfrow=c(1,1))
  text_3_lag_CIband=plot_simulatedCI(text_3_lag,probs=c(0.05,0.95),printFlag=F)
  plot(1:708,text_3_lag_CIband$mean,type="l",ylab="3-lag effect (texts)",xlab="Date",bty="n", cex.axis=2.5,cex.lab=2.5,ylim=c(-3,0.5))
  polygon(c(1:708,rev(1:708)),c(text_3_lag_CIband$upper,rev(text_3_lag_CIband$lower)),col="grey90",border="grey")
  points(1:708,text_3_lag_CIband$mean,type="l")
  abline(h=0,lty=3,lwd=2)
  legend("bottomright",legend=c("estimate", "95% CI"),
         pch = c(NA,15), lty=c(1,NA), col=c("black","grey"),
         bty = "n", # remove the bounder of the legend
         lwd = c(1,NA), pt.bg = c(NA,"grey90"),cex=2.5)
  dev.off()
  
  # text's 1-lag structural direct effect  (Figure 4c)
  pdf(file = paste(address_main,"lag1_controlled_direct_texts.pdf",sep=""),
      width = 10, # The width of the plot in inches
      height = 6) # The height of the plot in inches
  par(mar = c(2, 5, .1, .1))
  par(mfrow=c(1,1))
  text_1_lag_structural_direct_CIband = plot_simulatedCI(text_1_lag_structural_direct,probs=c(0.05,0.95),printFlag=F)
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
  pdf(file = paste(address_main,"step1_total_texts.pdf",sep=""),
      width = 10, # The width of the plot in inches
      height = 6) # The height of the plot in inches
  par(mar = c(2, 5, .1, .1))
  par(mfrow=c(1,1))
  text_1_step_CIband = plot_simulatedCI(text_1_step,probs=c(0.05,0.95),printFlag=F)
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
  pdf(file = paste(address_main,"step3_general_0101_texts.pdf",sep=""),
      width = 10, # The width of the plot in inches
      height = 6) # The height of the plot in inches
  par(mar = c(5, 5, .1, .1))
  par(mfrow=c(1,1))
  text_3_step_general_0101_CIband = plot_simulatedCI(text_3_step_general_0101,probs=c(0.05,0.95),printFlag=F)
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
  pdf(file = paste(address_main,"impulse_impact_texts.pdf",sep=""),
      width = 10, # The width of the plot in inches
      height = 6) # The height of the plot in inches
  par(mar = c(5, 5, .1, .1))
  par(mfrow=c(1,1))
  text_effects_vs_qlag_CIband=plot_simulatedCI(t(text_effects_vs_qlag),probs=c(0.05,0.95),printFlag=F)
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
  pdf(file = paste(address_main,"step_response_texts.pdf",sep=""),
      width = 10, # The width of the plot in inches
      height = 6) # The height of the plot in inches
  par(mar = c(5, 5, .1, .1))
  par(mfrow=c(1,1))
  text_effects_vs_qstep_CIband=plot_simulatedCI(t(text_effects_vs_qstep),probs=c(0.05,0.95),printFlag=F)
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

####################################################################################################
######         personalized intervention effect for texts (used in the Main Paper)             #####
####################################################################################################
qstep1001001 = simulate.counterfactual_path_singlet(t=600,tx=c(1,0,0,1,0,0,1,rep(0,10)),y_coeffi_table=y_coeffi_table_text,c_coeffi_table=c_coeffi_table_text,raw_data=raw_data_text,CI=F)
qstep0101010 = simulate.counterfactual_path_singlet(t=600,tx=c(0,1,0,1,0,1,0,rep(0,10)),y_coeffi_table=y_coeffi_table_text,c_coeffi_table=c_coeffi_table_text,raw_data=raw_data_text,CI=F)
qstep1110000 = simulate.counterfactual_path_singlet(t=600,tx=c(1,1,1,0,0,0,0,rep(0,10)),y_coeffi_table=y_coeffi_table_text,c_coeffi_table=c_coeffi_table_text,raw_data=raw_data_text,CI=F)
qstep0000000 = simulate.counterfactual_path_singlet(t=600,tx=c(0,0,0,0,0,0,0,rep(0,10)),y_coeffi_table=y_coeffi_table_text,c_coeffi_table=c_coeffi_table_text,raw_data=raw_data_text,CI=F)

# Figure 6b  (used in the Main Paper)
pdf(file = paste(address_main,"optimal_tx.pdf",sep=""),
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

#################################################################################################
######         test positivity assumption for texts (used in the Main Paper)                #####
#################################################################################################
positivity_calls=test.positivity(tx="keycontacts_call_totaldegree_merged",data=data[,c(1:7)],length=10)
positivity_texts=test.positivity(tx="keycontacts_text_reciprocity_degree_merged",data=data[,c(1:7)],length=10)

nrow(positivity_texts$details[[7]][rowSums(positivity_texts$details[[7]][,-8])==3,])

par(mar = c(5, 5, 2.5, 0.5))
barplot(height=positivity_calls$percentage,names=as.character(1:10),ylab="positivity % (calls)",xlab="exposure duration (days)",type="h",bty="n",cex.axis=1.5,cex.lab=2,las=1)
barplot(height=positivity_texts$percentage,names=as.character(1:10),ylab="positivity % (texts)",xlab="exposure duration (days)",type="h",bty="n",cex.axis=1.5,cex.lab=2,las=1)

