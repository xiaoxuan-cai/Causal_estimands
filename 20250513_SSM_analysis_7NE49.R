# File: 20251031_SSM_analysis_7NE49.R
# Date: 2025.10.31
# Author: Xiaoxuan Cai (contact nmcaixiaoxuan@gmail.com for questions)
# This file is created for the real data analysis of the paper "Causal estimands and identification of time-varying effects in non-stationary time series from N-of-1 mobile device data"
# Goal: Fit a time-varying SSM for the outcome variable only (not the covariates).
#       The purpose is to illustrate the heterogeneity in participant 7NE79's outcome SSM fitting.
#       No causal estimands are computed in this analysis; therefore, SSM modeling of covariates is not required.
#       This script is adapted from "20230810_5BT65_causalestimands.R" and includes only the sections up to the SSM modeling of the outcome variable.
# Note: 1. the methodology of the paper only applies to binary exposure and # obs on level 3 is too few for proper inference
#          therefore, we merge <keycontacts_call_totaldegree>(range 0~3) and <keycontacts_text_reciprocity_degree> (range 0~3)
#          to binary variables
#          we consider the following binary variable as exposure of interest:
#          <keycontacts_call_totaldegree> (binary) and <keycontacts_text_reciprocity_degree> (binary)
#       2. with no proper strategies for missing data, we use ignore case analysis for now with additional work to be developed to deal with missing data
#       3. we also update SSM fitting function to their new versions
#          old version "GitHub/SSMimpute/helper_multipleimputation.R" -> new version
#       4. keycontacts_text_totaldegree or keycontacts_text_reciprocity_degree? reciprocity shows real conversation than totaldegree
#       5. take a logit transformation on TAM_phone, to make it normally distributed from -infty to +infty
# Results: The fitted SSM outcome plot for participant 7NE49 is included in the Appendix.
source("/Users/xiaoxuancai/Documents/GitHub/mHealth_data_processing/Summary 1_packages.R")
source("/Users/xiaoxuancai/Documents/GitHub/mHealth_data_processing/Summary 2_helper functions.R")
source("/Users/xiaoxuancai/Documents/GitHub/Causal_estimands/helper_causal estimands.R")
source("/Users/xiaoxuancai/Documents/GitHub/SSMimpute2/helper_SSM_version3.R") 
source("/Users/xiaoxuancai/Documents/GitHub/SSMimpute/helper_multipleimputation.R") # include merge point used for SSM and plot functions
# Note on 7NE49
# MotoG3: 2016-02-12 19_00_00.csv -> < 2016-10-27 18_00_00.csv
# LGLS770: 2016-10-27 18_00_00.csv -> < 2018-07-25 18_00_00.csv
# LGLS755: 2018-07-25 18_00_00.csv -> end
load("/Users/xiaoxuancai/Desktop/Processed_data_linda/subject_7NE49/Data/processed/combined_all_7NE49.RData")
address_main = "/Users/xiaoxuancai/Dropbox/MHealthPsychSummerProject2020/Xiaoxuan_Cai/[Paper 1] Causal estimands for time series data/Graphs_main/"
address_appendix = "/Users/xiaoxuancai/Dropbox/MHealthPsychSummerProject2020/Xiaoxuan_Cai/[Paper 1] Causal estimands for time series data/Graphs_appendix/"

# -------------------------------------------------------------------- #
#             Organize and add variables needed for analysis           #
# -------------------------------------------------------------------- #
# change values for certain variables
TAM_phone_modified=data_7NE49$TAM_phone
TAM_phone_modified[which(TAM_phone_modified==1)]=0.9999
TAM_phone_modified[which(TAM_phone_modified==0)]=0.0001
table(data_7NE49$keycontacts_call_totaldegree) # merge level 2&3&4 to level 1 -> turn into a binary variable
keycontacts_call_totaldegree_binary=data_7NE49$keycontacts_call_totaldegree;keycontacts_call_totaldegree_binary[keycontacts_call_totaldegree_binary %in% c(2,3,4)]=1;table(keycontacts_call_totaldegree_binary)
table(data_7NE49$keycontacts_text_reciprocity_degree) # merge level 2&3 to level 1 -> turn into a binary variable
keycontacts_text_reciprocity_degree_binary=data_7NE49$keycontacts_text_reciprocity_degree;keycontacts_text_reciprocity_degree_binary[keycontacts_text_reciprocity_degree_binary>1]=1;table(keycontacts_text_reciprocity_degree_binary)
data=data.frame(Date=data_7NE49$Date,negative_total=data_7NE49$negative_total,
                keycontacts_call_totaldegree_binary,
                keycontacts_text_reciprocity_degree_binary,
                logit_TAM_phone=log(TAM_phone_modified/(1-TAM_phone_modified))) 
param=list(list(lagged_param=list(variables=colnames(data)[-1],param=rep(1,length(colnames(data)[-1])))))
data=add_variables_procedures(data,param);
data=data.frame(intercept=1,data);colnames(data)
{
  # carefully review the raw data and then there is a clear change of pattern 
  # -> in preparation of potential change points
  plot(data$negative_total,type="l") # seems to change around t=1000 
  plot(data_7NE49$keycontacts_call_totaldegree) # seems to change around t=600 and t=900 
  plot(data_7NE49$keycontacts_text_reciprocity_degree) # super clear pattern change between t=250 and t=900, zero text messages at all 
  plot(data_7NE49$TAM_phone) # super clear pattern change between t=600 and t=900
  plot(data$logit_TAM_phone_1) # super clear pattern change between t=600 and t=900
}

###############################################################################################################
######          SSM ignore modeling of covariates  (used in Appendix)            #####
###############################################################################################################
# This paper cannot handle missing data for both covariates and outcome
#   thus, we use SSM ignore temporarily
# outcome regression: negative_total~ intercept + negative_total_1 
#                                     + keycontacts_call_totaldegree_binary + keycontacts_call_totaldegree_binary_1
#                                     + keycontacts_text_reciprocity_degree_binary + keycontacts_text_reciprocity_degree_binary
#                                     + logit_TAM_phone_1
# ssm_y model (used in the Appendix)
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
  plot.KFS(ssm_y_init_out,range=10:1460)
  ssm_y_init_out$model$Q; ssm_y_init_out$model$H
  # check: intercept(1.090972e-06), negative_total_1 (6.2912e-05), 
  #        keycontacts_call_totaldegree_binary (0.004271464), keycontacts_call_totaldegree_binary_1 (9.194463e-08),
  #        keycontacts_text_reciprocity_degree_binary (0.001547042), keycontacts_text_reciprocity_degree_binary_1 (1.778678e-07)
  #        logit_TAM_phone_1 (3.060366e-08)
  
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
  # used on 2025.10.31
  # ss_param_y=list(inits=c(log(0.01),log(5)),m0=ssm_y_init_out$att[1460,],C0=diag(rep(10^3),7),
  #                 AR1_coeffi=NULL,rw_coeffi=c("intercept"),v_cp_param=NULL,
  #                 w_cp_param=list(list(variable="keycontacts_text_reciprocity_degree_binary",segments=2,changepoints=c(600),fixed_cpts=F),
  #                                 #list(variable="keycontacts_text_reciprocity_degree_binary_1",segments=2,changepoints=c(900),fixed_cpts=F),
  #                                 list(variable="logit_TAM_phone_1",segments=3,changepoints=c(400,600),fixed_cpts=F)),
  #                 max_iteration=100)
  # updated on 2025.11.05, after careful review of the raw data and the previously fitted ssm
  ss_param_y=list(inits=c(log(0.01),log(5)),m0=ssm_y_init_out$att[1460,],C0=diag(rep(10^3),7),
                  AR1_coeffi=NULL,rw_coeffi=c("intercept"),v_cp_param=NULL,
                  w_cp_param=list(list(variable="keycontacts_call_totaldegree_binary",segments=3,changepoints=c(600,1000),fixed_cpts=F),
                                  list(variable="keycontacts_text_reciprocity_degree_binary",segments=2,changepoints=c(1000),fixed_cpts=F),
                                  list(variable="logit_TAM_phone_1",segments=3,changepoints=c(400,600),fixed_cpts=F)),
                  max_iteration=100)
  ssm_y=run.SSM(data_ss=data_ss_ignore_y,formula=formula,ss_param_temp=ss_param_y,max_iteration=100,
                cpt_learning_param=list(cpt_method="mean",burnin=1/10,mergeband=20,convergence_cri=10),
                cpt_merge_option="separate",dlm_cpt_learning_option="filter",
                bandwidth=5,cpt_V =1,printFlag=T)
  ssm_y$estimated_cpts
  # $keycontacts_call_totaldegree_binary
  # [1]  579 1005
  # $keycontacts_text_reciprocity_degree_binary
  # [1] 1005
  # $logit_TAM_phone_1
  # [1] 423 629
  plot.dlm(ssm_y$out_filter,benchmark = rep(0,7),option="filtered_state",range=30:1460)
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
  #                                                        Estimate  Std.Error  99% CI          95% CI          90% CI
  # (Intercept)                                            5.866**     1.001  (3.288,8.445)  (3.905,7.828)   (4.22,7.513)
  # negative_total_1                                       0.220**     0.037  (0.125,0.315)  (0.148,0.292)    (0.16,0.28)
  # keycontacts_call_totaldegree_binary(period1)           0.129     0.673 (-1.605,1.863) (-1.191,1.449) (-0.979,1.236)
  # keycontacts_call_totaldegree_binary(period2)           0.981     0.617  (-0.608,2.57)  (-0.228,2.19) (-0.034,1.996)
  # keycontacts_call_totaldegree_binary(period3)          -0.676     0.503 (-1.972,0.619) (-1.662,0.309) (-1.504,0.151)
  # keycontacts_call_totaldegree_binary_1                  0.486     0.342 (-0.394,1.367) (-0.184,1.156) (-0.076,1.049)
  # keycontacts_text_reciprocity_degree_binary(period1)   -0.074     0.752 (-2.012,1.864) (-1.549,1.401) (-1.312,1.164)
  # keycontacts_text_reciprocity_degree_binary(period2)    0.772     0.545 (-0.632,2.175)  (-0.296,1.84) (-0.125,1.668)
  # keycontacts_text_reciprocity_degree_binary_1          -0.021     0.432 (-1.134,1.092) (-0.868,0.826)  (-0.732,0.69)
  # logit_TAM_phone_1(period1)                             0.587     0.704 (-1.227,2.401) (-0.793,1.967) (-0.571,1.745)
  # logit_TAM_phone_1(period2)                            -0.647     0.606 (-2.207,0.913)  (-1.834,0.54) (-1.644,0.349)
  # logit_TAM_phone_1(period3)                             0.191     0.222 (-0.381,0.763) (-0.244,0.626) (-0.174,0.556)
}

pdf(file = paste(address_appendix,"Rplot_both_7NE49.pdf",sep=""),  
    width = 15, height = 4)
par(mfrow=c(1,2),mar =c(2.5, 8, .5, 1.5))
# outcome
n=length(data$negative_total)
plot(1:n,data$negative_total,type="l",bty="n",xlab="Time (days)",ylab="Negative mood",cex.lab=1.5,cex.axis=2)
mtext("Subject 3 (Schizoaffective)",side = 2, adj = 0.5, padj=-5,cex=1.8)
text(x=10,y=26,label="c)",cex=2)
# time-invariant coefficient for outdegree of texts yesterday, beta_22
a=cpt.mean(ssm_y$out_filter$m[,5],penalty="Manual",Q=1,method="BinSeg")
plot(1:n,ssm_y$out_filter$m[,5],type="l",ylab=expression(paste(beta["21,t"]," (text connectivity)")),xlab="Time (days)" ,bty="n",cex.lab=1.5,cex.axis=2, ylim=c(-5,5))
sd=unlist(lapply(1:n,function(i){sqrt(diag(dlmSvd2var(ssm_y$out_filter$U.C[[i]],ssm_y$out_filter$D.C[i,]))[5])}))
polygon(c(1:n,rev(1:n)),c(ssm_y$out_filter$m[,5]+qnorm(0.975)*sd,rev(ssm_y$out_filter$m[,5]-qnorm(0.975)*sd)),col="grey90",border="grey")
points(1:n,ssm_y$out_filter$m[,5],type="l")
lines(c(1:n)[c(1,cpts(a)[1])],rep(param.est(a)$mean[1],2),col="red")
lines(c(1:n)[c(cpts(a)[1]+1,n)],rep(param.est(a)$mean[2],2),col="red")
abline(h=0,col="black",lty=2,lwd=4)
legend("bottomright",legend=c("estimate", "95% CI"),
       pch = c(NA,15), lty=c(1,NA), col=c("black","grey"),
       bty = "n", # remove the bounder of the legend
       lwd = c(1,NA), pt.bg = c(NA,"grey90"),cex=1.5)
text(x=20,y=4.9,label="d)",cex=2)
dev.off()
