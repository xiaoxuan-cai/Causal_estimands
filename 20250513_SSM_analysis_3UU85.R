# File: 20251031_SSM_analysis_3UU85.R
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
#       6. adjust criteria of selecting variables for random walk from extremely small values to <0.01 given the short follow up period
#            thus, random walk variables change from "intercept","text_outdegree_1","TAM_phone_ave_1" to "call_outdegree_1" (variance = 2.602547e-03)
#            no clue why {"intercept","text_outdegree_1","TAM_phone_ave_1"} are selected for random walk in the previous analysis
# Results: The fitted SSM outcome plot for participant 3UU85 is included in the Appendix.
source("/Users/xiaoxuancai/Documents/GitHub/mHealth_data_processing/Summary 1_packages.R")
source("/Users/xiaoxuancai/Documents/GitHub/mHealth_data_processing/Summary 2_helper functions.R")
source("/Users/xiaoxuancai/Documents/GitHub/Causal_estimands/helper_causal estimands.R")
source("/Users/xiaoxuancai/Documents/GitHub/SSMimpute2/helper_SSM_version3.R") 
source("/Users/xiaoxuancai/Documents/GitHub/SSMimpute/helper_multipleimputation.R") # include merge point used for SSM and plot functions
load("/Users/xiaoxuancai/Desktop/Processed_data_linda/subject_3UU85/Data/processed/combined_all_3UU85.RData")
address_main = "/Users/xiaoxuancai/Dropbox/MHealthPsychSummerProject2020/Xiaoxuan_Cai/[Paper 1] Causal estimands for time series data/Graphs_main/"
address_appendix = "/Users/xiaoxuancai/Dropbox/MHealthPsychSummerProject2020/Xiaoxuan_Cai/[Paper 1] Causal estimands for time series data/Graphs_appendix/"

# -------------------------------------------------------------------- #
#             Organize and add variables needed for analysis           #
# -------------------------------------------------------------------- #
# change values for certain variables
TAM_phone_modified=data_3UU85$TAM_phone
TAM_phone_modified[which(TAM_phone_modified==1)]=0.9999
TAM_phone_modified[which(TAM_phone_modified==0)]=0.0001
table(data_3UU85$keycontacts_call_totaldegree) # merge level 2&3 to level 1 -> turn into a binary variable
keycontacts_call_totaldegree_binary=data_3UU85$keycontacts_call_totaldegree;keycontacts_call_totaldegree_binary[keycontacts_call_totaldegree_binary>1]=1;table(keycontacts_call_totaldegree_binary)
table(data_3UU85$keycontacts_text_reciprocity_degree) # merge level 2 to level 1 -> turn into a binary variable
keycontacts_text_reciprocity_degree_binary=data_3UU85$keycontacts_text_reciprocity_degree;keycontacts_text_reciprocity_degree_binary[keycontacts_text_reciprocity_degree_binary>1]=1;table(keycontacts_text_reciprocity_degree_binary)
data=data.frame(Date=data_3UU85$Date,negative_total=data_3UU85$negative_total,
                keycontacts_call_totaldegree_binary,
                keycontacts_text_reciprocity_degree_binary,
                logit_TAM_phone=log(TAM_phone_modified/(1-TAM_phone_modified))) 
param=list(list(lagged_param=list(variables=colnames(data)[-1],param=rep(1,length(colnames(data)[-1])))))
data=add_variables_procedures(data,param);
data=data.frame(intercept=1,data);colnames(data)

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
  plot.KFS(ssm_y_init_out,range=30:144)
  ssm_y_init_out$model$Q; ssm_y_init_out$model$H
  # check: intercept(0.0002032787), negative_total_1 (1.608733e-07), 
  #        keycontacts_call_totaldegree_binary (0.06969114), keycontacts_call_totaldegree_binary_1 (0.005263435),
  #        keycontacts_text_reciprocity_degree_binary (2.220929e-05), keycontacts_text_reciprocity_degree_binary_1 (0.05931756)
  #        logit_TAM_phone_1 (0.0004104951)
  #       6. adjust criteria of selecting variables for random walk from extremely small values to <0.01 given the short follow up period
  #          thus, random walk variables are "keycontacts_call_totaldegree_binary ", "keycontacts_text_reciprocity_degree_binary_1" 
  
  
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
  ss_param_y=list(inits=c(log(5)),m0=ssm_y_init_out$att[144,],C0=diag(rep(10^3),7),
                     AR1_coeffi=NULL,rw_coeffi=NULL,v_cp_param=NULL,
                     w_cp_param=list(list(variable="keycontacts_call_totaldegree_binary",segments=2,changepoints=c(80),fixed_cpts=F),
                                     list(variable="keycontacts_text_reciprocity_degree_binary_1",segments=2,changepoints=c(80),fixed_cpts=F)),
                     max_iteration=100)
  par(mfrow=c(1,2))
  ssm_y=run.SSM(data_ss=data_ss_ignore_y,formula=formula,ss_param_temp=ss_param_y,max_iteration=100,
                cpt_learning_param=list(cpt_method="mean",burnin=1/10,mergeband=20,convergence_cri=10),
                cpt_merge_option="separate",dlm_cpt_learning_option="filter",
                bandwidth=5,cpt_V =1,printFlag=T)
  plot.dlm(ssm_y$out_filter,benchmark = rep(0,7),option="filtered_state")
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
  # (Intercept)                                              3.955**     1.312   (0.575,7.334)   (1.383,6.526)   (1.797,6.113)
  # negative_total_1                                        -0.037       0.140  (-0.399,0.324)  (-0.313,0.238)  (-0.268,0.194)
  # keycontacts_call_totaldegree_binary(period1)             2.792**     0.963   (0.311,5.273)    (0.904,4.68)   (1.208,4.377)
  # keycontacts_call_totaldegree_binary(period2)            -0.564       1.090  (-3.371,2.244)    (-2.7,1.572)  (-2.356,1.229)
  # keycontacts_call_totaldegree_binary_1                    2.246**     0.804   (0.175,4.316)    (0.67,3.821)   (0.924,3.568)
  # keycontacts_text_reciprocity_degree_binary               0.580       0.811    (-1.51,2.67)    (-1.01,2.17)  (-0.755,1.914)
  # keycontacts_text_reciprocity_degree_binary_1(period1)   -2.933**     1.068 (-5.686,-0.181) (-5.028,-0.839) (-4.691,-1.176)
  # keycontacts_text_reciprocity_degree_binary_1(period2)   -0.767       0.967  (-3.258,1.724)  (-2.662,1.129)  (-2.358,0.824)
  # logit_TAM_phone_1                                       -0.020       0.222  (-0.592,0.552)  (-0.455,0.415)  (-0.385,0.345)
}

pdf(file = paste(address_appendix,"Rplot_both_3UU85.pdf",sep=""),  
    width = 15, height = 4)
par(mfrow=c(1,2),mar =c(2.5, 8, .5, 1.5))
# outcome
n=length(data$negative_total)
plot(1:n,data$negative_total,type="l",bty="n",xlab="Time (days)",ylab="Negative mood",cex.lab=1.5,cex.axis=2)
mtext("Subject 2 (Bipolar I)",side = 2, adj = 0.5, padj=-5,cex=1.8)
text(x=2,y=16,label="a)",cex=2)
# time-invariant coefficient for outdegree of texts yesterday, beta_22
plot(1:n,ssm_y$out_filter$m[,5],type="l",ylab=expression(paste(beta["21,t"]," (text connectivity)")),xlab="Time (days)" ,bty="n",cex.lab=1.5,cex.axis=2, ylim=c(-4,7))
sd=unlist(lapply(1:n,function(i){sqrt(diag(dlmSvd2var(ssm_y$out_filter$U.C[[i]],ssm_y$out_filter$D.C[i,]))[5])}))
polygon(c(1:n,rev(1:n)),c(ssm_y$out_filter$m[,5]+qnorm(0.975)*sd,rev(ssm_y$out_filter$m[,5]-qnorm(0.975)*sd)),col="grey90",border="grey")
points(1:n,ssm_y$out_filter$m[,5],type="l")
abline(h=0,col="black",lty=2,lwd=4)
legend("bottomright",legend=c("estimate", "95% CI"),
       pch = c(NA,15), lty=c(1,NA), col=c("black","grey"),
       bty = "n", # remove the bounder of the legend
       lwd = c(1,NA), pt.bg = c(NA,"grey90"),cex=1.5)
text(x=2,y=6.9,label="b)",cex=2)
dev.off()
