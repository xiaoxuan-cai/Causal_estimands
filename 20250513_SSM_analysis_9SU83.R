# File: 20251031_SSM_analysis_9SU83.R
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
#       6. random walk variables change from {"intercept","TAM_phone_ave_1"} to "intercept"
#            period-stable variables change from "text_outdegree" to {"text_outdegree","text_outdegree_1","TAM_phone_ave"}
# Results: The fitted SSM outcome plot for participant 9SU83 is included in the Appendix.


# File: 20250513_MCEMSSM_analysis_9SU83.R
# Date: 2025.05.13
# Author: Xiaoxuan Cai (contact nmcaixiaoxuan@gmail.com for questions)
# Goal: 1) produce descriptive graphs for raw data
#       2) apply SSMimpute, MCEM_SSM and compare their performances on real data
# Note: compared to "/Users/xiaoxuancai/Documents/GitHub/SSMimpute/20221114_SSMimpute_analysis_9SU83.R"
#       SSM and SSMimpute versions have been updated to version 3, providing identical results to the previous versions.
#       analysis results have changed due to the 
#         1. switch from using "TAM_phone_ave" to "TAM_phone_ave_1"
#         2. random walk variables change from {"intercept","TAM_phone_ave_1"} to "intercept"
#            period-stable variables change from "text_outdegree" to {"text_outdegree","text_outdegree_1","TAM_phone_ave"}
source("/Users/xiaoxuancai/Documents/GitHub/mHealth_data_processing/Summary 1_packages.R")
source("/Users/xiaoxuancai/Documents/GitHub/mHealth_data_processing/Summary 2_helper functions.R")
source("/Users/xiaoxuancai/Documents/GitHub/Causal_estimands/helper_causal estimands.R")
source("/Users/xiaoxuancai/Documents/GitHub/SSMimpute2/helper_SSM_version3.R") 
source("/Users/xiaoxuancai/Documents/GitHub/SSMimpute/helper_multipleimputation.R") # include merge point used for SSM and plot functions
load("/Users/xiaoxuancai/Desktop/Processed_data_linda/subject_9SU83/Data/processed/combined_all_9SU83.RData")
address_main = "/Users/xiaoxuancai/Dropbox/MHealthPsychSummerProject2020/Xiaoxuan_Cai/[Paper 1] Causal estimands for time series data/Graphs_main/"
address_appendix = "/Users/xiaoxuancai/Dropbox/MHealthPsychSummerProject2020/Xiaoxuan_Cai/[Paper 1] Causal estimands for time series data/Graphs_appendix/"
data_9SU83=data_9SU83[data_9SU83$Date<"2016-11-14",] # information after 2016-11-14 is too sparse and no good to use

# -------------------------------------------------------------------- #
#             Organize and add variables needed for analysis           #
# -------------------------------------------------------------------- #
# change values for certain variables
TAM_phone_modified=data_9SU83$TAM_phone
TAM_phone_modified[which(TAM_phone_modified==1)]=0.9999
TAM_phone_modified[which(TAM_phone_modified==0)]=0.0001
table(data_9SU83$keycontacts_call_totaldegree) # merge level 2&3&4&5 to level 1 -> turn into a binary variable
keycontacts_call_totaldegree_binary=data_9SU83$keycontacts_call_totaldegree;keycontacts_call_totaldegree_binary[keycontacts_call_totaldegree_binary >1]=1;table(keycontacts_call_totaldegree_binary)
table(data_9SU83$keycontacts_text_reciprocity_degree) # merge level 2&3&4&5 to level 1 -> turn into a binary variable
keycontacts_text_reciprocity_degree_binary=data_9SU83$keycontacts_text_reciprocity_degree;keycontacts_text_reciprocity_degree_binary[keycontacts_text_reciprocity_degree_binary>1]=1;table(keycontacts_text_reciprocity_degree_binary)
data=data.frame(Date=data_9SU83$Date,negative_total=data_9SU83$negative_total,
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
  plot.KFS(ssm_y_init_out)
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
  ss_param_y=list(inits=c(log(0.01),log(3)),m0=ssm_y_init_out$att[163,],C0=diag(rep(10^3),7),
                  AR1_coeffi=NULL,rw_coeffi=c("negative_total_1"),v_cp_param=NULL,
                  w_cp_param=list(list(variable="keycontacts_call_totaldegree_binary",segments=2,changepoints=c(100),fixed_cpts=F),
                                  list(variable="keycontacts_call_totaldegree_binary_1",segments=2,changepoints=c(50),fixed_cpts=F)),
                  max_iteration=100)
  ss_param_y=list(inits=c(log(3)),m0=ssm_y_init_out$att[163,],C0=diag(rep(10^3),7),
                    AR1_coeffi=NULL,rw_coeffi=NULL,v_cp_param=NULL,w_cp_param=NULL,
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
  # (Intercept)                                            5.457     1.040  (2.778,8.136)  (3.418,7.495)  (3.746,7.167)
  # (Intercept)                                      -0.168     1.002 (-2.749,2.414) (-2.132,1.797) (-1.816,1.481)
  # negative_total_1                                     NA        NA        (NA,NA)        (NA,NA)        (NA,NA)
  # keycontacts_call_totaldegree_binary(period1)      0.010     0.524  (-1.34,1.359) (-1.017,1.036) (-0.852,0.871)
  # keycontacts_call_totaldegree_binary(period2)     -0.302     0.637  (-1.944,1.34) (-1.551,0.948)  (-1.35,0.747)
  # keycontacts_call_totaldegree_binary_1(period1)    0.560     0.635 (-1.077,2.196) (-0.685,1.804) (-0.485,1.604)
  # keycontacts_call_totaldegree_binary_1(period2)    3.171     0.500  (1.883,4.459)  (2.191,4.151)  (2.349,3.994)
  # keycontacts_text_reciprocity_degree_binary       -0.364     0.598 (-1.904,1.175) (-1.536,0.807) (-1.347,0.619)
  # keycontacts_text_reciprocity_degree_binary_1      0.258     0.712 (-1.575,2.091) (-1.136,1.653) (-0.912,1.429)
  # logit_TAM_phone_1                                -0.281     0.238 (-0.894,0.333) (-0.747,0.186) (-0.672,0.111)
}
pdf(file = paste(address_appendix,"Rplot_both_9SU83.pdf",sep=""),  
    width = 15, height = 4)
par(mfrow=c(1,2),mar =c(2.5, 8, .5, 1.5))
# outcome
n=length(data$negative_total)
plot(1:n,data$negative_total,type="l",bty="n",xlab="Time (days)",ylab="Negative mood",cex.lab=1.5,cex.axis=2)
mtext("Subject 4 (Bipolar I)",side = 2, adj = 0.5, padj=-5,cex=1.8)
text(x=10,y=7,label="e)",cex=2)
# time-invariant coefficient for outdegree of texts yesterday, beta_22
plot(1:n,ssm_y$out_filter$m[,5],type="l",ylab=expression(paste(beta["21,t"]," (text connectivity)")),xlab="Time (days)" ,bty="n",cex.lab=1.5,cex.axis=2, ylim=c(-4,3))
sd=unlist(lapply(1:n,function(i){sqrt(diag(dlmSvd2var(ssm_y$out_filter$U.C[[i]],ssm_y$out_filter$D.C[i,]))[5])}))
polygon(c(1:n,rev(1:n)),c(ssm_y$out_filter$m[,5]+qnorm(0.975)*sd,rev(ssm_y$out_filter$m[,5]-qnorm(0.975)*sd)),col="grey90",border="grey")
points(1:n,ssm_y$out_filter$m[,5],type="l")
abline(h=0,col="black",lty=2,lwd=4)
legend("bottomright",legend=c("estimate", "95% CI"),
       pch = c(NA,15), lty=c(1,NA), col=c("black","grey"),
       bty = "n", # remove the bounder of the legend
       lwd = c(1,NA), pt.bg = c(NA,"grey90"),cex=1.5)
text(x=10,y=2.9,label="f)",cex=2)
dev.off()

