# Set the working dir such that all files in the zip file are available
setwd('XXX/submission_code')
library(igraph)
source('SCREEN_code_for_submission.R')
source('twogroups_methods_for_submission.R')
source('simulation_functions.R')

# Load the cancer DEG data
cancer_deg_pvals = as.matrix(read.delim('cancer_deg_pvals.txt'))
# Load the HLA data
hla_pvals = as.matrix(read.delim('hla_pvals.txt'))

############################ Analyze the cancer DEG dataset #############################
pvals = cancer_deg_pvals
k_range = 2:25
lfdr_thr = 0.2
col = "SCREEN 20"
nH = 10000

# The results of each methods can be extracted from rep_analysis_results$all_lfdr_rep_scores
rep_analysis_results = analyze_pval_matrix_slim(pvals,ks=k_range,nH=nH,lfdr_method='znormix',use_power=T,standardNull=T)
# Run SCREEN, and SCREEN-ind
screen_results = SCREEN(pvals,ks=k_range,nH=nH,lfdr_method='znormix',use_power=T,standardNull=T)
screen_ind_results = SCREEN_ind(pvals,ks=k_range,lfdr_method='znormix',use_power=T,standardNull=T)

# Compare to the paper's results with a much lower nH value
load('geo_cancer_data_genepvals_analysis_results_use_power.RData')
est = znormix_est_genepvals_theoretical_null
colnames(est$all_lfdr_rep_scores)
paper_results = est$all_lfdr_rep_scores[,"Clust 20"]
table(paper_results<0.2,screen_results[,col]<0.2)
in_original_publication = paper_results<0.2
new_results = screen_results[,col]
l = list()
l[["Significant in lower nH"]] = new_results[in_original_publication]
l[["Not significant in lower nH"]] = new_results[!in_original_publication]
library(vioplot)
vioplot(l[[1]],l[[2]],names=names(l))

############################ Analyze the HLA dataset #############################
pvals = hla_pvals
k_range = 2:5
lfdr_thr = 0.2
col = "SCREEN 4"
nH = 10000

# The results of each methods can be extracted from rep_analysis_results$all_lfdr_rep_scores
rep_analysis_results = analyze_pval_matrix_slim(pvals,ks=k_range,nH=nH,lfdr_method='znormix',use_power=T,standardNull=T)
# Run SCREEN, and SCREEN-ind
screen_results = SCREEN(pvals,ks=k_range,nH=nH,lfdr_method='znormix',use_power=T,standardNull=T)
screen_ind_results = SCREEN_ind(pvals,ks=k_range,lfdr_method='znormix',use_power=T,standardNull=T)

############################ EM comment #############################
# As explained in the discussion, one can improve the algorithm
# by suggesting better initialization of the EM or a different algorithm to estimate the prior.
# If you wish to explore this direction, the only different in implementation is to change the function below.
# For example the, the method below set a random starting point instead of uniform as in the 
# original implementation of SCREEN.
estimate_prior<-function(f_1_mat,f_0_mat,eps = 1e-50,
			convergenceEps = 1e-6,H=NULL,p0 = NULL,maxiter=1200){
	# m is the number of datasets
	m = ncol(f_1_mat)
	if (is.null(H)){
		H = t(sapply(0:(2^m-1),function(x,m){as.integer(intToBits(x))[1:m]},m=m))
		if (m==1){H=t(H)}
	}
	if (is.null(p0)){
		unif_sample = runif(nrow(H))
		start_point_h_p = unif_sample/sum(unif_sample)
	}
	else{
		start_point_h_p = p0
	}
	start_point_h_p_H	 = H
	z_given_h_mat = get_z_given_h(start_point_h_p_H,f_0_mat,f_1_mat,eps)
	p_0 = start_point_h_p
	p_t = runEM(z_given_h_mat,p_0,convergenceEps,maxiter)
	return (cbind(start_point_h_p_H,p_t))
}
# For example, after changing the function above you can rerun the algorithm:
alt_em_screen_results = list()
for(j in 1:reps){
	alt_em_screen_results[[j]] = SCREEN(pvals,ks=k_range,nH=1024,lfdr_method='znormix',use_power=T,standardNull=T)
}
length(alt_screen_results)
alt_em_results = sapply(alt_em_screen_results,function(x,y)x[,y],y=col)

############################ How to run simulations #############################
# Example: simulate data of 4 clusters with r=0.8 and run all methods
ngenes = 5000
nsig_within = 150
beta1_params = c(1,1)
beta2_params = c(100,1)
beta3_params = beta2_params[2:1]
m=10
sigmaMat = matrix(0.8,m,m)
diag(sigmaMat)=1
simulation_obj = simulations_case2(m,ngenes,nsig_within,
		sigmaMat,use_beta=T,beta1_params,beta2_params,beta3_params,num_clusters=4)
pmat = simulation_obj$pmat
nonull_inds = simulation_obj$nonull_inds
normix_results = analyze_pval_matrix_slim(pmat,ks=2:5,nH=512,lfdr_method='znormix')
# get the Jaccard and FDP scores
scores_normix = get_pred_scores2(normix_results,pmat,nonull_inds,thr=0.2)

