# This is the code implementing the replicability analysis algorithms in the paper
# For Screen: igraph is needed
# To run the method: source("twogroups_methods_for_submission.R") 
#	- these are implementations of locfdr wrappers and the normix method
# For executables of SCREEN and SCREEN-ind, or an executable of all six tested methods see the code 
# at the bottom.

############## The EM algorithm, similar to Yekutieli and Heller #############
# z_given_h_mat - fixed, n x 2^m matrix, rows are genes
# pi_t_vec - latest estimation of the prior - a vector of size 2^m
# multiply each row of z_given_h by pi_t_vec
# normalize each row 
get_z_given_h<-function(H,F_0,F_1,eps=1e-50){
	p = exp(log(F_0+eps)%*%t(1-H) + log(F_1+eps)%*%t(H))
	return (p)
}
# For QA only
#get_z_given_h_slow<-function(H,F_0,F_1,eps=1e-50){
#	p = matrix(nrow=nrow(F_0),ncol=nrow(H))
#	for (i in 1:nrow(F_0)){
#		a1 = F_0[i,];a2 = F_1[i,]
#		for (j in 1:nrow(H)){
#			h = H[j,]
#			p[i,j] = prod((a1^(1-h))*(a2^h))
#		}
#	}
#	return (p)
#}
E_step<-function(z_given_h_mat,pi_t_vec){
	m = t(t(z_given_h_mat)*pi_t_vec)
	m = m / rowSums(m)
	return(m)
}
M_step <- function(h_given_z_mat){return (colMeans(h_given_z_mat))}
runEM<-function(z_given_h_mat,p_0,convergenceEps = 1e-6,maxiter=1200){
	p_t = p_0
	num_iters = 0
	while(num_iters <= maxiter){
		if (num_iters > 1 && num_iters %% 100 == 0){
			print(paste("completed EM iterations: ",num_iters))
		}
		num_iters = num_iters + 1
		e_step = E_step(z_given_h_mat,p_t)
		new_p = M_step(e_step)
		if (max(abs(new_p - p_t)) < convergenceEps){break}
		p_t = new_p
	}
	print(paste("**** completed all EM iterations: ",num_iters))
	return (p_t)
}
# A method to get the pi(h) values under independence assumption
# This can be used as a starting point for the EM
get_indep_assumption_pi_hat<-function(H,p0){
	pi = apply(H,1,function(x,p0){prod(pmin(p0^(1-x),(1-p0)^x))}, p0=p0)
	return (pi)
}
# If H is null create all configs
# If p0 is null, start with a uniform distribution
estimate_prior<-function(f_1_mat,f_0_mat,eps = 1e-50,
			convergenceEps = 1e-6,H=NULL,p0 = NULL,maxiter=1200){
	# m is the number of datasets
	m = ncol(f_1_mat)
	if (is.null(H)){
		H = t(sapply(0:(2^m-1),function(x,m){as.integer(intToBits(x))[1:m]},m=m))
		if (m==1){H=t(H)}
	}
	if (is.null(p0)){
		start_point_h_p = rep(1/nrow(H),nrow(H))
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
#####################################################################################
# A simple method for obtaining empirical configurations
# using marginal lfdr estimation
get_empirical_configs<-function(lfdrs,H=NULL,thrs=c(0.1,0.2),minCount=2){
	#### get the empirical dist for the given threshold
	lfdrs_h_p = c()
	for (thr in thrs){
		currlfdrs_H = lfdrs <=  thr; mode(currlfdrs_H) = 'numeric'
		lfdrs_h_vec = apply(currlfdrs_H,1,paste,collapse=';',sep='')
		tab = table(lfdrs_h_vec)
		tab = tab[tab>=minCount]
		lfdrs_h_p = union(lfdrs_h_p,names(tab))
	}
	
	#### add configs to the given H
	in_empirical = rep(T,length(lfdrs_h_p))
	names(in_empirical) = lfdrs_h_p
	if (length(H) > 0){
		H_vec = apply(H,1,paste,collapse=';',sep='')
		to_add = setdiff(H_vec,lfdrs_h_p)
		if (length(to_add)>0){
			lfdrs_h_p = union(lfdrs_h_p,to_add)
			in_empirical[to_add] = F
		}
	}
	in_empirical = in_empirical[lfdrs_h_p]
	####
	start_point_h_p_H = t(sapply(lfdrs_h_p,function(x)strsplit(x,split=';')[[1]]))
	mode(start_point_h_p_H) = 'numeric'
	return (cbind(start_point_h_p_H,in_empirical))
}
####################################################################################
######################## This is our new EM-like algorithm #########################
# Here we iteratively add a new study
# We keep:  (1) the top nH configs from the last iteration
#		(2) configs from the empirical dist
Space_constrained_estimation<-function(lfdrs,f_1_mat,f_0_mat,pi0Vec=NULL,
			nH=512,use_empirical=F,thrs=seq(0.1,0.2,0.1),...){
	nstart = round(log(nH,2))-1
	ind = min(nstart,ncol(f_1_mat))
	ests = estimate_prior(f_1_mat[,1:ind],f_0_mat[,1:ind],p0=pi0Vec,...)
	if (ind==ncol(f_1_mat)){return(list(H=ests,eps=0,sumEstimated=1))}
	currH = ests[,1:ind];sumEstimated=1;eps=0;currPi = ests[,ind+1]
	for (ind in (nstart+1):ncol(lfdrs)){
		newH = rbind(cbind(currH,rep(0,nrow(currH))),cbind(currH,rep(1,nrow(currH))))
		newP0 = c(currPi/2,currPi/2)
		if (use_empirical){
			emp = get_empirical_configs(lfdrs[,1:ind],newH,thrs)
			currNewH = emp[,1:ind]; in_empirical = emp[,ind+1]
			ests = estimate_prior(f_1_mat[,1:ind],f_0_mat[,1:ind],H=currNewH,p0=newP0,...)
		}
		else{
			ests = estimate_prior(f_1_mat[,1:ind],f_0_mat[,1:ind],H=newH,p0=newP0,...)
			in_empirical = rep(F,nrow(newH))
		}
		currH = ests[,1:ind];currPi = ests[,ind+1]
		currPi = currPi*sumEstimated
		currEps = sort(currPi,decreasing=T)[min(nH,length(currPi))]
		boolVec = currPi > currEps | in_empirical
		if (sum(boolVec) > nH && !use_empirical){boolVec = currPi > currEps}
		if (sum(boolVec) < length(boolVec)){eps = max(eps,max(currPi[!boolVec]))}
		currPi = currPi[boolVec];currH = currH[boolVec,]
		sumEstimated = sum(currPi)
		print (paste("num configs:", nrow(currH),"xi is:", sumEstimated))
	}
	return(list(H=cbind(currH,currPi),eps=eps,sumEstimated=sumEstimated))
}

get_solution_score<-function(f_1_mat,f_0_mat,pi_h,H,eps=1e-50){
	x1 = exp(log(f_1_mat+eps) %*% t(H))
	x0 = exp(log(f_0_mat+eps) %*% t(1-H))
	z_given_h_mat = x1*x0
	return (sum(z_given_h_mat %*% as.matrix(pi_h)))
}
####################################################################################
################## Below are DP-algorithms for analyzing a gene ####################
# The DP algorithm for analyzing a gene
# under independence assumption
# A[i,j] = the prob of < i 1's in the config vector of 1,...,i
# OUTPUT: sum over the last column in A
# gives prob of <= i 1's
# Therefore, when anlyzing a specific k, use this algorithm
# with k-1 to get genes with at least k non-nulls.
get_gene_analysis<-function(fdrs,k,returnMatrix=F){
	N = length(fdrs)
	A = matrix(1,nrow=k+1,ncol=N)
	A[1,] = exp(cumsum(log(fdrs)))
	if (k==0){return(A[1,length(fdrs)])}
	for (l in 2:nrow(A)){
		for (j in (l-1):length(fdrs)){
			if ((l-1) == j){A[l,j] = prod(1-fdrs[1:j]) ; next}
			A[l,j] = fdrs[j]*A[l,j-1] + (1-fdrs[j])*A[l-1,j-1]
		}
	}
	#print (A)
	if(returnMatrix){return(A)}
	return (sum(A[,length(fdrs)]))
}
# QA: should be 0.5
# get_gene_analysis(rep(0.5,5),2)

fishersMethod = function(x) pchisq(-2 * sum(log(x)),df=2*length(x),lower=FALSE)
# A general version of the DP algorithm
# let m1 = m[,1:ncol(m)/2];m2 = m[,((ncol(m)/2)+1):ncol(m)]
# m1 and m2 are vectors of the same length, defined above
# Use m1 for adding null info, and m2 for non-null
# Assumption: k <= (m/2)
# Keep in mind that for a specific k the computed
# table goes from 0 to k (including k).
get_gene_analysis2<-function(m,k,returnMatrix=F){
	nn = length(m)
	m1 = m[1:(nn/2)];m2 = m[((nn/2)+1):nn]
	N = length(m1)
	A = matrix(1,nrow=k+1,ncol=N)
	A[1,] = exp(cumsum(log(m1)))
	if (k==0){return(A[1,length(m1)])}
	for (l in 2:nrow(A)){
		for (j in (l-1):N){
			if ((l-1) == j){A[l,j] = prod(m2[1:j]) ; next}
			A[l,j] = m1[j]*A[l,j-1] + (m2[j])*A[l-1,j-1]
		}
	}
	#print (A)
	if(returnMatrix){return(A)}
	return (sum(A[,length(m1)]))
}
# Tests:
# both should be 0.5:
#get_gene_analysis2(rep(0.5,10),2)
#get_gene_analysis(rep(0.5,5),2)
# should be: 0.5^5 + 0.5^4*0.25*5 + choose(5,2)*0.25^2*0.5^3
#get_gene_analysis2(c(rep(0.5,5),rep(0.25,5)),2)

####################################################################################
# This is the basic method for obtaining an upper bound for the lfdr 
get_lfdr_upper_bound<-function(H_hat,pi_hat,epsilon,k,f_0_mat,f_1_mat){
	z_given_h_vals = get_z_given_h(H_hat,f_0_mat,f_1_mat)
	null_area_configs = rowSums(H_hat) < k
	term1 = c()
	if (any(null_area_configs)){
		H_hat_k = H_hat[null_area_configs,]
		if(is.null(dim(H_hat_k))){H_hat_k = t(as.matrix(H_hat[null_area_configs,]))}
		z_given_h_vals2 = get_z_given_h(H_hat_k,f_0_mat,f_1_mat)
		term1 = apply(z_given_h_vals2,1,function(x,y){x*y},y=pi_hat[null_area_configs]-epsilon)
		if (!is.null(dim(term1))){term1 = colSums(term1)}
	}
	term2 = apply(cbind(f_0_mat,f_1_mat),1,get_gene_analysis2,k=k-1)
	p_z_apprs = colSums(apply(z_given_h_vals,1,function(x,y){x*y},y=pi_hat))
	if (length(term1)>0){
		lfdr_upper_bound = (term1 + epsilon*term2)/p_z_apprs
	}
	else{
		lfdr_upper_bound = (epsilon*term2)/p_z_apprs
	}
	return (pmin(1,lfdr_upper_bound))
}

# The method above calculates an upper bound for the lfdr.
# Here, we calculate and return the entire table of probabilities.
# That is, for each k we calculate an upper bound for the probability
# that the gene has exactly k 1's (i.e., non-nulls).
get_prob_upper_bound_exact<-function(H_hat,pi_hat,epsilon,f_0_mat,f_1_mat){
	term2 = t(apply(cbind(f_0_mat,f_1_mat),1,
			function(x)get_gene_analysis2(x,length(x)/2,T)[,length(x)/2]))
	z_given_h_vals = get_z_given_h(H_hat,f_0_mat,f_1_mat)
	term1 = c()
	for (k in 0:ncol(f_0_mat)){
		null_area_configs = rowSums(H_hat) == k
		v = rep(0,nrow(f_1_mat))
		if (any(null_area_configs)){
			H_hat_k = H_hat[null_area_configs,]
			if(is.null(dim(H_hat_k))){H_hat_k = t(as.matrix(H_hat[null_area_configs,]))}
			z_given_h_vals2 = get_z_given_h(H_hat_k,f_0_mat,f_1_mat)
			v = apply(z_given_h_vals2,1,function(x,y){x*y},y=pi_hat[null_area_configs]-epsilon)
			if (!is.null(dim(v))){v = colSums(v)}
		}
		term1 = cbind(term1,v)
	}
	p_z_apprs = colSums(apply(z_given_h_vals,1,function(x,y){x*y},y=pi_hat))
	probs = c()
	for (j in 1:ncol(term2)){
		probs = cbind(probs,(term1[,j] + epsilon*term2[,j])/p_z_apprs)
		probs[,j] = pmin(probs[,j],1)
	}
	return(probs)
}
boostrtap_EM_UB<-function(pmat,f_1_mat,f_0_mat,reps=10){
	UBs = c()
	for (j in 1:reps){
		inds = sample(1:nrow(pmat),replace=T)
		pi_h = Space_constrained_estimation(pmat[inds,],
			f_1_mat[inds,],f_0_mat[inds,],convergenceEps = 1e-6,maxiter=1000,nH=512)
		H_hat = pi_h[[1]][,1:ncol(pmat)]
		pi_hat = pi_h[[1]][,ncol(pmat)+1]
		epsilon = pi_h$eps
		lfdrs_rep = get_lfdr_upper_bound(H_hat,pi_hat,epsilon,4,f_0_mat,f_1_mat)
		UBs = cbind(UBs,lfdrs_rep)
	}
	return (UBs)
}

# Here we use the EM process to obtain pairwise covariances
# between studies
getStudyPairwiseEstimation<-function(f_1_mat,f_0_mat,...){
	study_pair_ests = list()
	study_pair_ests[[1]] = list()
	for (i in 2:ncol(f_1_mat)){
		study_pair_ests[[i]] = list()
		for (j in 1:(i-1)){
			curr_est = estimate_prior(f_1_mat[,c(i,j)],f_0_mat[,c(i,j)],...)
			study_pair_ests[[j]][[i]] = curr_est[c(1,3,2,4),]
			study_pair_ests[[i]][[j]] = curr_est
		}
	}
	return (study_pair_ests)
}

get_study_pair_corr_matrix<-function(f_1_mat,f_0_mat,B=50,...){
	if (B<=1){inds = 1:nrow(f_0_mat)}
	corrs = c()
	for (j2 in 1:B){
		inds = sample(1:nrow(f_0_mat),replace=T)
		pairwise_scores = getStudyPairwiseEstimation(f_1_mat[inds,],f_0_mat[inds,],...)
		study_corrs = matrix(0,nrow=ncol(f_1_mat),ncol=ncol(f_1_mat))
		for (i in 2:ncol(study_corrs)){
			for (j in 1:(i-1)){
				curr_est = pairwise_scores[[i]][[j]]
				p1 = curr_est[4,3]+curr_est[2,3]
				p2 = curr_est[4,3]+curr_est[3,3]
				study_corrs[i,j] = (curr_est[4,3] - p1*p2) / sqrt(p1*(1-p1)*p2*(1-p2))
				study_corrs[j,i] = study_corrs[i,j]
			}
		}
		diag(study_corrs) = 1
		study_corrs[is.na(study_corrs)] = 0
		if (j2==1){corrs = study_corrs}
		else{ corrs = corrs + study_corrs}
		inds = sample(1:nrow(f_0_mat),replace=T)
	}
	corrs= corrs/B
	return (corrs)
}
get_study_pair_corr_matrix2<-function(f_1_mat,f_0_mat,B=50,func=median,...){
	if (B<=1){inds = 1:nrow(f_0_mat)}
	corrs = c()
	corr_lists = list()
	for (i in 2:ncol(study_corrs)){
		corr_lists[[i]] = list()
		for (j in 1:(i-1)){
			corr_lists[[i]][[j]] = rep(0,B)
		}
	}
	for (j2 in 1:B){
		inds = sample(1:nrow(f_0_mat),replace=T)
		pairwise_scores = getStudyPairwiseEstimation(f_1_mat[inds,],f_0_mat[inds,],...)
		study_corrs = matrix(0,nrow=ncol(f_1_mat),ncol=ncol(f_1_mat))
		for (i in 2:ncol(study_corrs)){
			for (j in 1:(i-1)){
				curr_est = pairwise_scores[[i]][[j]]
				p1 = curr_est[4,3]+curr_est[2,3]
				p2 = curr_est[4,3]+curr_est[3,3]
				corr_lists[[i]][[j]][j2] = (curr_est[4,3] - p1*p2) / sqrt(p1*(1-p1)*p2*(1-p2))
			}
		}
	}
	study_corrs = matrix(0,nrow=ncol(f_1_mat),ncol=ncol(f_1_mat))
	diag(study_corrs) = 1
	for (i in 2:ncol(study_corrs)){
		for (j in 1:(i-1)){
			study_corrs[i,j] = func(corr_lists[[i]][[j]])
			study_corrs[j,i] = study_corrs[i,j]
		}
	}
	return (study_corrs)
}


####################################################################################

clustering_based_lfdr_analysis<-function(lfdrs,clustering,f_0_mat,f_1_mat,...){
	# After the loop below this matrix should hold in each
	# column j, the exact probability that the number of nulls
	# of each gene is j. Therefore, it goes from 0 to m, and
	# has m+1 columns.
	final_probs = c()
	for (C in sort(unique(clustering))){
		currinds = which(clustering==C)
		if (length(currinds)==1){
			cluster_probs = cbind(lfdrs[,currinds],1-lfdrs[,currinds])
		}
		else{
			print("Performing restricted EM")
			pi_h = Space_constrained_estimation(lfdrs[,currinds],
					f_1_mat[,currinds],f_0_mat[,currinds],...)
			print ("Done")

			H_hat = pi_h[[1]][,1:length(currinds)]
			pi_hat = pi_h[[1]][,length(currinds)+1]
			epsilon = pi_h$eps
			cluster_probs = get_prob_upper_bound_exact(H_hat,pi_hat,epsilon,f_0_mat[,currinds],f_1_mat[,currinds])
		}
		# merge the scores
		if(length(final_probs)==0){final_probs = cluster_probs;next}
		newm = matrix(0,nrow=nrow(final_probs),ncol=ncol(final_probs)+ncol(cluster_probs)-1)
		for (j in 1:ncol(cluster_probs)){
			currk = j-1
			for (j2 in 1:ncol(final_probs)){
				currk2 = j2-1
				newk = currk+currk2
				newv = cluster_probs[,j]*final_probs[,j2]
				newm[,newk+1] = newm[,newk+1] + newv
			}
		}
		final_probs = newm
	}
	# We now apply cumsum. Column 1 will have the prob
	# that the number of nulls is at most 0. Using the same argument
	# column j will hold the prob of at most j-1 non-nulls, which
	# is exactly the lfdr value we need (for each column).
	final_probs = t(apply(final_probs,1,cumsum))
	for (j in 1:ncol(final_probs)){final_probs[,j] = pmin(1,final_probs[,j])}
	colnames(final_probs) = 0:(ncol(final_probs)-1)
	rownames(final_probs) = rownames(lfdrs)
	return(final_probs)
}

perform_clustering_based_lfdr_analysis<-function(lfdrs,f_0_mat,f_1_mat,...){
	study_corrs = get_study_pair_corr_matrix(f_1_mat,f_0_mat,B=50,convergenceEps=1e-4)
	clustering = run_leadingeigen_clustering(study_corrs)
	print ("########## Cluster sizes: ###########")
	print(table(clustering))
	return (clustering_based_lfdr_analysis(lfdrs,clustering,f_0_mat,f_1_mat,...))
}

# This is an old method used for calculating lfdrs given two parameters: delta and d,
# (see the writeup from May 2016)
get_lfdr_upper_bound_given_clustering<-function(pvals,lfdrs,clustering,delta,d,f_0_mat,f_1_mat,nH=512){
	cluster_lfdrs = c()
	for (C in sort(unique(clustering))){
		currinds = which(clustering==C)
		if (length(currinds)==1){
			cluster_lfdrs = cbind(cluster_lfdrs,lfdrs[,currinds])
			next
		}

		print("Performing restricted EM")
		pi_h = Space_constrained_estimation(pvals[,currinds],
				f_1_mat[,currinds],f_0_mat[,currinds],nH)
		print ("Done")

		H_hat = pi_h[[1]][,1:length(currinds)]
		pi_hat = pi_h[[1]][,length(currinds)+1]
		epsilon = pi_h$eps

		# Can we use the correction algorithm?
		H_hat = H_hat[pi_hat > epsilon,]
		pi_hat = pi_hat[pi_hat > epsilon]

		print ("Calculating the lfdr upper bound")
		k = ceiling(delta*length(currinds))
		curr_lfdr_upper_bound = get_lfdr_upper_bound(H_hat,pi_hat,epsilon,k,
			f_0_mat[,currinds],f_1_mat[,currinds])
		cluster_lfdrs = cbind(cluster_lfdrs,curr_lfdr_upper_bound)
	}
	colnames(cluster_lfdrs) = sort(unique(clustering))
	rownames(cluster_lfdrs) = rownames(pvals)
	final_lfdrs = apply(cluster_lfdrs,1,get_gene_analysis,k=d-1)
	return(final_lfdrs)
}

run_simple_clustering<-function(x,cor_thr=0.1){
	y = x
	x = x >= cor_thr
	mode(x) = 'numeric';diag(x)=0
	res = list()
	g = graph.adjacency(x,mode='undirected',weighted=T)
	return (clusters(g)$membership)
}

try({library(igraph);library(kernlab)})
run_leadingeigen_clustering<-function(x,cor_thr=0.1){
	x = x >= cor_thr
	mode(x) = 'numeric';diag(x)=1
	g = graph.adjacency(x,mode='undirected',weighted=T)
	return (cluster_infomap(g)$membership)
}


run_specc<-function(x,improvement_thr = 0.05,cor_thr=0.1){
	y = x
	x = x >= cor_thr
	mode(x) = 'numeric';diag(x)=0
	res = list()
	g = graph.adjacency(x,mode='undirected',weighted=T)
	return (clusters(g)$membership)
	
	# Old code
	if (length(E(g))==0){return(1:ncol(x))}
	ws = E(g)$weight
	bin_ws = ws; bin_ws[ws<0] = 0
	num_clusters = 1:ncol(x)
	modularities = c()
	for(k in num_clusters){
		res[[as.character(k)]] = 1:nrow(x)
		if (k < nrow(x)){try({res[[as.character(k)]] = kmeans(y,k)$cluster})}
		modularities[as.character(k)] = modularity(g,membership = res[[as.character(k)]])
	}
	# select the best
	worth_it = rep(T,length(modularities))
	for (i in 2:length(modularities)){
		worth_it[i] = modularities[i]-modularities[i-1] >= improvement_thr
	}
	best_k = num_clusters[max(which(worth_it))]
	return(res[[as.character(best_k)]])
}
####################################################################################
###  Executables ###

get_study_marginal_estimation<-function(pvals,lfdr_method='znormix',use_power=T,threegroups=T,...){
	if (lfdr_method == 'locfdr'){
		lfdr_objs = apply(pvals,2,get_f1_f0_locfdr,dfs=7,...)
		f_1_mat = sapply(lfdr_objs,function(x)x[[1]])
		f_1_mat[is.nan(f_1_mat)] = 0
		f_0_mat = sapply(lfdr_objs,function(x)x[[2]])
		lfdr_objs = apply(pvals,2,get_lfdrs_locfdr,dfs=7,...)
		pi0_study_vec = sapply(lfdr_objs,function(x) x$pi_0)
		lfdrs = sapply(lfdr_objs,function(x)x$fdr)
		pws = sapply(lfdr_objs,function(x)x$Efdr)[1,]
	}
	if (lfdr_method=='znormix'){
		lfdr_objs = apply(pvals,2,get_f1_f0_znormix,threegroups=threegroups,...)
		f_1_mat = sapply(lfdr_objs,function(x)x[[1]])
		f_1_mat[is.nan(f_1_mat)] = 0
		f_0_mat = sapply(lfdr_objs,function(x)x[[2]])
		lfdr_objs = apply(pvals,2,get_lfdrs_znormix,threegroups=threegroups,...)
		pi0_study_vec = sapply(lfdr_objs,function(x) x$params[1])
		lfdrs = sapply(lfdr_objs,function(x)x$fdr)
		pws = c()
		for (j in 1:ncol(f_0_mat)){
			curr_tdrs = sum(1-lfdrs[,j])
			curr_rejs = sum(lfdrs[,j]<0.3)
			pws[j] = 0.1
			if (curr_tdrs > 0){pws[j] = max(0.1,curr_rejs/curr_tdrs)}
		}
	}
	if (lfdr_method=='heavy_1_tail'){
		lfdr_objs = apply(pvals,2,run_two_groups_with_heavy_1_tail)
		f_1_mat = sapply(lfdr_objs,function(x)x$f_1)
		f_1_mat[is.nan(f_1_mat)] = 0
		f_0_mat = sapply(lfdr_objs,function(x)x$f_0)
		pi0_study_vec = sapply(lfdr_objs,function(x) x$pi_0)
		lfdrs = sapply(lfdr_objs,function(x)x$fdr)
		pws = c()
		for (j in 1:ncol(f_0_mat)){
			curr_tdrs = sum(1-lfdrs[,j])
			curr_rejs = sum(lfdrs[,j]<0.3)
			pws[j] = 0.1
			if (curr_tdrs > 0){pws[j] = max(0.1,curr_rejs/curr_tdrs)}
		}
	}
	if (lfdr_method=='repfdr_isotonic'){
		lfdr_objs = apply(pvals,2,run_isotonic_repfdr)
		f_1_mat = sapply(lfdr_objs,function(x)x[[2]])
		f_1_mat[is.nan(f_1_mat)] = 0
		f_0_mat = sapply(lfdr_objs,function(x)x[[3]])
		pi0_study_vec = sapply(lfdr_objs,function(x) x[[1]])
		lfdrs = sapply(lfdr_objs,function(x)x$fdr)
		pws = c()
		for (j in 1:ncol(f_0_mat)){
			curr_tdrs = sum(1-lfdrs[,j])
			curr_rejs = sum(lfdrs[,j]<0.3)
			pws[j] = 0.1
			if (curr_tdrs > 0){pws[j] = max(0.1,curr_rejs/curr_tdrs)}
		}
	}

	# These are used to shrink the effect of f_1
	if (use_power){
	# multiply by the power of the dataset
		for (j in 1:ncol(f_0_mat)){
			pw = pws[j];print (pw)
			f_1_mat[,j] = f_1_mat[,j]*pw
		}
	}
	return(list(lfdr_objs=lfdr_objs,f_1_mat=f_1_mat,f_0_mat=f_0_mat,pi0_study_vec=pi0_study_vec,pws=pws))
}

# In all methods below:
#	pvals: a p-value matrix
#	ks: a range of k-values for local fdr calculations
#	nH - number of allowed configurations
#	emEps: convergence parameter for the EM
#	lfdr_method: either "znormix" (normix) or "locfdr"
#	use_power: binary, specify shrinking the f1 densities
#	threegroups: binary, relevant for normix. TRUE: fit three Gaussians; FALSE - fit two
#	... - additional parameters for znormix or locfdr
analyze_pval_matrix_slim<-function(pvals,ks = 3:6,nH=1024,emEps =1e-6, lfdr_method='znormix',use_power=T,threegroups=T,...){
	print ("Analyzing each study")
	mar_est = get_study_marginal_estimation(pvals,lfdr_method=lfdr_method,use_power=use_power,threegroups=threegroups,...)
	lfdr_objs = mar_est$lfdr_objs
	f_1_mat = mar_est$f_1_mat
	f_0_mat = mar_est$f_0_mat
	pi0_study_vec = mar_est$pi0_study_vec
	print ("Done")

	# 1. Analysis under independence assumption
	print ("lfdr analysis under independence assumption")
	lfdrs = sapply(lfdr_objs,function(x)x$fdr)
	lfdrs_indep_assumption = c()
	for (k in ks){
		lfdrs_indep_assumption = cbind(lfdrs_indep_assumption,
			apply(lfdrs,1,get_gene_analysis,k=k-1))
	}
	colnames(lfdrs_indep_assumption) = paste("Indep",ks)
	print ("Done")

	# 2. EM-based analysis
	print("Performing restricted EM")
	pi_h = Space_constrained_estimation(pvals,f_1_mat,f_0_mat,nH=nH,convergenceEps=emEps)
	print ("Done")
	H_hat = pi_h[[1]][,1:ncol(pvals)]
	pi_hat = pi_h[[1]][,ncol(pvals)+1]
	epsilon = pi_h$eps
	print ("Calculating the lfdr upper bound")
	lfdr_upper_bound = c()
	for (k in ks){
		lfdr_upper_bound = cbind(lfdr_upper_bound,
			get_lfdr_upper_bound(H_hat,pi_hat,epsilon,k,f_0_mat,f_1_mat))
	}
	colnames(lfdr_upper_bound) = paste('EM-UB',ks)

	# 3. EM-Clustering-DP algorithm
	clust_analysis = perform_clustering_based_lfdr_analysis(lfdrs,f_0_mat,f_1_mat,nH=nH,convergenceEps=emEps)
	clust_analysis = clust_analysis[,ks]
	colnames(clust_analysis) = paste("Clust",ks)

	# Other, much more "naive analyses"
	# 1. Standard meta-analysis
	fisher_pvals = apply(pvals,1,fishersMethod)

	# 2. Expected number of non-nulls
	expected_datasets = rowSums(1-lfdrs)

	# 3. Naive p-vals analysis
	qmat = apply(pvals,2,p.adjust,method='fdr')
	num_q_datasets = rowSums(qmat <= 0.1)

	# Put all q-val / lfdr estimation in the same matrix
	all_lfdr_rep_scores = cbind(p.adjust(fisher_pvals,method='fdr'),lfdrs_indep_assumption)
	all_lfdr_rep_scores = cbind(all_lfdr_rep_scores,lfdr_upper_bound)
	all_lfdr_rep_scores = cbind(all_lfdr_rep_scores,clust_analysis)
	#all_lfdr_rep_scores = cbind(all_lfdr_rep_scores,lfdrs_hierarch)
	rownames(all_lfdr_rep_scores) = rownames(pvals)
	colnames(all_lfdr_rep_scores)[1] = "fisher_qval"

	# Put all "count-based" scores in another mat
	all_count_rep_scores = cbind(expected_datasets,num_q_datasets)
	colnames(all_count_rep_scores) = c("Expected count","Num q<=0.1")

	obj = list(
		lfdrs = lfdrs,f_0_mat = f_0_mat,f_1_mat = f_1_mat,
		pi0_study_vec = pi0_study_vec,
		H_hat = H_hat, pi_hat=pi_hat, epsilon=epsilon,
		all_lfdr_rep_scores=all_lfdr_rep_scores,
		all_count_rep_scores = all_count_rep_scores
	)
	return (obj)
}

SCREEN<-function(pvals,ks = 3:6,nH=1024,emEps =1e-6, lfdr_method='znormix',use_power=T,threegroups=T,...){
	print ("Analyzing each study")
	mar_est = get_study_marginal_estimation(pvals,lfdr_method=lfdr_method,use_power=use_power,threegroups=threegroups,...)
	lfdr_objs = mar_est$lfdr_objs
	f_1_mat = mar_est$f_1_mat
	f_0_mat = mar_est$f_0_mat
	pi0_study_vec = mar_est$pi0_study_vec
	print ("Done")
	print ("obtain local fdrs in each study")
	lfdrs = sapply(lfdr_objs,function(x)x$fdr)
	print ("Done")
	clust_analysis = perform_clustering_based_lfdr_analysis(lfdrs,f_0_mat,f_1_mat,nH=nH,convergenceEps=emEps)
	clust_analysis = clust_analysis[,ks]
	colnames(clust_analysis) = paste("SCREEN",ks)
	return (clust_analysis)
}

SCREEN_ind<-function(pvals,ks = 3:6,lfdr_method='znormix',use_power=T,threegroups=T,...){
	print ("Analyzing each study")
	mar_est = get_study_marginal_estimation(pvals,lfdr_method=lfdr_method,use_power=use_power,threegroups=threegroups,...)
	lfdr_objs = mar_est$lfdr_objs
	f_1_mat = mar_est$f_1_mat
	f_0_mat = mar_est$f_0_mat
	pi0_study_vec = mar_est$pi0_study_vec
	print ("Done")
	print ("lfdr analysis under independence assumption")
	lfdrs = sapply(lfdr_objs,function(x)x$fdr)
	lfdrs_indep_assumption = c()
	for (k in ks){
		lfdrs_indep_assumption = cbind(lfdrs_indep_assumption,
			apply(lfdrs,1,get_gene_analysis,k=k-1))
	}
	colnames(lfdrs_indep_assumption) = paste("SCREEN-ind",ks)
	print ("Done")
	return (lfdrs_indep_assumption)
}
