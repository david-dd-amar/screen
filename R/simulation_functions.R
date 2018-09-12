######### Functions ########
get_pred_scores<-function(rep_results,pmat,nonull_inds,sig_genes,thr=0.2){
	lfdrs_j = sum(rep_results$lfdrs<thr & nonull_inds==1)/sum(rep_results$lfdrs<thr | nonull_inds==1)
	lfdrs_FDR = sum(rep_results$lfdrs<thr & nonull_inds==0)/sum(rep_results$lfdrs<thr)
	gene2num = apply(nonull_inds,1,sum)
	rep_fdrs = rep_results$lfdr_rep_scores[,-1]
	ks = as.numeric(sapply(colnames(rep_fdrs),function(x)strsplit(x,split=' ')[[1]][2]))
	ks[is.na(ks)] = as.numeric(sapply(colnames(rep_fdrs)[is.na(ks)],function(x)strsplit(x,split='d=')[[1]][2]))
	fisherq = rep_results$lfdr_rep_scores[,1]
	fisher_rep_scores = c()
	method2score = c()
	# fisherq scores for the tested k values
	for (k in unique(ks)){
		j = lfdrs_vs_genenum_score(fisherq,gene2num,k,0.1)
		FDR = lfdrs_vs_genenum_score(fisherq,gene2num,k,0.1,"FDR")
		method2score = rbind(method2score,c(j,FDR))
		rownames(method2score)[nrow(method2score)] = paste("FisherQ",k)
	}
	for (i in 1:length(ks)){
		currk = ks[i]
		j = lfdrs_vs_genenum_score(rep_fdrs[,i],gene2num,currk,thr)
		FDR = lfdrs_vs_genenum_score(rep_fdrs[,i],gene2num,currk,thr,"FDR")
		method2score = rbind(method2score,c(j,FDR))
		rownames(method2score)[nrow(method2score)] = paste(colnames(rep_fdrs)[i],k)
	}
	colnames(method2score) = c("Jaccard","FDR")
	return(list(lfdrs_j=lfdrs_j,lfdrs_FDR=lfdrs_FDR,method2score=method2score))
}

get_pred_scores2<-function(rep_results,pmat,nonull_inds,thr=0.2){
	lfdrs_j = sum(rep_results$lfdrs<thr & nonull_inds==1)/sum(rep_results$lfdrs<thr | nonull_inds==1)
	lfdrs_FDR = sum(rep_results$lfdrs<thr & nonull_inds==0)/sum(rep_results$lfdrs<thr)
	gene2num = apply(nonull_inds,1,sum)
	rep_fdrs = rep_results$all_lfdr_rep_scores[,-1]
	ks = as.numeric(sapply(colnames(rep_fdrs),function(x)strsplit(x,split=' ')[[1]][2]))
	ks[is.na(ks)] = as.numeric(sapply(colnames(rep_fdrs)[is.na(ks)],function(x)strsplit(x,split='d=')[[1]][2]))
	methods1 = sapply(colnames(rep_fdrs),function(x)strsplit(x,split=' ')[[1]][1])

	fisherq = rep_results$all_lfdr_rep_scores[,1]
	jaccard_scores = c();fdp_scores = c()
	x1 = c();x2=c()
	# fisherq scores for the tested k values
	for (k in unique(ks)){
		x1[as.character(k)] = lfdrs_vs_genenum_score(fisherq,gene2num,k,0.1)
		x2[as.character(k)] = lfdrs_vs_genenum_score(fisherq,gene2num,k,0.1,"FDR")
	}
	jaccard_scores = rbind(jaccard_scores,x1);fdp_scores = rbind(fdp_scores,x2)
	rownames(jaccard_scores)[1] = "FisherQ"

	methods1 = unique(methods1)
	ks = unique(ks)
	for (j in 1:length(unique(methods1))){
		x1=c();x2=c()
		for (i in 1:length(unique(ks))){
			currk = ks[i]
			currname = paste(methods1[j],currk)
			print(currname)
			x1[as.character(currk)] = lfdrs_vs_genenum_score(rep_fdrs[,currname],gene2num,currk,thr)
			x2[as.character(currk)] = lfdrs_vs_genenum_score(rep_fdrs[,currname],gene2num,currk,thr,"FDR")
		}
		jaccard_scores = rbind(jaccard_scores,x1);fdp_scores = rbind(fdp_scores,x2)
		rownames(jaccard_scores)[nrow(jaccard_scores)] = methods1[j]
	}

	count_scores = rep_results$all_count_rep_scores
	methods2 = colnames(count_scores)
	for (j in 1:length(unique(methods2))){
		x1=c();x2=c()
		for (i in 1:length(unique(ks))){
			currk = ks[i]
			x1[as.character(currk)] = lfdrs_vs_genenum_score(-count_scores[,methods2[j]],gene2num,currk,-currk)
			x2[as.character(currk)] = lfdrs_vs_genenum_score(-count_scores[,methods2[j]],gene2num,currk,-currk,"FDR")
		}
		jaccard_scores = rbind(jaccard_scores,x1);fdp_scores = rbind(fdp_scores,x2)
		rownames(jaccard_scores)[nrow(jaccard_scores)] = methods2[j]
	}

	rownames(fdp_scores) = rownames(jaccard_scores)
	return(list(lfdrs_j=lfdrs_j,lfdrs_FDR=lfdrs_FDR,fdp_scores=fdp_scores,jaccard_scores=jaccard_scores))
}
# x = lfdrs for a given x
# y = the real gene2num
lfdrs_vs_genenum_score<-function(x,y,k,thr,method="jaccard"){
	g1 = x<=thr
	if (sum(g1)==0){return(0)}
	g2 = y>=k
	if (method=="jaccard"){
		return (sum(g1&g2)/sum(g1|g2))
	}
	if (method=="FDR"){
		return (sum(g1&!g2)/sum(g1))
	}
	return (0)
}
# x = lfdrs for a given x
# y = the real gene2num
lfdrs_vs_set_score<-function(x,setinds,k,thr,method="jaccard"){
	y = rep(0,length(x))
	y[setinds] = k
	return(lfdrs_vs_genenum_score(x,y,k,thr,method))
}	

simulations_case1<-function(m,ngenes,nsigs_across,nsig_studies,nsig_within,use_beta=T,
		params1,params2,params3){
	pmat = matrix(0,nrow=ngenes,ncol=m)
	rownames(pmat) = paste("gene",1:ngenes)
	colnames(pmat) = paste("study",1:m)
	nonull_inds = matrix(0,nrow=ngenes,ncol=m)
	rownames(nonull_inds) = paste("gene",1:ngenes)
	colnames(nonull_inds) = paste("study",1:m)
	for (j in 1:m){
		if (use_beta){
			ps = rbeta(ngenes,params1[1],params1[2])
		}
		else{
			ps = pnorm(rnorm(ngenes,params1[1],params1[2]))
		}
		names(ps) = rownames(pmat)
		other_genes = sample(rownames(pmat))[1:(2*nsig_within)]
		nonull_inds[other_genes,j] = 1
		if (use_beta){
			ps[other_genes[1:nsig_within]] = rbeta(nsig_within,params3[1],params3[2])
			ps[other_genes[(nsig_within+1):length(other_genes)]] = rbeta(nsig_within,params2[1],params2[2])
		}else{
			ps[other_genes[1:nsig_within]] = pnorm(rnorm(nsig_within,params3[1],params3[2]))
			ps[other_genes[(nsig_within+1):length(other_genes)]] = pnorm(rnorm(nsig_within,params2[1],params2[2]))
		}
		pmat[,j] = ps
	}
	sig_genes = 1:nsigs_across
	nonull_inds[sig_genes,1:nsig_studies] = 1
	for (s in sig_genes){
		curr_studies = sample(1:m)[1:nsig_studies]
		if (use_beta){
			pmat[s,curr_studies] = rbeta(nsig_studies,params3[1],params3[2])
		}else{
			pmat[s,curr_studies] = pnorm(rnorm(nsig_studies,params3[1],params3[2]))
		}
		nonull_inds[s,curr_studies] = 1
	}
	gene2num=apply(nonull_inds,1,sum)
	return(list(pmat=pmat,nonull_inds=nonull_inds,gene2num=gene2num))
}

library(MASS)
simulations_case2<-function(m,ngenes,nsig_within,sigmaMat,use_beta=T,
		params1,params2,params3,num_clusters=1){

	x = c();y=c()
	for (kk in 1:num_clusters){
		nonull_inds = mvrnorm(ngenes,mu=rep(0,m),Sigma=sigmaMat)
		nonull_inds = nonull_inds >= quantile(c(nonull_inds),1-(2*nsig_within)/ngenes)
		mode(nonull_inds) = "numeric"
		rownames(nonull_inds) = paste("gene",1:ngenes)
		colnames(nonull_inds) = paste("study",1:m)
		all_null_inds = nonull_inds==0

		pmat = matrix(0,nrow=ngenes,ncol=m)
		rownames(pmat) = paste("gene",1:ngenes)
		colnames(pmat) = paste("study",1:m)

		if (use_beta){
			pmat[all_null_inds] = rbeta(sum(all_null_inds),params1[1],params1[2])
		}else{
			pmat[all_null_inds] = pnorm(rnorm(sum(all_null_inds),params1[1],params1[2]))
		}
		for (j in 1:m){
			nonull_genes = sample(which(nonull_inds[,j] == 1))
			nn = length(nonull_genes)
			if (use_beta){
				pmat[nonull_genes[1:(nn/2)],j] = rbeta(nn/2,params3[1],params3[2])
				pmat[nonull_genes[(nn/2+1):nn],j] = rbeta(length((nn/2+1):nn),params2[1],params2[2])
			}else{
				pmat[nonull_genes[1:(nn/2)],j] = pnorm(rnorm(nn/2,params3[1],params3[2]))
				pmat[nonull_genes[(nn/2+1):nn],j] = pnorm(rnorm(length((nn/2+1):nn),params2[1],params2[2]))
			}
		}
		x = cbind(x,pmat)
		y = cbind(y,nonull_inds)
	}
	gene2num=apply(y,1,sum)
	return(list(pmat=x,nonull_inds=y,gene2num=gene2num))
}
