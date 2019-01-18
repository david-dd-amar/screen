source("https://raw.githubusercontent.com/david-dd-amar/screen/master/R/SCREEN.R")
source("https://raw.githubusercontent.com/david-dd-amar/screen/master/R/twogroups_methods.R")

# real data
scr_path = "/Users/David/Desktop/MoTrPAC/PA_database/screen_res/"
pvals_files = list.files(scr_path)
pvals_files = pvals_files[grepl("pvals.txt",pvals_files)]

par(mfrow=c(2,2))
pvals = read.delim(paste(scr_path,pvals_files[2],sep=""))
pvals = pvals[,1:8]
apply(pvals,2,hist)

pvals_cut = apply(pvals,2,function(x)as.numeric(cut(x,breaks = 200,ordered_result = T))/200)
colnames(pvals_cut) = colnames(pvals)
rownames(pvals_cut) = rownames(pvals)
hist(pvals_cut[,4])
pvals_cut2 = apply(pvals,2,function(x)pmax(x,0.0001))
colnames(pvals_cut2) = colnames(pvals)
rownames(pvals_cut2) = rownames(pvals)

marg_est = get_study_marginal_estimation(pvals_cut2,lfdr_method = "znormix",use_power = F,threegroups = F)
marg_est = get_study_marginal_estimation(pvals_cut2,lfdr_method = 'bum',use_power = F)
marg_est = get_study_marginal_estimation(pvals,lfdr_method = 'bum',use_power = F)
marg_est = get_study_marginal_estimation(pvals,lfdr_method = "znormix",use_power = T,threegroups = F)

par(mfrow=c(2,4))
for(j in 1:4){
  plot(marg_est$f_1_mat[,j],pvals[,j],pch=20,cex=0.5)
  plot(marg_est$f_0_mat[,j],pvals[,j],pch=20,cex=0.5)
}
par(mfrow=c(2,4))
for(j in 1:8){
  plot(marg_est$f_1_mat[,j]/marg_est$f_0_mat[,j],pvals[,j],pch=20,cex=0.5)
}

corrs = get_study_pair_corr_matrix(marg_est$f_1_mat,marg_est$f_0_mat,B = 5,convergenceEps=1e-3)
library(corrplot)
corrplot(corrs)
corrplot(corrs>0.2)

scr = SCREEN(pvals_cut2,corr_net_B = 5,corr_net_eps = 1e-3,lfdr_method = "bum",use_power = F)
scr = SCREEN(pvals,corr_net_B = 10,corr_net_eps = 1e-3,lfdr_method = "znormix",use_power = F)
colSums(scr < 0.2)
gs = rownames(pvals)[scr[,ncol(scr)-1]<0.2]
barplot(as.numeric(-log(pvals["1282",],10)))
xx = as.matrix(-log(pvals[gs,],10))
library(gplots)
heatmap.2(xx,trace = "none",scale = "none")
yy = xx>2
mode(yy) = "numeric"
heatmap.2(yy)
