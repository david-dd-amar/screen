# Beta-uniform estimation of two-groups using BioNet
library(BioNet)
run_bum<-function(ps){
  model = fitBumModel(ps,plot=F)
  beta0 = c(1,1)
  beta1 = c(model$a,1)
  pi0 = (model$lambda) + (1-model$lambda)*model$a
  f0 = dbeta(ps,beta0[1],beta0[2])
  f1 = dbeta(ps,beta1[1],beta1[2])
  fdr = pi0*f0 / (pi0*f0 + (1-pi0)*f1)
  params = c(pi0,beta0,beta1)
  names(params) = c("pi0","b0_shape1","b0_shape2","b1_shape1","b1_shape2")
  return(list(fdr = fdr, params=params,f_0=f0,f_1=f1))
}
get_f0_f1_bum<-function(ps,obj=NULL){
  if(is.null(obj)){obj = run_bum(ps)}
  return(obj)
}

######################### Wrapper methods for the locfdr package ########################
# nulltype is 1 or 2
run_locfdr<-function(zs,dfs=15,nulltype=1,...){
	try({
		obj = locfdr(zs,df=dfs,nulltype=nulltype)
		if (nulltype==1){p0_1 = obj$fp0[3,3]}
		if (nulltype==0){p0_1 = min(0.9999,obj$fp0[1,3])}
		if (nulltype>1){p0_1 = obj$fp0[5,3]}
		if (p0_1 < 1){obj$pi_0 = p0_1;return(obj)}
	})
	ntypes = c(1,0)
	original_ntype = nulltype
	ntypes = setdiff(ntypes,nulltype)
	for (nulltype in ntypes){
		print (nulltype)
		try({
			obj = locfdr(zs,df=dfs,nulltype=nulltype)
			if (nulltype==1){p0_1 = obj$fp0[3,3]}
			if (nulltype==0){p0_1 = obj$fp0[1,3]}
			if (nulltype>1){p0_1 = obj$fp0[5,3]}
			if (p0_1 < 1){obj$pi_0 = p0_1;return(obj)}
		})
	}
	# If both failed, scale the Zs
	zs = scale(zs)
	ntype = c(original_ntype,zs)
	for (nulltype in ntypes){
		print (nulltype)
		try({
			obj = locfdr(zs,df=dfs,nulltype=nulltype)
			if (nulltype==1){p0_1 = obj$fp0[3,3]}
			if (nulltype==0){p0_1 = obj$fp0[1,3]}
			if (nulltype>1){p0_1 = obj$fp0[5,3]}
			if (p0_1 < 1){obj$pi_0 = p0_1;return(obj)}
		})
	}
	return (NULL)
}

get_lfdrs_locfdr<-function(p,...){
	zs = qnorm(1-p);
	zs[is.na(zs)] = 0
	zs[zs<0 & is.infinite(zs)] = min(zs[!is.infinite(zs)])
	zs[zs>0 & is.infinite(zs)] = max(zs[!is.infinite(zs)])
	obj = run_locfdr(zs,...)
	return (obj)
}

get_f1_f0_locfdr<-function(p,...){
	zs = qnorm(1-p);
	zs[is.na(zs)] = 0
	zs[zs<0 & is.infinite(zs)] = min(zs[!is.infinite(zs)])
	zs[zs>0 & is.infinite(zs)] = max(zs[!is.infinite(zs)])
	#hist(zs,freq=F);curve(dnorm,add=T)
	obj = run_locfdr(zs,...)
	#obj = run_locfdr(zs)
	f0dist = obj$mat[,"f0"]/length(zs)
	names(f0dist) = obj$mat[,"x"]
	p0 = min(obj$pi_0,0.999)
	f1dist = (obj$mat[,"p1f1"])/(length(zs)*(1-p0))
	names(f1dist) = obj$mat[,"x"]
	#par(mfrow=c(1,2));plot(f0dist);plot(f1dist)
	# tests  - plot should result in a line
	#calculated_lfdrs = (f0dist*p0)/(f0dist*p0+f1dist*(1-p0))
	#plot(calculated_lfdrs,obj$mat[,"fdr"]);abline(0,1)
	binned_zs = sapply(zs,function(x,y)y[abs(x-y)==min(abs(x-y))],y=as.numeric(names(f1dist)))
	return (list(f_1 = f1dist[as.character(binned_zs)], 
		f_0 = f0dist[as.character(binned_zs)],pi_0=p0))
}

################################## Normix implementation #########################
get_lfdrs_znormix<-function(p,eps=1e-5,niter=2000,threegroups=T,...){
	p[is.nan(p)] = 1;	p[is.na(p)] = 1; newobj = list()
	if(threegroups){
		obj = modified_znormix_halfnorm(p,eps=eps,niter=niter,start.pi0=0.85,...)
	}else{
		obj = modified_znormix(p,eps=eps,niter=niter,start.pi0=0.85,...)
	}
	newobj[['fdr']] = attr(obj,"lfdr")
	newobj[['Fdr']] = attr(obj,"fdr")
	newobj[['params']] = obj[1:5]
	return (newobj)
}

#newobj[['fdr']] = rowMeans(sapply(objs,function(obj)attr(obj,"lfdr")))
#newobj[['params']] = rowMeans(sapply(objs,function(obj)obj[1:5]))
#ps = pmin(p,1-p)
#newobj[['fdr']] = correct_mono(ps,newobj[['fdr']])
correct_mono<-function(ps,x,func="min"){
	cutt_ps = cut(pmin(ps,1-ps),100,ordered_result=T)
	ls = levels(cutt_ps)
	for (i in 1:(length(ls)-1)){
		currl = ls[i]
		currinds = cutt_ps==currl
		curr_ps = ps[currinds];curr_x = x[currinds]
		nextinds = cutt_ps > currl
		next_x = x[nextinds]
		if (func=="min"){
			x[currinds] = pmin(x[currinds],min(next_x))
		}else{
			x[currinds] = pmax(x[currinds],max(next_x))
		}
	}
	return(x)
}
correct_mono2<-function(ps,x,func="min"){
	cutt_ps = cut(ps,100,ordered_result=T)
	ls = levels(cutt_ps)
	for (i in 1:(length(ls)-1)){
		currl = ls[i]
		currinds = cutt_ps==currl
		curr_ps = ps[currinds];curr_x = x[currinds]
		nextinds = cutt_ps > currl
		next_x = x[nextinds]
		x[currinds] = pmax(x[currinds],max(next_x))
	}
	return(x)
}

# If powerthr is null do not use the power of the dataset.
get_f1_f0_znormix<-function(p,eps=1e-5,niter=2000,threegroups=T,...){
	obj = get_lfdrs_znormix(p,eps=eps,niter=niter,threegroups=threegroups,...)
	params=obj$params
	zs = qnorm(1-p)
	zs[is.na(zs)] = 0
	zs[zs<0 & is.infinite(zs)] = min(zs[!is.infinite(zs)])
	zs[zs>0 & is.infinite(zs)] = max(zs[!is.infinite(zs)])
	print(params)
	f0dist = dnorm(zs,mean=params[2],sd=params[3])
	if(threegroups){
		f1dist = dnorm(abs(zs),mean=params[4],sd=params[5])
		# correct the monotonicity issue
		f1dist = correct_mono(pmin(p,1-p),f1dist,'max')
	}else{
		f1dist = dnorm(zs,mean=params[4],sd=params[5])
		# correct the monotonicity issue
		f1dist = correct_mono2(p,f1dist)
	}
	return (list(f_1 = f1dist, f_0 = f0dist))
}

modified_znormix_halfnorm <- function (p, start.pi0=0.9, 
		eps = 1e-05, niter = 2000, verbose = FALSE,
		correctMono=F,standardNull=F,sameSDs=F,sdval=1){
    z = as.matrix(qnorm(1 - p))
    z[is.infinite(z) & z < 0] = min(z[is.finite(z)])
    z[is.infinite(z) & z > 0] = max(z[is.finite(z)])
    absz = abs(z)
    G = length(z)
    stopifnot(G >= 4)
    if (start.pi0 <= 0) 
        start.pi0 = 0.001
    if (start.pi0 >= 1) 
        start.pi0 = 1 - 0.001
    zcut = quantile(absz, start.pi0)
    z0.idx = which(absz < zcut)
    last.par = c(start.pi0, 0, sd(drop(z[z0.idx])), 
    mean(abs(z[-z0.idx])), sd(drop(abs(z[-z0.idx]))))
    iter = 1
    new.par = last.par
    repeat {
        f0 = last.par[1] * 2*dnorm(absz, last.par[2], last.par[3])
        ppee = pmin(1 - 1e-06, pmax(1e-06, f0/(f0 + (1 - last.par[1]) * 
            dnorm(absz, last.par[4], last.par[5]))))
        new.par[1] = mean(ppee)
        sum.ppee = sum(ppee)
        new.par[4] = crossprod(absz, 1 - ppee)/(G - sum.ppee)
        new.par[5] = sqrt(crossprod((absz - new.par[4])^2, 1 - ppee)/(G - 
            sum.ppee))
        new.par[2] = crossprod(z, ppee)/sum.ppee
        new.par[3] = sqrt(crossprod((z - new.par[2])^2, ppee)/sum.ppee)
 	  new.par[2] = 0
	  new.par[3] = max(new.par[3],1)   
	  if(standardNull){new.par[3]=1}
	  if (sameSDs){new.par[3]=sdval;new.par[5]=sdval}
        if (iter >= niter || max(abs(new.par - last.par)) < eps) 
            break
        last.par = new.par
        iter = iter + 1
    }
    ord = order(ppee)
    fdr = numeric(G)
    fdr[ord] = cumsum(ppee[ord])/(1:G)

    # correct for the monotonicity issue when sd.z1 is too low
    if (correctMono){
	    for (j in 1:20){
			print (cor(fdr,pmin(p,1-p),method='spearman'))
			plot(fdr,pmin(p,1-p))
			if (cor(fdr,pmin(p,1-p),method='spearman') > 0.999999){break}
			new.par[5] = new.par[5]+0.01
			print ('in correction loop')
			f0 = new.par[1] * 2*dnorm(absz, new.par[2], new.par[3])
			ppee = pmin(1 - 1e-06, pmax(1e-06, f0/(f0 + (1 - new.par[1]) * 
            		dnorm(absz, new.par[4], new.par[5]))))
			new.par[1] = mean(ppee)
			ord = order(ppee);fdr = numeric(G)
			fdr[ord] = cumsum(ppee[ord])/(1:G)
	    }
    }

    names(new.par) = c("pi0", "mean.z0", "sd.z0", "mean.z1","sd.z1")
    attr(new.par, "converged") = iter < niter
    attr(new.par, "iter") = iter
    attr(new.par, "call") = match.call()
    attr(new.par, "lfdr") = ppee
    attr(new.par, "fdr") = fdr
    class(new.par) = "znormix"
    new.par
}

modified_znormix <- function (p, start.pi0=0.85, eps = 1e-03, niter = 1000, verbose = FALSE,min_mean_z1=1,theoretical_null=F){
    z = as.matrix(qnorm(1 - p))
    z[is.infinite(z) & z < 0] = min(z[is.finite(z)])
    z[is.infinite(z) & z > 0] = max(z[is.finite(z)])
    G = length(z)
    stopifnot(G >= 4)
    if (start.pi0 <= 0) 
        start.pi0 = 0.001
    if (start.pi0 >= 1) 
        start.pi0 = 1 - 0.001
    zcut = quantile(z, start.pi0)
    z0.idx = which(z < zcut)
    last.par = c(start.pi0, mean(z[z0.idx]), sd(drop(z[z0.idx])), 
    mean(z[-z0.idx]), sd(drop(z[-z0.idx])))
    iter = 1
    new.par = last.par
    repeat {
        f0 = last.par[1] * dnorm(z, last.par[2], last.par[3])
        ppee = pmin(1 - 1e-06, pmax(1e-06, f0/(f0 + (1 - last.par[1]) * 
            dnorm(z, last.par[4], last.par[5]))))
        new.par[1] = mean(ppee)
        sum.ppee = sum(ppee)
        new.par[4] = crossprod(z, 1 - ppee)/(G - sum.ppee)
        new.par[5] = sqrt(crossprod((z - new.par[4])^2, 1 - ppee)/(G - 
            sum.ppee))
       
        new.par[2] = crossprod(z, ppee)/sum.ppee
        new.par[3] = sqrt(crossprod((z - new.par[2])^2, ppee)/sum.ppee)
        if (abs(new.par[2]) > abs(new.par[4])) {
            tmp = new.par[2]
            new.par[2] = new.par[4]
            new.par[4] = tmp
            tmp = new.par[3]
            new.par[3] = new.par[5]
            new.par[5] = tmp
        }
	  new.par[2] = 0
	  new.par[3] = max(new.par[3],1)
	  if(theoretical_null){
	    new.par[3] = 1
	  }
	  new.par[4] = max(min_mean_z1,new.par[4])
        if (isTRUE(verbose)) 
            cat("iter", iter, "\tparameters=", new.par, "\tmax.diff=", 
                max(abs(new.par - last.par)), fill = TRUE)
        if (iter >= niter || max(abs(new.par - last.par)) < eps) 
            break
        last.par = new.par
        iter = iter + 1
    }
    ord = order(ppee)
    fdr = numeric(G)
    fdr[ord] = cumsum(ppee[ord])/(1:G)
    names(new.par) = c("pi0", "mean.z0", "sd.z0", "mean.z1", 
        "sd.z1")
    attr(new.par, "converged") = iter < niter
    attr(new.par, "iter") = iter
    attr(new.par, "call") = match.call()
    attr(new.par, "lfdr") = ppee
    attr(new.par, "fdr") = fdr
    class(new.par) = "znormix"
    new.par
}

#ps = c(runif(10000),rbeta(500,1,1000))
#ps = c(runif(10000),rbeta(500,1,1000),rep(1,10000))
#est = modified_znormix(ps)
#est = modified_znormix_halfnorm(ps)
#plot(ps,attr(est,"lfdr"))
#plot(ps,attr(est,"fdr"))

# Dani's new model
# Estimate pi.0 f.0 and f.1 for each study with smoothed isotonic regression
est.Pi	<- function(Z.vec,Z.bin,plt = F,smooth.par = 0.9)
{		
	k		<- length(Z.bin)
	z.med	<- (Z.bin[-k] + Z.bin[-1])/2
	f.0		<- diff(pnorm(Z.bin)) / diff(Z.bin)
	bz.vec	<- apply(outer(Z.bin,Z.vec,"<"),2,sum) 
	f.z		<- table(bz.vec) / (5000*diff(Z.bin))
	q.seq	<- seq(0.05,0.95,length = k-1)

	pred.smth.spl	<- predict(smooth.spline(q.seq,log(f.z/f.0),spar = smooth.par))$y
	pred.inc		<- isoreg(q.seq,pred.smth.spl)$yf
	f.z.inc			<- f.0 * exp(pred.inc)
	f.1				<- f.z.inc - f.0*exp(pred.inc[1])
	
	pi.0	<- exp(pred.inc[1])
	pb.0	<- diff(pnorm(Z.bin)) / (pnorm(Z.bin[k]) - pnorm(Z.bin[1]))
	pb.1	<- f.1	* diff(Z.bin) / sum(f.1	* diff(Z.bin) )

	if(plt)
	{
		plot(q.seq,c(log(f.z/f.0)))
		lines(q.seq,pred.smth.spl,col=3,lty=2)
		lines(q.seq,pred.inc,col=3)
		abline(h = 0)

		plot(z.med,f.0,type = "b",col=2)
		lines(z.med,f.z,col=3,lty = 2)
		lines(z.med,f.z.inc,col=3,lty = 1)

		plot(z.med,f.z.inc,col=1,lty = 1,type = "l")
		lines(z.med,f.0*exp(pred.inc[1]),col=2)
		lines(z.med,f.1,col=3)
	}
	return(list(pi.0 = pi.0,pb.0 = pb.0,pb.1 = pb.1,bz.vec = bz.vec))
}
# Estimation of Pi for all studies with repfdr
est.Pi.mat	<- function(Z.mat,Z.bin)
{		
	Pi.1	<- est.Pi(Z.mat[,1],Z.bin,plt = F,smooth.par = 0.75)
	Pi.2	<- est.Pi(Z.mat[,2],Z.bin,plt = F,smooth.par = 0.75)
	Pi.3	<- est.Pi(Z.mat[,3],Z.bin,plt = F,smooth.par = 0.75)

	bz.smp     		<- cbind(Pi.1$bz.vec,Pi.2$bz.vec,Pi.3$bz.vec)
	pbz.smp			<- array(dim=c(3,length(Z.bin) - 1,2))
	pbz.smp[,,1]	<- rbind(Pi.1$pb.0,Pi.2$pb.0,Pi.3$pb.0)
	pbz.smp[,,2]	<- rbind(Pi.1$pb.1,Pi.2$pb.1,Pi.3$pb.1)
	output.repfdr 	<- repfdr(pbz.smp, bz.smp, "meta-analysis",control = em.control(pi.initial = c(0.86,rep(0.02,7)), tol = 1e-10,verbose = F))

	return(list(Pi = output.repfdr$Pi[,4],binned.z = bz.smp,prob.binned.z = pbz.smp))
}
run_isotonic_repfdr<-function(pval_vector){
	Z.mat = qnorm(1-pval_vector)
	Z.mn	<- floor(min(Z.mat)*10)/10
	Z.mx	<- ceiling(max(Z.mat)*10)/10
	Z.bin	<- c(Z.mn, qnorm(seq(0.05,0.98,length = 25)),Z.mx)
	est = est.Pi(Z.mat,Z.bin,plt = T,smooth.par = 0.9)
	pi_0 = est[[1]]
	f0_dist = est[[2]]
	f1_dist = est[[3]]
	zbins = est[[4]]
	f_0 = f0_dist[zbins]
	f_1 = f1_dist[zbins]
	fdr = (f_0 * pi_0) / (f_0 * pi_0 + (1-pi_0)*f_1)
	return(list(pi_0=pi_0,f_0=f_0,f_1=f_1,fdr=fdr))
}

