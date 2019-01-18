args = commandArgs(trailingOnly=TRUE)

if(length(args) < 4){
  print("usage: ")
  print("1. A path to a matrix of p-values")
  print("2. 1 or 0, 1 means studies are not assumed independent, 0 otherwise")
  print("3. The name of the method for two-groups estimation: locfdr, znormix, or bum")
  print("4. A name of an output path")
  print("5. (OPTIONAL): path=p, where p is a path for a dir with the required R packages:")
  print("igraph, locfdr,BioNet,devtools")
  print("6. (OPTIONAL): nH=n, where n is the maximal number of allowed configurations for SCREN")
  print("7. (OPTIONAL): emEps=e, where e is the convergence parameter for the EM (default=1e-6)")
  print("8. (OPTIONAL): minP=x, where x is the minimal allowed p-value, default is 0")
  q("no")
}

pvals = as.matrix(read.delim(args[1],row.names = 1,header=T,check.names = F))
pvals[is.na(pvals)] = 0.5
model_study_dep = args[2]=="1"
lfdr_method = args[3]
out_path = args[4]
libs = c("igraph","BiocGenerics","graph","kernlab","RBGL","locfdr","BioNet")

emEps = 1e-6
nH=10000
rpath=NULL
use_power=T
minP = 0

if(length(args > 4)){
  optional_parms = args[5:length(args)]
  if(any(grepl("path=",optional_parms))){
    ind = which(grepl("path=",optional_parms))
    rpath = gsub("path=",replacement = "",x=optional_parms[ind])
  }
  if(any(grepl("nH=",optional_parms))){
    ind = which(grepl("nH=",optional_parms))
    nH = as.numeric(gsub("nH=",replacement = "",x=optional_parms[ind]))
  }
  if(any(grepl("emEps=",optional_parms))){
    ind = which(grepl("emEps=",optional_parms))
    emEps = as.numeric(gsub("emEps=",replacement = "",x=optional_parms[ind]))
  }
  if(any(grepl("minP=",optional_parms))){
    ind = which(grepl("minP=",optional_parms))
    minP = as.numeric(gsub("minP=",replacement = "",x=optional_parms[ind]))
    use_power=F
  }
}

nH = min(nH,2^(ncol(pvals)+1))

if(!is.null(rpath)){
  sapply(libs,library,lib.loc=rpath,character.only = T)
}
if(is.null(rpath)){
  sapply(libs,library,character.only = T)
}

source("https://raw.githubusercontent.com/david-dd-amar/screen/master/R/SCREEN.R")
source("https://raw.githubusercontent.com/david-dd-amar/screen/master/R/twogroups_methods.R")

ks = 2:ncol(pvals)
use_power = use_power && (lfdr_method != "bum")

print("Running SCREEN with the following parameters:")
print(paste("input_matrix:",args[1]))
print(paste("Model study dependence:",model_study_dep))
print(paste("Two groups method:",lfdr_method))
print(paste("Convergence epsilon:",emEps))
print(paste("nH:",nH))

if(model_study_dep){
  res = SCREEN(pvals,ks=ks,lfdr_method = lfdr_method,emEps=emEps,nH=nH,use_power = use_power)
}
if(!model_study_dep){
  res = SCREEN_ind(pvals,ks=ks,lfdr_method = lfdr_method,use_power = use_power)
}

res = cbind(rownames(pvals),res)
write.table(res,sep="\t",row.names=F,col.names = T,quote=F,file=out_path)





