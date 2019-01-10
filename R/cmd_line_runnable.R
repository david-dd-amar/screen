args = commandArgs(trailingOnly=TRUE)

if(length(args) < 4){
  print("usage: ")
  print("1. A path to a matrix of p-values")
  print("2. 1 or 0, 1 means studies are not assumed independent, 0 otherwise")
  print("3. The name of the method for two-groups estimation: locfdr, znormix, or bum")
  print("4. A name of an output path")
  print("5. (OPTIONAL): a path for a dir with the required R packages:")
  print("igraph, locfdr,BioNet")
  q("no")
}

pvals = read.delim(args[1],row.names = 1,header=T)
is_ind = args[2]=="0"
lfdr_method = args[3]
out_path = args[4]
libs = c("igraph","locfdr","BioNet")

if(length(args)==5){
  rpath = args[5]
  sapply(libs,library,lib.loc=rpath,character.only = T)
}
if(length(args)<5){
  sapply(libs,library,character.only = T)
}
