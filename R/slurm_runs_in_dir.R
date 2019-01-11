
get_sh_prefix_one_node_specify_cpu_and_mem<-function(err="",log="",Ncpu=4,mem_size=32000,time="24:00:00"){
  return(
    c(
      "#!/bin/bash",
      "#",
      paste("#SBATCH --time=",time,sep=""),
      "#SBATCH --partition=euan,mrivas,normal,owners",
      "#SBATCH --nodes=1",
      paste("#SBATCH -c",Ncpu),
      paste("#SBATCH --mem=",mem_size,sep=""),
      paste("#SBATCH --error",err),
      paste("#SBATCH --out",log),
      "#SBATCH -x sh-113-15",
      "",
      "module load R"
    )
  )
}
print_sh_file<-function(path,prefix,cmd){
  cmd = c(prefix,cmd)
  write.table(file=path,t(t(cmd)),row.names = F,col.names = F,quote = F)
}
run_sbatch_rscript_command<-function(cmd,out_path,name,batch_script_func=get_sh_default_prefix,...){
  err_path = paste(out_path,name,".err",sep="")
  log_path = paste(out_path,name,".log",sep="")
  curr_sh_file = paste(out_path,name,".sh",sep="")
  batch_script_prefix = batch_script_func(err_path,log_path,...)
  print_sh_file(curr_sh_file,batch_script_prefix,cmd)
  system(paste("sbatch",curr_sh_file))
}

curr_dir = "/home/users/davidama/motrpac_metaanalysis/"
setwd(curr_dir)
pvals_files = list.files()
pvals_files = pvals_files[grepl("pvals.txt",pvals_files)]
for (ff in pvals_files){
  for(m in c("bum","znormix")){
    currname = gsub("_pvals.txt","",ff)
    currname = paste(currname,"_",m,sep="")
    outfile = paste(curr_dir,currname,".txt",sep="")
    curr_cmd = paste("Rscript /home/users/davidama/repos/screen/R/cmd_line_runnable.R",
                     ff,"1",m,outfile,"path=/home/users/davidama/R/packages/ emEps=1e-5 nH=10000")
    if(m=="bum"){
      curr_cmd = paste("Rscript /home/users/davidama/repos/screen/R/cmd_line_runnable.R ",
                       ff,"1",m,outfile,"path=/home/users/davidama/R/packages/ emEps=1e-3 nH=10000")
    }
    run_sbatch_rscript_command(curr_cmd,curr_dir,currname,get_sh_prefix_one_node_specify_cpu_and_mem)
  }
}




