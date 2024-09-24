#' @title sci_fromfastq
#' @description From fastq to gene expression matrix for SCI data.
#' @param group, a character string. Two options: sudo or docker, depending to which group the user belongs
#' @param folder, a character string indicating the path of the working folder containing the input files
#' @param sample.name, a character string indicating the name of the experiment
#' @param UMI.cutoff, integer indicating the minimum number of UMI per cell to consider the cell valid
#'
#' @author Maria Luisa Ratto, marialuisa.ratto [at] unito [dot] it, UNITO
#'
#' @return a gene expression matrix + QC plots and stats
#'
#' @examples
#'\dontrun{
#'
#'sci_fromfastq(
#'  group="docker",
#'  folder="/20tb/ratto/catcheR/tomatrix/", 
#'  sample.name="H001AS8", UMI.cutoff=500)
#'
#' @export

sci_fromfastq <- function(
  group=c("docker","sudo"),
  folder, sample.name, UMI.cutoff){
  
  #running time 1
  ptm <- proc.time()
  #setting the folder as working folder
  if (!file.exists(folder)){
    cat(paste("\nIt seems that the ",folder, " folder does not exist\n"))
    return(2)
  }
  
  #storing the position of the home folder
  home <- getwd()
  setwd(folder)
  #initialize status
  #print("echo 0 > ExitStatusFile")
  
  #testing if docker is running
  test <- dockerTest()
  if(!test){
    cat("\nERROR: Docker seems not to be installed in your system\n")
    print("echo 10 > ExitStatusFile")
    setwd(home)
    return(10)
  }
  
  #check  if scratch folder exist
  #  if (!file.exists(scratch.folder)){
  #    cat(paste("\nIt seems that the ",scratch.folder, " folder does not exist\n"))
  #    system("echo 3 > ExitStatusFile 2>&1")
  setwd(folder)
  #   return(3)
  # }
  #  tmp.folder <- gsub(":","-",gsub(" ","-",date()))
  scrat_tmp.folder=folder
  writeLines(scrat_tmp.folder,paste(folder,"/tempFolderID", sep=""))
  # cat("\ncreating a folder in scratch folder\n")
  #dir.create(file.path(scrat_tmp.folder))
  #preprocess matrix and copying files
  if (!dir.exists(paste(folder,"/fastq",sep=""))){
    cat(paste("\n It Seems that fastaq read1 file is not in ",folder,"\n"))
    print("echo 3 > ExitStatusFile 2>&1 &")
    setwd(folder)
    return(3)
  }
  
  #executing the docker job
  run_in_docker(
    image_name = "docker.io/repbioinfo/sci_tomatrix_genome",
    volumes = list(
      c(folder, "/data/scratch/")
    ),
    additional_arguments = c(
      "/home/tomatrix.sh",
      #"/data/scratch/",
      sample.name,
      UMI.cutoff
    )
  )
  
  run_in_docker(
    image_name = "docker.io/repbioinfo/monocle",
    volumes = list(
      c(folder, "/data/scratch/")
    ),
    additional_arguments = c(
      "Rscript /home/endpipeline.R"
    )
  )
  
  #waiting for the end of the container work
  # if(resultRun==0){
  #   cat("\nData filtering is finished\n")
  # }
  
  #saving log and removing docker container
  #container.id <- readLines(paste(folder,"/dockerID", sep=""), warn = FALSE)
  #system(paste("docker logs ", substr(container.id,1,12), " >& ",folder,"/", substr(container.id,1,12),".log", sep=""))
  system(paste("docker logs ", substr(container.id,1,12), " > ",folder,"/", substr(container.id,1,12),".log 2>&1", sep=""))
  system(paste("docker rm ", container.id, sep=""))
  
  
  #removing temporary folder
  cat("\n\nRemoving the temporary file ....\n")
  # system(paste("rm -R ",scrat_tmp.folder))
  #file.remove(paste0(folder,"out.info"))
  #file.remove(paste0(folder,"dockerID"))
  #file.remove(paste0(folder,"tempFolderID"))
  #system(paste("cp ",paste(path.package(package="rCASC"),"containers/containers.txt",sep="/")," ",data.folder, sep=""))
  setwd(home)
} 

