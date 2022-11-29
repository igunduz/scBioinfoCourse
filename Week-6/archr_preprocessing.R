#Download the data
if (!dir.exists("data")) {
  dir.create("data")
}
files <- c("GSM4138888_scATAC_BMMC_D5T1.fragments.tsv.gz",
  "wgeGSM4138890_scATAC_CD34_D7T1.fragments.tsv.gz",
  "GSM4138891_scATAC_CD34_D8T1.fragments.tsv.gz")

urls <- c("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM4138888&format=file&file=GSM4138888%5FscATAC%5FBMMC%5FD5T1%2Efragments%2Etsv%2Egz",
          "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM4138890&format=file&file=GSM4138890%5FscATAC%5FCD34%5FD7T1%2Efragments%2Etsv%2Egz",
          "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM4138891&format=file&file=GSM4138891%5FscATAC%5FCD34%5FD8T1%2Efragments%2Etsv%2Egz")
names(urls) <- files

if (!all(file.exists(files))) {
  for (file in files) {
    download.file(urls[[file]], destfile = paste0("data/",file,".tar"))  # Method = curl to avoid corruption
    untar(paste0("data/",file,".tar", exdir = paste0("data/",file)))
  }
}

#call libraries
library(ArchR)
library(dplyr)

#Link to ArchR file (the most updated version,the book is not updated yet)
#https://www.archrproject.com/index.html

addArchRGenome("hg19") # set the reference genome
addArchRThreads(threads = cores) # set the cores
#NOTE!! Parallization won't work in Windows!!

#list the files
atac_files <- list.files("data",pattern = "*.tsv.gz", full.names = TRUE)

#name the files
names(atac_files) <- c("scATAC_BMMC_D5T1","scATAC_CD34_D7T1", "scATAC_CD34_D8T1")

#create the arrow files
ArrowFiles <- createArrowFiles(
  inputFiles = atac_files,
  sampleNames = names(atac_files), #TSS enrichment scores and nucleosome info
  minTSS = 4, #Dont set this too high because you can always increase later
  minFrags  = 1000, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE)

#show the files
ArrowFiles

#identify doublets
doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, 
  knnMethod = "UMAP", 
  LSIMethod = 1
)

#create an ArchR project
project <- ArchR::ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "data/project/",
  copyArrows = TRUE, #save arrow files as arrow
  showLogo = FALSE  
)

#Once you created an ArchR project, you can load it via this function
#project <- ArchR::loadArchRProject(path = "data/project/",showLogo = F)


project

#How many cells does your project include?
#What is the median TSS-value and the median of the number of fragments?
#What are the dimensions of your dataset?

#To save ArchR project
saveArchRProject(ArchRProj = project, outputDirectory = "data/project/", load = FALSE)