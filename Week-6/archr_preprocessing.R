#Download the data
if (!dir.exists("data")) {
  dir.create("data")
}
files <- c("GSM4138888_scATAC_BMMC_D5T1.fragments.tsv.gz",
  "GSM4138890_scATAC_CD34_D7T1.fragments.tsv.gz",
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

addArchRGenome("hg19") # set the reference genome
addArchRThreads(threads = cores) # set the cores

#NOTE!! Parallization won't work in Windows!!


