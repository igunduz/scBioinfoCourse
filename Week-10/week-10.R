library(Seurat)
library(dplyr)
library(SCDC)

set.seed(12)

# read in image file
image1 <-  Read10X_Image("Section_1/spatial/")

#read the data
spdata <- Load10X_Spatial(data.dir = "Section_1/",
                          filename = "V1_Mouse_Brain_Sagittal_Posterior_filtered_feature_bc_matrix.h5",
                          slice = "slice1",image = image1)
#set the ident
spdata <- SetIdent(spdata, value = "slice1")


#read slice 2
image2 <-  Read10X_Image("Section_2/spatial/")
spdata2 <- Load10X_Spatial(data.dir =  "Section_2/",
                           filename = "V1_Mouse_Brain_Sagittal_Posterior_Section_2_filtered_feature_bc_matrix.h5",
                           slice = "slice2",image = image2)
spdata2 <- SetIdent(spdata2, value = "slice2")


#access the images
spdata@images

#access the expression data
spdata@assays

#add the percentage of mitochondrial features
spdata$percent.mito <- PercentageFeatureSet(spdata,pattern = "^mt-")

#overlaying the per-spot values for gene count
SpatialFeaturePlot(spdata, "nCount_Spatial")

#overlaying the per-spot values for mitochondrial  features
SpatialFeaturePlot(spdata, "percent.mito")

#pick random gene to visualize
SpatialFeaturePlot(spdata, features = c("Hpca", "Ttr"),
                   pt.size.factor = 1, alpha = c(0.1, 1))

#plot violin plots
VlnPlot(spdata, features = c("nCount_Spatial", "nFeature_Spatial","percent.mito"),
        pt.size = 0.1, ncol = 2) +
  NoLegend()

#filtering based on cut-off determined by violin plots
spdata <- subset(spdata, nFeature_Spatial >= 200)

#INTERMEDIATE STEPS ARE SKIPPED HERE!!!

# maxSize for PrepSCTIntegration
options(future.globals.maxSize = 2000 * 1024^2)

#THÝS MAY BE RELEVANT FOR SOME OF YOU

#Merge the dataset
experiment.merged <- merge(spdata, spdata2)
