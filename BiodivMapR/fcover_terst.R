#Rscript

#######################################
##    fraction of vegetation cover   ##
#######################################

#####Packages : RStoolbox
#               mapview

#####Load arguments

library(RStoolbox)
library(mapview)

#####Import data
#setwd("/home/pndb/Scriptr_MJ/data/L1A/L1A/localisation/shape/")

#table1 <- read.table("data_FrenchBBS.tabular", sep = "\t", dec = ".", header = TRUE)
normlas.file = system.file("extdata", "lidar_example.laz", package="leafR")
#colnames(table1) <- c("site", "annee", "espece", "abond")

#####Your analysis
  library(raster)
  library(caret)
  ## Create fake input images
  data(rlogo)
  lsat <- rlogo
  agg.level <- 9
  modis <- aggregate(lsat, agg.level)
  
  ## Perform classification
  lc      <- unsuperClass(lsat, nClass=2)
  
  ## Calculate the true cover, which is of course only possible in this example, 
  ## because the fake corse resolution imagery is exactly res(lsat)*9
  trueCover <- aggregate(lc$map, agg.level, fun = function(x, ...){sum(x == 1, ...)/sum(!is.na(x))})
  
  ## Run with randomForest and support vector machine (radial basis kernel)
  ## Of course the SVM is handicapped in this example due to poor tuning (tuneLength)
  par(mfrow=c(2,3))
    fc <- fCover(
      classImage = lc$map ,
      predImage = modis,
      classes=1,
      nSample = 50,
      number = 5,
      tuneLength=2
    )        
mapview(fc)

  
