#Rscript

###########################################
##    Mapping alpha and beta diversity   ##
###########################################

#####Packages : stars
#               utils
#               biodivmapr
#               raster
#               sf
#               mapview
#               leafpop
#               RColorBrewer
#               labdsv
#               rgdal
#               ggplot2
#               gridExtra
remotes::install_github('jbferet/biodivMapR')
#####Load arguments

args <- commandArgs(trailingOnly = TRUE)

#####Import the S2 data

if (length(args) < 1) {
    stop("This tool needs at least 1 argument")
}else{
    data_raster <- args[1]
    rasterheader <- args[2]
    text_compo <- args[3]
    plots_zip <- args[4]
    type <- as.character(args[5])
}

################################################################################
##              DEFINE PARAMETERS FOR DATASET TO BE PROCESSED                 ##
################################################################################
# expected to be in ENVI HDR  

Input_Image_File <- file.path(getwd(), data_raster, fsep = "/")
Input_Header_File <- file.path(getwd(), rasterheader, fsep = "/")
# path for the Mask raster corresponding to image to process
# expected to be in ENVI HDR format, 1 band, integer 8bits
# expected values in the raster: 0 = masked, 1 = selected
# set to FALSE if no mask available
Input_Mask_File <- FALSE

# relative or absolute path for the Directory where results will be stored
# For each image processed, a subdirectory will be created after its name
Output_Dir  = 'RESULTS'

# SPATIAL RESOLUTION
# resolution of spatial units for alpha and beta diversity maps (in pixels), relative to original image
# if Res.Map = 10 for images with 10 m spatial resolution, then spatial units will be 10 pixels x 10m = 100m x 100m surfaces
# rule of thumb: spatial units between 0.25 and 4 ha usually match with ground data
# too small window_size results in low number of pixels per spatial unit, hence limited range of variation of diversity in the image
window_size <- 10

# PCA FILTERING: Set to TRUE if you want second filtering based on PCA outliers to be processed. Slower
FilterPCA <- TRUE

# type of PCA:
# PCA: no rescaling of the data
# SPCA: rescaling of the data
TypePCA <-'SPCA'


################################################################################
##                    DEFINE PARAMETERS FOR METHOD                            ##
################################################################################
nbCPU <- 4
MaxRAM <- 0.5
nbclusters <- 50

################################################################################
##                              PROCESS IMAGE                                 ##
################################################################################
# 1- Filter data in order to discard non vegetated / shaded / cloudy pixels
NDVI_Thresh <- 0.5
Blue_Thresh <- 500
NIR_Thresh  <- 1500
Continuum_Removal <- TRUE

print("PERFORM RADIOMETRIC FILTERING")
ImPathShade <- biodivMapR::perform_radiometric_filtering(
  Image_Path = Input_Image_File, Mask_Path = Input_Mask_File, Output_Dir = Output_Dir, NDVI_Thresh = NDVI_Thresh, 
  Blue_Thresh = Blue_Thresh, NIR_Thresh = NIR_Thresh)

print("PERFORM PCA ON RASTER")
PCA_Output <- biodivMapR::perform_PCA(Input_Image_File = Input_Image_File, Input_Mask_File = Input_Mask_File,
                          Output_Dir = Output_Dir, TypePCA = TypePCA, FilterPCA = FilterPCA, nbCPU = nbCPU, MaxRAM = MaxRAM)

PCA_Files <- PCA_Output$PCA_Files
Pix_Per_Partition <- PCA_Output$Pix_Per_Partition
nb_partitions <- PCA_Output$nb_partitions
Input_Mask_File <- PCA_Output$ImPathShade
PCA_model <- PCA_Output$PCA_model
SpectralFilter <- PCA_Output$SpectralFilter
# path for the updated mask
Input_Mask_File <- PCA_Output$MaskPath

# 3- Select principal components from the PCA raster
# Select components from the PCA/SPCA/MNF raster
Image_Name <- tools::file_path_sans_ext(basename(Input_Image_File))
Output_Dir_Full <- file.path(Output_Dir, Image_Name, TypePCA, "PCA")
data_components <- read.table(text_compo, sep = "\t", dec = ".", fill = TRUE, encoding = "UTF-8")
write.table(data_components, paste0(Output_Dir_Full, "/Selected_Components.txt"))

Sel_PC <-  file.path(Output_Dir_Full, "Selected_Components.txt")

################################################################################
##                      MAP ALPHA AND BETA DIVERSITY                          ##
################################################################################
print("MAP SPECTRAL SPECIES")

Kmeans_info <- biodivMapR::map_spectral_species(Input_Image_File = Input_Image_File, Output_Dir = Output_Dir, PCA_Files = PCA_Files, Input_Mask_File = Input_Mask_File, Pix_Per_Partition = Pix_Per_Partition, nb_partitions = nb_partitions, nbCPU = nbCPU, MaxRAM = MaxRAM, nbclusters = nbclusters, TypePCA = TypePCA)

if (type == "alpha" | type == "all") {
  print("MAP ALPHA DIVERSITY")
  Index_Alpha <- c('Shannon')
  alpha_div <- biodivMapR::map_alpha_div(Input_Image_File = Input_Image_File, Output_Dir = Output_Dir, TypePCA = TypePCA, window_size = window_size, nbCPU = nbCPU, MaxRAM = MaxRAM, Index_Alpha = Index_Alpha, nbclusters = nbclusters)

alpha_path <- file.path(Output_Dir, Image_Name, TypePCA, "ALPHA", "Shannon_10_Fullres.zip")
dir.create("alpha_fold")
unzip(alpha_path, exdir = "alpha_fold")
alpha_data <- list.files("alpha_fold")
alpha_raster <- raster::raster(alpha_data)
get_alpha <- raster::rasterToPoints(alpha_raster, spatial = T)
stop(head(get_alpha))
#get_alpha <- grep(pattern = paste0(window_size, '_Fullres.zip'), x = alpha_path)

#alpha_plot <- rasterVis::levelplot(alpha_data,layout=c(0,1,1), main="alpha")
}
#Create a random raster layer

#alpha_map <- ggplot2::ggplot(data = rdf)+
 # ggplot2::geom_raster(mapping = ggplot2::aes(x = x, y = y))
#mapper_alpha_div <- ggplot2::ggsave("Alpha_diversity.png", alpha_map, scale = 0.65, width = 12, height = 9, units = "in", dpi = 200, limitsize = TRUE)
#mapper_alpha <- plot(list.files(file.path("alpha_beta")))
#save(mapper_alpha, file = "Alpha_map.png")

##### Trouver le moyen de mapview le truc la 

if (type == "beta" | type == "all") {
  print("MAP BETA DIVERSITY")
  beta_div <- biodivMapR::map_beta_div(Input_Image_File = Input_Image_File, Output_Dir = Output_Dir, TypePCA = TypePCA, window_size = window_size, nb_partitions = nb_partitions, nbCPU = nbCPU, MaxRAM = MaxRAM, nbclusters = nbclusters)
} 


################################################################################
##          COMPUTE ALPHA AND BETA DIVERSITY FROM FIELD PLOTS                 ##
################################################################################
## read selected features from dimensionality reduction 
Selected_Features <- read.table(Sel_PC)[[1]]
## path for selected components

if (type == "funct" | type == "all") {
mapper <- biodivMapR::map_functional_div(Original_Image_File = Input_Image_File, Functional_File = PCA_Files,  Selected_Features = Selected_Features, Output_Dir = Output_Dir, window_size = window_size, nbCPU = nbCPU, MaxRAM = MaxRAM, TypePCA = TypePCA)
}
# location of the directory where shapefiles used for validation are saved
dir.create("VectorDir")
unzip(plots_zip, exdir = "VectorDir")

#Path_Vector <- list.files("VectorDir")

#VectorDir <- destunz
# list vector data
Path_Vector <- biodivMapR::list_shp("VectorDir")
Name_Vector <- tools::file_path_sans_ext(basename(Path_Vector))


# location of the spectral species raster needed for validation
Path_SpectralSpecies <- Kmeans_info$SpectralSpecies
# get diversity indicators corresponding to shapefiles (no partitioning of spectral dibversity based on field plots so far...)

Biodiv_Indicators <- biodivMapR::diversity_from_plots(Raster_SpectralSpecies = Path_SpectralSpecies, Plots = Path_Vector, nbclusters = nbclusters, Raster_Functional = PCA_Files, Selected_Features = Selected_Features)

Shannon_RS <- c(Biodiv_Indicators$Shannon)[[1]]
FRic <- c(Biodiv_Indicators$FunctionalDiversity$FRic)
FEve <- c(Biodiv_Indicators$FunctionalDiversity$FEve)
FDiv <- c(Biodiv_Indicators$FunctionalDiversity$FDiv)
# if no name for plots
Biodiv_Indicators$Name_Plot = seq(1, length(Biodiv_Indicators$Shannon[[1]]), by = 1)


####################################################
# write RS indicators
####################################################
# write a table for Shannon index
if (type == "alpha" | type == "comparison" | type == "all") {
write.table(Shannon_RS, file = "ShannonIndex.tabular", sep = "\t", dec = ".", na = " ", row.names = Biodiv_Indicators$Name_Plot, col.names = c("Shannon_Index"), quote = FALSE)
}

if (type == "funct" | type == "comparison" | type == "all") {            
# write a table for all spectral diversity indices corresponding to alpha diversity
Results <- data.frame(Name_Vector, Biodiv_Indicators$Richness, Biodiv_Indicators$Fisher,
                      Biodiv_Indicators$Shannon, Biodiv_Indicators$Simpson,
                      Biodiv_Indicators$FunctionalDiversity$FRic,
                      Biodiv_Indicators$FunctionalDiversity$FEve,
                      Biodiv_Indicators$FunctionalDiversity$FDiv)

names(Results)  = c("ID_Plot", "Species_Richness", "Fisher", "Shannon", "Simpson", "FRic", "FEve", "FDiv")
write.table(Results, file = "AlphaDiversity.tabular", sep = "\t", dec = ".", na = " ", row.names = F, col.names = T, quote = FALSE)
}

if (type == "beta" | type == "comparison" | type == "all") {
# write a table for Bray Curtis dissimilarity
BC_mean <- Biodiv_Indicators$BCdiss
colnames(BC_mean) <- rownames(BC_mean) <- Biodiv_Indicators$Name_Plot
write.table(BC_mean, file = "BrayCurtis.tabular", sep = "\t", dec = ".", na = " ", row.names = F, col.names = T, quote = FALSE)
}

if (type == "comparison" | type == "all") {
####################################################
# illustrate results
####################################################
# apply ordination using PCoA (same as done for map_beta_div)

MatBCdist <- as.dist(BC_mean, diag = FALSE, upper = FALSE)
BetaPCO <- labdsv::pco(MatBCdist, k = 3)

# assign a type of vegetation to each plot, assuming that the type of vegetation 
# is defined by the name of the shapefile         
 
nbSamples <- shpName <- c()
for (i in 1:length(Path_Vector)){
  shp <- Path_Vector[i]
  nbSamples[i] <- length(rgdal::readOGR(shp,verbose = FALSE))
  shpName[i] <- tools::file_path_sans_ext(basename(shp))
}

Type_Vegetation = c()
for (i in 1: length(nbSamples)){
  for (j in 1:nbSamples[i]){
    Type_Vegetation = c(Type_Vegetation,shpName[i])
  }
}

#data frame including a selection of alpha diversity metrics and beta diversity expressed as coordinates in the PCoA space
Results <- data.frame('vgtype' = Type_Vegetation, 'pco1' = BetaPCO$points[,1], 'pco2' = BetaPCO$points[,2], 'pco3' = BetaPCO$points[,3], 'shannon' = Shannon_RS, 'FRic' = FRic, 'FEve' = FEve, 'FDiv' = FDiv)
                      
#plot field data in the PCoA space, with size corresponding to shannon index
g1 <- ggplot2::ggplot(Results, ggplot2::aes(x = pco1, y = pco2, color = vgtype, size = shannon)) + ggplot2::geom_point(alpha = 0.6) + ggplot2::scale_color_manual(values = c("#e6140a", "#e6d214", "#e68214", "#145ae6"))

g2 <- ggplot2::ggplot(Results, ggplot2::aes(x = pco1, y = pco3, color = vgtype, size = shannon)) + ggplot2::geom_point(alpha = 0.6) + ggplot2::scale_color_manual(values = c("#e6140a", "#e6d214", "#e68214", "#145ae6"))

g3 <- ggplot2::ggplot(Results, ggplot2::aes(x = pco2, y = pco3, color = vgtype, size = shannon)) + ggplot2::geom_point(alpha = 0.6) + ggplot2::scale_color_manual(values = c("#e6140a", "#e6d214", "#e68214", "#145ae6"))                     

#extract legend                 
get_legend <- function(a.gplot) {
    tmp <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)
}

legend <- get_legend(g3)
gAll <- gridExtra::grid.arrange(gridExtra::arrangeGrob(g1 + ggplot2::theme(legend.position = "none"), g2 + ggplot2::theme(legend.position = "none"), g3 + ggplot2::theme(legend.position = "none"), nrow=1), legend, nrow = 2, heights = c(3, 2)) 


filename <- ggplot2::ggsave("BetaDiversity_PcoA1_vs_PcoA2_vs_PcoA3.png", gAll, scale = 0.65, width = 12, height = 9, units = "in", dpi = 200, limitsize = TRUE)

filename
}
