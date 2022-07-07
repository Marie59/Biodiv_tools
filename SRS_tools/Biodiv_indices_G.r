#Rscript

###########################################
##    Mapping alpha and beta diversity   ##
###########################################

#####Packages :  raster,
#                rgdal,
#                sp,
#                rasterdiv,
#                ggplot2,

#####Load arguments

args <- commandArgs(trailingOnly = TRUE)

#####Import the S2 data

if (length(args) < 1) {
    stop("This tool needs at least 1 argument")
}else{
    data_raster <- args[1]
    rasterheader <- args[2]
    alpha <- as.numeric(args[3])
}

########################################################################
##                   COMPUTE BIODIVERSITY INDICES                     ##
########################################################################
# Read raster
copNDVI <- raster::raster(data_raster)

copNDVIlr <- raster::reclassify(copNDVI, cbind(252, 255, NA), right=TRUE)

#Resample using raster::aggregate and a linear factor of 10
copNDVIlr <- raster::aggregate(copNDVIlr, fact = 50)
#Set float numbers as integers to further speed up the calculation
storage.mode(copNDVIlr[]) = "integer"

# Convert raster to SpatialPointsDataFrame
r_pts <- raster::rasterToPoints(copNDVIlr, spatial = T)

# reproject sp object
geo.prj <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0" 
r_pts <- sp::spTransform(r_pts, sp::CRS(geo.prj)) 


# Assign coordinates to @data slot, display first 6 rows of data.frame
r_pts@data <- data.frame(r_pts@data, longitude = sp::coordinates(r_pts)[, 1],
                         latitude = sp::coordinates(r_pts)[, 2])   

##Potting indices 
spectrale_indices <- function(data, indice_choice) {
  graph_indices <- ggplot2::ggplot(data) +
  ggplot2::geom_point(ggplot2::aes_string(x = data[, 2], y = data[, 3], color = data[, indice_choice]), shape = "square", size = 2) + ggplot2::scale_colour_gradient(low = "blue", high = "orange", na.value = "grey50") +
  ggplot2::xlab("Longitude") + ggplot2::ylab("Latitude") +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1), plot.title = ggplot2::element_text(color = "black", size = 12, face = "bold")) + ggplot2::ggtitle(indice_choice)
  
ggplot2::ggsave(paste0(indice_choice, ".png"), graph_indices, width = 12, height = 10, units = "cm")
write.table(data, file = "BiodivIndex.tabular", sep = "\t", dec = ".", na = " ", row.names = F, col.names = T, quote = FALSE)
}


#Shannon's Diversity
sha <- rasterdiv::Shannon(copNDVIlr, window = 9, np = 1)
sha_df <- data.frame(raster::rasterToPoints(sha, spatial = T))
sha_name <- "Shannon"
r_pts@data[, sha_name] <- sha_df[, 1]

#Renyi's Index
ren <- rasterdiv::Renyi(copNDVIlr, window = 9, alpha = alpha, np = 1)
ren_df <- data.frame(raster::rasterToPoints(ren[[1]]))
ren_name <- "Renyi"
r_pts@data[, ren_name] <- ren_df[, 3]

#Berger-Parker's Index
ber <- rasterdiv::BergerParker(copNDVIlr, window = 9, np = 1)
ber_df <- data.frame(raster::rasterToPoints(ber, spatial = T))
ber_name <- "Berger-Parker"
r_pts@data[, ber_name] <- ber_df[, 1]

#Pielou's Evenness
#pie <- rasterdiv::Pielou(copNDVIlr, window = 9, np = 1)
#pie_df <- data.frame(raster::rasterToPoints(pie, spatial = T))
#stop(length(pie_df[, 1]), "et", length(r_pts@data[, 1]))
#pie_name <- "Pielou"
#r_pts@data[, pie_name] <- pie_df[, 1]

#Hill's numbers
hil <- rasterdiv::Hill(copNDVIlr, window = 9, alpha = alpha, np = 1)
hil_df <- data.frame(raster::rasterToPoints(hil[[1]]))
hil_name <- "Hill"
r_pts@data[, hil_name] <- hil_df[, 3]

#Parametric Rao's quadratic entropy with alpha ranging from 1 to 5
prao <- rasterdiv::paRao(copNDVIlr, window = 9, alpha = alpha, dist_m = "euclidean", np = 1)
prao_df <- data.frame(raster::rasterToPoints(prao$window.9[[1]]))
prao_name <- "Prao"
r_pts@data[, prao_name] <- prao_df[, 3]

#Cumulative Residual Entropy
cre <- rasterdiv::CRE(copNDVIlr, window = 9, np = 1)
cre_df <- data.frame(raster::rasterToPoints(cre, spatial = T))
cre_name <- "CRE"
r_pts@data[, cre_name] <- cre_df[, 1]

## Plotting all the graph and writing a tabular
list_indice <- list("Shannon", "Renyi", "Berger-Parker", "Hill", "Prao", "CRE")
for (indice in list_indice) {
  spectrale_indices(data = r_pts@data, indice_choice = indice)
}

write.table(r_pts@data, file = "BiodivIndex.tabular", sep = "\t", dec = ".", na = " ", row.names = F, col.names = T, quote = FALSE)



