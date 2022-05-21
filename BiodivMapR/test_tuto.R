#Rscript
library(dplyr)
library(readr)
library(raster)
library(sf)

##########################
##    GIS with R test   ##
##########################
#####Source :https://bookdown.org/michael_bcalles/gis-crash-course-in-r/data.html

#####Packages : sf
#               raster
#               crsuggest
#               osmdata
#               cancensus
#               readr
#               dplyr

#####Load arguments

###gEOGRAPHIC VECTOR AND RASTER DATA
#Point

a_single_point <- st_point(x = c(1, 3))

attributes(a_single_point)

#create five sfg point objects
point1 <- st_point(x = c(1,3)) 
point2 <- st_point(x = c(2, 4))
point3 <- st_point(x = c(3, 3))
point4 <- st_point(x = c(4, 3))
point5 <- st_point(x = c(5, 2))

points <- st_sfc(point1,point2,point3,point4,point5) # create a single sfc points object

points_wgs <- st_sfc(point1,point2,point3,point4,point5,crs=4326) # create a single sfc points object, and define the CRS

attributes(points)
plot(points_wgs,col="red",pch=16) # map the sfc object using the plot function

points_attribute_data <- data.frame(transport_mode = c("Bicycle","Pedestrian","Motor Vehicle","Motor Vehicle","Motor Vehicle"))

points_sf <- st_sf(points_attribute_data,geometry=points)

plot(points_sf,pch=16)

#LINES
a_single_line_matrix <- rbind(c(1,1), 
                              c(2, 4), 
                              c(3, 3), 
                              c(4, 3), 
                              c(5,2))# create a matrix with 

a_single_line <-  st_linestring(a_single_line_matrix)
plot(a_single_line,col="darkblue")

line1 <- st_linestring(rbind(c(1,1), 
                             c(2,4), 
                             c(3,3), 
                             c(4,3), 
                             c(5,2)))

line2 <- st_linestring(rbind(c(3,3), 
                             c(3,5), 
                             c(4,5)))

line3 <- st_linestring(rbind(c(2,4), 
                             c(-1,4), 
                             c(-2,-2)))

lines <- st_sfc(line1,line2,line3,crs=4326)
plot(lines)

line_attribute_data <- data.frame(road_name = c("Vinmore Avenue","Williams Road","Empress Avenue"), 
                                  speed_limit = c(30,50,40))

lines_sf <- st_sf(line_attribute_data,geometry=lines)

plot(lines_sf[1])

#POLYGONS
a_polygon_vertices_matrix <- rbind(c(1,1), c(2, 4), c(3, 3), c(4, 3), c(5,2),c(1,1))
a_polygon_list = list(a_polygon_vertices_matrix)

a_polygon <-  st_polygon(a_polygon_list)

lake_vertices_matrix <- rbind(c(1.5, 1.5),c(1.5,1.75),c(1.75,1.75),c(1.75,1.5),c(1.5,1.5))
lake_vertices_matrix

a_single_polygon_with_a_hole_list = list(a_polygon_vertices_matrix,lake_vertices_matrix)
a_single_polygon_with_a_hole_list

a_single_polygon_with_a_hole <-  st_polygon(a_single_polygon_with_a_hole_list)

park1 <- a_single_polygon_with_a_hole
park2 <- st_polygon(list(
  rbind(
    c(6,6),c(8,7),c(11,9), c(10,6),c(8,5),c(6,6)
  )
))

park_attributes <- data.frame(park_name = c("Gilmore","Dixon"))

parks_sf <- st_sf(park_attributes, geometry = st_sfc(park1,park2,crs=4326))

plot(parks_sf)

#RASTER
van_raster <- raster(resolution=1000, #specify spatial resolution
                     xmn=483691, xmx=498329, #specify extent from east to west in metres
                     ymn=5449535, ymx=5462381, #specify extent from north to south in metres
                     crs=26910) #set projection

number_cells <- ncell(van_raster) #the number of grid cells in the RasterLayer object

cell_values <- runif(number_cells,min = -1,max=1)

values(van_raster) <- cell_values # assign each cell a value from the random generation

plot(van_raster)

#READ SPATIAL DATA
icbc_crash <- read_csv("/home/pndb/Scriptr_MJ/GIS_with_R/data/ICBC_cyclist_crashes.csv") # use the file location of where you have stored your data
str(icbc_crash)#look at data frame

is.na(icbc_crash$Latitude) %>% summary() #check for missing values in latitude
is.na(icbc_crash$Longitude) %>% summary() #check for missing values in longitude

icbc_crash <- icbc_crash %>% 
  filter(!is.na(Latitude)) # remove the records which have missing values

is.na(icbc_crash$Latitude) %>% summary()
is.na(icbc_crash$Longitude) %>% summary()

icbc_crash_sf <- st_as_sf(icbc_crash, coords = c("Longitude","Latitude"), crs = 4326)

str(icbc_crash_sf)
plot(icbc_crash_sf["Cyclist Flag"])
#st_crs(icbc_crash_sf) <- 4326

#possible_crs <- suggest_crs(icbc_crash_sf,units = "m",limit = 20)

#possible_crs <- tibble( y = icbc_crash_sf)#View the the top 10 suggestions
#possible_crs

icbc_crash_sf_proj <- st_transform(icbc_crash_sf,crs=3005)

st_crs(icbc_crash_sf_proj)

icbc_crash_van <- icbc_crash_sf %>% 
  filter(`Municipality Name`=="VANCOUVER")

#possible_crs_van <- suggest_crs(icbc_crash_van)

icbc_crash_van_proj <- st_transform(icbc_crash_van,crs=26910)

plot(icbc_crash_van_proj["Cyclist Flag"])


van_boundary <- read_sf("/home/pndb/Scriptr_MJ/GIS_with_R/data/CoV_boundary.shp") # use the file location of where you have stored your data
van_boundary
plot(van_boundary["Type"])

library(osmdata)

van_bb <- getbb("Vancouver BC")

#define 
opq_van <- opq(bbox = van_bb)

osm_van <- add_osm_feature(opq = opq_van,
                           key = "highway" 
)

osm_van_sf <- osmdata_sf(osm_van)

osm_van_network <- osm_van_sf$osm_lines

#osm_van_network <- st_transform(osm_van_network,crs=26910) # project data
st_crs(osm_van_network) <- 26910

plot(osm_van_network["highway"],key.pos = 1)

library(cancensus)

options(cancensus.api_key = "CensusMapper_42c0d61ae9af82d716ad8bab03d031eb")

find_census_vectors("Main mode of Commute",#search term
                    type = "total", # pick one of all, total, male or female in terms of counts
                    dataset = "CA16", #specify the 2016 census
                    query_type = "semantic" 
                    
)
list_census_regions("CA16")

#list_census_regions("CA16") %>% 
#  filter(name=='Vancouver') #isolate the census regions with the name "Vancouver"

vancouver_da_2016 <- get_census(dataset='CA16', #specify the census of interest
                                regions=list(CSD="5915022"),#specify the region of interest
                                vectors=c("v_CA16_5792","v_CA16_5795","v_CA16_5798","
                                    v_CA16_5801","v_CA16_5804","v_CA16_5807",
                                          "v_CA16_5810"), # specify the attribute data of interest
                                labels = "short",# here we specify short variabhle names, otherwise they get very long
                                level='DA', #specify the census geography
                                geo_format = "sf",#specify the geographic data format as output
                                
                                )

vancouver_da_2016
plot(vancouver_da_2016["v_CA16_5807"])

#DYNAMIC MAPS
library(mapview)
#mapview(vancouver_da_2016,zcol = c("Population")) + 
#  mapview(osm_van_network)


###GIS ANALYSIS TOOL KIT
#ATTRIBUTE QUERIES

find_census_vectors("Median total income", tyoe = "total", dataset ="CA16", query_type = "semantic")
mv_census <- get_census(dataset='CA16', #specify the census of interest
                        regions=list(CSD="5915"),#specify the region of interest
                        vectors=c("v_CA16_2397"), # specify the attribute data of interest
                        labels = "short",# here we specify short variabhle names, otherwise they get very long
                        level='DA', #specify the census geography
                        geo_format = "sf", #specify the geographic data format as output
                        quiet = TRUE #turn off download progress messages (set to FALSE to see messages)
                        
)


cov_high_income <- mv_census %>% 
  filter(`Region Name`=="Vancouver" & v_CA16_2397 >175000)

mapview(mv_census) + mapview(cov_high_income,col.regions = "red")


streets <- van_bb %>%
  opq()%>%
  add_osm_feature(key = "highway") %>%
  osmdata_sf()

streets_sf <- streets$osm_lines

bicycle <- streets_sf %>% 
  filter(highway=="cycleway" | (highway=="path"&bicycle=="yes") | (highway=="path"&bicycle=="designated") | (highway=="footway"&bicycle=="yes") | (highway=="footway"&bicycle=="designated"))


mapview(streets_sf) + mapview(bicycle,color = "red") #marche pas meme erreur que d'habitude comprends pas

#SPATIAL QUERIES
van_boundary <- read_sf("/home/pndb/Scriptr_MJ/GIS_with_R/data/CoV_boundary.shp") # use the file location of where you have stored your data

mapview(mv_census) + mapview(van_boundary,col.regions = "red")

#st_intersects(mv_census,van_boundary,sparse = FALSE)

crs_to_use <- st_crs(van_boundary) # store crs information from the van_boundary dataset
crs_epsg <- crs_to_use$epsg # store the epsg code from the van_boundary dataset

mv_census_proj <- st_transform(mv_census,crs=crs_epsg)

st_crs(van_boundary)==st_crs(mv_census_proj)

st_intersects_output <- st_intersects(mv_census_proj,van_boundary,sparse = FALSE)

head(st_intersects_output)

census_intersects <- filter(mv_census_proj,st_intersects(x = mv_census_proj, y = van_boundary, sparse = FALSE))

mapview(van_boundary,col.region = "red") + mapview(census_intersects) 

mv_centroids <- st_centroid(mv_census_proj)

mapview(van_boundary,col.regions="red") + mapview(mv_census_proj) + mapview(mv_centroids,col.regions="black") 

census_intersects_cntrd <- filter(mv_census_proj,st_intersects(x = st_centroid(mv_census_proj), y = van_boundary, sparse = FALSE))

mapview(van_boundary,col.region = "red") + mapview(census_intersects_cntrd) 

#ATTRIBUTE JOINS
library(readxl)
deprivation_da <- read_xlsx("/home/pndb/Scriptr_MJ/GIS_with_R/data/bc_scores_quintiles_EN.xlsx",sheet = 1)

vancouver_da_dep <- left_join(census_intersects_cntrd,deprivation_da, by = c("GeoUID"="PRCDDA"))

mapview(vancouver_da_dep, zcol = "Residential instability Scores")

#BUFFERs
icbc_crash <- read_csv("/home/pndb/Scriptr_MJ/GIS_with_R/data/ICBC_cyclist_crashes.csv") %>%# use the file location of where you have stored your data
  filter(!is.na(Latitude)) %>% # remove the recrods which have missing values
  st_as_sf(.,coords = c("Longitude","Latitude"),crs = 4326) %>% 
  mutate(uniqid = row_number()) #assign unique identifier to each crash record

str(icbc_crash)#look at data frame

icbc_crash_sf <- st_as_sf(icbc_crash, coords = c("Longitude","Latitude"), crs = 4326)

str(icbc_crash_sf)
plot(icbc_crash_sf["Cyclist Flag"])
#st_crs(icbc_crash_sf) <- 4326

#possible_crs <- suggest_crs(icbc_crash_sf,units = "m",limit = 20)

#possible_crs <- tibble( y = icbc_crash_sf)#View the the top 10 suggestions
#possible_crs

icbc_crash_sf_proj <- st_transform(icbc_crash_sf,crs=3005)

st_crs(icbc_crash_sf_proj)

icbc_crash_van <- icbc_crash_sf %>% 
  filter(`Municipality Name`=="VANCOUVER")

#possible_crs_van <- suggest_crs(icbc_crash_van)

icbc_crash_van_proj <- st_transform(icbc_crash_van,crs=26910)

mapview(icbc_crash_van_proj)

st_crs(icbc_crash_van_proj)$units # my map projection uses metres as a unit - 

crash_buffer_icbc <- st_buffer(icbc_crash_van_proj,dist = 150)

mapview(crash_buffer_icbc)
#bicycle_buffer <- st_buffer(bicycle,dist=20)

st_crs(bicycle) <- st_crs(icbc_crash_van)$epsg

bicycle_buffer <- st_buffer(bicycle,dist=20)

mapview(bicycle_buffer)

#UNION
crash_buffer_union <- st_union(crash_buffer_icbc)

mapview(crash_buffer_union)

#INTERSECTION
intersect_van_census_crash_buffer <- st_intersection(census_intersects_cntrd,crash_buffer_icbc)

intersect_van_census_crash_buffer

plot(intersect_van_census_crash_buffer["v_CA16_2397"])

#average_income_by_crash_buffer <- intersect_van_census_crash_buffer %>% 
#  group_by(uniqid) %>% 
#  summarise(average_income = mean(v_CA16_2397,na.rm=TRUE))
#Uniqid pb pas compris
#plot(average_income_by_crash_buffer["average_income"])

#SPATIAL JOINS
icbc_crash_van_join_da <- st_join(icbc_crash_van_proj,vancouver_da_dep)
head(icbc_crash_van_join_da)

no_icbc_crashes_by_da <-   icbc_crash_van_join_da %>% 
  st_drop_geometry() %>%
  group_by(GeoUID) %>% 
  summarise(total_crashes = sum(`Crash Count`))

vancouver_da_dep_crashes <- left_join(vancouver_da_dep,no_icbc_crashes_by_da,by="GeoUID")

head(vancouver_da_dep_crashes)

library(tidyr)
vancouver_da_dep_crashes <- vancouver_da_dep_crashes %>% 
  mutate(total_crashes = replace_na(total_crashes,0))
mapview(vancouver_da_dep_crashes,zcol = "total_crashes") + mapview(icbc_crash_van)

st_snap_points <- function(x, y, namevar, max_dist = 1000) {
  
  # this evaluates the length of the data
  if (inherits(x, "sf")) n = nrow(x)
  if (inherits(x, "sfc")) n = length(x)
  
  # this part: 
  # 1. loops through every piece of data (every point)
  # 2. snaps a point to the nearest line geometries
  # 3. calculates the distance from point to line geometries
  # 4. retains only the shortest distances and generates a point at that intersection
  out = do.call("c",
                lapply(seq(n), function(i) {
                  nrst = st_nearest_points(st_geometry(x)[i], y)
                  nrst_len = st_length(nrst)
                  nrst_mn = which.min(nrst_len)
                  if (as.vector(nrst_len[nrst_mn]) > max_dist) return(st_geometry(x)[i])
                  return(st_cast(nrst[nrst_mn], "POINT")[2])
                })
  )
  # this part converts the data to a dataframe and adds a named column of your choice
  out_xy <- st_coordinates(out) %>% as.data.frame()
  out_xy <- out_xy %>% 
    mutate({{namevar}} := x[[namevar]]) %>% 
    st_as_sf(coords=c("X","Y"), crs=st_crs(x), remove=FALSE)
  
  return(out_xy)
}

icbc_snp <- st_snap_points(x = icbc_crash_van,#crash van 
                          y = bicycle, #network infrasructure data
                          max_dist = 20,
                          namevar = "uniqid") #this only returns the geometry - doesn't preserve attributes

icbc_snp <-  left_join(icbc_snp,
                       icbc_crash_van %>% 
                         st_drop_geometry(),
                       by = "uniqid")#return the attributes to the snapped locations


mapview(icbc_snp) + mapview(bicycle)

icbc_snp_bicycle_inf <- st_join(icbc_snp,
                                bicycle["osm_id"], #select only the column "osm_id"
                                join = st_is_within_distance, dist = 0.001) #for each point assign the LIXID that it falls on by indicating a very small distance tolerance (here it is 1mm)


icbc_snp_bicycle_inf %>% 
  select(uniqid,`Crash Count`,osm_id) # look at the results


no_crashes_by_osm_id <- icbc_snp_bicycle_inf %>% 
  group_by(osm_id) %>% 
  summarise(n_crashes = sum(`Crash Count`)) %>%
  st_drop_geometry() # drop the geometry from this - we only need the uniqid for the bicycle infrastrucutre data and the counts



bicycle_infra_crash_count <- left_join(bicycle,no_crashes_by_osm_id , by = "osm_id") %>%
  mutate(n_crashes = replace_na(n_crashes,0))

mapview(bicycle_infra_crash_count,zcol="n_crashes") + mapview(icbc_snp)


bicycle_cut <- bicycle %>% 
  mutate(length_m = st_length(.))

bicycle_cut


mapview(bicycle_cut) + mapview(vancouver_da_dep)

cut_join_da <- st_join(bicycle_cut, vancouver_da_dep,join = st_is_within_distance,dist=5)

length_by_da <- cut_join_da %>% 
  st_drop_geometry() %>%
  group_by(GeoUID) %>%
  summarise(total_length_km = sum(as.numeric(length_m))/1000)

head(length_by_da)

vancouver_da_sub <- left_join(vancouver_da_dep,length_by_da) %>%
  mutate(total_length_km = replace_na(total_length_km,0))

vancouver_da_sub

mapview(vancouver_da_sub,zcol = "total_length_km") + mapview(bicycle)

vancouver_da_sub <- vancouver_da_sub %>% 
  mutate(total_area_km2 = as.numeric(st_area(.)/1000/1000),
         total_bike_dens = total_length_km/ total_area_km2                     
  )

mapview(vancouver_da_sub,zcol = "total_bike_dens")  + mapview(bicycle)
  