

setwd("F:\\Black_mypassport_1.8TB\\Smith_vanLeeuwen_work\\Babst_T-lidar\\tlidar_workspace")



##This is a 'meta-package' that lists all of the R packages for lidar and forest analysis. 
#devtools::install_github("atkinsjeff/ForestAnalysisInR") 
#library(ForestAnalysisInR)
#library(shiny)
#launchRFA()


##install TreeLS package that is not on CRAN
#remotes::install_github('tiagodc/TreeLS')

##install latest version of lidr tools (bug fixes) from github. Latest version is not on CRAN yet. 
#remotes::install_github('Jean-Romain/lidr')


library(raster)
library(lidR)
library(TreeLS)
library(rgdal) 
library(sf) 
library(rgeos) 
library(sp)
#lidR_parallelism

set_lidr_threads(30)
get_lidr_threads()

##Bring las point cloud into Rstudio
#points = readLAS("TLS_SantaRita_UBplot_pointcloud_041922_aligned_CC.las")

##Bring FARO Scene georeferenced point cloud (ASCII) into R. 
points = readTLS("TLS_SantaRita_UBplot_pointcloud_041922_utm.xyz")

points = readTLS("UB_tree1.xyz")


##Voxelize to reduce point cloud density to 1 pnt/3 cubic cm
points_voxel = voxelize_points(points, 0.03)
plot(points_voxel, color ='Z')


##Delete the points on periphery of the cloud that are outside of the 'plot'
##This clipping procedure will be based on user defined plot center coordinate (GPS) and radius of plot
#The plot center coordinates and radii are currently found in the file 'Site_overview_SantaRitaMts_October2021.xls

center_coord_easting_x = 513637.25

center_coord_northin_y = 3507116.26
  
points_voxel_clip = clip_circle(points_voxel, center_coord_easting_x, center_coord_northin_y, 30)

plot(points_voxel_clip, color = 'Z')

##Separate trees from ground using lidr cloth simulation filter
points_classify = classify_ground(points_voxel_clip, algorithm = csf(sloop_smooth = TRUE, class_threshold = 0.1, cloth_resolution =  0.5, rigidness = 2))

#points_classify = classify_ground(points_voxel, algorithm = csf(sloop_smooth = TRUE, class_threshold = 0.1, cloth_resolution =  0.5, rigidness = 2))

plot(points_classify, color = 'Classification')

##Make a point cloud file for just ground points
ground_points = filter_poi(points_classify, Classification == 2)
plot(ground_points)

#Make a point cloud file for just non-ground (trees) points
tree_points = filter_poi(points_classify, Classification != 2)
plot(tree_points)

##Create a grided digital terrain model from ground points
DTM = grid_terrain(ground_points, res = 0.5, algorithm = knnidw(k = 10, p = 2))

##Calculate vertical distance of each tree point above the DTM. This creates a vegetation height point cloud.
normalized_height = normalize_height(las = tree_points, algorithm = DTM)

##Remove points that are below 0
AGL_clean = filter_poi(normalized_height, Z > 0.5)
plot(AGL_clean, color = 'Z')
##Export normalized point cloud to .las
#writeLAS(AGL_clean, "F:\\Black_mypassport_1.8TB\\Smith_vanLeeuwen_work\\Babst_T-lidar\\tlidar_workspace\\UB_AGL.las")


##Create canopy height model raster. I don't think this is totally necessary. 
#CHM = rasterize_canopy(AGL_clean, res = 0.05)




##Identify stems with 'TreeLS'

##Identify tree occurrences for a normalized point cloud. It uses a Hough Transform which is a circle search
##Merge is a distance between stems. If two occurrences are less than the specified distance, they will be merged. 
##I don't completely understand how this algorithm works and what the parameters mean yet. It does seem to work, though. 

#tree_map2 = treeMap(AGL_clean, method = map.hough(min_h = 1, max_h = 3, h_step = 0.5, pixel_size = 0.025, 
                                                # max_d = 0.85, min_density = 0.1, min_votes = 3), merge = 0.2, positions_only = FALSE)

#The eigen method seems to find more tree occurrences than the hough method
tree_map_eigen = treeMap(AGL_clean, method = map.eigen.knn(max_curvature = 0.1,
                                                           max_verticality = 10,
                                                           max_mean_dist = 0.1,
                                                           max_d = 0.5,
                                                           min_h = 1.5,
                                                           max_h = 3), merge = 0.2, positions_only = FALSE)


##Plot the voxelized point cloud with the identified tree stems
#AGL_clean_topped = filter_poi(AGL_clean, Z < 5)
x = plot(AGL_clean, color = "Z")
add_treeMap(x, tree_map_eigen, color = 'yellow', size = 3)



##Count the number of unique trees
unique_trees = unique(tree_map_eigen$TreeID)
length(unique_trees)




##Create 2D point layer showing the locations of the tree stems
##These point locations will be compared with tree coordinates collected in the field
xymap = treeMap.positions(tree_map_eigen)

## Export 2D tree stem location points in a CSV
#write.csv(xymap,"F:\\Black_mypassport_1.8TB\\Smith_vanLeeuwen_work\\Babst_T-lidar\\tlidar_workspace\\UB_stems_eigen_29april2022.csv", row.names = TRUE)

 

##Identify all of the points that belong to each tree occurrence found in 'treeMap'
tree_points = treePoints(AGL_clean, tree_map_eigen, trp.voronoi())

##Plot AGL_clean with tree points and tree IDs
x = plot(AGL_clean, color = "Z")
add_treePoints(x, tree_points, size=4)
add_treeIDs(x, tree_points, color='yellow', cex=2)

##create single tree point cloud .las
single_tree_1 = filter_poi(tree_points, TreeID == '58'|TreeID == '32')
plot(single_tree_1, color = 'Z')


##Identify stem points in each segmented tree. The 'h_step' parameter sets the height of each cylinder
stem_points_hough1 = stemPoints(tree_points, method = stm.hough(h_step = 0.15,
                                                               max_d = 1,
                                                               h_base = c(1, 2.5),
                                                               pixel_size = 0.05,
                                                               min_density = 0.1,
                                                               min_votes = 3))

##plot AGL_clean point cloud with stem points
x = plot(AGL_clean, color = "Z")
add_stemPoints(x, stem_points_hough1, color='red', size=8)

#plot only stem points
stem_points_hough_isolated1 = filter_poi(stem_points_hough1, Stem == TRUE)
plot(stem_points_hough_isolated1 , color = "Segment")

##Identify stem points for a single tree
stem_points_single = stemPoints(single_tree_1, method = stm.hough(h_step = 0.15,
                                                                max_d = 1,
                                                                h_base = c(1, 2.5),
                                                                pixel_size = 0.05,
                                                                min_density = 0.1,
                                                                min_votes = 3))


#Plot single tree stem points with segmented tree points
x = plot(single_tree_1, color = "Z")
add_stemPoints(x, stem_points_single, color='red', size=8)



#plot only single tree stem points
stem_points_single_tree = filter_poi(stem_points_single, Stem == TRUE)
plot(stem_points_single_tree , color = "Segment")


##different method to identify stemPoints using eigen.knn method. 
##stem_points_eigen = stemPoints(tree_points, method = stm.eigen.knn(h_step = 0.5,
                                                                #   max_curvature = 0.1,
                                                                  # max_verticality = 10,
                                                                 #  voxel_spacing = 0.025,
                                                                  # max_d = 0.5,
                                                                  # votes_weight = 0.2,
                                                                  # v3d = FALSE))



writeLAS(stem_points_hough_isolated1, "F:\\Black_mypassport_1.8TB\\Smith_vanLeeuwen_work\\Babst_T-lidar\\tlidar_workspace\\stem_points_hough_isolated1.las")



##fit many cylinders to stem points
stem_seg = stemSegmentation(stem_points_single_tree, sgt.ransac.circle(tol = 0.1, n = 10, 
                                                                      conf = 0.95, inliers = 0.9))

##plot stem segments, WATCH the Cylinders stacked into tree stems
x = plot(stem_points_single_tree)
add_stemSegments(x, stem_seg, color='blue')
#add_treeIDs(x, tree_points, color='yellow', cex=1)

##Alternative way to watch the cylinders stack into a tree stem
tlsPlot(stem_seg)


##Calculate the cylinder volume for each stem cylinder, then add a column in the 'stem_seg' table showing the calculation
stem_seg_2 = cbind(stem_seg, stem_seg$Radius^2*3.1415*0.25)
colnames(stem_seg_2) = c('TreeID', 'Segment', 'X', 'Y', 'Radius', 'Error', 'AvgHeight', 'N', 'cylinder_volume')



##### Sum the cylinder volumes for each tree
tree_volume = aggregate(x= stem_seg_2$cylinder_volume,
          by= list(stem_seg$TreeID),
          FUN=sum, na.rm=TRUE)
colnames(tree_volume) = c('TreeID', 'Cubic Meters')


##This will create a data.frame table that lists every tree and the radius at DBH (1.3 m AGL)
inv = tlsInventory(stem_points_hough_isolated1, dh = 1.3, dw = 0.5, hp = 1, 
                                    d_method = shapeFit(shape = "circle", 
                                    algorithm = "ransac", n = 15, n_best = 20))
inv_2 = cbind(inv, inv$Radius*2)
colnames(inv_2) = c('TreeID', 'X', 'Y', 'Radius', 'Error', 'H', 'h_radius', 'DBH')

##plot the tree inventory. Graphic shows height of each tree and DBH
tlsPlot(inv)




writeLAS(stem___points, "F:\\Black_mypassport_1.8TB\\Smith_vanLeeuwen_work\\Babst_T-lidar\\tlidar_workspace\\stem___points.las")




#########The following code will probably not be included. Most are inferior methods to segment trees#################

##Individual tree detection 
#Almost all of the algorithms I have found in R use a local maximum filter to identify tree tops

#Identify individual trees. Tree tops are identified by a user defined moving window looking for local maximums across CHM
tree_tops = locate_trees(AGL_clean, lmf(ws = 4, shape = "circular",hmin=4))

#Identified tree tops are the starting points for a region grow routine that adds new pixels to the tree based on some user defined thresholds
#th_seed = pixel is added to tree if its height is greater than user defined proportion multiplied by local max height
#For example, if a tree top is 10 m high, and the parameter is set to 0.25, then the pixel needs to be at least 2.5 m high to be added.
#th_cr = similar to th_seed, except instead of using local max height, it uses mean height of existing region. 
#Tip: to grow the regions bigger, make th_seed and th_cr have small values (e.g. < 0.25)
treeseg_dalponte= lastrees(AGL_clean, dalponte2016(chm = CHM, treetops = tree_tops, th_tree = 1.5, th_seed = 0.20,
                                                      th_cr = 0.20, max_cr = 35, ID = "treeID"))

#Calculate a convex hull and draw a polygon for each tree
metric = tree_metrics(treeseg_dalponte, .stdtreemetrics)
hulls4  = tree_hulls(treeseg_dalponte)
hulls4@data = dplyr::left_join(hulls4@data, metric@data)

#Export the convex hull polygons into a shapefile 
writeOGR(hulls4, dsn = "E:\\Jemez_drone_Oct2020\\Casa_angelica_mid_canyon", layer = "mid_canyon_treeseg", driver = "ESRI Shapefile", overwrite_layer = TRUE)

#Export the treetop points into a shapefile. This can be used to understand how trees were identified. Not totally necessary in the workflow. 
writeOGR(treetops, dsn = "D:\\Black_mypassport_1.8TB\\Jemez_TNC\\GIS_data", layer = "treeseg_dalponte3_tops", driver = "ESRI Shapefile", overwrite_layer = TRUE)


##Identify individual trees using the algorithm by Li et al. 2012
#This algorithm segments one tree at a time. It starts by finding the highest point in the normalized point cloud
#and assumes it is the top of a tree. Using a top-down classification approach, it moves down to points below the highest
#point and classifies them as part of the same tree or not. 
#dt1 and dt2 are Euclidian distance (2D) thresholds from the point to be classified and the points that are already classified as 
#the tree. If the distance is greater than the threshold, then the point is classified as not part of the tree. Conversely,
#if the point is within the threshold distance, then it is classified as part of the tree. There are two thresholds (dt1 & dt2)to lend
#flexibility to adjust the threshold distance depending on the height of the point to be classified. This is recognition that the tops of trees
#may be further spaced apart than the middle of the trees (works well for pine tree species). 

#Zu is the point height parameter you set which triggers the use of dt1 or dt2. A point > Zu will use dt2. A point <= to Zu will use dt1.

treeseg_Li = segment_trees(AGL_clean, li2012(dt1 = 3, dt2 = 2, R = 0, 
                                         Zu = 10, hmin = 4, speed_up = 5), attribute = "treeID")


treeseg_silva = segment_trees(AGL_clean, silva2016(chm = CHM, treetops = locate_trees2, max_cr_factor = 0.6, exclusion = 0.3, ID = "treeID"), attribute = "treeID")

treeseg_watershed = segment_trees(AGL_clean, watershed(chm = CHM, th_tree = 2, tolerance = 1), attribute = "treeID")

#Count the nmber of unique trees
unique_trees = unique(tree_map2$TreeID)
length(unique_trees)



plot(locate_trees)


writeLAS(treeseg, "F:\\Black_mypassport_1.8TB\\Smith_vanLeeuwen_work\\Babst_T-lidar\\tlidar_workspace\\treeseg1.las")

writeLAS(AGL_clean, "F:\\Black_mypassport_1.8TB\\Smith_vanLeeuwen_work\\Babst_T-lidar\\tlidar_workspace\\normalized_heights.las")



