#' Generate fake data
#' 
#' Generate a dataframe of locations and labels for use with testing the 
#' error identification and correction algorithm
#' @param iso The all-caps, 3-letter ISO code of a country
#' @param n_communities The number of communities
#' @param n_locations The number of individual locations
#' @param max_distance The maximum distance in meters a location can be from
#' its community's centroid
#' @param sd_size The standard deviation in the number of locations per community
#' @return A dataframe with one row for each location (appoximately, but 
#' most often not exactly \code{n_locations} rows) and columns indicating the
#' longitude, latitude and community of each point

generate_fake_data <- function(iso = 'KEN',
                               n_communities = 15,
                               n_locations = 1000,
                               max_distance = 10000,
                               sd_size = 15){
  
  require(sp)
  require(dplyr)
  require(raster)
  require(geosphere)
  
  # Get a shapefile
  map <- getData(name = 'GADM', country = iso, level = 3)
  # From the shapefile, get vectors of possible x and y coordinates
  xs <- seq(bbox(map)[1,1],
            bbox(map)[1,2],
            length = n_locations * 100)
  ys <- seq(bbox(map)[2,1],
            bbox(map)[2,2],
            length = n_locations * 100)
  # Create community centroids
  centroids <- data_frame(id = 1:n_communities,
                          x = rep(NA, n_communities),
                          y = rep(NA, n_communities))
  for (i in 1:nrow(centroids)){
    message(paste0('Creating a centroid for community number ', i))
    ok <- FALSE
    # Make sure we get a point in the map
    while(!ok){
      # Create a random point
      random_point <- data_frame(x = sample(xs, 1),
                                 y = sample(ys, 1))
      # Make spatial
      random_point_spatial <- random_point
      coordinates(random_point_spatial) <- ~x+y  
      proj4string(random_point_spatial) <- proj4string(map)
      # Determine if it's in the map
      ok <- as.logical(!is.na(over(random_point_spatial, polygons(map))))
    }
    # Now that it's ok, add that point to the centroids dataframe
    centroids$x[i] <- random_point$x
    centroids$y[i] <- random_point$y
  }
  # Make centroids spatial
  centroids$lng <- centroids$x
  centroids$lat <- centroids$y
  centroids <- data.frame(centroids)
  coordinates(centroids) <- ~x+y
  proj4string(centroids) <- proj4string(map)
  
  # Create actual points
  average_community_size <- round(n_locations / n_communities)
  locations_per_community <- round(rnorm(n = n_communities, 
                                         mean = average_community_size, 
                                         sd = sd_size))
  which_negatives <- which(locations_per_community < 1)
  negatives <- length(which(locations_per_community < 1))
  remainder <- n_locations - 
    sum(locations_per_community[locations_per_community > 0])
  n_negatives <- length(which_negatives)
  if(n_negatives > 0){
    locations_per_community[which_negatives] <- 
      round(remainder / n_negatives)
  }
  results_list <- list()
  counter <- 0
  for (i in 1:nrow(centroids)){
    message('Generating random points for centroid ', i)
    # Extract this centroid
    this_centroid <- centroids[i,]
    this_n <- locations_per_community[i]
    enough <- FALSE
    master_points <- this_centroid[0,]
    while(!enough){
      # Create a random set of points
      these_points <- data.frame(id = centroids$id[i],
                                 x = sample(xs, n_locations * 100),
                                 y = sample(ys, n_locations * 100))
      coordinates(these_points) <- ~x+y
      proj4string(these_points) <- proj4string(map)
      # Get distance from centroid to each location
      distances <- as.numeric(distm(x = this_centroid,
                                    y = these_points,
                                    fun = distHaversine))
      # Keep only those below the max distance
      these_points <- these_points[which(distances <= max_distance),]
      master_points <- rbind(these_points, master_points)
      if(length(master_points) >= this_n){
        master_points <- master_points[1:this_n,]
        enough <- TRUE
      }
    }
    # Add points to the master list
    results_list[[i]] <- master_points
  }
  # Combine all the results and return
  results <- do.call('rbind', results_list)
  results@data$longitude <- coordinates(results)[,1]
  results@data$latitude <- coordinates(results)[,2]
  return(results)
}
