if(require("MetricGraph") == FALSE) {install.packages("MetricGraph")}
library(MetricGraph)
if(require("rSPDE") == FALSE) {install.packages("rSPDE")}
library(rSPDE)
if(require("INLA") == FALSE) {install.packages("INLA")}
library(INLA)
if(require("inlabru") == FALSE) {install.packages("inlabru")}
library(inlabru)
if(require("sf") == FALSE) {install.packages("sf")}
library(sf)
# if(require("ggplot2") == FALSE) {install.packages("ggplot2")}
# library(ggplot2) # package for plotting
if(require("tidyverse") == FALSE) {install.packages("tidyverse")}
library(tidyverse)

## Model parameters----
sigma <- 1
range <- 1000
nu <- 0.5
## Likelihood parameters ----
noise_lines <- 0.25
noise_points <- 0.01
#Setting the number of replicates----
nweeks <- 5 #1,  5, 25
nsim <- 1
#Import the network that your graph will be based off of
#See vignettes of MetricGraph for other graph creation methods
# https://davidbolin.github.io/MetricGraph/articles/metric_graph.html
path <-  ""
folder <- "geometries"
file_road <- "road.gpkg"

# Must be downloaded/stored in your working directory (following file can be found at
# https://github.com/karinalilleborge/LinesMetricGraph)
source("metric_graph.R")
road_links <- sf::st_read(paste0(path, folder,file_road))

#Construct graph
graphs <- graph_components$new(
  edges = road_links$geom,
  longlat = T,
  perform_merges = T,
  vertex_unit = "m",
  length_unit = "m"
)
graph <- graphs$get_largest()
#Free up some memory
rm(road_links)
rm(graphs)
# remove vertices of order 2
graph$prune_vertices()
# build the mesh
graph$build_mesh(h = 70) # meter


## Operator for simulating the truth ----
rspde.order <- 2
op <- matern.operators(
  nu = nu, range = range, sigma = sigma,
  parameterization = "matern",
  m = rspde.order, graph = graph
)

## Latent field parameters and covariate ----
beta0 <- 1
beta <- 1
cov_op <- matern.operators(
  nu = 1.5, range = 6000, sigma = 3,
  parameterization = "matern",
  m = rspde.order, graph = graph
)
#Optional seed
#set.seed(310125)
cov <- simulate(cov_op)
#Standardizing the covariate to ensure zero-mean
cov <- (cov-mean(cov))/sd(as.vector(cov))

#Free up memory again
rm(cov_op)

## Model ----
rspde_model <- rspde.metric_graph(graph, nu = nu,
                                  parameterization = "matern",
                                  prior.range.nominal = 700,
                                  prior.std.dev.nominal = 1)

# Make four paths through the graph ----
#Can be replaced with other st_LINESTRING objects
line1_sf <- sf::st_read(paste0(path, folder, "lines1.gpkg"))
line2_sf <- sf::st_read(paste0(path, folder, "lines2.gpkg"))
line3_sf <- sf::st_read(paste0(path, folder, "lines3.gpkg"))
line4_sf <- sf::st_read(paste0(path, folder, "lines4.gpkg"))

#Get the lengths of the lines, for scaling of observation noise
scale1 <- sf::st_length(line1_sf)
scale2 <- sf::st_length(line2_sf)
scale3 <- sf::st_length(line3_sf)
scale4 <- sf::st_length(line4_sf)
scales <- c(scale1, scale2, scale3, scale4)

#Get the geometry
line1_sf <- sf::st_geometry(line1_sf)
line2_sf <- sf::st_geometry(line2_sf)
line3_sf <- sf::st_geometry(line3_sf)
line4_sf <- sf::st_geometry(line4_sf)

#Converting to valid paths
path1 <- geom_path_to_path_MGG(line1_sf, graph)
path2 <- geom_path_to_path_MGG(line2_sf, graph)
path3 <- geom_path_to_path_MGG(line3_sf, graph)
path4 <- geom_path_to_path_MGG(line4_sf, graph)

#Construct the full sampler
paths <- list()
i <- 0
for (id in unique(path1$ID)) {
  i <- i + 1
  # check the length of computed path here
  paths[[i]] <- path1[path1$ID == id, c("paths")]
}
for (id in unique(path2$ID)) {
  i <- i + 1
  paths[[i]] <- path2[path2$ID == id, c("paths")]
}
for(id in unique(path3$ID)){
  i = i +1
  paths[[i]] <- path3[path3$ID==id, c("paths")]
}
for (id in unique(path4$ID)) {
  i <- i + 1
  paths[[i]] <- path4[path4$ID == id, c("paths")]
}

#Free up memory
rm(path1)
rm(path2)
rm(path3)
rm(path4)

#Determine mid points for the wrong support model (WSM)
midpoints1_sf <- sf::st_centroid(line1_sf)
midpoints2_sf <- sf::st_centroid(line2_sf)
midpoints3_sf <- sf::st_centroid(line3_sf)
midpoints4_sf <- sf::st_centroid(line4_sf)


#Free up memory
rm(line1_sf)
rm(line2_sf)
rm(line3_sf)
rm(line4_sf)

# Sampler ----
test_sampler <- tibble::tibble(x = paths, weight = rep(1, length(paths)))
# Doing the integration ----
agg <- bru_mapper_aggregate(rescale = TRUE)
#Finding integration points
ips <- fm_CSM(graph, test_sampler)
#Repeat for the replicates
ips_full <- do.call("rbind", replicate(nweeks, ips, simplify = FALSE))
nblocks <- length(unique(ips$.block))
#Update the $.block for new replicates
if(nweeks>1){
  for(i in 2:nweeks){
    ips_full$.block[((i-1)*nrow(ips)+1):(i*nrow(ips))] = ips_full$.block[((i-1)*nrow(ips)+1):(i*nrow(ips))] + nblocks*(i-1)
  }
}
#Construct bru_mapper for aggregation
agg_full <- bru_mapper_aggregate(rescale = TRUE, n_block =length(unique(ips_full$.block)))


# Point stations -----
points_sf <- sf::st_read(paste0(path, "case study data/","point_obs_locations.gpkg"))
#remove two stations outside the graph
points_sf <- points_sf[-c(3:4), ]
stations <- fm_bary(graph, points_sf, MGG=TRUE)
# replicates
point_loc <- do.call("rbind", replicate(nweeks,
                                        stations,
                                        simplify = FALSE))
#free up memory
rm(points_sf)

#replicates lines
replicates_l <- data.frame(".block" = ips_full$.block,
                          "repl" = rep(seq_len(nweeks), each=nrow(ips)))
# add column specifying the replicate ID
ips_full$repl <- rep(seq_len(nweeks), each=nrow(ips))
ips2 <- ips_full
colnames(ips2) <- c("x", "weight", ".block", "repl")
# replicate points
replicates_p <- rep(seq_len(nweeks), each=nrow(stations))

# Determine the observation locations for each obs:
# Midpoints, graph coordinates:
obs_loc <- rbind(fm_bary(graph, midpoints1_sf),
                 fm_bary(graph, midpoints2_sf),
                 fm_bary(graph, midpoints3_sf),
                 fm_bary(graph, midpoints4_sf))
#WSM
replicates_l_WSM <- rep(1:nweeks, each = nrow(obs_loc))
obs_loc <- do.call("rbind", replicate(nweeks, obs_loc, simplify = FALSE))

#formulas
form_lines <- y ~ ibm_eval(agg,
                           input = list(block = .block, weights = weight),
                           state = beta0 + covariate + spde
)
form_points <- y ~  beta0 + covariate + spde
form_WSM_lines <- y ~ beta0 + ibm_eval(agg, input = list(block = .block, weights = weight), state = covariate) + spde
form_WSM_points <- y ~ beta0 + covariate + spde

#covariate (global environment)
cov1 <- as.vector(cov)

#priors (can be set here)
prec.prior.lines <- list(theta = list(prior = "loggamma", param = c(1, 5e-05)),
                         initial = 4, fixed = FALSE)
prec.prior.points <- list(theta = list(prior = "loggamma", param = c(1, 5e-05)),
                          initial = 4, fixed = FALSE)

#### CRPS: -----
mean_crps <- function(mean, sd, obs){
  z = as.vector((obs - mean) / sd)
  return( mean(sd * ( z*(2*pnorm(z)-1) + 2*dnorm(z) - 1/sqrt(pi) ) ))
}
#Optional seed
#set.seed(010125)
for(j in 1:nsim){
  u <- simulate(op, nsim = nweeks)
  ### Construct the data ----
  y_lines <- matrix(nrow = length(unique(ips$.block)), ncol = nweeks)
  y_points <- matrix(nrow = nrow(stations), ncol = nweeks)
  for(i in seq_len(nweeks)){
    y_lines[,i] <- rnorm(nrow(test_sampler), sd = noise_lines/scales) + beta0 +
      beta * with(ips, ibm_eval(agg,
                                input = list(block = .block, weights = weight),
                                state = fm_evaluate(graph, loc = x, field = cov)
      )) +
      with(ips, ibm_eval(agg,
                         input = list(block = .block, weights = weight),
                         state = fm_evaluate(graph, loc = x, field = as.vector(u[,i]))
      ))
    y_points[,i] <- rnorm(nrow(stations), sd = noise_points) + beta0 +
      beta * fm_evaluate(graph, cov, loc = stations) +
      fm_evaluate(graph, u[,i], loc = stations)
  }
  # Intergration model ----
  ### Likelihood----
  lik_lines <- bru_obs(
    formula = form_lines,
    #control.family = list(hyper=prec.prior.lines),
    response_data = data.frame(y = c(y_lines)),
    scale = rep(scales**2, nweeks),
    data = ips2,
    allow_combine = TRUE
  )

  lik_points <- bru_obs(
    formula = form_points,
    #control.family = list(hyper=prec.prior.points),
    response_data = data.frame(y = c(y_points)),
    data = tibble::tibble(x = point_loc,
                          repl = replicates_p),
    allow_combine = TRUE
  )
  ### Fitting the model ----
  res_CSM <- bru(
    components = y ~ beta0(1) +
      covariate(fm_evaluate(graph, cov1, loc = x), model = "linear") +
      spde(x, model = rspde_model, mapper = bru_mapper(graph, n_rep = 1), replicate=repl),
    lik_lines, lik_points, options = list(verbose = TRUE, bru_verbose = 4, safe = FALSE)
  ) #specify runs
  spde_res_CSM <- rspde.result(res_CSM, "spde", rspde_model)
  rm(lik_points)
  rm(lik_lines)
  # Simple model ----
  lik_WSM_lines <- bru_obs(form_WSM_lines,
                              response_data = data.frame(y = c(y_lines)),
                              scale = rep(scales**2, nweeks),
                              data = list(
                                loc = obs_loc,
                                loc_cov = ips2$x, weight = ips2$weight, .block = ips2$.block,
                                repl = replicates_l_WSM
                              ),
                              allow_combine = TRUE
  )
  lik_WSM_points <- bru_obs(form_WSM_points,
                               response_data = tibble::tibble(y = c(y_points)),
                               data = list(
                                 loc = point_loc,
                                 loc_cov = point_loc,
                                 repl = replicates_p
                               ),
                               allow_combine = TRUE
  )
  ### Fitting the model ----
  res_WSM <- bru(
    components = y ~ beta0(1) +
      covariate(fm_evaluate(graph, cov1, loc = loc_cov), model = "linear") +
      spde(loc, model = rspde_model, mapper = bru_mapper(graph, n_rep = 1), replicate = repl),
    lik_WSM_lines, lik_WSM_points
  )
  spde_res_WSM <- rspde.result(res_WSM, "spde", rspde_model)
  rm(lik_WSM_lines)
  rm(lik_WSM_points)
  #*Store the results.*
  #Predict in all mesh locations
  pred_mesh <- as_MGG(graph$mesh$VtE)
  for(i in 1:nweeks){
    #Compute the truth
    truth <- beta0 + beta*cov + u[,i]
    pred_CSM <- predict(res_CSM,
                        newdata = list(x = pred_mesh,
                                       cov1 = as.vector(cov),
                                       repl=i),
                        formula = ~ beta0 + covariate + spde) #
    pred_WSM <- predict(res_WSM,
                        newdata = list(loc = pred_mesh,
                                       loc_cov = pred_mesh,
                                       cov1 = as.vector(cov),
                                       repl=i),
                        formula = ~ beta0 + covariate + spde)
    ## MSE (RMSE)
    mse_WSM <- mean((pred_WSM$mean-(truth))**2)
    mse_CSM <- mean((pred_CSM$mean-(truth))**2)
    ## CRPS (approximate)
    crps_WSM <- mean_crps(pred_WSM$mean, pred_WSM$sd, truth)
    crps_CSM <- mean_crps(pred_CSM$mean, pred_CSM$sd, truth)
    # coverage
    coverage_WSM <- mean( 1* (truth >= pred_WSM$mean - 1.96 * pred_WSM$sd &
                                truth <= pred_WSM$mean + 1.96 * pred_WSM$sd))
    coverage_CSM <- mean(1* (truth >= pred_WSM$mean - 1.96 * pred_WSM$sd &
                               truth <= pred_WSM$mean + 1.96 * pred_WSM$sd))
    scoring_rules <- data.frame(MSE = c(mse_WSM,mse_CSM),
                                CRPS=c(crps_WSM,crps_CSM),
                                coverage = c(coverage_WSM, coverage_CSM),
                                row.names = c("WSM", "CSM"))
    #*Store the results for this week*
  }
}

