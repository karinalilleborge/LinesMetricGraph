library(MetricGraph)
library(INLA)
library(inlabru)
library(rSPDE)
library(sf)
#source("line_support_MG.R")
devtools::load_all()
# Edges of a simple graph
edge1 <- rbind(c(0,0),c(1,0))
edge2 <- rbind(c(0,0),c(0,1))
edge3 <- rbind(c(0,1),c(-1,1))
theta <- seq(from=pi,to=3*pi/2,length.out = 20)
edge4 <- cbind(sin(theta),1+ cos(theta))
edge5 <- rbind(c(0,0), c(1,0))
edge6 <- rbind(c(1,0), c(1,1))
edge7 <- rbind(c(1,1), c(2,1))
edge8 <- rbind(c(0,1), c(1,1))
edges = list(edge1, edge2, edge3, edge4, edge5, edge6, edge7, edge8)
# Graph construction
graph <- MetricGraph::metric_graph$new(edges = edges)
# Build mesh
graph$build_mesh(h = 0.01)
#graph <- fm_as_MG(graph, MGG=FALSE)
# Make model construction, we fix smoothness here
rspde_model <- rspde.metric_graph(graph, nu=0.5)
# Define a valid path object from line segment to graph
line1 <- st_sfc(st_linestring(matrix(c(0.5,1.1,1,1),ncol=2)))
line2 <- st_sfc(st_linestring(matrix(c(1.3,2,1,1),ncol=2)))
#lines <- st_sfc(list(line1,line2))
path1 <- geom_path_to_path_MGG(line1, graph)
path2 <- geom_path_to_path_MGG(line2, graph)

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
for (id in unique(path1$ID)) {
  i <- i + 1
  # check the length of computed path here
  paths[[i]] <- path1[path1$ID == id, c("paths")]
}
for (id in unique(path2$ID)) {
  i <- i + 1
  paths[[i]] <- path2[path2$ID == id, c("paths")]
}

# Construct sampler using path and weight
sampler <- tibble::tibble(x = paths,
                          weight = rep(1,length(paths)))
# Mapper (here we choose aggregate and rescale for integral observations)
agg <- bru_mapper_aggregate(rescale =  TRUE, n_block = nrow(sampler))
# Find integration points in the mesh
ips <- fm_int(graph, sampler)
#ips2 <- ips
#ips2$.block <- ips2$.block + 2
#ips <- rbind(ips, ips2)


#simulate a true field u and its observations
rspde.order <- 2
nu <- 0.5
op <- matern.operators(
  nu = nu, range = 1.5, sigma = 3,
  parameterization = "matern",
  m = rspde.order, graph = graph
)
u <- simulate(op, nsim = 1)
y <- rnorm(nrow(sampler), sd = 0.5) + 5 +
  with(ips, ibm_eval(agg,
                     input = list(block = .block, weights = weight),
                     state = fm_evaluate(graph, loc = x, field = as.vector(u))
  ))


# Formula
formula <- y ~ ibm_eval(agg,
                        input = list(block = .block, weights = weight),
                        state = Intercept + spde)
# Observation data for inlabru
obs <- bru_obs(formula = formula,
               response_data = data.frame(y=y),
               data = ips,
               allow_combine = TRUE)

# inlabru call
bru_res <- bru(components = y ~ Intercept(1) +
                 spde(x, model = rspde_model, mapper = bru_mapper(graph)),
               obs)
# inlabru summary
summary(bru_res)
# rSPDE parameter summary
summary(rspde.result(bru_res, "spde", rspde_model))
