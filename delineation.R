#!/usr/bin/env Rscript
# Description: Computes urban delineations
# Author: Clement Gorin, gorin@gate.cnrs.fr
# Date: January 2021

# Packages
suppressMessages(if(!require('pacman')) install.packages('pacman', repos = 'https://cloud.r-project.org/'))
pacman::p_load(compiler, fst, future.apply, imager, memuse, optparse, raster, rgdal, tictoc)
options(warn = -1)

# Parameters
arguments <- list(
  make_option('--density',   type = 'character', default = '',  help = 'Path to input raster.'),
  make_option('--unlivable', type = 'character', default = '',  help = 'Path to unlivable raster.'),
  make_option('--outdir',    type = 'character', default = '',  help = 'Path to output directory.'),
  make_option('--tmpdir',    type = 'character', default = '',  help = 'Path to temporary directory.'),
  make_option('--nboots',    type = 'integer',   default = 100, help = 'Number of bootstraps (default 100).'),
  make_option('--bandwidth', type = 'integer',   default = 15,  help = 'Kernel bandwidth in pixels (default 15).'),
  make_option('--quantile',  type = 'integer',   default = 95,  help = 'Quantile for threshold (default 95).'),
  make_option('--replace',   type = 'integer',   default = 1,   help = 'Bootstrap with replacement (default 1).'),
  make_option('--joinall',   type = 'integer',   default = 1,   help = 'Joins distant delineations:\n 0: None\n 1: Joins delineations one pixel appart (default)\n 2: 1 and joins delineations up to two non-buildable pixels appart'),
  make_option('--joinunl',   type = 'integer',   default = 2,   help = 'Joins distant delineations:\n 0: None\n 1: Joins delineations one pixel appart (default)\n 2: 1 and joins delineations up to two non-buildable pixels appart'),
  make_option('--filter',    type = 'integer',   default = 1,   help = 'Removes delineations that are smaller than specified value (default 1)'),
  make_option('--workers',   type = 'integer',   default = -1,  help = 'Available cores (default -1 for all cores)'),
  make_option('--memory',    type = 'integer',   default = -1,  help = 'Available memory in GB (default -1 for all memory)'),
  make_option('--seed',      type = 'integer',   default = 1,   help = 'Bootstrap seed (default 1)')
)

params <- parse_args(OptionParser(usage = 'Computes urban delineations', option_list = arguments))
params <- subset(params, names(params) != 'help')

# Checks directories
if(!dir.exists(params$outdir)) dir.create(params$outdir)
if(!dir.exists(params$tmpdir)) dir.create(params$tmpdir)
unlink(dir(params$tmpdir, full.names = T), recursive = T)

# Checks parameters
tests <- list(
  density   = file.exists, 
  unlivable = file.exists, 
  outdir    = file.exists, 
  tmpdir    = file.exists,
  nboots    = function(.) is_greater_than(., 0),
  bandwidth = function(.) is_greater_than(., 1),
  quantile  = function(.) is_in(., 0:100),
  replace   = function(.) is_in(., 0:1),
  joinall   = function(.) is_weakly_greater_than(., 0),
  joinunl   = function(.) is_weakly_greater_than(., 0),
  filter    = function(.) is_weakly_greater_than(., 0),
  workers   = function(.) equals(., -1) | is_greater_than(., 0),
  memory    = function(.) equals(., -1) | is_greater_than(., 0),
  seed      = is.integer)

tests <- mapply(function(test, param) do.call(test, list(param)), tests, params)
if(sum(tests) != length(tests)) stop('Wrong argument(s): ', paste(names(tests)[which(!tests)], collapse = ', '), '\n')
rm(arguments, tests)

# Checks workers and memory
maxcor <- availableCores() - 1
maxmem <- Sys.meminfo()$totalram@size
params$workers <- ifelse(params$workers == -1, maxcor, min(params$workers, maxcor))
params$memory  <- ifelse(params$memory  == -1, maxmem, min(params$memory,  maxmem))
rm(maxcor, maxmem)

# Prints parameters
cat('\nComputes urban delineations (version 22/01/04)\n')
cat('\nParameters:', sprintf('- %-10s= %s', names(params), unlist(params)), sep = '\n')

# Sets up workers
cat('\nOperations:\n- Optimising computations')
plan(multiprocess, workers = params$workers, gc = T)

# Functions ---------------------------------------------------------------

# Displays cimg
display <- cmpfun(function(cimg, discrete = F) {
  ebimg <- EBImage::as.Image(cimg)
  if(discrete) ebimg <- EBImage::colorLabels(ebimg, normalize = T)
  EBImage::display(ebimg, method = 'raster')
})

# Reads raster as cimg
read_foo <- cmpfun(function(file, livable = NULL) {
  file <- raster(file)
  file <- as.cimg(file, maxpixels = ncell(file))
  if(!is.null(livable)) {
    file <- replace(file, !livable, 0)
  }
  return(file)
})

# Creates file names
filename_foo <- cmpfun(function(label, params) {
  filename <- gsub('.tif$', '', basename(params$density))
  filename <- file.path(params$outdir, paste0(filename, 'd', params$bandwidth, 'b', params$nboots, '_', label, '.tif'))
  return(filename)
})

# Writes cimg as raster
write_foo <- cmpfun(function(delineation, label, reference, params) {
  delineation <- setValues(reference, c(delineation))
  filename    <- filename_foo(label, params)
  writeRaster(delineation, filename, NAflag = 0, overwrite = T)  
})

# Optimises computations
optimise_foo <- cmpfun(function(livable, params) {
  bootsize   <- 8 * sum(livable) / 1024^3
  sliceindex <- ceiling(bootsize * params$nboots / (0.1 * params$memory))
  sliceindex <- round(seq(0, sum(livable), length.out = sliceindex + 1))
  usedisk    <- as.integer(length(sliceindex) > 2)
  params     <- modifyList(params, list(sliceindex = sliceindex, usedisk = usedisk))
  return(params)
})

# Computes bi-squared kernel
kernel_foo <- cmpfun(function(bandwidth) {
  size     <- ifelse(bandwidth %% 2 == 0, bandwidth + 1, bandwidth)
  kernel   <- matrix(0, size, size)
  centre   <- ceiling(size / 2)
  distance <- sqrt((col(kernel) - centre)^2 + (row(kernel) - centre)^2)
  distance <- ifelse(distance <=  bandwidth / 2, (1 - (distance / bandwidth * 2)^2)^2, 0)
  kernel   <- as.cimg(distance / sum(distance))
  return(kernel)
})

# Computes single bootstrap
bootstrap_foo <- cmpfun(function(density, livable, kernel, params) {
  bootstrap <- replace(density, livable, sample(density[livable], replace = params$replace))
  bootstrap <- convolve(bootstrap, kernel)
  bootstrap <- bootstrap[livable]
  return(bootstrap)
})

# Computes multiple bootstraps
bootstraps_foo <- cmpfun(function(density, livable, kernel, params) {
  if(params$usedisk) {
    bootstraps <- future_sapply(1:params$nboots, function(.) {
      bootstrap <- bootstrap_foo(density, livable, kernel, params)
      file      <- tempfile('boot_', params$tmpdir, '.fst')
      write_fst(data.frame(bootstrap), file, compress = 0)
      gc()
      return(file)
    }, future.seed = params$seed)
  } else {
    bootstraps <- future_replicate(params$nboots, bootstrap_foo(density, livable, kernel, params), future.seed = params$seed)
  }
  return(bootstraps)
})

# Computes quantiles
quantile_foo <- cmpfun(function(bootstraps, params) {
  index <- params$nboots * params$quantile / 100
  bootstraps <- matrix(bootstraps[order(row(bootstraps), bootstraps)], ncol = ncol(bootstraps), byrow = T)
  bootstraps <- rowMeans(bootstraps[, c(floor(index), ceiling(index))])
  return(bootstraps)
})

# Computes threshold
threshold_foo <- cmpfun(function(bootstraps, livable, params) {
  if(params$usedisk) {
    threshold <- lapply(head(seq(params$sliceindex), -1), function(i) {
      slice <- future_sapply(bootstraps, function(file) as.matrix(read_fst(file, from = params$sliceindex[i] + 1, to = params$sliceindex[i + 1])))
      slice <- quantile_foo(slice, params)
      gc()
      return(slice)
    })
    threshold <- do.call(c, threshold)
  } else {
    threshold <- quantile_foo(bootstraps, params)
  }
  threshold <- replace(as.cimg(livable), livable, threshold)
  return(threshold)
})

# Computes delineations
delineation_foo <- cmpfun(function(density, livable, threshold, kernel) {
  delineation <- convolve(density, kernel)
  delineation <- as.cimg(delineation > threshold)
  delineation <- replace(delineation, !livable, 0)
  return(delineation)
})

# Computes ranks
rank_foo <- cmpfun(function(identifier) {
    ids <- subset(identifier, identifier > 0)
    rnk <- rank(-tabulate(ids), ties = "first")
    rnk <- factor(ids, seq(max(ids)), rnk)
    rnk <- as.integer(levels(rnk))[rnk]
    identifier <- replace(identifier, identifier > 0, rnk)
    return(identifier)
})

# Computes identifiers
identifier_foo <- cmpfun(function(delineation, livable, params) {
  identifier <- delineation
  if(params$joinall > 0) {
    identifier <- identifier | (fill(delineation, params$joinall + 1) & !delineation) 
  }
  if(params$joinunl > 0) {
    identifier <- identifier | (fill(delineation, params$joinunl + 1) & !livable) 
  }
  identifier <- label(identifier, high_connectivity = T) + 1
  identifier <- replace(identifier, !delineation, 0)
  if(params$filter > 0) {
    filtered   <- seq(max(identifier))[tabulate(identifier) <= params$filter]
    identifier <- replace(identifier, identifier %in% filtered, 0)
  }
  identifier <- rank_foo(identifier)
  return(identifier)
})

# Computations ------------------------------------------------------------

tic('Runtime')

# Cities data
cat('\n- Loading data')
reference <- raster(raster(params$density))
livable   <- read_foo(params$unlivable)
livable   <- replace(livable, px.na(livable), 1) == 0
density   <- read_foo(params$density, livable)
kernel    <- kernel_foo(params$bandwidth)

# Cities computations
cat("\n- Computing cities:")
params      <- optimise_foo(livable, params)                       ; cat(" bootstraps...")
bootstraps  <- bootstraps_foo(density, livable, kernel, params)    ; cat(" thresholds...")
threshold   <- threshold_foo(bootstraps, livable, params)          ; cat(" delineations...")
delineation <- delineation_foo(density, livable, threshold, kernel); cat(" identifiers")
cities      <- identifier_foo(delineation, livable, params)
write_foo(threshold, 'ut', reference, params)                                 
write_foo(cities, 'ur', reference, params)
unlink(dir(params$tmpdir, full.names = T), recursive = T)
rm(bootstraps, threshold, delineation)

# Cores data
cat("\n- Computing cores: ")
livable <- imsub(cities > 0)
density <- replace(density, !livable, 0)

# Core computations
params      <- optimise_foo(livable, params)                       ; cat(" bootstraps...")
bootstraps  <- bootstraps_foo(density, livable, kernel, params)    ; cat(" thresholds...")
threshold   <- threshold_foo(bootstraps, livable, params)          ; cat(" delineations...")
delineation <- delineation_foo(density, livable, threshold, kernel); cat(" identifiers")
cores <- replace(cities, !delineation, 0)
write_foo(threshold, 'ct', reference, params)
write_foo(cores, 'co', reference, params)
unlink(params$tmpdir, recursive = T)
rm(bootstraps, threshold, delineation)

# Cities with cores computations
cat('\n- Computing cities with cores\n\n')
citycores <- replace(cities, !(cities %in% unique(subset(cities, cores > 0))), 0)
write_foo(citycores, 'cc', reference, params)

toc()

# # Testing -----------------------------------------------------------------
# setwd('~/Dropbox/research/delineation_public')
# params <- modifyList(params, list(
#   density   = 'input/density.tif',
#   unlivable = 'input/unlivable.tif',
#   outdir    = 'output',
#   tmpdir    = '~/Desktop/tmp',
#   nboots    = 1000))
# cat(paste0(c("Rscript", "~/Dropbox/research/delineation_public/delineation.R", paste0("--", names(params), "=", params)), collapse  =  " "))
# -------------------------------------------------------------------------