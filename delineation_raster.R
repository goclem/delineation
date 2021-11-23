#!/usr/bin/env Rscript
# Description: Computes urban delineations from raster
# Author: Clement Gorin
# Contact: gorinclem@gmail.com
# Version: 2021.11.23

# Packages
suppressMessages(if(!require('pacman')) install.packages('pacman', repos = 'https://cloud.r-project.org/'))
pacman::p_load(compiler, fst, future.apply, imager, matrixStats, memuse, optparse, raster, rgdal, tictoc)
options(warn = -1, future.globals.maxSize = (2048 * 1024^2))

# Parameters
arguments <- list(
  make_option('--density',   type = 'character', default = '',  help = 'Path to density raster'),
  make_option('--unlivable', type = 'character', default = '',  help = 'Path to unlivable raster'),
  make_option('--outdir',    type = 'character', default = '',  help = 'Path to output directory'),
  make_option('--tmpdir',    type = 'character', default = '',  help = 'Path to temporary directory'),
  make_option('--nboots',    type = 'integer',   default = 100, help = 'Number of bootstrapped counterfactual densities (default 100)'),
  make_option('--niter',     type = 'integer',   default = 1,   help = 'Number of iteration for the cores (default 1)'),
  make_option('--bandwidth', type = 'integer',   default = 15,  help = 'Kernel bandwidth in pixels (default 15)'),
  make_option('--quantile',  type = 'integer',   default = 95,  help = 'Quantile for threshold (default 95)'),
  make_option('--replace',   type = 'integer',   default = 1,   help = 'Bootstrap with replacement (default 1)'),
  make_option('--joinall',   type = 'integer',   default = 0,   help = 'Joins delineations up to [value] pixels apart (default 0)'),
  make_option('--joinunl',   type = 'integer',   default = 0,   help = 'Joins delineations up to [value] unlivable pixels apart (default 0)'),
  make_option('--filter',    type = 'integer',   default = 0,   help = 'Keeps delineations that are larger than [value] pixels (default 0)'),
  make_option('--workers',   type = 'integer',   default = -1,  help = 'Number of available cores (default -1 for all cores)'),
  make_option('--memory',    type = 'integer',   default = -1,  help = 'Available memory in GB (default -1 for all memory)'),
  make_option('--seed',      type = 'integer',   default = 1,   help = 'Bootstrap seed (default 1)')
)

params <- parse_args(OptionParser(usage = 'Computes urban delineations from raster', option_list = arguments))
params <- subset(params, names(params) != 'help')

# (!) Testing only --------------------------------------------------------
# setwd('~/github/delineation')
# # Interactive
# params <- modifyList(params, list(
#   density   = 'input/raster/density.tif',
#   unlivable = 'input/raster/unlivable.tif',
#   outdir    = 'output/raster',
#   tmpdir    = 'temporary'
#   ))
# # Command line
# cat(paste(c('Rscript', 'delineation_raster.R', paste0('--', names(params), '=', params)), collapse =  ' '))
# -------------------------------------------------------------------------

# Checks directories
if(!dir.exists(params$outdir)) dir.create(params$outdir)
if(!dir.exists(params$tmpdir)) dir.create(params$tmpdir)
unlink(dir(params$tmpdir, pattern = '^boot_.*\\.fst$', full.names = T), recursive = T)

# Checks parameters
tests <- list(
  density   = file.exists,
  unlivable = file.exists,
  outdir    = file.exists,
  tmpdir    = file.exists,
  nboots    = function(.) is_greater_than(., 0),
  niter     = function(.) is_greater_than(., 0),
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
maxmem <- round(Sys.meminfo()$totalram@size)
params$workers <- ifelse(params$workers == -1, maxcor, min(params$workers, maxcor))
params$memory  <- ifelse(params$memory  == -1, maxmem, min(params$memory,  maxmem))
rm(maxcor, maxmem)

# Prints parameters
cat('\nComputes urban delineations from raster (version 2021.11.23)\n')
cat('\nParameters:', sprintf('- %-10s= %s', names(params), unlist(params)), sep = '\n')

# Sets up workers
cat('\nOperations:')
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
  image <- raster(file)
  image <- as.cimg(image, maxpixels = ncell(image))
  if(!is.null(livable)) {
    image <- pad(image, (nrow(livable) - nrow(image)), 'xy')
    image <- replace(image, !livable, 0)
  }
  return(image)
})

# Creates file names
filename_foo <- cmpfun(function(label, params) {
  filename <- sub('\\.csv$|\\.tif$', '', basename(params$density))
  filename <- file.path(params$outdir, paste0(filename, 'd', params$bandwidth, 'b', params$nboots, '_', label, '.tif'))
  return(filename)
})

# Writes cimg as raster
write_foo <- cmpfun(function(image, label, params, navalue = -1) {
  container <- raster(params$density) < 0
  image     <- crop.borders(image, nPix = ((nrow(image) - ncol(container)) / 2))
  image     <- setValues(container, c(image))
  image     <- mask(image, container)
  filename  <- filename_foo(label, params)
  writeRaster(image, filename, NAflag = navalue, overwrite = T)  
})

# Optimises computations
optimise_foo <- cmpfun(function(livable, params, memshr = 0.1) {
  bootsize   <- 8 * sum(livable) / 1024^3
  sliceindex <- ceiling(bootsize * params$nboots / (memshr * params$memory))
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
  distance <- ifelse(distance <= bandwidth / 2, (1 - (distance / bandwidth * 2)^2)^2, 0)
  kernel   <- as.cimg(distance / sum(distance))
  return(kernel)
})

# Computes single bootstrap
bootstrap_foo <- cmpfun(function(density, livable, kernel, params) {
  bootstrap <- replace(density, livable, sample(density[livable], replace = params$replace))
  bootstrap <- convolve(bootstrap, kernel)
  bootstrap <- subset(bootstrap, livable)
  return(bootstrap)
})

# Computes single bootstrap for cores
bootstrap_cores_foo <- cmpfun(function(density1, livable1, density2, livable2, kernel, params) {
  bootstrap1 <- replace(density1, livable1, sample(density1[livable1], replace = params$replace))
  bootstrap2 <- replace(density2, livable2, sample(density2[livable2], replace = params$replace))
  bootstrap  <- bootstrap1 + bootstrap2
  bootstrap  <- convolve(bootstrap, kernel)
  bootstrap  <- subset(bootstrap, livable)
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

# Computes multiple bootstraps for cores
bootstraps_cores_foo <- cmpfun(function(density1, livable1, density2, livable2, kernel, params) {
  if(params$usedisk) {
    bootstraps <- future_sapply(1:params$nboots, function(.) {
      bootstrap <- bootstrap_cores_foo(density1, livable1, density2, livable2, kernel, params)
      file      <- tempfile('boot_', params$tmpdir, '.fst')
      write_fst(data.frame(bootstrap), file, compress = 0)
      gc()
      return(file)
    }, future.seed = params$seed)
  } else {
    bootstraps <- future_replicate(params$nboots, bootstrap_cores_foo(density1, livable1, density2, livable2, kernel, params), future.seed = params$seed)
  }
  return(bootstraps)
})

# Computes threshold
threshold_foo <- cmpfun(function(bootstraps, livable, params) {
  if(params$usedisk) {
    threshold <- lapply(head(seq(params$sliceindex), -1), function(i) {
      slice <- future_sapply(bootstraps, function(file) as.matrix(read_fst(file, from = params$sliceindex[i] + 1, to = params$sliceindex[i + 1])))
      slice <- rowQuantiles(slice, probs = (params$quantile / 100))
      gc()
      return(slice)
    })
    threshold <- do.call(c, threshold)
  } else {
    threshold <- rowQuantiles(bootstraps, probs = (params$quantile / 100))
  } 
  threshold <- replace(as.cimg(livable), livable, threshold)
  return(threshold)
})

# Computes delineations
delineation_foo <- cmpfun(function(density, livable, threshold, kernel) {
  delineation <- convolve(density, kernel)
  delineation <- replace(delineation, !livable, 0)
  delineation <- as.cimg(delineation > threshold)
  return(delineation)
})

# Computes ranks
rank_foo <- cmpfun(function(identifier) {
  values     <- subset(identifier, identifier > 0)
  position   <- rank(-tabulate(values), ties = 'first')
  position   <- position[match(values, seq(max(values)))]
  identifier <- replace(identifier, identifier > 0, position)
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

# Urban areas -------------------------------------------------------------

tic('Runtime')

# Urban areas data
cat('\n- Computing urban areas: data...')
kernel  <- kernel_foo(params$bandwidth)
livable <- read_foo(params$unlivable)
livable <- replace(livable, px.na(livable), 1) == 0
livable <- pad(livable, nPix = params$bandwidth, 'xy')
density <- read_foo(params$density, livable)

# Urban areas computations
params      <- optimise_foo(livable, params)                       ; cat(' bootstraps...')
bootstraps  <- bootstraps_foo(density, livable, kernel, params)    ; cat(' thresholds...')
threshold   <- threshold_foo(bootstraps, livable, params)          ; cat(' delineations...')
delineation <- delineation_foo(density, livable, threshold, kernel); cat(' identifiers...')
urban       <- identifier_foo(delineation, livable, params)        ; cat(' saving')

# Writes urban areas
write_foo(threshold, 'ut', params)
write_foo(urban, 'ur', params, navalue = 0)
if(params$usedisk) unlink(bootstraps, recursive = T)
rm(bootstraps, threshold, delineation)

# Urban cores -------------------------------------------------------------

# Reference delineation
previous <- urban

for(iter in seq(params$niter)) {
  # First iteration is not indexed
  corelab <- ifelse(iter == 1, '', as.character(iter))
 
  # Urban cores data
  cat('\n- Computing urban cores: data...')
  livable1 <- previous > 0
  livable2 <- convolve(livable1, kernel) > 0 
  livable2 <- replace(livable2, !livable | livable1, F)
  livable  <- livable1 | livable2
  density1 <- replace(density, !livable1, 0)
  density2 <- replace(density, !livable2, 0)
  
  # Cores computations
  params      <- optimise_foo(livable, params)                                                ; cat(' bootstraps...')
  bootstraps  <- bootstraps_cores_foo(density1, livable1, density2, livable2, kernel, params) ; cat(' thresholds...')
  threshold   <- threshold_foo(bootstraps, livable, params)                                   ; cat(' delineations...')
  delineation <- delineation_foo(density, livable, threshold, kernel)                         ; cat(' identifiers...')
  delineation <- replace(delineation, previous == 0, 0)
  cores       <- identifier_foo(delineation, livable, params)                                 ; cat(' saving')
  
  # Writes cores
  write_foo(threshold, sprintf('ct%s', corelab), params)
  write_foo(cores, sprintf('co%s', corelab), params, navalue = 0)
  if(params$usedisk) unlink(bootstraps, recursive = T)  
  rm(bootstraps, threshold, delineation, density1, density2, livable1, livable2)
  
  # Urban areas with cores
  cat('\n- Computing urban areas with cores\n\n')
  urbancores <- urban %in% unique(subset(urban, cores > 0))
  urbancores <- replace(urban, !urbancores, 0)
  write_foo(urbancores, sprintf('cc%s', corelab), params, navalue = 0)
  rm(corelab)
  
  # Updates reference
  previous <- cores
}

toc()