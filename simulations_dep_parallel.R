# Simulation study (under mixing)

require("cubature") # adaptive multidimensional integration (e.g., Cuhre)
require("doFuture") # registers future-based backends for foreach (%dopar%)
require("fs") # cross-platform filesystem ops (paths, files, dirs)
require("future.batchtools") # runs futures via batchtools (Slurm templates/resources)
require("ksm") # SPD kernel density tools (kdens_symmat, integrate_spd, bandwidth_optim)
require("parallel") # base R parallelism (PSOCK clusters, detectCores)

# Define the list of libraries to load on each cluster node

libraries_to_load <- c(
  "cubature",
  "doFuture",
  "fs",
  "future.batchtools",
  "parallel",
  "ksm" # ADDED
)

# Define the list of variables/functions to export to the worker nodes

vars_to_export <- c(
  "combo",
  "cores_per_node",
  "criteria",
  "ISE",
  "IAE",
  "kernels",
  "libraries_to_load",
  "logISE",
  "Mmod",
  "K",
  "models",
  "nobs",
  "path",
  "replication",
  "resources_list",
  "RR",
  "setup_parallel_cluster",
  "simu_rWAR",
  "simu_fdens_WAR",
  "Smod",
  "timing"
)

# Sets up a parallel cluster, loads necessary libraries, and exports required variables globally

setup_parallel_cluster <- function() {
  num_cores <<- detectCores()
  cl <<- makeCluster(num_cores)

  # Export the list of libraries to the worker nodes
  clusterExport(cl, varlist = "libraries_to_load")

  # Load necessary libraries on each cluster node
  invisible(clusterEvalQ(cl, {
    lapply(libraries_to_load, library, character.only = TRUE)
  }))

  # Export all necessary objects, functions, and parameters to the worker nodes
  clusterExport(cl, varlist = vars_to_export)

  return(cl) # Return the cluster object
}

# Initialize all variables in the list as NULL except vars_to_export and setup_parallel_cluster

invisible(
  lapply(
    vars_to_export[
      !(vars_to_export %in% c("vars_to_export", "setup_parallel_cluster"))
    ],
    function(x) assign(x, NULL, envir = .GlobalEnv)
  )
)

#' Simulation study serial dependence simulation from WAR(1)
#' @param n sample size
#' @param Mmod model for autoregressive coefficients matrix
#' @param Smod covariance matrix of VAR components
#' @param K degrees of freedom
#' @return a cube of dimension 2 by 2 by \code{n}.
simu_rWAR <- function(
  n,
  Mmod = c("M1", "M2", "M3"),
  Smod = c("S1", "S2", "S3"),
  K = 2L
) {
  Mmod <- match.arg(Mmod)
  Smod <- match.arg(Smod)
  n <- as.integer(n)
  K <- as.integer(K)
  stopifnot(n >= 1, K >= 2) # otherwise matrix not invertible / degenerate
  if (Mmod == "M1") {
    M <- matrix(
      c(0.9, 0, 1, 0),
      byrow = TRUE,
      nrow = 2,
      ncol = 2
    )
  } else if (Mmod == "M2") {
    M <- matrix(
      c(0.3, -0.3, -0.3, 0.3),
      byrow = TRUE,
      nrow = 2,
      ncol = 2
    )
  } else {
    M <- diag(rep(0.5, 2))
  }
  if (Smod == "S1") {
    S <- diag(2)
  } else if (Smod == "S2") {
    S <- matrix(0.5, 2, 2) +
      diag(rep(0.5, 2))
  } else {
    S <- matrix(
      c(1, 0.9, 0.9, 1),
      byrow = TRUE,
      nrow = 2,
      ncol = 2
    )
  }
  ksm::rWAR(
    n = n,
    M = M,
    S = S,
    K = K
  )
}
#' Simulation study asymptotic density for dependent
#' @param x array of positive definite matrices of size 2 by 2 by \code{n}.
#' @param Mmod model for autoregressive coefficients matrix
#' @param Smod covariance matrix of VAR components
#' @param K degrees of freedom
#' @return a vector of density functions
simu_fdens_WAR <- function(
  x,
  Mmod = c("M1", "M2", "M3"),
  Smod = c("S1", "S2", "S3"),
  K = 2L
) {
  Mmod <- match.arg(Mmod)
  Smod <- match.arg(Smod)
  n <- as.integer(n)
  K <- as.integer(K)
  stopifnot(n >= 1, K >= 2) # otherwise matrix not invertible / degenerate
  if (Mmod == "M1") {
    M <- matrix(
      c(0.9, 0, 1, 0),
      byrow = TRUE,
      nrow = 2,
      ncol = 2
    )
  } else if (Mmod == "M2") {
    M <- matrix(
      c(0.3, -0.3, -0.3, 0.3),
      byrow = TRUE,
      nrow = 2,
      ncol = 2
    )
  } else {
    M <- diag(rep(0.5, 2))
  }
  if (Smod == "S1") {
    S <- diag(2)
  } else if (Smod == "S2") {
    S <- matrix(0.5, 2, 2) +
      diag(rep(0.5, 2))
  } else {
    S <- matrix(
      c(1, 0.9, 0.9, 1),
      byrow = TRUE,
      nrow = 2,
      ncol = 2
    )
  }

  Sigma_inf <- ksm::Riccati(M = M, S = S)$solution
  c(ksm::dWishart(x, df = K, S = Sigma_inf, log = FALSE))
}

#' Integrated squared error of kernel density estimator for symmetric matrices
#'
#' Given a sample \code{x} from a target density, and the optimal bandwidth, compute the integrated squared error via numerical integration.
#' @param x sample of symmetric matrix observations from which to build the kernel density kernel
#' @param bandwidth double for the bandwidth of the kernel
#' @param fdens density of the target function from which points are drawn.
#' @param kernel string, one of \code{Wishart}, \code{smlnorm} (log-Gaussian) or \code{smnorm} (Gaussian).
#' @param tol double for tolerance of numerical integral
#' @param lb lower bound for integration range of eigenvalues
#' @param ub upper bound for integration range of eigenvalues
#' @param return double with value of the integrated squared error
#' @param ... additional parameters, currently ignored

logISE <- function(
  x,
  bandwidth,
  fdens,
  kernel = c("Wishart", "smlnorm", "smnorm"),
  tol = 1e-3,
  lb = 1e-5,
  ub = Inf,
  neval = 1e6L,
  method = c("suave", "hcubature"),
  ...
) {
  method <- match.arg(
    method[1],
    choices = c("hcubature", "pcubature", "cuhre", "divonne", "suave", "vegas")
  )
  kernel <- match.arg(kernel)
  dim <- ncol(x)
  stopifnot(
    dim %in% c(2L, 3L),
    length(lb) == 1L,
    length(ub) == 1L,
    lb >= 0,
    ub > lb,
    tol > 0,
    is.numeric(bandwidth) & length(bandwidth) == 1L
  )
  integrand <- function(vars, kernel) {
    # Construct SPD matrix S
    if (length(vars) == 3L) {
      S <- ksm::rotation_scaling(vars[1], vars[2:3])
      jacobian_value <- abs(vars[2] - vars[3]) / 4
    } else {
      S <- ksm::rotation_scaling(vars[1:3], vars[4:6])
      jacobian_value <- abs(vars[4] - vars[5]) *
        abs(vars[4] - vars[6]) *
        abs(vars[5] - vars[6]) *
        sin(vars[2]) /
        24
    }
    S <- array(S, dim = c(nrow(S), ncol(S), 1))

    # Compute the squared difference between the kernel and the target density
    diff_squared <- (ksm::kdens_symmat(
      x = S,
      xs = x,
      b = bandwidth,
      kernel = kernel,
      log = FALSE
    ) -
      fdens(S))^2

    # Return the product of diff_squared and the Jacobian factor
    return(diff_squared * jacobian_value)
  }
  if (dim == 2L) {
    # Define the limits of integration
    lower_limit <- c(0, rep(lb, length.out = dim))
    upper_limit <- c(2 * pi, rep(ub, length.out = dim))
  } else if (dim == 3L) {
    lower_limit <- c(rep(0, 3), rep(lb, length.out = dim))
    upper_limit <- c(2 * pi, pi, 2 * pi, rep(ub, length.out = dim))
  }
  result <- cubature::cubintegrate(
    f = integrand,
    lower = lower_limit,
    upper = upper_limit,
    method = method,
    relTol = 1e-3,
    maxEval = neval,
    kernel = kernel
  )

  # print(result)
  # Return the result of the integral
  return(log(result$integral))
}

##################
## Set the path ##
##################

path <- getwd()
# setwd(path)

##############################
## Parallelization on nodes ##
##############################

cores_per_node <- 64 # number of cores for each node in the super-computer

resources_list <- list(
  cpus_per_task = cores_per_node,
  mem = "240G",
  walltime = "0:05:00",
  nodes = 1
  # Omit 'partition' to let SLURM choose
)

##############################
## Hyper-parameters         ##
##############################

RR <- 1:128 # replications
nobs <- c(125L)
Mmod <- c("M1")
Smod <- c("S1")
kernels <- c("Wishart", "smlnorm")
K <- 4

###############
## Main code ##
###############

.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))

# Disable the check for random number generation misuse in doFuture
options(doFuture.rng.onMisuse = "ignore")

# Register the doFuture parallel backend
doFuture::registerDoFuture()

# Tweak the batchtools_slurm with the custom template and resources
myslurm <- future::tweak(
  batchtools_slurm,
  template = "batchtools.slurm.dependent.tmpl",
  resources = resources_list
)

# Set the plan for future
future::plan(list(myslurm, multisession))

# Create empty data frames to store the results
raw_results <- data.frame(
  nobs = integer(),
  Mmod = integer(),
  Smod = integer(),
  kernel = character(),
  logISE = numeric(),
  bandwidth = numeric(),
  timing = numeric()
)

# Capture the start time
start_time <- Sys.time()

# Split RR into node-sized blocks
RR_blocks <- split(RR, ceiling(seq_along(RR) / cores_per_node))

# For each over blocks (one Slurm job per block); inside each, parallelize reps across node cores
res <- foreach::foreach(
  block = RR_blocks,
  .combine = "rbind",
  .export = vars_to_export,
  .packages = libraries_to_load
) %dopar%
  {
    # spin up a per-node PSOCK cluster
    cl <- setup_parallel_cluster()

    block_results <- data.frame(
      nobs = integer(),
      Mmod = integer(),
      Smod = integer(),
      kernel = character(),
      logISE = numeric(),
      bandwidth = numeric(),
      timing = numeric()
    )

    # each r in this block runs in parallel across the node's cores
    r_list <- parallel::parLapply(cl, X = block, fun = function(r_local) {
      set.seed(r_local)

      out_rows_all <- vector(
        "list",
        length(nobs) * length(Mmod) * length(Smod)
      )
      idx_out <- 0L

      for (i in seq_along(nobs)) {
        for (j in seq_along(Mmod)) {
          for (k in seq_along(Smod)) {
            local_rows <- vector("list", length(kernels))

            xs <- simu_rWAR(
              n = nobs[i],
              Mmod = Mmod[j],
              Smod = Smod[k],
              K = K # degrees of freedom must be at least d
            )
            for (l in seq_along(kernels)) {
              start_time <- Sys.time()
              band <- ksm:::bandwidth_optim(
                x = xs,
                criterion = "lscv",
                h = ceiling(nobs[i]^0.25), # adjustment for serial dependence
                kernel = kernels[l]
              )

              log_ISE <- try(
                logISE(
                  x = xs,
                  bandwidth = band,
                  fdens = function(x) {
                    simu_fdens_WAR(x = x, Mmod = Mmod[j], Smod = Smod[k], K = K)
                  },
                  kernel = kernels[k],
                  lb = 1e-5, # avoid problems with matrices nearly non PSD, which are caught by sensitive Cpp routine
                  ub = Inf,
                  tol = 1e-3,
                  method = "hcubature"
                ),
                silent = FALSE
              )
              if (inherits(log_ISE, "try-error")) {
                log_ISE <- NA
              }
              timing <- as.numeric(difftime(
                Sys.time(),
                start_time,
                units = "secs"
              ))

              local_rows[[l]] <- data.frame(
                nobs = nobs[i],
                Mmod = j,
                Smod = k,
                kernel = kernels[l],
                logISE = log_ISE,
                bandwidth = band,
                timing = timing
              )
            }

            idx_out <- idx_out + 1L
            out_rows_all[[idx_out]] <- do.call(rbind, local_rows)
          }
        }
      }

      do.call(rbind, out_rows_all)
    })

    block_results <- rbind(block_results, do.call(rbind, r_list))

    # tear down the per-node cluster
    try(parallel::stopCluster(cl), silent = TRUE)

    block_results
  }

# Combine results from all nodes
raw_results <- res

# Stop parallel execution
future::plan(sequential)

# Calculate the duration in minutes
elapsed_time_minutes <- as.numeric(difftime(
  Sys.time(),
  start_time,
  units = "mins"
))
print(paste("Elapsed time:", round(elapsed_time_minutes, 2), "minutes"))

# Save the raw results to a CSV file
raw_output_file <- file.path(path, "raw_ISE_dep.csv")
write.csv(raw_results, raw_output_file, row.names = FALSE)

print("Raw results saved to raw_ISE_dep.csv")

#########################
## Process the results ##
#########################

# Build a summary table by (nobs, Mmod, Smod, kernel)

# Guard against accidental factor coercion
raw_results$kernel <- as.character(raw_results$kernel)

# Split by grouping keys
.grp <- interaction(
  raw_results$nobs,
  raw_results$Mmod,
  raw_results$Smod,
  raw_results$kernel,
  drop = TRUE
)

summ_list <- lapply(split(raw_results, .grp), function(df) {
  data.frame(
    nobs = df$nobs[1],
    Mmod = df$Mmod[1],
    Smod = df$Smod[1],
    kernel = df$kernel[1],
    # ISE stats
    mean_ISE = mean(df$logISE, na.rm = TRUE),
    sd_ISE = sd(df$logISE, na.rm = TRUE),
    median_ISE = median(df$logISE, na.rm = TRUE),
    IQR_ISE = IQR(df$logISE, na.rm = TRUE),
    mean_time_min = mean(df$timing, na.rm = TRUE) / 60,
    stringsAsFactors = FALSE
  )
})

summary_results <- do.call(rbind, summ_list)
row.names(summary_results) <- NULL

# Save the summary to CSV next to the raw file
summary_output_file <- file.path(path, "summary_dep.csv")
write.csv(summary_results, summary_output_file, row.names = FALSE)

print("Summary results saved to summary_dep.csv")
