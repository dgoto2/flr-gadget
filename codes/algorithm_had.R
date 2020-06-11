###############################################################################
# Testing the FLR/mse-Gadget MSE framework with the simple haddock model
# updated: 11-jun-2020
###############################################################################

# Install required packages (if not already)
remotes::install_github("REDUS-IMR/gadget", ref = "gadgetr")
#remotes::install_github("flr/mse")

# Load packages
library(mse)
library(dplyr)
library(FLa4a)
library(FLash)
library(FLAssess)
library(ggplotFL)
library(FLBRP)
library(FLCore)
library(MASS)
library(FLSAM)
library(filelock)
#library(profvis) # Performance measurement

#=====================================================================
# Set up an MSE loop
#=====================================================================
runOneTimeline <- function(iterSim, saveRaw) {
	
	# Set seed
	set.seed(0)

	# Create directories for the Gadget OM and helper functions
	codeDir <- paste0(homeDir, "codes")
	paramFileDir <- paste0(homeDir, "paramfiles")

	# Load GadgetR locally
	library(gadgetr)

	# Load helper functions
	source(paste0(codeDir, "/gadget-fls.R"), local = T)
	source(paste0(codeDir, "/assessment.R"), local = T)

	# Load the stock and files
	# Set a directory for the Gadget OM 
	setwd(paste0(homeDir, "/models/", modelName))

	# Load Gadget parameters
	paramsfile <- "refinputfile"
	gadget(c("-s", "-main", "main", "-i", paramsfile))

	# Initialize simulation
	initSim()

	# Set simulation time parameters
	firstYear <- 1978
	projYear <- 2000 
	finalYear <- 2020 
	
	# Set up a stock (or stocks) for Gadget 
	stockList <- c("had")
	had.fleets <- c("comm", "survey", "future")
	had.stocks <- c("had")
	had.stocks.mature <- c("had")
	had.surveys <- c("survey")
	had.forecasts <- c("future")
	
	# Set stock parameters for the OM 
	had.params <- list(stockStep = 2, minage = 1, maxage = 10, minfbar = 2, maxfbar = 8, 
	                   startf = 0.56, endf = 0.65, areas = c(1), m1 = c(0.2), m2 = NULL)
	# m2 = NULL means we calculate m2 from gadget result, m2=0 means we use only residual mortality (m1)
	# StockStep: in which step we should observe the stock number
	
	# Set quarterly allocation of TAC 
	had.forecasts.tac.proportion <- c(0.232, 0.351, 0.298, 0.119)
	
	# Set HCR function parameters (in this example, ICES HCR)
	had.hcr.params <- list(method = ices.hcr, args = list(blim = 100000, bsafe = 200000, fmin = 0.05, ftrg = 0.25))
	
	# Set recruitment parameters
	# If read csv (as a data frame), it will apply the values accordingly
	# If a constant value, it will apply the value as mux
	# If NULL, it will leave the recruitment params as it is
	#had.recruit.params <- read.csv(paste0(paramFileDir, "/muxfactors_cod.csv"))
	had.recruit.params <- 90000000 # constant recruit number
	#had.recruit.params <- NULL ## recruit params
	
	# Set assessment model parameters (truePlusNoise, SCAA, or SAM)
	# Select an assessment model
	had.assessment <- "SAM" ## "truePlusNoise", "SCAA" or "SAM"
	
	# If truePlusNoise is chosen, the noise using the residual params will be applied to stock only
	# Set noise parameters
	noise_age <- data.frame(age = c (1:10), mean = array(0.001283, dim = c(10,1)), sd = array(0.053986, dim = c(10,1))) # generate simple noise
	had.residual.params.catch <- noise_age
	had.residual.params.index <- noise_age
	had.residual.params.stock <- noise_age ########################???
	had.residual.params.mean.stock <- NULL
	had.residual.params.vcratios.stock <- NULL
	# If SCAA is chosen, the noise will be applied to both catch and index
	# If no error is applied
	#had.residual.params <- NULL
	had.noteating.forecast <- FALSE

	# Set Gadget output 
	gadgetOut <- list()
  # performance measurements
  #test <- profvis({
	# Run Gadget until the projection year - 1 for conditioning the OM
	gadgetOut <- runUntil(projYear - 1)
  #})
	#print(test)
  #browser()
	
	# Set up a projection loop
	prepareStock  <- function(stockNameGl) {
    
	  #------------------------------------------------------------------------------
	  # OM object
	  #------------------------------------------------------------------------------
	  # Use Gadget output to condition the OM 
		gadget.ret <- gadgetOut[[stockNameGl]]
		stk <- gadget.ret$stk
		idx <- FLIndices(a = gadget.ret$idx)

    # Set MSE simualtion parameters
		it <- 1                     # iterations
		fy <- finalYear             # final year
		y0 <- range(stk)["minyear"] # initial data year
		dy <- range(stk)["maxyear"] # final data year
		iy <- projYear              # initial year of projection (also intermediate)
		ny <- fy - iy + 1           # number of years to project from intial year
		nsqy <- 3                   # number of years to compute status quo metrics
		vy <- ac(iy:fy)             # vector of years to be projected
		management_lag <- 1         # For ICES HCR

    # Use short-term forecasts to create an initial stock object
		stk <- stf(stk, fy-dy, nsqy, nsqy)

		# Set fleet parameters
		#fb <- mseCtrl(method = hyperstability.fb, args = list(beta = 0.8))
    
		# Create an OM object 
		om <- FLom(stock = stk)#, fleetBehaviour = fb)
		#save(om, it, fy, y0, dy, iy, ny, nsqy, vy, fit, file = "om.RData")
		
		#-----------------------------------------------------------------------------
		# OEM object
		#-----------------------------------------------------------------------------
		# Set survey indices
		idx <- FLIndices(a = gadget.ret$idx)
		stk <- stock(om)
		stk0 <- stk

		idcs <- FLIndices()		# Use all indices
		
		for (i in 1:length(idx)){
		  
			# Estimate survey catchability from the a4a fit (w/o simulation) 
      # the index is based on 01 January abundances
			lst <- mcf(list(idx[[i]]@index, stock.n(stk0)))
			idx.lq <- log(lst[[1]]/lst[[2]])
			idx.qmu <- idx.qsig <- stock.n(iter(stk, 1)) ## set a data frame for indices
			idx.qmu[] <- yearMeans(idx.lq) ## use constant catchability
			idx.qsig[] <- sqrt(yearVars(idx.lq))
			idx.q <- FLQuant(NA, dimnames = dimnames(stock.n(stk)))
			
			# Estimate survey catchability based on lognormal distribution with mean and sd calculated above
			idx.q <- rlnorm(it, idx.qmu, idx.qsig)
			#idx.q[, ac(y0:iy)] <- idx.q[,ac(y0:iy)]
			idx_temp <- idx.q * stock.n(stk)
			
			## Generate an initial index
			idx_temp <- FLIndex(index = idx_temp, index.q = idx.q)
			range(idx_temp)[c("startf", "endf")] <- c(0, 0)
			idcs[[i]] <- idx_temp
		}
		names(idcs) <- names(idx)
		#idx <- FLIndices(a=idcs$a)
		idxDev <- lapply(idcs, index.q)
		names(idxDev) <- "index.q"
		stkDev <- FLQuant()
		dev <- list(idx = idxDev, stk = stkDev)
		obs <- list(idx = idcs[1], stk = stk)
		
		## Set deviances for catch.n ################################################################???
		#catch.dev <- log(catch.n(stk))
		#catch.dev <- catch.dev-iterMeans(catch.dev)
		#Sig <- apply(catch.dev[,ac(y0:dy),1,1,,drop=TRUE], 3, function(x) cov(t(x)))
		#Sig <- apply(Sig, 1, mean)
		#Sig <- matrix(Sig, ncol=dim(catch.dev)[1])
		#catch.dev[,ac(vy)][] <- t(mvrnorm(it * length(vy), rep(0, nrow(Sig)), Sig))
		#catch.dev <- exp(catch.dev)
		
		# Create an OEM object
		#oem <- FLoem(method = sampling.oem, args = list(oe = "index"), observations = obs, deviances = dev)
		oem <- FLoem()
		#save(oem, file = "oem.RData")

		# Set implementation error parameters
		#iem <- FLiem(method=noise.iem, args=list(fun="rlnorm", mean=0, sd=0.1, multiplicative=TRUE))
		iem <- FLiem()

		#-----------------------------------------------------------------------
		# Management procedure (MP) object
		#-----------------------------------------------------------------------
    # Set MP parameters
		mpPars <- list(seed = 1234, fy = fy, y0 = y0, dy = dy, iy = iy, management_lag = management_lag, nsqy = nsqy, it = it)

		# Set simulation scenarios
		
		# Tell stocks to stop eating (if requested) #####################################################???
		if(eval(parse(text=paste0(stockNameGl, ".noteating.forecast")))){
			stockCat <- eval(parse(text = paste0(stockNameGl, ".stocks")))
			tmp <- lapply(stockCat, stopEating)
			print("Stocks stop eating now")
			print(tmp)
		}

		## Set stock-specific HCR parameters
		hcrParams <- eval(parse(text=paste0(stockNameGl, ".hcr.params")))

		## Set assessment model parameters
		saParam <- eval(parse(text = paste0(stockNameGl, ".assessment")))
		if(saParam == "truePlusNoise") saMethod <- truePlusNoise.sa
		else if(saParam == "SCAA") saMethod <- sca.sa
		else if(saParam == "SAM") saMethod <- sam.sa

    # Set MP object parameters
		ctrl <- list(hcr = mseCtrl(method = hcrParams[["method"]], args = hcrParams[["args"]]), 
		             isys = mseCtrl(method = tac.is), est = mseCtrl(method = saMethod))

		# Set scenario names
		scenarioName <- paste0(stockNameGl, ".", "iter", iterSim)

		return(list(opModel = om, indices = idx, obsModel = oem, impModel = NULL, ctrl.mp = ctrl, 
		            mpPars = mpPars, scenario = scenarioName, tracking = NULL))
	}

	# Load helper functions
	source(paste0(codeDir, "/gadget-fwd.R"), local = T)
	source(paste0(codeDir, "/mp-methods-gadget.R"), local = T)

	# 
	inputPre <- lapply(stockList, prepareStock) ################################???
	names(inputPre) <- stockList
	res <- mp.gadget(inputPre)

	#return(list(mseResults = res, gadgetResults = gadgetOut))
	return(list(mseResults = res))
}


#======================================================================
# Run MSE simulations
#======================================================================
# Enable below to run directly from R shell
combIndex <- 1
iterIndex <- 1

# Set global variables
homeDir <- paste0(getwd(), "/")
modelName <- "had"
saveAllRawData <- FALSE

# Read effort combination
#print(paste("combination no.", combIndex, "iteration", iterIndex))
#fComb <- read.csv(paste0(homeDir, "paramfiles/effort_combination.csv"))

# Run with combination and iterIndex
resultFinal <- runOneTimeline(iterIndex, saveAllRawData)

# Name an output file
outFileName <- paste0(homeDir, "/results-combination", combIndex, ".rds")

# Use lock to prevent race condition when combining results
lck <- lock(paste0(outFileName, ".lock"))

# Check if any output files exist
if(file.exists(outFileName)) {
  # Load existing results
	allResults <- readRDS(outFileName)
} else {
  # First result
	allResults <- list()
}

# Combine all results
allResults[[iterIndex]] <- resultFinal

# Plot outputs 
stk.plot <- plot(FLStocks(stk.om = allResults[[iterIndex]]$mseResults$had$mse@stock, 
                          stk.mp = allResults[[iterIndex]]$mseResults$had$sa.result$stk0)) + 
  theme(legend.position = "top") + geom_vline(aes(xintercept = 2000))
stk.plot

# Save combination info
write.table(1, file = (paste0(outFileName,".info.txt")))

# Save it back
saveRDS(allResults, file = outFileName)

# Unlock the file
unlock(lck)
