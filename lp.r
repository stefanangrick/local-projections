#' Calculate and plot local projection-based impulse responses
#'
#' Calculates impulse responses from local projections based on single-equation
#' least squares models.
#'
#' @param data Data frame containing variables to be analysed.
#' @param y.var Character string. Specifies response variable.
#' @param shock.var Character vector. Specifies one or more shocks. These will
#'        be included in the equation as contemporaneous terms. Note that if you
#'        include y.var here, you will receive an error message for a 'near
#'        perfect fit' when horizon 0 is calculated, since a shock of one
#'        variable to itself is one.
#' @param shock.var.endo Character vector. Specifies one or more shocks which
#'        are presumed to be endogenous. For each endogenous shock variable, an
#'        equal or greater number of instruments (iv.var) is required.
#' @param iv.var Character vector. Specifies one or more instruments for
#'        endogenous shock variables (shock.var.endo). Instrumental variables
#'        disable automatic lag order selection and stargazer output of
#'        regression summaries.
#' @param contr.var Character vector. Specifies control variables. These make up
#'        the lag terms in the regression models. If no control variables are
#'        specified explicitly, the variables from the data frame not contained
#'        in y.var and shock.var are automatically added as controls. Note that
#'        if you want to include lags of the response variable and shock
#'        variables, this must be explicitly stated here. Also note that if you
#'        want to include contemporaneous control terms, add them to shock.var
#'        instead.
#' @param state.var Character vector. Specifies binary state variable used to
#'        calculate state-dependent responses. This must be a binary integer
#'        dummy containing ones (1) and zeros (0).
#' @param cum.y Logical. If TRUE, the cumulative response of the response
#'        variable from time t to t+h will be calculated.
#' @param const Logical. If TRUE, a constant will be included. Default is TRUE.
#' @param trend Logical. If TRUE, a linear trend will be included. Default is
#'        FALSE.
#' @param trend.sq Logical. If TRUE, a squared trend term will be inculuded.
#'        Default is FALSE.
#' @param trend.cu Logical. If TRUE, a cubic trend will be included. Default is
#'        FALSE.
#' @param trend.qu Logical. If TRUE, a quartic trend will be included. Default
#'        is FALSE.
#' @param n.lag Integer. The number of lags for control variables (default = 2).
#' @param infocrit Character string. Specifies the information criterion to be
#'        used to automatically determine lag order at each horizon. One of
#'        "AIC", "BIC", or "HQC". Criteria will check the optimal lag order up
#'        to n.lag.
#' @param h Integer. Specifies the number of impulse response horizons.
#' @param ci Numeric. Three numbers between 0 and 1 that specify the confidence
#'           interval for error bands. Defaults are 0.95 (95%), 0.68 (68%), and
#'           0.38 (38%).
#' @param lw Numeric. Specifies line width for plots.
#' @param overplot Logical. If TRUE, if state-dependent responses will be
#'        plotted on top of each other. Requires specification of a state
#'        variable.
#' @param output Character or file path. Specifies location where summary of
#'        all regression models (one per horizon) is to be saved.
#' @param debug Logical. Specifies if the model parameters should be printed.
#'        For debugging.
#'
#' @keywords local projection, impulse response, time series
#'
#' @return Returns a list containing the coefficients (coefs),
#'         heteroskedasticity and autocorrelation-consistent (HAC) standard
#'         errors (se), a list of regression models (lmm) for each impulse
#'         response horizon, and a list of ggplot2 plot objects containing the
#'         impulse response graphs (g).
#'
#' @export

lpirf <- function(data, y.var, shock.var = NULL, shock.var.endo = NULL,
                  iv.var = NULL, contr.var = NULL, state.var = NULL,
                  cum.y = FALSE, const = TRUE, trend = FALSE, trend.sq = FALSE,
                  trend.cu = FALSE, trend.qu = FALSE, n.lag = 2,
                  infocrit = NULL, hor = 12, ci = c(0.95, 0.68, 0.38),
                  lw = 1, overplot = FALSE, output = FALSE, debug = FALSE) {
  # Check parameter values
  y.var <- y.var[1]
  
  # Horizons specify how far out we're going, starting after horizon 0
  hor <- hor + 1
  
  # Check if endogenous variables are given for IV regression
  if (!is.null(iv.var)) {
    if (is.null(shock.var.endo)) {
      stop(paste0("Need to specify at least one endogenous shock variable for",
                  " IV regression."))
    }
  }
  
  # Check endogenous variables
  if (!is.null(shock.var.endo)) {
    # Add shock.var.endo to shock.var, if it is not in that list already
    shock.var <- union(shock.var, shock.var.endo)
    
    # We need at least one IV per endogenous regressor
    if (length(iv.var) < length(shock.var.endo)) {
      stop("Need at least one IV per endogenous regressor.")
    }
  } else {
    # Remove IV vars if no endogenous regressors are given
    iv.var <- NULL
  }
  
  # Make sure we have at least one shock variable
  if (is.null(shock.var)) {
    stop("Need at least one shock variable.")
  }
  
  # If no control variables given, use the columns in data frame not specified
  # in y.var, shock.var, or iv.var
  if (is.null(contr.var)) {
    contr.var <-
      colnames(data)[!(colnames(data) %in% unique(c(y.var, shock.var, iv.var)))]
  }
  
  # Overplotting only works with a state vector
  if ((overplot == TRUE) & (is.null(state.var))) {
    message("Overplotting only works with a state vector. Setting overplot = FALSE.")
    overplot <- FALSE
  }
  
  # Make sure state vector is binary
  if (!is.null(state.var)) {
    if (!all(data[, state.var] %in% c(0, 1))) {
      stop("State variable must consist of ones (1) and zeros (0) only.")
    }
  }
  
  # Check confidence intervals given
  if (any(c((ci <= 0), (ci >= 1)))) {
    stop("Confidence intervals must lie between 0 and 1.")
  }
  
  # Convert probabilities to normal distribution standard errors
  cse <- qnorm(1 - (1-ci)/2)
  
  # Pad up cse in case second and third confidence band are missing
  cse <- sort(cse, decreasing = TRUE)
  while(length(cse) < 3) {
    cse <- c(cse, cse[length(cse)] * 0.5)
  }
  
  # Make sure column names contain no dots
  if (length(grep("\\.", colnames(data))) > 0) {
    stop("No dots in column names allowed.")
  }
  
  # Change variable name to indicate cumulation
  if (cum.y) {
    colnames(data)[which(colnames(data) == y.var)] <- paste0("cum.", y.var)
    contr.var[which(contr.var == y.var)] <- paste0("cum.", y.var)
    y.var <- paste0("cum.", y.var)
  }
  
  # Check for NA values in data frame
  if(any(is.na(data[, unique(c(y.var, shock.var, iv.var, contr.var))]))) {
    message("NA values present in data set.")
  }
  
  # Print debug information if requested
  if (debug) {
    message(paste("y.var:          ", y.var))
    message(paste("shock.var:      ", paste(shock.var, collapse = ", ")))
    message(paste("shock.var.endo: ", paste(shock.var.endo, collapse = ", ")))
    message(paste("iv.var:         ", paste(iv.var, collapse = ", ")))
    message(paste("contr.var:      ", paste(contr.var, collapse = ", ")))
    message(paste("state.var:      ", paste(state.var, collapse = ", ")))
    message(paste("trend:          ", trend))
    message(paste("n.lag:          ", n.lag))
  }
  
  # Get data frame in format we can use
  dp       <- .getdata(data, n.lag, y.var, shock.var, iv.var, contr.var,
                       state.var, const, trend, trend.sq, trend.cu, trend.qu)
  coefs    <- dp$rhs[0, ]
  se       <- dp$rhs[0, ]
  lmm      <- list()
  lmmse    <- list()
  df       <- list()
  lagorder <- rep(n.lag, hor)
  aicstat  <- rep(NA, hor)
  bicstat  <- rep(NA, hor)
  hqcstat  <- rep(NA, hor)
  
  # Iterate over horizons and calculate model for each horizon
  for (i in 1:hor) {
    # Shift y forward to horizon i
    dd <- .shiftdata(dp, i)
    
    # If cumulating y is requested
    if ((i >= 2) & (cum.y)) {
      dd$lhs <- dd$lhs + dd.last$lhs[1:(nrow(dd.last$lhs)-1),]
    }
    
    # Bind everything into one data frame for convenience
    ds <- cbind(dd$lhs, dd$rhs)
    
    # Select lag length optimally if requested
    if (!is.null(infocrit)) {
      if (is.null(iv.var)) {
        dl <- ds
        ic <- data.frame("AIC" = rep(NA, 6),
                         "BIC" = rep(NA, 6),
                         "HQC" = rep(NA, 6))
        
        for(l in n.lag:1) {
          lml          <- lm(paste0("y.", y.var, " ~ 0 + ."), data = dl)
          ic[l, "AIC"] <- AIC(lml)
          ic[l, "BIC"] <- BIC(lml)
          ic[l, "HQC"] <- .HQC(lml)
          dl           <- dl[, -grep(paste0(".l", l), colnames(dl))]
        }
        
        new.lag     <- which.min(ic[[infocrit]])
        lagorder[i] <- new.lag
        
        ds <- .getdata(data, new.lag, y.var, shock.var, iv.var, contr.var,
                       state.var, const, trend, trend.sq, trend.cu, trend.qu)
        ds <- .shiftdata(ds, i)
        # If cumulating y is requested
        if ((i >= 2) & (cum.y)) {
          ds$lhs <- ds$lhs + dd.last$lhs[1:(nrow(dd.last$lhs)-1),]
        }
        ds <- cbind(ds$lhs, ds$rhs)
      } else {
        message(paste0("Automatic lag order selection not compatible with IV.",
                       " Defaulting to lag order n.lag = ", n.lag))
        infocrit <- NULL
      }
    }
    
    # Estimate model for this forecast horizon
    if (is.null(iv.var)) {
      # Standard lm() regression
      lmp <- lm(paste0("y.", y.var, " ~ 0 + ."), data = ds)
      
      aicstat[i] <- AIC(lmp)
      bicstat[i] <- BIC(lmp)
      hqcstat[i] <- .HQC(lmp)
      
    } else {
      # Instrumental variable ivreg() regression
      rhs1 <- colnames(ds)[(!colnames(ds) %in% c(paste0("y.", y.var),
                                                 iv.var))]
      rhs2 <- colnames(ds)[(!colnames(ds) %in% c(paste0("y.", y.var),
                                                 shock.var.endo))]
      ivreg.formula <- as.formula(paste0("y.", y.var, " ~ 0 + ",
                                         paste(rhs1, collapse = " + "), " | ",
                                         paste(rhs2, collapse = " + ")))
      
      # Note: stargazer output does not work when calling AER::ivreg()
      #lmp <- AER::ivreg(ivreg.formula, data = ds)
      lmp <- ivreg(ivreg.formula, data = ds)
    }  
    
    # lag setting follows Ramey (2017)
    nwse       <- sqrt(diag(sandwich::NeweyWest(lmp, lag = i,
                                                prewhite = FALSE)))
    lmm[[i]]   <- lmp
    lmmse[[i]] <- nwse
    df[[i]]    <- ds
    
    # Store coefficients and HAC standard errors
    coefs     <- dplyr::bind_rows(coefs, lmp$coefficients)
    se        <- dplyr::bind_rows(se, nwse)
    
    # Store current dd as dd.last
    dd.last   <- dd
  }
  
  # Save output to html file if requested
  if (!is.null(output)) {
    capture.output(stargazer::stargazer(
      lmm, type = "html", out = output, se = lmmse,
      dep.var.labels = as.character(y.var), no.space = TRUE, digits = 3,
      digits.extra = 0, omit.stat = c("ser","f"),
      star.cutoffs = sort(1 - ci, decreasing = TRUE),
      add.lines = list(c("Horizon", 0:(hor-1)),
                       c("Lag order", lagorder),
                       c("AIC", round(aicstat, 2)),
                       c("BIC", round(bicstat, 2)),
                       c("HQC", round(hqcstat, 2)))) 
    )
  }
  
  if (!is.null(state.var)) {
    shock.varA <- paste0("A.", shock.var)
    shock.varB <- paste0("B.", shock.var)
    shock.var  <- c(shock.varA, shock.varB)
  }
  
  # Plot
  irfplots <- list()
  
  if (overplot == FALSE) {
    for (r in 1:length(shock.var)) {
      local({
        r <- r
        
        zz <- data.frame(
          idx   = 0:(hor-1),
          lp    = coefs[, shock.var[r]],
          se    = se[, shock.var[r]],
          ci1up = coefs[, shock.var[r]] + (se[, shock.var[r]]*cse[1]),
          ci1lo = coefs[, shock.var[r]] - (se[, shock.var[r]]*cse[1]),
          ci2up = coefs[, shock.var[r]] + (se[, shock.var[r]]*cse[2]),
          ci2lo = coefs[, shock.var[r]] - (se[, shock.var[r]]*cse[2]),
          ci3up = coefs[, shock.var[r]] + (se[, shock.var[r]]*cse[3]),
          ci3lo = coefs[, shock.var[r]] - (se[, shock.var[r]]*cse[3]))
        
        g <- ggplot2::ggplot(zz) +
          geom_hline(aes(yintercept = 0), colour = "#000000", size = 0.5*lw) +
          geom_ribbon(aes(x = idx, ymin = ci1up, ymax = ci1lo), alpha = 0.2) +
          geom_ribbon(aes(x = idx, ymin = ci2up, ymax = ci2lo), alpha = 0.2) +
          geom_ribbon(aes(x = idx, ymin = ci3up, ymax = ci3lo), alpha = 0.2) +
          geom_line(aes(x = idx, y = lp), colour = "#000000", size = lw) +
          theme_bw() +
          labs(title = paste0(shock.var[r], " -> ", y.var),
               x = paste(paste(paste0(ci * 100, "%"), collapse = ", "),
                         "confidence bands")) +
          theme(axis.title.y = element_blank(), axis.title.x = element_blank())
        
        irfplots[[r]] <<- g
      })      
    }
  } else {
    for (r in 1:length(shock.varA)) {
      local({
        r <- r
        
        zz <- data.frame(
          idx    = 0:(hor-1),
          lpA    = coefs[, shock.varA[r]],
          seA    = se[, shock.varA[r]],
          ci1upA = coefs[, shock.varA[r]] + (se[, shock.varA[r]]*cse[1]),
          ci1loA = coefs[, shock.varA[r]] - (se[, shock.varA[r]]*cse[1]),
          ci2upA = coefs[, shock.varA[r]] + (se[, shock.varA[r]]*cse[2]),
          ci2loA = coefs[, shock.varA[r]] - (se[, shock.varA[r]]*cse[2]),
          ci3upA = coefs[, shock.varA[r]] + (se[, shock.varA[r]]*cse[3]),
          ci3loA = coefs[, shock.varA[r]] - (se[, shock.varA[r]]*cse[3]),
          lpB    = coefs[, shock.varB[r]],
          seB    = se[, shock.varB[r]],
          ci1upB = coefs[, shock.varB[r]] + (se[, shock.varB[r]]*cse[1]),
          ci1loB = coefs[, shock.varB[r]] - (se[, shock.varB[r]]*cse[1]),
          ci2upB = coefs[, shock.varB[r]] + (se[, shock.varB[r]]*cse[2]),
          ci2loB = coefs[, shock.varB[r]] - (se[, shock.varB[r]]*cse[2]),
          ci3upB = coefs[, shock.varB[r]] + (se[, shock.varB[r]]*cse[3]),
          ci3loB = coefs[, shock.varB[r]] - (se[, shock.varB[r]]*cse[3]))
        
        g <- ggplot2::ggplot(zz) +
          geom_hline(aes(yintercept = 0), colour = "#000000", size = 0.5*lw) +
          geom_line(aes(x = idx, y = ci3upB), colour = "#0000ff", size = 0.5*lw,
                    linetype = 5) +
          geom_line(aes(x = idx, y = ci3loB), colour = "#0000ff", size = 0.5*lw,
                    linetype = 5) +
          geom_line(aes(x = idx, y = ci2upB), colour = "#0000ff", size = 0.5*lw,
                    linetype = 2) +
          geom_line(aes(x = idx, y = ci2loB), colour = "#0000ff", size = 0.5*lw,
                    linetype = 2) +
          geom_line(aes(x = idx, y = ci1upB), colour = "#0000ff", size = 0.5*lw,
                    linetype = 3) +
          geom_line(aes(x = idx, y = ci1loB), colour = "#0000ff", size = 0.5*lw,
                    linetype = 3) +
          geom_line(aes(x = idx, y = lpB), colour = "#0000ff", size = lw) +
          geom_line(aes(x = idx, y = ci3upA), colour = "#ff0000", size = 0.5*lw,
                    linetype = 5) +
          geom_line(aes(x = idx, y = ci3loA), colour = "#ff0000", size = 0.5*lw,
                    linetype = 5) +
          geom_line(aes(x = idx, y = ci2upA), colour = "#ff0000", size = 0.5*lw,
                    linetype = 2) +
          geom_line(aes(x = idx, y = ci2loA), colour = "#ff0000", size = 0.5*lw,
                    linetype = 2) +
          geom_line(aes(x = idx, y = ci1upA), colour = "#ff0000", size = 0.5*lw,
                    linetype = 3) +
          geom_line(aes(x = idx, y = ci1loA), colour = "#ff0000", size = 0.5*lw,
                    linetype = 3) +
          geom_line(aes(x = idx, y = lpA), colour = "#ff0000", size = lw) +
          theme_bw() +
          labs(title = paste0(substring(shock.var[r], 3), " -> ", y.var),
               x = paste(paste(paste0(ci * 100, "%"), collapse = ", "),
                         "confidence bands")) +
          theme(axis.title.y = element_blank(), axis.title.x = element_blank())
        
        g
        
        irfplots[[r]] <<- g
      })
    }
  }
  
  return(list(coefs = coefs, se = se, lm_models = lmm, df = df, g = irfplots))
}

#' Lag matrix (internal)
#'
#' Internal function to obtain a 'lag matrix' containing contemporaneous and
#' lagged terms of the variables in a time-series matrix.

.lagmatrix <- function(x, k, only.lags = FALSE) {
  xn   <- embed(x, k+1)
  
  # Add identifiers to lags
  cn <- colnames(x)
  for (i in 1:k) {
    tn <- paste(colnames(x), ".l", i, sep = "")
    cn <- c(cn, tn)
  }
  colnames(xn) <- cn
  
  # Cut out contemporaneous terms if requested
  if (only.lags) {
    xn <- xn[, -(1:ncol(x))]
  }
  
  return(xn)
}

#' Prepare projection data (internal)
#'
#' Internal function to prepare data frame with lags and deterministic
#' components.

.getdata <- function(data, n.lag, y.var, shock.var, iv.var, contr.var,
                     state.var, const, trend, trend.sq, trend.cu, trend.qu) {
  
  data.orig   <- data
  data.lagged <- .lagmatrix(as.matrix(data), n.lag)
  
  # Left-hand side
  lhs           <- data.frame(data.lagged[, y.var])  # Response variable
  colnames(lhs) <- paste0("y.", y.var)
  lhs           <- lhs
  
  # Right-hand side
  rhs <- lhs[, 0]
  
  # Add constant if requested
  if (const) {
    rhs <- cbind(rhs, constant = rep(1, nrow(data.lagged)))
  }
  
  # Add linear trend if requested
  if (trend) {
    rhs <- cbind(rhs, tln = seq(1, nrow(data.lagged), by = 1))
  }
  
  # Add squared trend if requested
  if (trend.sq) {
    rhs <- cbind(rhs, tsq = seq(1, nrow(data.lagged), by = 1)^2)
  }
  
  # Add cubic trend if requested
  if (trend.cu) {
    rhs <- cbind(rhs, tcu = seq(1, nrow(data.lagged), by = 1)^3)
  }
  
  # Add quartic trend if requested
  if (trend.qu) {
    rhs <- cbind(rhs, tqu = seq(1, nrow(data.lagged), by = 1)^4)
  }
  
  # Shock variables
  shock.var.data <- data.frame(data.lagged[, c(shock.var, iv.var)])
  colnames(shock.var.data) <- c(shock.var, iv.var)
  
  # Control variables, only keep lagged
  if (length(contr.var) > 0) {
    contr.var.data <- .lagmatrix(as.matrix(data[, contr.var]), n.lag)
    contr.var.data <- contr.var.data[, -which(colnames(contr.var.data) %in%
                                                contr.var)]
  }
  
  rhs <- cbind(rhs, shock.var.data, contr.var.data)
  
  # Split right-hand side by state if requested
  if (!is.null(state.var)) {
    rhs0           <- (!as.logical(data[-(1:n.lag), state.var])) * rhs
    colnames(rhs0) <- paste0("A.", colnames(rhs0))
    
    rhs1           <- as.logical(data[-(1:n.lag), state.var]) * rhs
    colnames(rhs1) <- paste0("B.", colnames(rhs1))
    
    rhs       <- cbind(rhs0, rhs1)
  }
  
  return(list(lhs = lhs, rhs = rhs, data.orig = data.orig))
}

#' Shift data set forward (internal)
#'
#' Internal function to shift y (left-hand side) forward to i-th impulse
#' response horizon and adjust regressors (right-hand side) for lost
#' observations.

.shiftdata <- function(dp, i, cum.y) {
  lhs           <- dp$lhs[i:nrow(dp$lhs),]
  lhs           <- data.frame(lhs)
  colnames(lhs) <- colnames(dp$lhs)
  
  rhs           <- dp$rhs[1:(nrow(dp$rhs)-i+1),]
  
  return(list(lhs = lhs, rhs = rhs))
}

#' Get Hannan-Quinn information criterion (internal)
#' 
#' Internal function to obtain Hannan-Quinn information criterion for an R model
#' object. Formula adapted from BIC() method "BIC.default".

.HQC <- function (object) {
  lls <- logLik(object)
  nos <- attr(lls, "nobs")
  hqc <- -2 * as.numeric(lls) + 2 * log(log(nos)) * attr(lls, "df")
  
  return(hqc)
}
