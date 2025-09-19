require(bezier)

controlPointsBySlope <- function(P0, P3, slope_start, slope_end, alpha, beta) {
	# (Helper function)
	# Find control points P1 and P2
	#	when given P0, a start slope vector, and segment (P0,P1) magnitude alpha
	#	and  given P3, an end slope vector,  and segment (P2,P3) magnitude beta
	#
	# Arguments:
	# 	P0 = starting control point
	# 	P3 = ending control point for cubic Bezier
	# 	slope_start = slope vector from the starting point
	# 	slope_end = slope vector to the ending point
	# 	alpha = magnitude of segment from P0 to P1
	# 	beta = magnitude of segment from P2 to P3
	#
	# Returns:
	#   A set of control points that define a cubic Bezier curve

	P1 <- P0 + alpha * slope_start	# May have to be (in case of 1-dim vectors being mis-cast) :  P1 <- as.numeric(P0) + alpha * as.numeric(slope_start)
	P2 <- P3 - beta * slope_end		# May have to be (in case of 1-dim vectors being mis-cast) :  P2 <- as.numeric(P3) - beta * as.numeric(slope_end)
	return(rbind(P0, P1, P2, P3))   # row-bind into matrix and return
}

tForX <- function(X_goal, P, tol=1e-8, max_iter=20) {
	# (Helper function)
	# Estimates the set of t parameters (T) that will give
	# 	the desired X values via Newton's method (Newton-Raphson)
	# Assumes that the Bezier curve is monotonic on X
	# (Which it will be if P1x and P2x are within [P0x, P3x])
	#
	# Arguments:
	# 	X_goal = vector of x points to find (t)s for
	# 	P = matrix of control points	(for cubic Bezier, 4 x (n-dim) )
	# 		Note that, because of matrix notation, P0,P1,P2,P3 are P1,P2,P3,P4 below
	# 	tol = tolerance (how close to the root should we get?)
	# 	max_iter = max iterations (when should we give up on root finding?)
	#
	# Returns:
	#   the approximate t values (T) that produce X_goal
	#
	
	
	# Get X co-ordinates of control points (first column = x)
	Px <- P[,1]

	# T0 -- initial guesses by linear interpolation (proprotional from endpoints)
	T <- (X_goal - Px[1]) / (Px[4] - Px[1])

	for (i in 1:max_iter) {
	
		# x(T)	-- what are the x values for our T values?
		X_T <- 	(1-T)^3*Px[1] +
				3*(1-T)^2*T*Px[2] +
				3*(1-T)*T^2*Px[3] +
				T^3*Px[4]

		# dx/dt	-- derivative of the above, per Bezier (1986)
		dXdT <- 3*(1-T)^2*(Px[2]-Px[1]) +
				6*(1-T)*T*(Px[3]-Px[2]) +
				3*T^2*(Px[4]-Px[3])

		# The function we are trying to solve for is:
		# f(t) = x(t) - X_goal = 0
		f_T <- X_T - X_goal

		# Core Newton's method step
		T_next <- T - f_T/dXdT
		
		# if any t values in T_next are out of [0,1] range
		# clamp to [0,1] using parallel min / max
		T_next <- pmin(pmax(T_next,0),1)
		
		# Return vector of (t)s if all within tolerance
		if (max(abs(f_T)) < tol) return(T_next)
		
		# otherwise, iterate
		T <- T_next
	}

	warning(sprintf("Newton's Method did not converge for X=%s and control points: %s",
                toString(X_goal), paste(P, collapse=", ")))
	return(T)
}

bezRSS <- function(par, m, P0, P3, slope_start, slope_end, max_iter) {
	# (Helper function)
	# Calculates the Residual Sum of Squares for a
	#	restricted cubic Bezier from a set of parameters (from optim()),
	#	a data matrix, start and end points,
	#	and slope vectors from (P0,P1) and (P2,P3)
	#
	# Arguments:
	#	par = parameters to be fit, passed by optim(), contains:
	#		alpha = magnitude of segment from P0 to P1
	# 		beta = magnitude of segment from P2 to P3
	# 	P0 = starting control point
	# 	P3 = ending control point for cubic Bezier
	# 	slope_start = slope vector from the starting point
	# 	slope_end = slope vector to the ending point
	#
	# Returns:
	#   a single numeric value
	#
	
	# Get alpha and beta values from parameters (from optim())
	alpha <- par[1]
	beta <- par[2]
	
	# Get estimated control points to pass to Bezier
	P <- controlPointsBySlope(P0, P3, slope_start, slope_end, alpha, beta)
	
	# Estimate Bezier T values from X values (presumed first column of matrix)
	T <- tForX(X_goal = m[, 1], P = P, tol=1e-8, max_iter=max_iter)
	
	# Return (x,y,...) values from T values along Bezier
	predicted <- bezier(T, P)
	
	# Calculate RSS over the rightmost matrix column (presumed dependent variable)
	d <- ncol(m)
	sum((m[, d] - predicted[, d])^2)
}


viewParamMap <- function(a0, b0, m, P0, P3, slope_start, slope_end){
	# Heatmap of alpha,beta values around starting values (a0,b0)
	# Takes other parameters common to bezRSS
	
	alphas <- seq(a0*0.1, a0*10, length=30)
	betas  <- seq(b0*0.1, b0*10, length=30)
	rss_fit <- matrix(NA, nrow=length(alphas), ncol=length(betas))
	
	for(i in seq_along(alphas)){
		for(j in seq_along(betas)){
			r <- bezRSS(c(alphas[i], betas[j]), m, P0, P3, slope_start, slope_end)
			rss_fit[i,j] <- sum(r^2)
		}
	}
	
	image(alphas, betas, log(rss_fit), xlab="alpha", ylab="beta", main="log(Resid. Sum of Squares)")
	contour(alphas, betas, log(rss_fit), add=TRUE)
}


restrictedCubicBezierFit <- function(m, ends_from_data = FALSE,
	P0 = NULL, P3 = NULL, slope_start, slope_end,
	na_fill = FALSE, max_iter = 20) {
	#
	# Modified from "bezierCurveFit.R" in packge "bezier" by Aaron Olsen
	#	(ALL CAPS comments are his / his code)
	# 
	# Fits a restricted cubic Bezier curve
	#	(where P0 (starting) and P3 (ending) points
	#	 and slopes from (P0,P1) and (P2,P3) are given)
	# Assumes a monotonic Bezier
	#	(which it will be if P1 and P2 are within [P0,P3])
	#
	# Arguments:
	#	m = matrix of datapoints (can be in multiple dimensions -- each column is a dimension)
	# 		Will use first and last rows as control points P0 and P3 (start and end points)
	#	ends_from_data = (Boolean)
	#		(if true) Get starting and ending points from data matrix
	#	P0, P3 = (if ends_from_data not specified)
	#		Starting and ending points, each as n-dim (x,y,...) vectors
	# 	slope_start = Slope vector from starting point
	#	slope_end   = Slope vector from ending point
	#	na_fill = (Boolean) interpolate missing values in data matrix
	#	max_iter = maximum number of iterations to use when fitting
	#
	# Returns:
	#   A fit of a cubic Bezier curve to the data, including:
	#	fit$par = alpha,beta parameters (magnitude of (P0,P1) and (P2,P3) segments)
	#	fit$value = Residual Sum of Squares
	#	fit$Points = set of control points for best fit cubic Bezier
	#	and information on convergence
	#


	####                              ####
	#### Pre-process passed arguments ####
	####                              ####
	
	if(na_fill){

		# FIND NA VALUES
		is_na <- is.na(m[, 1])
	
		# IF NA VALUES, FILL
		i <- 1
		while(i < nrow(m)){
	
			# SKIP IF NOT NA
			if(!is_na[i]){i <- i + 1;next}
	
			# FIND NEXT NON-NA VALUE
			next_non_na <- i + which(!is_na[(i+1):nrow(m)])[1]
	
			# FIND POINTS ON LINE BETWEEN NON-NA VALUES
			m_fill <- matrix(m[next_non_na, ] - m[i-1, ], nrow=next_non_na-i+2, ncol=ncol(m), byrow=TRUE)*seq(0, 1, length=next_non_na-i+2) + matrix(m[i-1, ], nrow=next_non_na-i+2, ncol=ncol(m), byrow=TRUE)

			# ENTER VALUES INTO M MATRIX
			m[(i-1):next_non_na, ] <- m_fill
	
			# WHERE TO START SEARCH FOR NEXT NA
			i <- next_non_na
		}
	}

	# CONVERT INPUT PARAMETER POINTS TO MATRIX
	if(!is.matrix(m)) m <- as.matrix(m)

	# Provide warning message if (P0 or P3) and ends_from_data are provided
	if(ends_from_data && (!any(is.null(P0)) || !any(is.null(P3)))){
		warning("Warning: fix.end.start specified; taking P0, P3 from data. Ignoring user specified P0, P3.")
	}
	
	# Set P0 and P3 from matrix of datapoints
	if(ends_from_data){
		P0 <- as.vector(m[1, ])			# first line
		P3 <- as.vector(m[nrow(m), ])	# last line
	}
	
	# Rescale slope vectors on [-1,1]
	slope_start <- slope_start / max(abs(slope_start))
	slope_end <- slope_end / max(abs(slope_end))



	####            ####
	#### Procedures ####
	####            ####
	
	# Fit magnitude of P0,P1 (alpha) and P2,P3 (beta) control point segments
	# using optim() and the Residual Sum of Squares for points off the spline
	# (optim() defaults to minimization)
	fit <- optim(
		par = c(alpha=0.1, beta=0.1),
		fn = bezRSS,
		m = m, P0 = P0, P3 = P3,
		slope_start = slope_start, slope_end = slope_end,
		max_iter = max_iter,
		method = "L-BFGS-B",
		lower = c(0,0)   # optional: constrain >0
	)
	
	# Save Control Point Estimates to fit
	parameter_estimate <- fit$par
	fit$Points <- rbind(
		P0,
		P0 + parameter_estimate['alpha'] * slope_start,	# P1
		P3 - parameter_estimate['beta'] * slope_end,	# P2
		P3
	)
	
	## Save additional fit metrics
	
	# Get points on fit Bezier
	d <- ncol(m)	# rightmost matrix column (presumed dependent variable) (= # of dimensions)
	n <- nrow(m)
	T <- tForX(X_goal = m[, 1], P = fit$Points, tol=1e-8, max_iter=max_iter)	# T values for X values from data
	predicted <- bezier(T, fit$Points)
	
	# Total sum of squares (for calculations)
	tss <- sum((m[, d] - mean(m[, d]))^2)
	
	# Save fit metrics to "fit"
	fit$residuals 	<- m[, d] - predicted[, d]
	fit$rss 		<- sum(fit$residuals^2)
	fit$rmse 		<- sqrt(mean(fit$residuals^2))
	fit$r2			<- 1 - fit$rss / tss
	fit$r2_adj		<- 1 - ((fit$rss/(n - 2))/(tss/(n - 1)))	# 2 free paramters (alpha,beta)
	
	
	# Final check for monotonicity
	monotonic <- 
		apply(fit$Points, 2, function(col) {
			all(col[1] <= col[2:3] & col[2:3] <= col[4]) || 
			all(col[1] >= col[2:3] & col[2:3] >= col[4])
		}) |> all()
	if(!monotonic){
		warning("Warning: Central control points exceed bounds of start/end points. May not be monotonic, convergence may not be accurate.")
	}
	
	# Return fit
	class(fit) <- "Restricted Cubic Bezier Fit"
	return(fit)
}