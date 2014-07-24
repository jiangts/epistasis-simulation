#This function adjusts the m12/m21 null so we can
#use it as a symmetric distribution


symmetrize.null <- function(null.dist, emp.dist, plot.dist = FALSE, plot.title = NULL){
	
	if(is.null(plot.title)){
		plot.title <- "m12.m21.Distribution.Normalization.pdf"
		}
	
	if(plot.dist){
		layout.mat <- matrix(c(1,2,3,5,4,5), ncol = 2, byrow = TRUE)	
		pdf(plot.title, width = 8, height = 8)
		# quartz(width = 8, height = 8)
		layout(layout.mat, heights = c(2,1))
		# layout.show(5)
		
		null.dens <- density(null.dist)
		emp.dens <- density(emp.dist)
		xmin <- min(c(null.dens$x, emp.dens$x)); xmax <- max(c(null.dens$x, emp.dens$x))
		ymin <- min(c(null.dens$y, emp.dens$y)); ymax <- max(c(null.dens$y, emp.dens$y))	
		plot(density(null.dist), col = "blue", type = "l", lwd = 3, main = "Distribution of Standardized m12/m21", xlim = c(xmin, xmax), ylim = c(ymin, ymax), cex.main = 2, cex.axis = 2)
		points(emp.dens, col = "red", type = "l", lwd = 3)
		legend("topleft", legend = c("Null", "Obs."), col = c("blue", "red"), lty = 1, cex = 2, lwd = 3)
		}

	
	min.val <- min(null.dist, emp.dist)
	max.val <- max(null.dist, emp.dist)
	mid.val <- (max.val + min.val)/2
	min.val <- (1.1*(min.val - mid.val)) + mid.val
	max.val <- (1.1*(max.val - mid.val)) + mid.val
		
	#generate a smooth cumulative 
	#distribution function
	#for the null distribution
	null.smooth.fun <- smooth.cdf(null.dist, min.val = min.val, max.val = max.val)

	if(plot.dist){
		plot(null.smooth.fun, xlim = c(min.val, max.val), col = "blue", lwd = 3, main = "Null CDF", cex.axis = 2, cex.main = 2)
		}
	
	#apply it to both the null and empirical distributions
	smooth.null <- null.smooth.fun(null.dist)
	norm.null <- qnorm(smooth.null)
	
	smooth.emp <- null.smooth.fun(emp.dist)
	norm.emp <- qnorm(smooth.emp)
	
	if(plot.dist){
		par(mar = c(4,4,3,2))
		hist(smooth.null, main = "Null Fit to CDF", freq = FALSE)
		hist(smooth.emp, main = "Empirical Fit to CDF", freq = FALSE)	
	
		par(mar = c(5,4,4,2))
		xmin <- min(c(density(norm.null)$x, density(norm.emp)$x)); xmax <- max(c(density(norm.null)$x, density(norm.emp)$x))
		ymin <- min(c(density(norm.null)$y, density(norm.emp)$y)); ymax <- max(c(density(norm.null)$y, density(norm.emp)$y))
		plot(density(norm.null), col = "blue", xlim = c(xmin, xmax), lwd = 3, main = "Normalized Distributions", ylim = c(ymin, ymax), cex.axis = 2, cex.main = 2)
		points(density(norm.emp), col = "red", type = "l", lwd = 3)
		dev.off()
		}
	
	smoothed.dists <- list(norm.null, norm.emp)
	names(smoothed.dists) <- c("normalized.null", "normalized.empirical")
	return(smoothed.dists)	
	
}