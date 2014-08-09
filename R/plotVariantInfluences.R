plotVariantInfluences <-
function(data.obj, p.or.q = 0.05, all.markers = FALSE, standardize = TRUE, not.tested.col = "lightgray", show.marker.labels = FALSE, show.chr = TRUE, label.chr = TRUE, scale.effects = c("log10", "sqrt", "none"), pheno.width = 11, phenotype.labels = NULL){
	
	if(length(grep("n", scale.effects)) > 0){
		scale.effects <- "none"
		}
	if(length(scale.effects) == 1){
		if(scale.effects != "log10" & scale.effects != "sqrt" & scale.effects != "none"){
			stop("scale.effects must be 'log10', 'sqrt' or 'none.'")
			}
		}
	
	var.influences <- data.obj$var.to.var.p.val
		
	pheno.inf <- data.obj$max.var.to.pheno.influence
	if(is.null(phenotype.labels)){
		pheno.names <- names(data.obj$max.var.to.pheno.influence)
		}else{
		pheno.names <- phenotype.labels
		if(length(pheno.names) != length(names(data.obj$max.var.to.pheno.influence))){
			stop("I am detecting the wrong number of phenotype labels for the phenotypes present.")
			}
		}
	num.pheno <- length(pheno.names)
	
	if(not.tested.col == TRUE){
		not.tested.col = "lightgray"
		}
	
	if(is.null(var.influences)){
		stop("calc.p() must be run to calculate variant-to-variant influences.")
		}

	if(is.null(pheno.inf)){
		stop("direct.influence() must be run to calculate variant-to-trait influences.")
		}

	
	if(all.markers){
		unique.markers <- colnames(data.obj$geno)
		}else{
		unique.markers <- unique(c(as.vector(var.influences[,"Source"]), as.vector(var.influences[,"Target"]), pheno.inf[[1]][,"marker"]))
		}
		
	unique.marker.locale <- which(colnames(data.obj$geno) %in% unique.markers)
	#get coordinates of the chromosome boundaries
	if(show.chr){
		chromosomes <- data.obj$chromosome[sort(unique.marker.locale)]
		u_chr <- unique(chromosomes[which(!is.na(chromosomes))])
		chr.boundaries <- apply(matrix(u_chr, ncol = 1), 1, function(x) max(which(chromosomes == x))) + 0.5
		chr.boundaries <- c(0, chr.boundaries)
		if(label.chr){
			chr.names <- unique(chromosomes)
			}else{
			chr.names <- NULL
			}
		}else{
		chr.boundaries <- NULL
		chr.names <- NULL
		}
	
		
	marker.locale <- match(unique.markers, colnames(data.obj$geno))
	sorted.markers <- unique.markers[order(marker.locale)]
	
	var.influence.mat <- matrix(NA, nrow = length(unique.markers), ncol = length(unique.markers))
	var.pval.mat <- matrix(NA, nrow = length(unique.markers), ncol = length(unique.markers))
	colnames(var.influence.mat) <- rownames(var.influence.mat) <- colnames(var.pval.mat) <- rownames(var.pval.mat) <- sorted.markers

	pheno.influence.mat <- matrix(NA, nrow = length(unique.markers), ncol = num.pheno)
	pheno.pval.mat <- matrix(NA, nrow = length(unique.markers), ncol = num.pheno)
	colnames(pheno.influence.mat) <- colnames(pheno.pval.mat) <- pheno.names
	rownames(pheno.influence.mat) <- rownames(pheno.pval.mat) <- sorted.markers
	
	
	#fill the variant-to-variant matrix with test statistics with sources in rows and targets in columns
	var.sig.col <- as.vector(na.omit(match(c("qval", "lfdr", "p.adjusted"), colnames(var.influences))))
	for(i in 1:length(var.influences[,1])){
		if(standardize){
			var.influence.mat[as.character(var.influences[i,"Source"]), as.character(var.influences[i,"Target"])] <- as.numeric(as.vector(var.influences[i,"Effect"]))/as.numeric(as.vector(var.influences[i,"SE"]))
			}else{
			var.influence.mat[as.character(var.influences[i,"Source"]), as.character(var.influences[i,"Target"])] <- as.numeric(as.vector(var.influences[i,"Effect"]))		
				}
		var.pval.mat[as.character(var.influences[i,"Source"]), as.character(var.influences[i,"Target"])] <- as.numeric(as.vector(var.influences[i,var.sig.col]))
		}
	
	#fill the variant-to-phenotype matrix with test statistics 
	#(still with sources in rows and targets in columns)
	#use phenotypes or eigentraits based on user input
	pheno.sig.col <- as.vector(na.omit(match(c("qval", "lfdr", "p.adjusted"), colnames(pheno.inf[[1]]))))
	for(i in 1:length(unique.markers)){
			for(j in 1:length(pheno.names)){
				marker.locale <- which(pheno.inf[[j]][,"marker"] == unique.markers[i])
				if(length(marker.locale) > 0){
					if(standardize){	
						pheno.influence.mat[as.character(unique.markers[i]), pheno.names[j]] <- pheno.inf[[j]][marker.locale, "t.stat"]
						}else{
						pheno.influence.mat[as.character(unique.markers[i]), pheno.names[j]] <- pheno.inf[[j]][marker.locale, "coef"]
						}
				pheno.pval.mat[as.character(unique.markers[i]), pheno.names[j]] <- pheno.inf[[j]][marker.locale, pheno.sig.col]
				}else{
					pheno.influence.mat[as.character(unique.markers[i]), pheno.names[j]] <- NA
					pheno.pval.mat[unique.markers[i], pheno.names[j]] <- NA
					}
			}
		}
	
	
	#expand the phenotype influence matrix to give it more visual weight in the plot
	expanded.pheno.mat <- matrix(NA, nrow = dim(pheno.influence.mat)[1], ncol = dim(pheno.influence.mat)[2]*pheno.width)
	expanded.pheno.pval.mat <- matrix(NA, nrow = dim(pheno.influence.mat)[1], ncol = dim(pheno.influence.mat)[2]*pheno.width)
	expanded.pheno.names <- rep("", dim(pheno.influence.mat)[2]*pheno.width)
	
	start.col = 1
	for(i in 1:dim(pheno.influence.mat)[2]){
		expanded.pheno.mat[,start.col:(start.col+pheno.width-1)] <- pheno.influence.mat[,i]
		expanded.pheno.pval.mat[,start.col:(start.col+pheno.width-1)] <- pheno.pval.mat[,i]
		expanded.pheno.names[median(start.col:(start.col+pheno.width-1))] <- colnames(pheno.pval.mat)[i]
		start.col = start.col + pheno.width
		}

		
	full.inf.mat <- cbind(var.influence.mat, expanded.pheno.mat)
	full.pval.mat <- cbind(var.pval.mat, expanded.pheno.pval.mat)
	
	full.inf.mat <- apply(full.inf.mat, 2, as.numeric)
	full.pval.mat <- apply(full.pval.mat, 2, as.numeric)
	
	marker.locale <- which(colnames(data.obj$geno) %in% sorted.markers)
	rownames(full.inf.mat) <- data.obj$marker.names[marker.locale]
	colnames(full.inf.mat) <- c(data.obj$marker.names[marker.locale], expanded.pheno.names)
	
	#get the coordinates for all pairs not tested
	# not.tested.locale <- which(is.na(rotate.mat(full.inf.mat)), arr.ind = TRUE)
	not.tested.locale <- which(is.na(rotate.mat(full.inf.mat)), arr.ind = TRUE)
	
	if(not.tested.col == FALSE || is.na(not.tested.col)){
		not.tested.locale <- NULL
		}
	
	#take out any values that aren't significant by the cutoff
	full.inf.mat[which(full.pval.mat > p.or.q)] <- NA
	main <- "Variant Influences"
	
	if(scale.effects == "log10"){
		neg.locale <- which(full.inf.mat < 0)
		scaled.effects <- log10(abs(full.inf.mat))
		scaled.effects[neg.locale] <- scaled.effects[neg.locale]*-1	
		full.inf.mat <- scaled.effects
		main <- "log10 Variant Influences"
		}
	if(scale.effects == "sqrt"){
		neg.locale <- which(full.inf.mat < 0)
		scaled.effects <- sqrt(abs(full.inf.mat))
		scaled.effects[neg.locale] <- scaled.effects[neg.locale]*-1	
		full.inf.mat <- scaled.effects
		main <- "Square Root of Variant Influences"
		}
	
		
	if(length(which(na.omit(as.vector(full.inf.mat)) != 0)) == 0){
			plot.new()
			plot.window(xlim = c(0,1), ylim = c(0,1))
			text(0.5, 0.5, "No Significant Interactions")
		}else{
		data.obj$full.adjacency <- full.inf.mat
		myImagePlot(full.inf.mat, min.x = (max(abs(full.inf.mat), na.rm = TRUE)*-1), max.x = max(abs(full.inf.mat), na.rm = TRUE), main = main, xlab = "Target", ylab = "Source", mark.coords = not.tested.locale, mark.col = not.tested.col, show.labels = show.marker.labels, chromosome.coordinates = chr.boundaries, chr.names = chr.names, show.pheno.labels = TRUE)
		#add phenotype names
		if(!is.null(not.tested.locale)){
			legend("topright", legend = "not testable", col = not.tested.col, pch = 16)
			}
		}
	
	invisible(full.inf.mat)

	
	}
