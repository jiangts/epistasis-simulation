#We need to make sure all markers end up in the blocks!
#This is some code to start playing around with a 
#combination 2D peak-finding plus main-effect
#method for determining groups of markers that
#should be condensed. 

linkage.blocks.cormat <- function(data.obj, p.or.q = 0.05, effect.drop = 1.5, collapse.linked.markers = TRUE, threshold.power = 2, plot.results = TRUE, verbose = FALSE){

	# plot.results = TRUE
	# require(igraph)

	var.influences <- data.obj$var.to.var.p.val
	main.effects <- data.obj$max.var.to.pheno.influence
	data.obj$network.p.or.q <- p.or.q
	
	if(length(var.influences) == 0){
		stop("calc.p() must be run to calculate variant-to-variant influences.")
		}

	
	all.chr <- data.obj$chromosome
	u_chr <- unique(all.chr)
	source.names <- data.obj$marker.names[match(var.influences[,1], colnames(data.obj$geno))]
	target.names <- data.obj$marker.names[match(var.influences[,2], colnames(data.obj$geno))]
	
	ordered.effects <- lapply(main.effects, function(x) x[order(x[,1]),])
	main.marker.names <- lapply(ordered.effects, function(x) data.obj$marker.names[match(x[,1], colnames(data.obj$geno))])
	
	link.blocks <- vector(mode = "list", length = 1)
	var.sig.col <- as.vector(na.omit(match(c("qval", "lfdr", "p.adjusted"), colnames(var.influences))))
	pheno.sig.col <- as.vector(na.omit(match(c("qval", "lfdr", "p.adjusted"), colnames(main.effects[[1]]))))

	markers.used <- colnames(data.obj$geno.for.pairscan)
	marker.locale <- match(markers.used, colnames(data.obj$geno))
	all.marker.chr <- data.obj$chromosome[marker.locale]
	all.marker.names <- data.obj$marker.names[marker.locale]

	#if we are not collapsing the markers, just assign each marker to a block
	if(!collapse.linked.markers){		
		
		for(i in 1:length(markers.used)){
			link.blocks[[i]] <- markers.used[i]
			names(link.blocks)[i] <- all.marker.names[i]
			}
		
		}else{
	
		if(plot.results){
			pdf.name <- paste("Link.Block.Matrices.p.", p.or.q, ".soft.thresh.", threshold.power, ".effect.drop.", effect.drop, ".pdf", sep = "")		
			pdf(pdf.name, width = 40, height = 8)
			}
		num.blocks <- 1
		
		for(ch in 1:length(u_chr)){
			if(verbose){report.progress(ch, length(u_chr))}

			if(u_chr[ch] == 0){
				marker.locale <- which(all.marker.chr == 0)
				marker.num <- colnames(data.obj$geno.for.pairscan)[marker.locale]
				marker.name <- all.marker.names[marker.locale]
				for(m in 1:length(marker.locale)){
					link.blocks[[num.blocks]] <- marker.num[m]
					names(link.blocks)[num.blocks] <- marker.name[m]
					num.blocks <- num.blocks + 1
					}
				}else{

				chr.blocks <- 1
				
				chr.marker.locale <- which(all.marker.chr == u_chr[ch])
				chr.markers <- all.marker.names[chr.marker.locale]
				
				if(length(chr.markers) > 0){	
					source.locale <- which(source.names %in% chr.markers)
					target.locale <- which(target.names %in% chr.markers)
					
					source.mat <- matrix(NA, ncol = length(markers.used), nrow = length(chr.markers))
					target.mat <- matrix(NA, ncol = length(markers.used), nrow = length(chr.markers))
					max.mat <- matrix(NA, ncol = length(markers.used), nrow = length(chr.markers))
					source.pval <- matrix(NA, ncol = length(markers.used), nrow = length(chr.markers)) 
					target.pval <- matrix(NA, ncol = length(markers.used), nrow = length(chr.markers)) 
					min.pval <- matrix(NA, ncol = length(markers.used), nrow = length(chr.markers)) 
					rownames(source.mat) <- rownames(target.mat) <- rownames(max.mat) <- rownames(source.pval)  <- rownames(target.pval) <- rownames(min.pval) <- chr.markers
					colnames(source.mat) <- colnames(target.mat) <- colnames(max.mat) <- colnames(source.pval) <- colnames(target.pval) <- colnames(min.pval) <- all.marker.names
					
					for(s in source.locale){
						source.mat[as.character(source.names[s]), as.character(target.names[s])] <- var.influences[s,"Effect"]/var.influences[s,"SE"]
						source.pval[as.character(source.names[s]), as.character(target.names[s])] <- var.influences[s, var.sig.col]
						}
					
					for(s in target.locale){
						target.mat[as.character(target.names[s]), as.character(source.names[s])] <- var.influences[s,"Effect"]/var.influences[s,"SE"]
						target.pval[as.character(target.names[s]), as.character(source.names[s])] <- var.influences[s, var.sig.col]
						}
					
					
					#make a matrix containing the maximum interaction effect 
					#and one containing the minimum p value for each interaction
					for(i in 1:dim(source.mat)[1]){
						for(j in 1:dim(source.mat)[2]){
							max.mat[i,j] <- max(abs(source.mat[i,j]), abs(target.mat[i,j]))
							min.pval[i,j] <- min(abs(source.pval[i,j]), abs(target.pval[i,j]))
							}
						}
					
					
					overall.max <- max(abs(var.influences[,"Effect"]/var.influences[,"SE"]))
				
					#also get the main effects for each marker
					main.effect.mat <- matrix(NA, ncol = dim(data.obj$pheno)[2], nrow = length(chr.markers))
					main.pval <- matrix(NA, ncol = dim(data.obj$pheno)[2], nrow = length(chr.markers))
					colnames(main.effect.mat) <- colnames(main.pval) <- names(main.effects)
					rownames(main.effect.mat) <- rownames(main.pval) <- chr.markers
					for(ph in 1:length(main.effects)){
						main.effect.mat[,ph] <- ordered.effects[[ph]][which(main.marker.names[[ph]] %in% chr.markers),"|t.stat|"]
						main.pval[,ph] <- ordered.effects[[ph]][which(main.marker.names[[ph]] %in% chr.markers), pheno.sig.col]
						}
				
				
					#combine the interaction and main effect matrices
					combined.mat <- cbind(main.effect.mat, max.mat)
					combined.pval <- cbind(main.pval, min.pval)
					
					#do a hard threshold just to see where the peaks are
					non.sig.vals <- which(min.pval > p.or.q)
					thresh.int <- max.mat
					thresh.int[non.sig.vals] <- NA
					
					non.sig.vals <- which(main.pval > p.or.q)
					thresh.main <- main.effect.mat
					thresh.main[non.sig.vals] <- NA
	
					#threshold the matrix using the effect drop peak finder
					test.interactions <- get.peaks.2d(m.mat = max.mat, p.mat = min.pval, p.or.q, effect.drop) 
					# image(thresh.int);quartz();image(test.interactions)
					
					test.main <- get.peaks.1d(m.mat = main.effect.mat, p.mat = main.pval, p.or.q, effect.drop)
					# image(thresh.main); quartz(); image(test.main)
					
					combined.mat <- cbind(main.effect.mat, max.mat)
					combined.thresh.mat <- cbind(test.main, test.interactions)
					combined.sig.mat <- cbind(thresh.main, thresh.int)
					# image(combined.mat); quartz(); image(combined.thresh.mat); quartz(); image(combined.sig.mat)
					
					mult.mat <- combined.thresh.mat; mult.mat[which(is.na(mult.mat))] <- 0
					
					zero.row <- which(apply(mult.mat, 1, sum) == 0)
					if(length(zero.row) == dim(mult.mat)[1]){
						next()
						}
					if(length(zero.row) > 0){
						mult.mat <- mult.mat[-zero.row,,drop=FALSE]
						}
			
					mult.mat <- t(apply(mult.mat, 1, function(x) x/sqrt(sum(x^2))))
					cor.mat <- mult.mat %*% t(mult.mat)
					
					thresh.cor.mat <- cor.mat^threshold.power
					
					# image(cor.mat)
					# image(cor.mat, zlim = c(0,1))
					# image(thresh.cor.mat, zlim = c(0,1))
					sim.mat <- thresh.cor.mat; sim.mat[which(sim.mat < 0)] <- 0
					diag(sim.mat) <- 0
					net <- graph.adjacency(sim.mat, mode = "undirected", weighted = TRUE)
		
					comm <- fastgreedy.community(net)$membership
	
					comm.num <- 1
					adj.comm <- consec.pairs(comm)
					cm.changes <- which(!apply(adj.comm, 1, function(x) x[1] == x[2])) #find everywhere the community number changes
					if(length(cm.changes) == 0){ #if there are no changes, put the whole chromosome into the block
						marker.names <- chr.markers
						marker.num <- colnames(data.obj$geno)[match(marker.names, data.obj$marker.names)]
						link.blocks[[num.blocks]] <- marker.num
						names(link.blocks)[num.blocks] <- paste("Chr", u_chr[ch], "_", chr.blocks, sep = "")
						num.blocks <- num.blocks + 1
						}else{ #otherwise, step through the communities and add each one as a block
							for(cm in 1:(length(cm.changes)+1)){
								if(chr.blocks == 1){
									marker.names <- V(net)$name[1:cm.changes[cm]]
									}
								if(cm > length(cm.changes)){
									marker.names <- V(net)$name[(cm.changes[(cm-1)]+1):length(comm)]
									}
									
								if(cm <= length(cm.changes) && chr.blocks > 1){
									marker.names <- V(net)$name[(cm.changes[cm-1]+1):cm.changes[cm]]
									}
									
								marker.num <- colnames(data.obj$geno)[match(marker.names, data.obj$marker.names)]
								link.blocks[[num.blocks]] <- marker.num
								names(link.blocks)[num.blocks] <- paste("Chr", u_chr[ch], "_", chr.blocks, sep = "")
								num.blocks <- num.blocks + 1
								chr.blocks <- chr.blocks + 1
								}
							}
					
						if(plot.results){
						
						#This function adds axes to the images of the effect matrices
						add.axes <- function(){
							axis(1, at = seq(0,1,1/(length(data.obj$marker.names)+length(main.effects)-1)), labels = c(names(main.effects), data.obj$marker.names), las = 2, cex = 0.7)
						if(length(chr.markers) > 1){
							axis(2, at = seq(0,1,1/(length(chr.markers)-1)), labels = chr.markers, las = 2, cex = 0.7)
							}else{
							axis(2, at = 0.5, labels = chr.markers, las = 2, cex = 0.7)	
							}
						}
							
						layout(matrix(c(1), nrow = 1))
						zmin <- min(combined.mat, na.rm = TRUE); zmax <- max(combined.mat, na.rm = TRUE)
						my.palette <- colorRampPalette(c("lightblue2", "green4"),space = "rgb")
		
						image(rotate.mat(combined.mat), col = my.palette(50), axes = FALSE, main = paste("Combined Interactions and Main Effects Chr", u_chr[ch]), zlim = c(zmin, zmax))
						add.axes()
	
						image(rotate.mat(combined.sig.mat), col = my.palette(50), axes = FALSE, main = paste("Significant Interactions and Main Effects Chr", u_chr[ch]), zlim = c(zmin, zmax))
						add.axes()	
					
						image(rotate.mat(combined.thresh.mat), col = my.palette(50), axes = FALSE, main = paste("Thresholded Interactions and Main Effects Chr", u_chr[ch]), zlim = c(zmin, zmax))
						add.axes()
							
						layout(matrix(c(1,2,0), nrow = 1))
						image(1:dim(cor.mat)[1], 1:dim(cor.mat)[2], cor.mat, main = paste("Marker Correlation Chr", u_chr[ch]), xlim = c(0,(dim(cor.mat)[1]+1)), ylim = c(0,(dim(cor.mat)[1]+1)), col = my.palette(50))
						
					
						image(1:dim(thresh.cor.mat)[1], 1:dim(thresh.cor.mat)[2], thresh.cor.mat, main = paste("Marker Correlation Thresholded With Borders Chr", u_chr[ch]), xlim = c(0,(dim(thresh.cor.mat)[1]+1)), ylim = c(0,(dim(thresh.cor.mat)[1]+1)), col = my.palette(50))				
						block.names <- sapply(strsplit(names(link.blocks), "_"), function(x) x[1])
						chr.names <- sapply(strsplit(block.names, "Chr"), function(x) x[2])
						final.block.locale <- which(chr.names == u_chr[ch])
						start.block = 0.5
						#outline each block
						for(b in 1:length(final.block.locale)){
							end.block <- start.block + length(link.blocks[[final.block.locale[b]]])
							segments(x0 = start.block, y0 = start.block, x1 = start.block, y1 = end.block, lwd = 3)
							segments(x0 = start.block, y0 = start.block, x1 = end.block, y1 = start.block, lwd = 3)
							segments(x0 = end.block, y0 = start.block, x1 = end.block, y1 = end.block, lwd = 3)
							segments(x0 = start.block, y0 = end.block, x1 = end.block, y1 = end.block, lwd = 3)						
							start.block <- end.block
							}
	
	
						# diag(thresh.cor.mat) <- 0
						# hist(rowSums(thresh.cor.mat), main = "Degree Distribution of Thresholded Matrix", xlab = "Degree")
							
						} #end plotting routine
					} #end case for when the chromosome is not 0
				} #end case for chr.markers having a positive length
			} #end looping over chromosomes
		if(plot.results){dev.off()}
		} #end routine for case when collapse.markers is TRUE
		


	if(collapse.linked.markers){
		data.obj$linkage.blocks.collapsed <- link.blocks
		}else{
		data.obj$linkage.blocks.full <- link.blocks
		}

	return(data.obj)

	}