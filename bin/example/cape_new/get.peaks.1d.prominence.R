#This function finds peaks in a matrix
#by column based on p values and an effect 
#drop


get.peaks.1d.prominence <- function(m.mat, p.mat, p.or.q){

	# require(igraph)

	#restrict final peaks to only significant values.
	#these are the minimum values included for each
	#peak. This function extends peaks by prominence.
	final.peaks <- matrix(NA, nrow = dim(m.mat)[1], ncol = dim(m.mat)[2])
	colnames(final.peaks) <- colnames(m.mat)
	rownames(final.peaks) <- rownames(m.mat)
	
	sig.locale <- which(p.mat <= p.or.q)
	
	if(length(sig.locale) == 0){
		return(final.peaks)
		}

	final.peaks[sig.locale] <- m.mat[sig.locale]

	#===========================================================================
	#internal functions
	#===========================================================================

	#This function generates a lattice graph from a given matrix
	mat.lat.graph <- function(lat.mat){
		mat.lat <- graph.lattice(dimvector = dim(lat.mat))
		V(mat.lat)$vals <- unlist(apply(lat.mat, 2, list))
		
		dummy.mat <- matrix(1, nrow = dim(lat.mat)[1], ncol = dim(lat.mat)[2])
		all.ind <- which(dummy.mat == 1, arr.ind = TRUE)
		v.names <- apply(all.ind, 1, function(x) paste(rownames(lat.mat)[x[1]], colnames(lat.mat)[x[2]], collapse = ","))
		V(mat.lat)$name <- v.names
		V(mat.lat)$rows <- all.ind[,1]
		V(mat.lat)$cols <- all.ind[,2]
		V(mat.lat)$colors <- 1:vcount(mat.lat)
		return(mat.lat)
		}
		
	#This function removes vertices in a lattice graph based on
	#based on given values
	filter.lat <- function(lat.graph, vals = V(lat.graph)$vals, filter.val, filter.out = c("gr", "gre", "l", "le")){
		if(filter.out == "gr"){
			bad.verts <- which(vals > filter.val)			
			}
		if(filter.out == "gre"){
			bad.verts <- which(vals >= filter.val)
			}
		if(filter.out == "l"){
			bad.verts <- which(vals < filter.val)			
			}
		if(filter.out == "le"){
			bad.verts <- which(vals <= filter.val)
			}

		if(length(bad.verts) > 0){
			lat.graph <- delete.vertices(lat.graph, bad.verts)
			}
		na.verts <- which(is.na(V(lat.graph)$vals))
		if(length(na.verts) > 0){
			lat.graph <- delete.vertices(lat.graph, na.verts)
			}
		return(lat.graph)
		}
	
	
				
		#============================================================================
		#start at the row with the peak, and find the extent of the 
		#1d peak in this row. Fan out to the rest of the rows, each
		#time only looking in the region that previously contained
		#a peak		
		get.peak.prom <- function(vals, pvals, start.pt, plot.peaks = FALSE){

			peak.bin.mat <- rep(NA, length(vals))
			ind <- 1:length(vals)

			consec.vals <- consec.pairs(vals)
			consec.ind <- consec.pairs(ind)
			diffs <- apply(consec.vals, 1, function(x) x[2] - x[1])

			#find the peak
			peak.locale <- start.pt
			#and all significant values
			sig.locale <- which(pvals <= p.or.q)
				
			#find where the peak starts to go up
			upswing1 <- which(diffs > 0)
			if(length(upswing1) > 0){
				peak.start <- min(upswing1)
				}else{
				peak.start <- 1	
				}
			#if there are significant values beyond
			#the peak start, extend the start to the
			#first significant value
			if(length(sig.locale) > 0){
				min.sig <- min(sig.locale)
				if(min.sig < peak.start){
					peak.start <- min.sig
					}
				}
			
			#find the first position after the peak or
			#after the last significant value where it
			#starts to go up again
			
			if(length(sig.locale) > 0){
				max.sig <- max(sig.locale)
				}else{
				max.sig <- peak.locale
				}
			
			after.peak.ind <- c(max.sig:length(diffs))
			after.peak.vals <- diffs[after.peak.ind]
			upswing2 <- which(after.peak.vals > 0)
			if(length(upswing2) == 0){
				peak.end <- max(after.peak.ind)+1
				}else{
				peak.end <- after.peak.ind[min(upswing2)]
				}
				
			peak.bin.mat[peak.start:peak.end] <- vals[peak.start:peak.end]
						
			if(plot.peaks){
				quartz();plot(vals, type = "l")
				arrows(peak.start, y0 = (vals[peak.start]+max(vals)*0.1), y1 = vals[peak.start], col = "red", lwd = 3)
				arrows(peak.end, y0 = (vals[peak.end]+max(vals)*0.1), y1 = vals[peak.end], col = "red", lwd = 3)
				}

			return(peak.bin.mat)		
			}	

	
		#===========================================================================
		# end of internal functions
		#===========================================================================			
	
	#work on each phenotype individually
	for(i in 1:dim(m.mat)[2]){
		all.vals <- m.mat[,i]
		pvals <- p.mat[,i]
		sig.val.locale <- which(pvals <= p.or.q)
		
		#make a lattice graph out of the p values
		mat.lat <- mat.lat.graph(lat.mat = m.mat[,i,drop=FALSE])
	
		#delete the non-significant vertices
		p.lat.vals <- unlist(apply(p.mat[,i,drop=FALSE], 2, list))
		mat.lat <- filter.lat(lat.graph = mat.lat, vals = p.lat.vals, filter.val = p.or.q, filter.out = "gr")
		
		if(vcount(mat.lat) > 0){
			comm <- igraph::clusters(mat.lat)$membership
			comm.num <- max(comm)
			# plot(mat.lat, vertex.label = NA, vertex.color = comm+1)
			
			#each community represents a peak
			#for each one, find the node with 
			#the maximum effect size			
			for(cl in 1:comm.num){
				comm.locale <- which(comm == cl)
				vertex.names <- V(mat.lat)$name[comm.locale]
				marker.names <- unlist(lapply(strsplit(vertex.names, " "), function(x) x[1]))
				vertex.vals <- V(mat.lat)$vals[comm.locale]
				max.locale <- which(vertex.vals == max(vertex.vals))
				peak.locale <- which(rownames(m.mat) == marker.names[max.locale])
				
				peak.prom <- get.peak.prom(vals = all.vals, pvals = pvals, start.pt = peak.locale, plot.peaks = FALSE)
				prom.locale <- which(!is.na(peak.prom))
								
				final.peaks[prom.locale,i] <- peak.prom[prom.locale]
				}
			}
		}
	
	# quartz();image(rotate.mat(m.mat)); quartz();image(rotate.mat(final.peaks))
		
	return(final.peaks)
	}
	
