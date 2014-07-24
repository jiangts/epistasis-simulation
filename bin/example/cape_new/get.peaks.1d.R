#This function finds peaks in a matrix
#by column based on p values and an effect 
#drop


get.peaks.1d <- function(m.mat, p.mat, p.or.q, effect.drop){

	# require(igraph)

	final.peaks <- matrix(NA, nrow = dim(m.mat)[1], ncol = dim(m.mat)[2])
	colnames(final.peaks) <- colnames(m.mat)
	rownames(final.peaks) <- rownames(m.mat)

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
	
	
	#This function checks all edges of a matrix
	#to see if any values are still above the minimum effect
	#size. The orig.check vector 
	check.borders <- function(peak.mat, min.effect, border.check, step.num){
		#check the minimum and maximum rows and columns
		#to see which are below the minimum value
		#if all the values in one of the border vectors
		#is completely 0, change that vector's value
		#in border.check to TRUE
		
		min.row.check <- which(peak.mat[1,] < min.effect)
		if(length(min.row.check) == dim(peak.mat)[1] && border.check[1,1] == 0){
			border.check[1,1] <- step.num
			}
		max.row.check <- which(peak.mat[dim(peak.mat)[1],] < min.effect)
		if(length(max.row.check) == dim(peak.mat)[1] && border.check[1,2] == 0){
			border.check[1,2] <- step.num
			}
		min.col.check <- which(peak.mat[,1] < min.effect)
		if(length(min.col.check) == dim(peak.mat)[1] && border.check[2,1] == 0){
			border.check[2,1] <- step.num
			}
		max.col.check <- which(peak.mat[,dim(peak.mat)[2]] < min.effect)
		if(length(max.col.check) == dim(peak.mat)[1] && border.check[2,2] == 0){
			border.check[2,2] <- step.num
			}
		return(border.check)
		}
		#===========================================================================
		# end of internal functions
		#===========================================================================		

	sig.locale <- which(p.mat <= p.or.q)
	
	if(length(sig.locale) == 0){
		return(final.peaks)
		}
	
	for(i in 1:dim(m.mat)[2]){
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
				vertex.vals <- V(mat.lat)$vals[comm.locale]
				max.vert <- vertex.names[which(vertex.vals == max(vertex.vals))][1]
				vertex.pair <- strsplit(max.vert, " ")[[1]]
				
				#find the location of the maximum node in the effects matrix 
				row.num <- which(rownames(m.mat) == vertex.pair[1])
				col.num <- which(colnames(m.mat) == vertex.pair[2])
				
				#sweep out from the peak node and expand the peak 
				new.peak <- m.mat[row.num,col.num]
				min.row <- row.num; max.row <- row.num
				min.col <- col.num; max.col <- col.num
				
				min.effect <- min(vertex.vals) - effect.drop
				border.check <- matrix(0, nrow = 1, ncol = 2)
				num.borders.to.find <- length(which(border.check == 0))
				step.num = 1
				while(num.borders.to.find > 0){
					#add edges to the peak region at a time and check each border
					#to see if we've gone beyond the peak
					if(min.row == 1){border.check[1,1] <- step.num}; if(border.check[1,1] == 0){min.row <- min.row - 1}
					if(max.row == dim(m.mat)[1]){border.check[1,2] <- step.num}; if(border.check[1,2] == 0){max.row = max.row + 1}
					new.peak <- m.mat[min.row:max.row, min.col:max.col,drop=FALSE]
					# image(new.peak)
					# cat(min.row, max.row, "\n")
					border.check <- check.borders(new.peak, min.effect, border.check, step.num = step.num)
					num.borders.to.find <- length(which(border.check == 0))
					step.num <- step.num + 1
					}
				
				#put the values that are above the minimum effect size into the 
				low.vals <- which(new.peak < min.effect)
				new.peak[low.vals] <- NA
				col.pos <- which(colnames(m.mat) %in% colnames(new.peak))
				row.pos <- which(rownames(m.mat) %in% rownames(new.peak))
				
				final.peaks[row.pos,col.pos] <- new.peak
				
				}
			}
		}
	
		
	return(final.peaks)
	}
	
