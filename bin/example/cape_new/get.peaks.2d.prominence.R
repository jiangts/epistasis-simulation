#This function finds peaks in a matrix
#based on p values and an effect drop

get.peaks.2d.prominence <- function(m.mat, p.mat, p.or.q){

	# require(igraph)

	#final peaks contains at least all the significant values
	final.peaks <- matrix(NA, nrow = dim(m.mat)[1], ncol = dim(m.mat)[2])
	colnames(final.peaks) <- colnames(m.mat)
	rownames(final.peaks) <- rownames(m.mat)
	sig.locale <- which(p.mat <= p.or.q)
	if(length(sig.locale) == 0){
		return(final.peaks)			
		}
		

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
	

	#This function expands min/max row/col
	expand.range <- function(min.row, max.row, min.col, max.col, border.check, new.peak){
		
			#find out which row/col we are still expanding
			#based on the dimensions of the matrix
			expanding <- list(0, 0, 0, 0)
			names(expanding) <- names(border.check)
			
			not.expanding.val <- NULL
			#record the number of steps from the point we started
			#from.
			peak.coord <- which(new.peak == max(new.peak, na.rm = TRUE), arr.ind = TRUE)
			dim.mat <- unlist(lapply(border.check, length))[c(2,1)]
			if(min.row == 1){
				border.check$top[which(border.check$top == 0)] <- 1
				expanding$top <- not.expanding.val
				}
			if(max.row == dim(m.mat)[1]){
				border.check$bottom[which(border.check$bottom == 0)] <- 1
				expanding$bottom <- not.expanding.val
				}
			if(min.col == 1){
				border.check$left[which(border.check$left == 0)] <- 1
				expanding$left <- not.expanding.val
				}
			if(max.col == dim(m.mat)[2]){
				border.check$right[which(border.check$right == 0)] <- 1
				expanding$right <- not.expanding.val
				}

			
			#now look at whether we've encountered complete borders
			#If we have a complete border on the top, we don't want
			#to keep expanding it. If we are still adding rows to
			#the top, we need to expand the left and right border.checks
			#accordingly (a 0 on the top edge)
			if(length(which(border.check$top == 0)) > 0){ #if we still have at least a partial border
				min.row <- min.row - 1
				}else{
				expanding$top <- not.expanding.val #otherwise set the expanding value to NULL	
				}
			if(length(which(border.check$bottom == 0)) > 0){
				max.row = max.row + 1
				}else{
				expanding$bottom <- not.expanding.val
				}
			if(length(which(border.check$left == 0)) > 0){
				min.col <- min.col - 1
				}else{
				expanding$left <- not.expanding.val	
				}
			if(length(which(border.check$right == 0)) > 0){
				max.col = max.col + 1
				}else{
				expanding$right <- not.expanding.val	
				}
			# cat(min.row, max.row, min.col, max.col)
			
			#add the appropriate indexing terms to all sides. If the top
			#is expanding, add 0's to the right-hand side of the left and
			#right borders
			#if the top is not expanding, we stop adding anything to the
			#right-hand side of the left and right borders, but if these
			#two borders are still expanding, we can still add 0's to the
			#left and right of the top.
			border.check$top <- c(expanding$left, border.check$top, expanding$right)
			border.check$bottom <- c(expanding$left, border.check$bottom, expanding$right)
			border.check$left <- c(expanding$bottom, border.check$left, expanding$top)
			border.check$right <- c(expanding$bottom, border.check$right, expanding$top)			
	
		results <- list(min.row, max.row, min.col, max.col, border.check)
		names(results) <- c("min.row", "max.row", "min.col", "max.col", "border.check")
		return(results)
		}

	#This function checks all edges of a matrix
	#to see if any values are still above the minimum effect
	#size. The orig.check vector 
	compare.borders <- function(peak.mat, border.check){
		#check each border row and compare it with the last border
		#row. If any of the values in the new row are greater than
		#their neighbors in the previous row, we add a stop value 
		#to border.check.
		#border check is a list with one element for each border.
		#Each element is a vector the same length as the
		#current row. The step number is inserted when the edge
		#of the peak has been reached in that direction and we 
		#no longer need to compare the values. Numbers are 
		#carried from one step to the next until all values 
		#in the border check are > 0.
		
		for(i in 1:length(border.check)){
			old.check <- border.check[[i]]
			
			#if we haven't found a complete border yet...
			#if we are looking at the top or bottom, use 
			#columns, if we are looking at the left or right
			#borders, use rows
			
			
			if(length(which(border.check[[i]] == 0)) > 0){
				if(i == 1 || i == 4){
					new.check <- rep(0, dim(peak.mat)[2])
					}else{
					new.check <- rep(0, dim(peak.mat)[1])	
					}
				
				#carry any non-0's from before forward
				new.check[c(which(old.check != 0), which(is.na(old.check)))] <- old.check[c(which(old.check != 0), which(is.na(old.check)))]
	
				#get the values from each outer row/col
				#and the next row/col in
				if(i == 1){
					#top
					new.vals <- peak.mat[1,]
					old.vals <- peak.mat[2,]
					}
				if(i == 2){ 
					#left
					new.vals <- peak.mat[,1]
					old.vals <- peak.mat[,2]
					}
				if(i == 3){ 
					#right
					new.vals <- peak.mat[,dim(peak.mat)[2]]
					old.vals <- peak.mat[,(dim(peak.mat)[2]-1)]
					}
				if(i == 4){ 
					#bottom
					new.vals <- peak.mat[dim(peak.mat)[1],]
					old.vals <- peak.mat[(dim(peak.mat)[1]-1),]
					}
							
				#compare the values that line up with each other
				#subtract the new values from the old values
				border.comparisons <- apply(matrix(c(old.vals, new.vals), ncol = 2, byrow = FALSE), 1, function(x) x[1] - x[2])
				
				#replace the 0s with the appropriate number of steps from the peak
				peak.coord <- which(peak.mat == max(peak.mat, na.rm = TRUE), arr.ind = TRUE)
				dim.mat <- dim(peak.mat)
				
				new.check[intersect(which(new.check == 0), which(border.comparisons < 0))] <- 1
				
				border.check[[i]] <- new.check
				}#end case for continued expanding border
			} #end looping through borders
		
	
		return(border.check)
		}


		#this function takes in a region with a 2d peak and uses an 
		#iterative 1d peak-finding process to narrow the region to 
		#only the prominence of the 2d peak.
		find.peak <- function(peak.mat, start.pt, pmat){
				
			#============================================================================
			#start at the row with the peak, and find the extent of the 
			#1d peak in this row. Fan out to the rest of the rows, each
			#time only looking in the region that previously contained
			#a peak		
			get.1d.peaks <- function(row.col.num, peak.bin.mat, row.or.col = c("row", "col"), plot.peaks = FALSE){
				if(row.or.col == "row"){
					ind <- 1:dim(peak.mat)[2]
					}else{
					ind <- 1:dim(peak.mat)[1]	
					}

				last.row.max <- peak.mat[start.pt[1],start.pt[2]] #keep track of the peak maximum, stop when the maximum goes up between rows.cols
				row.col.diff <- 0
				pvals <- 0
				r = 1
				while(r <= length(row.col.num) && length(ind) > 1 && (row.col.diff <= 0 || length(which(pvals <= p.or.q)) > 0)){
					if(row.or.col == "row"){
						vals <- peak.mat[row.col.num[r],ind]
						pvals <- pmat[row.col.num[r], ind]
						}else{
						vals <- peak.mat[ind, row.col.num[r]]	
						pvals <- pmat[ind, row.col.num[r]]
						}
					
					if(length(which(!is.na(vals))) == 0){
						r <- length(row.col.num)
						next()
						}	
						
					consec.vals <- consec.pairs(vals)
					consec.ind <- consec.pairs(ind)
					diffs <- apply(consec.vals, 1, function(x) x[2] - x[1])
	
					#find the peak
					peak.locale <- which(vals == max(vals, na.rm = TRUE))
					row.col.diff <- max(vals, na.rm = TRUE) - last.row.max
					last.row.max <- max(vals, na.rm = TRUE)
					
					#find where the peak starts to go up
					upswing1 <- which(diffs > 0)
					if(length(upswing1) > 0){
						peak.start <- min(upswing1)
						}else{
						peak.start <- 1	
						}
					#also check to see if there are significant interactions
					#in this row. If they extend beyond upswing1, expand the
					#region to include all significant interactions
					sig.locale <- which(pvals <= p.or.q)
					if(length(sig.locale) == 0){
						sig.locale <- peak.start	
						}

					if(min(sig.locale) < peak.start){
						peak.start <- min(sig.locale)
						}
					
					#find the first position after the maximum of
					#the last significant value and the peak location 
					#where the vals start to go up again
					if(length(sig.locale) > 0){
						max.sig.locale <- max(sig.locale)
						}else{
						max.sig.locale <- 1	
						}
					max.pos <- max(c(max.sig.locale, peak.locale))
					after.peak.ind <- c(max.pos:length(diffs))
					after.peak.vals <- diffs[after.peak.ind]
					upswing2 <- which(after.peak.vals > 0)
					if(length(upswing2) == 0){
						peak.end <- max(after.peak.ind)+1
						}else{
						peak.end <- after.peak.ind[min(upswing2)]
						}
					
					ind <- ind[peak.start:peak.end]
					# print(ind)
					if(row.or.col == "row"){
						peak.bin.mat[row.col.num[r],ind] <- 1
						# peak.bin.mat[row.col.num[r],peak.start:peak.end] <- 1
						}else{
						peak.bin.mat[ind, row.col.num[r]] <- 1
						# peak.bin.mat[peak.start:peak.end, row.col.num[r]] <- 1	
						}
					r = r + 1
					
					if(plot.peaks){
						quartz();plot(vals, type = "l")
						arrows(peak.start, y0 = (vals[peak.start]+max(vals)*0.1), y1 = vals[peak.start], col = "red", lwd = 3)
						arrows(peak.end, y0 = (vals[peak.end]+max(vals)*0.1), y1 = vals[peak.end], col = "red", lwd = 3)
						}
					}

				if(plot.peaks){
					mat.to.plot <- peak.mat
					mat.to.plot[which(peak.bin.mat == 0)] <- 0
					quartz();plot.image.with.text(mat.to.plot)
					}
				return(peak.bin.mat)		
				}				
		#============================================================================			
						
				#find the coordinates of the overall peak
				# peak.coord <- which(peak.mat == max(peak.mat, na.rm = TRUE), arr.ind = TRUE)
				peak.coord = start.pt

				horiz.peaks <- matrix(0, nrow = dim(peak.mat)[1], ncol = dim(peak.mat)[2])
				rows.above <- peak.coord[1]:1
				rows.below <- peak.coord[1]:dim(peak.mat)[1] #always start with the peak row itself
				horiz.peaks <- get.1d.peaks(row.col.num = rows.above, peak.bin.mat = horiz.peaks, row.or.col = "row", plot.peaks = TRUE)
				horiz.peaks <- get.1d.peaks(row.col.num = rows.below, peak.bin.mat = horiz.peaks, row.or.col = "row", plot.peaks = TRUE)

				vert.peaks <- matrix(0, nrow = dim(peak.mat)[1], ncol = dim(peak.mat)[2])	
				col.left <- peak.coord[2]:1
				col.right <- peak.coord[2]:dim(peak.mat)[2]			
				vert.peaks <- get.1d.peaks(row.col.num = col.left, peak.bin.mat = vert.peaks, row.or.col = "col", plot.peaks = TRUE)
				vert.peaks <- get.1d.peaks(row.col.num = col.right, peak.bin.mat = vert.peaks, row.or.col = "col", plot.peaks = TRUE)
				
				peak.intersect <- horiz.peaks*vert.peaks
				#peak.intersect <- horiz.peaks+vert.peaks
				
				# image(rotate.mat(peak.intersect))
				
				delete.vals <- which(peak.intersect == 0)
				final.mat <- peak.mat
				final.mat[delete.vals] <- NA
				# image(rotate.mat(final.mat))
				return(final.mat)
			}
		
		
		#This function takes in border check, and returns a logical 
		#value indicating whether all edges of the peak have been
		#found.
		found.all.borders <- function(border.check){
			borders.incomplete <- unlist(lapply(border.check, function(x) length(which(x == 0))))
			if(length(which(borders.incomplete == 0)) == length(border.check)){
				return(TRUE)
				}else{
				return(FALSE)
				}			
			}
	
		
		
		#===========================================================================
		# end of internal functions
		#===========================================================================		

	
	#make a lattice graph out of the p values
	mat.lat <- mat.lat.graph(lat.mat = m.mat)

	#delete the non-significant vertices
	p.lat.vals <- unlist(apply(p.mat, 2, list))
	mat.lat <- filter.lat(lat.graph = mat.lat, vals = p.lat.vals, filter.val = p.or.q, filter.out = "gr")
	
	comm <- igraph::clusters(mat.lat)$membership
	comm.num <- max(comm)
	# plot(mat.lat, vertex.label = comm, vertex.color = comm+1)
	
	
	#each community represents a peak
	#for each one, find the node with 
	#the maximum effect size
	
	for(cl in 1:comm.num){
		# print(cl)
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
		
		#initialize border.check
		# border.check <- list(TRUE, TRUE, TRUE, TRUE)
		border.check <- list(0, 0, 0, 0)
		names(border.check) <- c("top", "left", "right", "bottom")
		found.borders <- found.all.borders(border.check)
		while(!found.borders){
			#add edges to the peak region at a time and check each border
			#to see if we've gone beyond the peak
			expanded <- expand.range(min.row, max.row, min.col, max.col, border.check, new.peak)
			min.row = expanded$min.row; max.row = expanded$max.row; min.col = expanded$min.col; max.col = expanded$max.col
			border.check <- expanded$border.check				
			# cat(min.row, max.row, min.col, max.col, "\n")
			# border.check
			
			new.peak <- m.mat[min.row:max.row, min.col:max.col,drop=FALSE]
			# dim(new.peak)
			# unlist(lapply(border.check, length))
			
			# quartz();image(new.peak, main = c(min.row, max.row, min.col, max.col))
			border.check <- compare.borders(peak.mat = new.peak, border.check)
			# border.check
			found.borders <- found.all.borders(border.check)
			}

		# quartz();plot.image.with.text(new.peak)

		#delete the markers beyond borders using the coordinates in 
		#border.check.
		#for each column across the top, for example, border.check
		#indicates how many rows up from the peak marker the border
		#exists.
		#specify the start point in case a peak shares a region with
		#a taller peak
		start.pt <- c(which(rownames(new.peak) == vertex.pair[1]), which(colnames(new.peak) == vertex.pair[2]))
		#also include the p values from this region. We want to include
		#at minimum all significant positions
		new.pmat <- p.mat[which(rownames(p.mat) %in% rownames(new.peak)), which(colnames(p.mat) %in% colnames(new.peak))]
		trimmed.peak <- find.peak(peak.mat = new.peak, start.pt = start.pt, pmat = new.pmat)
		# quartz();image(rotate.mat(new.peak));quartz();image(rotate.mat(trimmed.peak))
			
		#find out where these markers are and add them to final.peaks
		col.pos <- which(colnames(m.mat) %in% colnames(trimmed.peak))
		row.pos <- which(rownames(m.mat) %in% rownames(trimmed.peak))

		final.peaks[row.pos,col.pos] <- trimmed.peak
		}
	
	#make sure all the significant values are in the final matrix
	final.peaks[sig.locale] <- m.mat[sig.locale]
		
	return(final.peaks)
	}
	
