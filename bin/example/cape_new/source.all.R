source.all <-
function(){

	
	files <- matrix(list.files(pattern = "\\.R"), ncol = 1)
	
	for(i in 1:length(files)){
    print(files[i])
    if(files[i] != "source.all.R")
    {
      source(files[i])
      
    }
		}

	
	
	
	}
