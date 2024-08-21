library(data.tree)
library(VLMC)
library(stats)

#used in tDistance
vlmcTreePaths= function(vlmcObj)
{
	dendrogram= as.dendrogram(vlmcObj)
	#plot(dendrogram)
	allPaths=vector(mode="character")
	for(i in 1:length(names(dendrogram)))
	{
		string= capture.output(dendrogram[i])
		entries= which(grepl("$", string, fixed=TRUE))
		paths= string[entries]
		for( j in 1: length(paths))
		{
			paths[j]= gsub("$", "", paths[j], fixed=TRUE)
			paths[j]= gsub("`", "", paths[j], fixed=TRUE)
			paths[j]= gsub("NA", "", paths[j], fixed=TRUE)
		}
		allPaths=append(allPaths, paths)
	}
	return(allPaths)	
}

reverse = function(string)
{
  string_split = strsplit(as.character(string), split = "")
  reversed_split = string_split[[1]][nchar(string):1]
  paste(reversed_split, collapse = "")
} 



pDistance= function(VLMC1, VLMC2, textOutput1, textOutput2)
{
	#textOutput1 and textOutput2 should be the text output files scan(paste("Data/Poetry_Decoding/Poetry_",n,"_output", ".txt", sep="")) 
	Pa= vlmcTreePaths(VLMC1)
	#x11()
	Pb= vlmcTreePaths(VLMC2)
	
	diffAB= setdiff(Pa, Pb) #elements in Pa but not Pb
	diffBA= setdiff(Pb, Pa) #elements in Pb but not Pa
	union= intersect(Pa, Pb)
	distance=0
	if(length(diffAB)>0){
		for(i in 1:length(diffAB))
		{
			distance=distance+ 1/nchar(diffAB[i])
		}
	}
	if(length(diffBA)>0){
		for(i in 1:length(diffBA))
		{
		distance= distance+ 1/nchar(diffBA[i])
		}
	}
	
	#find distance between transition probabilities
	states=c("0","1","2","3","4")
	pDist=0
	text1=paste(textOutput1, collapse="")
	text2=paste(textOutput2, collapse="")
	
	for(i in 1: length(union))
	{
	
		sequence = reverse(union[i])
		occurences1 = gregexpr(sequence,text1)[[1]] #location where the path occurs 
		occurences2 = gregexpr(sequence,text2)[[1]]
		next1 = occurences1+nchar(sequence) #index of next term (for probabilities)
		next2 = occurences2+nchar(sequence)
		
		#check if pattern is at end of file (edge case)
		if(next1[length(next1)]> nchar(text1))
			next1=next1[-length(next1)]
			
		if(next2[length(next2)]> nchar(text2))
			next2=next2[-length(next2)]
			
				
		nextVal1 = as.vector(strsplit(text1, split=NULL)[[1]])[next1] #finds actual next value after pattern 
		nextVal2 = as.vector(strsplit(text2, split=NULL)[[1]])[next2]
		counts1= rep(0, length(states))
		counts2= rep(0, length(states))
		
		for(j in 1:length(nextVal1))
		{
			for(k in 1: length(states))
			{
				if(nextVal1[j]==states[k])
				{
					counts1[k]= counts1[k]+1
				}
			}
		}
		
		for(j in 1:length(nextVal2))
		{
			for(k in 1: length(states))
			{
				if(nextVal2[j]==states[k])
				{
					counts2[k]= counts2[k]+1
				}
			}
		}
		
		pDist=pDist + sum(abs((counts1/sum(counts1))-(counts2/sum(counts2))))
	}
	alpha= .5
	distance = alpha*distance + (1-alpha)*pDist
	
		return(distance)
}



##### TESTING 

numUS= 45
numUK= 45
vlmcs = list()

for (n in 1:numUS)
	{
	    file <-  paste("Data/US_Decoding/US_",n,"_output", ".txt", sep="")
		data_US = scan(file)
		vlmcs[[n]] = list(0,vlmc(dts=data_US), data_US)
	}

	for (n in 1:numUK)
	{
		file <-  paste("Data/Poetry_Decoding/Poetry_",n,"_output", ".txt", sep="")
	    #file <-  paste("Data/UK_Decoding/UK_",n,"_output", ".txt", sep="")
		data_UK = scan(file)
		vlmcs[[n+numUS]] = list(1,vlmc(dts=data_UK), data_UK)
	}



pDistance(vlmcs[[1]][[2]], vlmcs[[2]][[2]],vlmcs[[1]][[3]], vlmcs[[2]][[3]])







