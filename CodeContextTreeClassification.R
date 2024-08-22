library(data.tree)
library(VLMC)
library(stats)
library(gtools)
library(cluster)
library(factoextra)
library(mclust)
library(infotheo) ## for mutualinformation()
library(fossil) ## for rand.index()



strReverse <- function(x)
        sapply(lapply(strsplit(x, NULL), rev), paste, collapse="")


# 
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

# -- distance based only on the structure
tDistance= function(VLMC1, VLMC2)
{
	Pa= vlmcTreePaths(VLMC1)
	#x11()
	Pb= vlmcTreePaths(VLMC2)
	
	diffAB= setdiff(Pa, Pb) #elements in Pa but not Pb
	diffBA= setdiff(Pb, Pa) #elements in Pb but not Pa
	distance=0
	if (length(diffAB) > 0)
	{
		for(i in 1:length(diffAB))
		{
			distance=distance+ 1/nchar(diffAB[i])
		}
	}
	if (length(diffBA) > 0)
	{
	for(i in 1:length(diffBA))
	{
		distance= distance+ 1/nchar(diffBA[i])
	}
	}
	
	return(distance)
}

# distance based on the structure and the probabilities - need the data (textOutput1 and 2) to estimate the probabilities 
pDistance <- function(VLMC1, VLMC2, textOutput1, textOutput2)
{
	Pa = vlmcTreePaths(VLMC1)
	Pb = vlmcTreePaths(VLMC2)	
	union = intersect(Pa, Pb)
	
	#find distance between transition probabilities
	states = as.character(sort(unique(c(textOutput1, textOutput2))))
	pDist = 0
	text1 = paste(textOutput1, collapse="")
	text2 = paste(textOutput2, collapse="")
	
	if (length(union) > 0) ## one of the trees is just the root -> order 0
	for(i in 1: length(union))
	{
	
		sequence = strReverse(union[i])
		occurences1 = gregexpr(sequence,substr(text1,1,(nchar(text1)-1)))[[1]] #location where the path occurs  - the substr avoids the last observation (edge case there is no "next")
		occurences2 = gregexpr(sequence,substr(text2,1,(nchar(text2)-1)))[[1]]
		next1 = occurences1+nchar(sequence) #index of next term (for probabilities)
		next2 = occurences2+nchar(sequence)
		
			
				
		nextVal1 = as.vector(strsplit(text1, split=NULL)[[1]])[next1] #finds actual next value after pattern 
		nextVal2 = as.vector(strsplit(text2, split=NULL)[[1]])[next2]
		counts1 = table(nextVal1)
		counts2 = table(nextVal2)		
		
		if (!identical(names(table(nextVal1)), names(table(nextVal2)))) ## not all alphabet is in next symbols
		{	
			counts1_aux = rep(0,5)
			counts2_aux = rep(0,5)
			elements = c("0","1","2","3","4")
			for (j in 1:length(elements))
			{
				counts1_aux[j]  = counts1[elements[j]]
				counts2_aux[j]  = counts2[elements[j]]				
			}

			counts1_aux[is.na(counts1_aux)] = 0
			counts1 = counts1_aux
			counts2_aux[is.na(counts2_aux)] = 0
			counts2 = counts2_aux
		}
		
		pDist=pDist + sum(abs((counts1/sum(counts1))-(counts2/sum(counts2))))
	}

	return(pDist)
}






	

	
KNN <- function(data,classifyEntry, n.neighbors)
{ 
	#data is list with membership (0,1) and VLMC object
	#classifyEntry is a vlmcs object

 	k = n.neighbors   #how many nearest neighbors
	vlmcNEW = classifyEntry
		
	distances = 0  #will be vector of distances of vlmc to all other vlmc
	for (j in 1:length(data)) 
		{
			distances[j] = tDistance(classifyEntry, data[[j]][[2]])
		}
	
		majority = 0 # will add/subtract, look for sign of result to determine vote result
		order_dist = order(distances)[1:k] # orders which vlmc objects had least distance
		for (j in order_dist)
			majority = majority + ifelse(data[[j]][[1]] == 1,1,-1)
		#vote is 1 for label 1 and -1 for label 0
		return(sign(majority))
}


## classic
knn <- function(data_set_train, data_set_test, true_labels_train, n.neighbors)
{
	result = 0
	for (i in 1:length(data_set_test[,1]))
	{
		dist_to_train_data = rowSums(abs(sweep(data_set_train, 2, data_set_test[i,])))
		which_train = order(dist_to_train_data)[1:n.neighbors]
		majority = 0
		for (j in which_train)
			majority = majority + ifelse(true_labels_train[j] == 1,1,-1)
		
		if (majority == 0)
			result[i] = rbinom(1,1,.5)
		else
			result[i] = ifelse(majority > 0, 1, 0)
	}
	return(result)
}



KNN_dist_matrix <- function(fulldistanceMatrix, train_sample_index, test_sample_index, all_lables, n.neighbors)
{ 
	result = 0
	for (i in 1:length(test_sample_index))
	{
		distances = fulldistanceMatrix[,test_sample_index[i]]
		dist_to_train_data = distances[train_sample_index]
		order_dist = order(dist_to_train_data)[1:n.neighbors]
		which_train = train_sample_index[order_dist]

		result[i] = as.numeric(majorityVote(all_lables[which_train])$majority)
	}

	## returns a vector of the classification for each test data
	return(result)
}


	
	
Mismatch <- function(mm, true_clusters, K)
{
    sigma = permutations(n=K,r=K,v=1:K)

    Miss = length(which( true_clusters != mm))  ## for permutation 1, 2,... K

    mm_aux = mm
    for (ind in 2:dim(sigma)[1])
    {

       for (j in 1:K)
          mm_aux[which(mm == j)] = sigma[ind,j]
    
       Miss[ind] =  length(which( true_clusters != mm_aux))
    
    }

       
    return(min(Miss))
}
	







##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################
## Example
##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################


dataA = c(1, 3, 3, 3, 0, 0, 2, 1, 3, 3, 3, 2, 1, 0, 0, 3, 0, 4, 3, 2, 1, 2, 1, 0, 0, 4, 2, 1, 0, 3, 3, 0, 3, 3, 3, 3, 0, 0, 2, 1, 3, 0, 0, 3, 3, 3, 2, 1, 0, 0, 3, 0, 3, 3, 0, 3, 0, 3, 3, 0, 3, 0, 3, 0, 3, 4, 3, 2, 1, 0, 0, 3, 2, 1, 2, 1, 0, 3, 3, 3, 3, 2, 1, 3, 3, 0, 2, 1, 4, 3, 0, 3, 3, 3, 3, 4, 3, 0, 3, 3)
draw(vlmc(dataA))
vlmcA = vlmc(dataA)


dataB = c(3, 1, 3, 1, 3, 3, 3, 0, 2, 1, 2, 1, 3, 3, 3, 3, 3, 4, 3, 0, 3, 3, 0, 3, 3, 3, 0, 0, 3, 0, 3, 0, 3, 2, 1, 0, 0, 3, 3, 3, 0, 2, 0, 1, 0, 2, 1, 0, 0, 3, 2, 0, 1, 2, 0, 1, 0, 4, 3, 0, 3, 0, 2, 1, 4, 3, 1, 1, 0, 3, 3, 3, 3, 2, 1, 3, 3, 3, 3, 3, 0, 3, 3, 2, 1, 3, 2, 1, 3, 0, 3, 0, 0, 3, 4, 3, 2, 1, 4, 3, 0, 3, 3, 0, 4, 3, 0, 0, 3, 0, 3, 3, 2, 1, 0, 3, 2, 1, 3, 3, 3, 2, 1, 0, 3, 3, 4, 2, 1, 3, 4, 3, 2, 1, 3, 2, 0, 1, 3, 2)
draw(vlmc(dataB))
vlmcB = vlmc(dataB)

n.time_series = 500 ## how many observations of the chain we are generating
all_results = list() ## each element will hold all results of the simulation for a specific time series n

n = 80
n.neighbors = 7
vlmcs = list()
set.seed(1)

	rateKNN = matrix(0,4,1)
	rateKNNtrain = matrix(0,4,1)
	rate_Kmedoids = matrix(0,4,1)
	rate_Kmedoids_randindex = matrix(0,4,1)
	rate_Kmedoids_mutualinfo = matrix(0,4,1)
	rateClassicKNN = 0
	rateClassicKNN_train = 0
	rate_Kmeans = 0
	rate_Kmeans_rand_index = 0
	rate_Kmeans_mutual_info = 0

		cat("\n\n -- Scenario 1 ----- example:", " --- n.time_series: ",n.time_series, " ---------")
	
		##############################################################################################
		## simulate
		##############################################################################################
		for (index in 1:(n/2))
		{
			simulated_data = simulate(vlmcA, n.time_series)
			vlmcs[[index]] = list(0,vlmc(simulated_data, threshold.gen=10, alpha=0.05), as.numeric(simulated_data))
		}
		for (index in (n/2+1):n)
		{
			simulated_data = simulate(vlmcB, n.time_series)
			vlmcs[[index]] = list(1,vlmc(simulated_data, threshold.gen=10, alpha=0.05), as.numeric(simulated_data))
		}	
		
		vlmcs = vlmcs[sample(n)]
		
		##############################################################################################
		## tDistance is the distance based on the contexts only
		tdistanceMatrix = matrix(0, n, n)
		for(index in 2:n) #first entry a_11 = 0
		{
			for(j in 1:index)
			tdistanceMatrix[index,j]= tDistance(vlmcs[[index]][[2]], vlmcs[[j]][[2]])
		}
		tfulldistanceMatrix = tdistanceMatrix + t(tdistanceMatrix)
	
		##############################################################################################
		## pDistance is the distance based on the probabilities of the reduced subgraph
		pdistanceMatrix = matrix(0, n, n)
		for(index in 2:n) #first entry a_11 = 0
		{
			for(j in 1:index)
			pdistanceMatrix[index,j]= pDistance(vlmcs[[index]][[2]], vlmcs[[j]][[2]], vlmcs[[index]][[3]], vlmcs[[j]][[3]])
			
		}
		pdistanceMatrix[is.na(pdistanceMatrix)] = 0
		pfulldistanceMatrix = pdistanceMatrix + t(pdistanceMatrix)
	
		

		K = 2	
	
		all_lables = do.call("rbind", lapply(vlmcs, "[[", 1))
		
		
		
		
		##############################################################################################
		##############################################################################################
		## alpha homogeneity individually - choose alpha based on D(A,B)/alpha_divisor - CLUSTERING
		##############################################################################################
		##############################################################################################
		alpha_divisor = length(unique(dataA))^2/4 
		homogeneity = 0
		for (ind.alpha_divisor in 1:length(alpha_divisor))
		{
			alpha = pmin(tfulldistanceMatrix/alpha_divisor[ind.alpha_divisor],1)
			fulldistanceMatrix = alpha*tfulldistanceMatrix + (1-alpha)*pfulldistanceMatrix
			distanceMatrix = alpha*tdistanceMatrix + (1-alpha)*pdistanceMatrix
			k_medoids = pam(as.dist(distanceMatrix),k = K, medoids = sample(1:n,K))
			#### compute homogeneity of clusters
			homogeneity[ind.alpha_divisor] = 0
			for (ind.k in 1:K)
				homogeneity[ind.alpha_divisor] = homogeneity[ind.alpha_divisor] + sum(fulldistanceMatrix[k_medoids$medoids[ind.k],which(k_medoids$clustering == ind.k)])
		}
	
		alpha_divisor = alpha_divisor[which.min(homogeneity)]
			
		alpha = pmin(tfulldistanceMatrix/alpha_divisor,1)

		fulldistanceMatrix_indv_alpha = alpha*tfulldistanceMatrix + (1-alpha)*pfulldistanceMatrix
		distanceMatrix_indv_alpha = alpha*tdistanceMatrix + (1-alpha)*pdistanceMatrix
		##############################################################################################
		##############################################################################################

		##############################################################################################
		##############################################################################################
		## alpha homogeneity fixed - CLUSTERING
		##############################################################################################
		##############################################################################################
		alphas = seq(0.01, 0.99, length = 20)
		homogeneity = 0
		for (ind.alphas in 1:length(alphas))
		{
			alpha = alphas[ind.alphas]
			fulldistanceMatrix = alpha*tfulldistanceMatrix + (1-alpha)*pfulldistanceMatrix
			distanceMatrix = alpha*tdistanceMatrix + (1-alpha)*pdistanceMatrix
			k_medoids = pam(as.dist(distanceMatrix),k = K, medoids = sample(1:n,K))
			#### compute homogeneity of clusters
			homogeneity[ind.alphas] = 0
			for (ind.k in 1:K)
				homogeneity[ind.alphas] = homogeneity[ind.alphas] + sum(fulldistanceMatrix[k_medoids$medoids[ind.k],which(k_medoids$clustering == ind.k)])
		}
	
		alpha = alphas[which.min(homogeneity)]
			
		fulldistanceMatrix_alpha_fixed = alpha*tfulldistanceMatrix + (1-alpha)*pfulldistanceMatrix
		distanceMatrix_alpha_fixed = alpha*tdistanceMatrix + (1-alpha)*pdistanceMatrix
		##############################################################################################
		##############################################################################################


		##############################################################################################
		##############################################################################################
		## alpha homogeneity individually - choose alpha based on D(A,B)/alpha_divisor - CLASSIFICATION
		##############################################################################################
		##############################################################################################
		alpha_divisor = length(unique(dataA))^2/4 
		homogeneity = 0
		for (ind.alpha_divisor in 1:length(alpha_divisor))
		{
			alpha = pmin(tfulldistanceMatrix/alpha_divisor[ind.alpha_divisor],1)
			fulldistanceMatrix = alpha*tfulldistanceMatrix + (1-alpha)*pfulldistanceMatrix
			distanceMatrix = alpha*tdistanceMatrix + (1-alpha)*pdistanceMatrix
			#### compute homogeneity of clusters
			homogeneity[ind.alpha_divisor] = 0
			for (ind.k in 1:K)
			{
				## find the medoids of the groups from the classification labels
				med = pam(as.dist(distanceMatrix[which(all_lables == (ind.k-1)),which(all_lables == (ind.k-1))]),k = 1)
				homogeneity[ind.alpha_divisor] = homogeneity[ind.alpha_divisor] + sum(fulldistanceMatrix[med$medoids,which(all_lables == (ind.k-1))])
			}
		}
	
		alpha_divisor = alpha_divisor[which.min(homogeneity)]
			
		alpha = pmin(tfulldistanceMatrix/alpha_divisor,1)

		fulldistanceMatrix_indv_alpha_classif = alpha*tfulldistanceMatrix + (1-alpha)*pfulldistanceMatrix
		distanceMatrix_indv_alpha_classif = alpha*tdistanceMatrix + (1-alpha)*pdistanceMatrix
		##############################################################################################
		##############################################################################################

		##############################################################################################
		##############################################################################################
		## alpha homogeneity individually - choose alpha based on D(A,B)/alpha_divisor - CLASSIFICATION
		##############################################################################################
		##############################################################################################
		alphas = seq(0.01, 0.99, length = 20)
		homogeneity = 0
		for (ind.alphas in 1:length(alphas))
		{
			alpha = alphas[ind.alphas]
			fulldistanceMatrix = alpha*tfulldistanceMatrix + (1-alpha)*pfulldistanceMatrix
			distanceMatrix = alpha*tdistanceMatrix + (1-alpha)*pdistanceMatrix
			#### compute homogeneity of groups of classified objects
			homogeneity[ind.alphas] = 0
			for (ind.k in 1:K)
			{
				## find the medoids of the groups from the classification lables
				med = pam(as.dist(distanceMatrix[which(all_lables == (ind.k-1)),which(all_lables == (ind.k-1))]),k = 1)
				homogeneity[ind.alphas] = homogeneity[ind.alphas] + sum(fulldistanceMatrix[med$medoids,which(all_lables == (ind.k-1))])				
			}
		}
	
		alpha = alphas[which.min(homogeneity)]
			
		fulldistanceMatrix_alpha_fixed_classif = alpha*tfulldistanceMatrix + (1-alpha)*pfulldistanceMatrix
		distanceMatrix_alpha_fixed_classif = alpha*tdistanceMatrix + (1-alpha)*pdistanceMatrix
		##############################################################################################
		##############################################################################################








		###########################################################################################
		## K-medoids
		###########################################################################################
	
		n.init_guess = 100
		for (index.init_guess in 1:n.init_guess) ## do 100 initializations for random initial guess
		{
			# Partitioning (clustering) of the data into k clusters “around medoids”, a more robust version of K-means.
			k_medoids = pam(x = as.dist(tdistanceMatrix),k = K, medoids = sample(1:n,K))
			rate_Kmedoids[1,1] = rate_Kmedoids[1,1] + (n - Mismatch(k_medoids$clustering, (as.numeric(all_lables)+1), K = K))/n
			rate_Kmedoids_randindex[1,1] = rate_Kmedoids_randindex[1,1] + rand.index(k_medoids$clustering, (as.numeric(all_lables)+1))
			rate_Kmedoids_mutualinfo[1,1] = rate_Kmedoids_mutualinfo[1,1] + mutinformation(k_medoids$clustering, (as.numeric(all_lables)+1))

			k_medoids = pam(x = as.dist(pdistanceMatrix),k = K, medoids = sample(1:n,K))
			rate_Kmedoids[2,1] = rate_Kmedoids[2,1] + (n - Mismatch(k_medoids$clustering, (as.numeric(all_lables)+1), K = K))/n
			rate_Kmedoids_randindex[2,1] = rate_Kmedoids_randindex[2,1] + rand.index(k_medoids$clustering, (as.numeric(all_lables)+1))
			rate_Kmedoids_mutualinfo[2,1] = rate_Kmedoids_mutualinfo[2,1] + mutinformation(k_medoids$clustering, (as.numeric(all_lables)+1))


			k_medoids = pam(as.dist(distanceMatrix_indv_alpha),k = K, medoids = sample(1:n,K))
			rate_Kmedoids[3,1] = rate_Kmedoids[3,1] + (n - Mismatch(k_medoids$clustering, (as.numeric(all_lables)+1), K = K))/n
			rate_Kmedoids_randindex[3,1] = rate_Kmedoids_randindex[3,1] + rand.index(k_medoids$clustering, (as.numeric(all_lables)+1))
			rate_Kmedoids_mutualinfo[3,1] = rate_Kmedoids_mutualinfo[3,1] + mutinformation(k_medoids$clustering, (as.numeric(all_lables)+1))


			k_medoids = pam(as.dist(distanceMatrix_alpha_fixed),k = K, medoids = sample(1:n,K))
			rate_Kmedoids[4,1] = rate_Kmedoids[4,1] + (n - Mismatch(k_medoids$clustering, (as.numeric(all_lables)+1), K = K))/n
			rate_Kmedoids_randindex[4,1] = rate_Kmedoids_randindex[4,1] + rand.index(k_medoids$clustering, (as.numeric(all_lables)+1))
			rate_Kmedoids_mutualinfo[4,1] = rate_Kmedoids_mutualinfo[4,1] + mutinformation(k_medoids$clustering, (as.numeric(all_lables)+1))
		}
		rate_Kmedoids[1,1] = rate_Kmedoids[1,1]/n.init_guess
		rate_Kmedoids[2,1] = rate_Kmedoids[2,1]/n.init_guess
		rate_Kmedoids[3,1] = rate_Kmedoids[3,1]/n.init_guess
		rate_Kmedoids[4,1] = rate_Kmedoids[4,1]/n.init_guess
		rate_Kmedoids_randindex[1,1] = rate_Kmedoids_randindex[1,1]/n.init_guess
		rate_Kmedoids_randindex[2,1] = rate_Kmedoids_randindex[2,1]/n.init_guess
		rate_Kmedoids_randindex[3,1] = rate_Kmedoids_randindex[3,1]/n.init_guess
		rate_Kmedoids_randindex[4,1] = rate_Kmedoids_randindex[4,1]/n.init_guess
		rate_Kmedoids_mutualinfo[1,1] = rate_Kmedoids_mutualinfo[1,1]/n.init_guess
		rate_Kmedoids_mutualinfo[2,1] = rate_Kmedoids_mutualinfo[2,1]/n.init_guess
		rate_Kmedoids_mutualinfo[3,1] = rate_Kmedoids_mutualinfo[3,1]/n.init_guess
		rate_Kmedoids_mutualinfo[4,1] = rate_Kmedoids_mutualinfo[4,1]/n.init_guess

			cat("\n\n K-medoids rate:", (rate_Kmedoids[,1]))
			cat("\n\n K-medoids rand.index rate:", (rate_Kmedoids_randindex[,1]))
			cat("\n\n K-medoids mutualinfo rate:", (rate_Kmedoids_mutualinfo[,1]))

	
	
	
	
	
	
	
	
		##########################################################################################################
		## Proposed KNN - train data = all data
		##########################################################################################################
	
		classify = KNN_dist_matrix(tfulldistanceMatrix, 1:n, 1:n, all_lables = all_lables, n.neighbors)
		rateKNNtrain[1,1] = sum(classify == all_lables)/n
		classify = KNN_dist_matrix(pfulldistanceMatrix, 1:n, 1:n, all_lables = all_lables, n.neighbors)
		rateKNNtrain[2,1] = sum(classify == all_lables)/n
		classify = KNN_dist_matrix(fulldistanceMatrix_indv_alpha_classif, 1:n, 1:n, all_lables = all_lables, n.neighbors)
		rateKNNtrain[3,1] = sum(classify == all_lables)/n
		classify = KNN_dist_matrix(fulldistanceMatrix_alpha_fixed_classif, 1:n, 1:n, all_lables = all_lables, n.neighbors)
		rateKNNtrain[4,1] = sum(classify == all_lables)/n
	
		##########################################################################################################
		## Proposed KNN - test data - 500 random splits
		##########################################################################################################
	
		n.random_splits = 100		
		for (index.split in 1:n.random_splits)
		{	
			train_sample_index = sample(1:n,n*0.70)
			test_sample_index = (1:n)[-train_sample_index]
			
			train_vlmcs = vlmcs[train_sample_index]
			test_vlmcs = vlmcs[test_sample_index]
			true_labels_train = do.call("rbind", lapply(train_vlmcs, "[[", 1))
			true_labels_test = do.call("rbind", lapply(test_vlmcs, "[[", 1))
	
			classify = KNN_dist_matrix(tfulldistanceMatrix, train_sample_index, test_sample_index, all_lables = all_lables, n.neighbors)
			rateKNN[1,1] = rateKNN[1,1] + sum(classify == true_labels_test)/length(test_sample_index)
			classify = KNN_dist_matrix(pfulldistanceMatrix, train_sample_index, test_sample_index, all_lables = all_lables, n.neighbors)
			rateKNN[2,1] = rateKNN[2,1] + sum(classify == true_labels_test)/length(test_sample_index)
			classify = KNN_dist_matrix(fulldistanceMatrix_indv_alpha_classif, train_sample_index, test_sample_index, all_lables = all_lables, n.neighbors)
			rateKNN[3,1] = rateKNN[3,1] + sum(classify == true_labels_test)/length(test_sample_index)
			classify = KNN_dist_matrix(fulldistanceMatrix_alpha_fixed_classif, train_sample_index, test_sample_index, all_lables = all_lables, n.neighbors)
			rateKNN[4,1] = rateKNN[4,1] + sum(classify == true_labels_test)/length(test_sample_index)
	
		}
		rateKNN[1,1] = rateKNN[1,1]/n.random_splits
		rateKNN[2,1] = rateKNN[2,1]/n.random_splits
		rateKNN[3,1] = rateKNN[3,1]/n.random_splits
		rateKNN[4,1] = rateKNN[4,1]/n.random_splits

	
	
	
	
	
	
	
	
	
	
	
		##########################################################################################################
		## Classic Methods
		##########################################################################################################


		data_set_all = matrix(0, n, n.time_series)
		for (index in 1:n)
			data_set_all[index,] = vlmcs[[index]][[3]]
		classic_kmeans = kmeans(data_set_all, centers = K)
		rate_Kmeans[1] = (n - Mismatch(classic_kmeans$cluster, (as.numeric(all_lables)+1), K = K))/n
		rate_Kmeans_rand_index[1] = rand.index(classic_kmeans$cluster, (as.numeric(all_lables)+1))
		rate_Kmeans_mutual_info[1] = mutinformation(classic_kmeans$cluster, (as.numeric(all_lables)+1))
			
		classify = knn(data_set_all, data_set_all, all_lables, n.neighbors)
		rateClassicKNN_train[1] = sum(classify == all_lables)/n

		n.random_splits = 100	
		rateClassicKNN[1] = 0	
		for (index.split in 1:n.random_splits)
		{
			train_sample_index = sample(1:n,n*0.70)
			test_sample_index = (1:n)[-train_sample_index]
			
			train_vlmcs = vlmcs[train_sample_index]
			test_vlmcs = vlmcs[test_sample_index]
			true_labels_train = do.call("rbind", lapply(train_vlmcs, "[[", 1))
			true_labels_test = do.call("rbind", lapply(test_vlmcs, "[[", 1))

			data_set_train = matrix(0, length(train_vlmcs), n.time_series)
			for (index in 1:length(train_vlmcs))
				data_set_train[index,] = train_vlmcs[[index]][[3]]
			data_set_test = matrix(0, length(test_vlmcs), n.time_series)
			for (index in 1:length(test_vlmcs))
				data_set_test[index,] = test_vlmcs[[index]][[3]]
				
			classify = knn(data_set_train, data_set_test, true_labels_train, n.neighbors)
			rateClassicKNN[1] = rateClassicKNN[1] + sum(classify == true_labels_test)/length(test_sample_index)
		}
		rateClassicKNN[1] = rateClassicKNN[1]/n.random_splits
	
	
	
	
			cat("\n\n Kmeans:", mean(rate_Kmeans[1]))
			cat("\n\n Kmeans rand.index:", mean(rate_Kmeans_rand_index[1]))
			cat("\n\n Kmeans mutual.info:", mean(rate_Kmeans_mutual_info[1]))
			cat("\n\n KNN train rate:", (rateKNNtrain[,1]))
			cat("\n\n KNN test rate:", (rateKNN[,1]))
			cat("\n\n KNN Classic train rate:", mean(rateClassicKNN_train[1]))
			cat("\n\n KNN Classic test rate:", mean(rateClassicKNN[1]))
		
			
	
	

	
		
	
	
	


	
	

	
	

	
	

	
	

	
	

	
	

	
	

	
	

	
	

	
	

