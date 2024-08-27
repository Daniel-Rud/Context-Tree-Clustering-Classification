# Introduction

In this document, we provide code for the simulations described in the
manuscript *Context Tree Clustering and Classification*. We will first
run through an walk-through version of a single simulation scenario. At
the end of the manuscript, we provide the code used to run all
simulation scenarios and the results of the classification and
clustering for each simulation scenario.

## Simulation Walkthrough

We first introduce the context tree structural dissimilarity measures.

If you would like to run the results chunk by chunk, please refer to the
*Simulation-Walkthrough.Rmd* file.

## Relevant libraries

    if (!require("data.tree")) install.packages("data.tree", dependencies = TRUE); library("data.tree")

    ## Loading required package: data.tree

    if (!require("VLMC")) install.packages("VLMC", dependencies = TRUE); library("VLMC")

    ## Loading required package: VLMC

    if (!require("stats")) install.packages("stats", dependencies = TRUE); library("stats")
    if (!require("gtools")) install.packages("gtools", dependencies = TRUE); library("gtools")

    ## Loading required package: gtools

    if (!require("cluster")) install.packages("cluster", dependencies = TRUE); library("cluster")

    ## Loading required package: cluster

    if (!require("factoextra")) install.packages("factoextra", dependencies = TRUE); library("factoextra")

    ## Loading required package: factoextra

    ## Loading required package: ggplot2

    ## Welcome! Want to learn more? See two factoextra-related books at https://goo.gl/ve3WBa

    if (!require("mclust")) install.packages("mclust", dependencies = TRUE); library("mclust")

    ## Loading required package: mclust

    ## Package 'mclust' version 6.1.1
    ## Type 'citation("mclust")' for citing this R package in publications.

    if (!require("infotheo")) install.packages("infotheo", dependencies = TRUE); library("infotheo")

    ## Loading required package: infotheo

    ## 
    ## Attaching package: 'infotheo'

    ## The following object is masked from 'package:VLMC':
    ## 
    ##     entropy

    if (!require("fossil")) install.packages("fossil", dependencies = TRUE); library("fossil")

    ## Loading required package: fossil

    ## Loading required package: sp

    ## Loading required package: maps

    ## 
    ## Attaching package: 'maps'

    ## The following object is masked from 'package:mclust':
    ## 
    ##     map

    ## The following object is masked from 'package:cluster':
    ## 
    ##     votes.repub

    ## Loading required package: shapefiles

    ## Loading required package: foreign

    ## 
    ## Attaching package: 'shapefiles'

    ## The following objects are masked from 'package:foreign':
    ## 
    ##     read.dbf, write.dbf

    if (!require("future")) install.packages("future", dependencies = TRUE); library("future")

    ## Loading required package: future

    if (!require("future.apply")) install.packages("future.apply", dependencies = TRUE); library("future.apply")

    ## Loading required package: future.apply

# Loading Neccessary Functions

We first need some relevant helper functions.

    # reverse a string
    strReverse <- function(x)
      sapply(lapply(strsplit(x, NULL), rev), paste, collapse="")


    # generate all state space paths from a VLMC tree object 
    vlmcTreePaths= function(vlmcObj)
    {
      allPaths=vector(mode="character")
      if(vlmcObj$size["ord.MC"] != 0)
      {
        # if the VLMC is NOT order 0
        dendrogram= as.dendrogram(vlmcObj)
        #plot(dendrogram)
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
      }
      return(allPaths)  
    }

Now, we define the structural distance tDistance, the probability based
dissimilarity measure pDistance, and the

    # structural distance measure 
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

    # distance based on the structure and the probabilities - 
    # need the data (textOutput1 and 2) to estimate the probabilities 
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

Some functions for K nearest neighbors (KNN).

    ## classical KNN 
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



    KNN_dist_matrix <- function(fulldistanceMatrix, train_sample_index, test_sample_index, all_labels, n.neighbors)
    {
      # Performs KNN classification given distance matrix
      
      result = 0
      for (i in 1:length(test_sample_index))
      {
        distances = fulldistanceMatrix[,test_sample_index[i]]
        dist_to_train_data = distances[train_sample_index]
        order_dist = order(dist_to_train_data)[1:n.neighbors]
        which_train = train_sample_index[order_dist]
        
        result[i] = as.numeric(majorityVote(all_labels[which_train])$majority)
      }
      
      ## returns a vector of the classification for each test data
      return(result)
    }

We now define a function to measure cluster accuracy.

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

## Simulation Study – Scenario 1

First, we define two generating VLMCs from two seperate state sequences.
vlmcA and vlmcB serve as the two *population* context trees from which
we will simulate state sequences from.

    dataA = c(1, 3, 3, 3, 0, 0, 2, 1, 3, 3, 3, 2, 1, 0, 0, 3, 0, 4, 3, 2, 1, 2, 1, 0, 0, 4, 2, 1, 0, 3, 3, 0, 3, 3, 3, 3, 0, 0, 2, 1, 3, 0, 0, 3, 3, 3, 2, 1, 0, 0, 3, 0, 3, 3, 0, 3, 0, 3, 3, 0, 3, 0, 3, 0, 3, 4, 3, 2, 1, 0, 0, 3, 2, 1, 2, 1, 0, 3, 3, 3, 3, 2, 1, 3, 3, 0, 2, 1, 4, 3, 0, 3, 3, 3, 3, 4, 3, 0, 3, 3)
    draw(vlmc(dataA))

    ## [x]-( 27 12 12 43 5| 99)
    ##  +--[2]-( 0 12 0 0 0| 12) <25.32>-T
    ##  '--[3]-( 14 0 6 20 2| 42) <5.52>-T

    vlmcA = vlmc(dataA)


    dataB = c(3, 1, 3, 1, 3, 3, 3, 0, 2, 1, 2, 1, 3, 3, 3, 3, 3, 4, 3, 0, 3, 3, 0, 3, 3, 3, 0, 0, 3, 0, 3, 0, 3, 2, 1, 0, 0, 3, 3, 3, 0, 2, 0, 1, 0, 2, 1, 0, 0, 3, 2, 0, 1, 2, 0, 1, 0, 4, 3, 0, 3, 0, 2, 1, 4, 3, 1, 1, 0, 3, 3, 3, 3, 2, 1, 3, 3, 3, 3, 3, 0, 3, 3, 2, 1, 3, 2, 1, 3, 0, 3, 0, 0, 3, 4, 3, 2, 1, 4, 3, 0, 3, 3, 0, 4, 3, 0, 0, 3, 0, 3, 3, 2, 1, 0, 3, 2, 1, 3, 3, 3, 2, 1, 0, 3, 3, 4, 2, 1, 3, 4, 3, 2, 1, 3, 2, 0, 1, 3, 2)
    draw(vlmc(dataB))

    ## [x]-( 32 22 19 58 8| 139)
    ##  +--[0]-( 5 4 4 17 2| 32) <1.03>
    ##  |   '--[2]-( 0 4 0 0 0| 4) <8.32>-T
    ##  +--[2]-( 4 14 0 0 0| 18) <22.15>-T
    ##  '--[4]-( 0 0 1 7 0| 8) <5.09>-T

    vlmcB = vlmc(dataB)

The simulations in Table 1 iterate over the length of the state sequence
sampled from each of the 40 vlmcA and 40 vlmcB context trees. We
showcase a single iteration, where the length of the state sequence is
500.

    n.time_series = 500 ## how many observations of the chain we are generating
    all_results = list() ## each element will hold all results of the simulation for a specific time series n

    n = 80 # total number of VLMC observations -- 40 + 40 
    n.neighbors = 7 # number of neighbors for classification 
    vlmcs = list()

    # for reproducability
    set.seed(2024)

    # number of unique classes K
    K = 2   

Now, we begin the simulations. First, we simulate the VLMC state
sequences from vlmcA and vlmcB.

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

    # shuffle the observations 
    vlmcs = vlmcs[sample(n)]

    # true labels of shuffled VLMCs
    all_labels = do.call("rbind", lapply(vlmcs, "[[", 1))

After simulating the state sequences, we can compute both the structural
tDistance and the probability based pDistance to generate distance
matrices that will be used by the K-medoids algorithm. The computations
for *D* and *Δ* are performed here.

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

### Clustering Specific Distances

Here, we compute the distance
*D*<sub>*α*<sub>|*χ*|</sub></sub><sup>\*</sup>. Note that the for loop
over alpha divisor allows the user to investigate a set of alpha divisor
values. However, only the alpha divisor |*χ*|<sup>2</sup>/4 is used.

    ##############################################################################################
    ##############################################################################################
    ## alpha homogeneity individually - choose alpha based on D(A,B)/alpha_divisor - CLUSTERING
    ##############################################################################################
    ##############################################################################################
    alpha_divisor = length(unique(dataA))^2/4 

    alpha = pmin(tfulldistanceMatrix/alpha_divisor,1)

    fulldistanceMatrix_indv_alpha = alpha*tfulldistanceMatrix + (1-alpha)*pfulldistanceMatrix
    distanceMatrix_indv_alpha = alpha*tdistanceMatrix + (1-alpha)*pdistanceMatrix

Here, the distance *D*<sub>*α*<sub>*W**C**S**S*</sub></sub><sup>\*</sup>
is computed.

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

### Classification Specific Distances

Computation of *D*<sub>*α*<sub>|*χ*|</sub></sub><sup>\*</sup> for
classification. Note that one could pre-specify a sequence of alpha
divisor values, where the alpha divisor that minimizes the sum of the
distances to the medoids would be selected. However, only the alpha
divisor |*χ*|<sup>2</sup>/4 is used.

    alpha_divisor = length(unique(dataA))^2/4 
    alpha = pmin(tfulldistanceMatrix/alpha_divisor,1)

    fulldistanceMatrix_indv_alpha_classif = alpha*tfulldistanceMatrix + (1-alpha)*pfulldistanceMatrix
    distanceMatrix_indv_alpha_classif = alpha*tdistanceMatrix + (1-alpha)*pdistanceMatrix

Computation of *D*<sub>*α*<sub>*m**e**d*</sub></sub><sup>\*</sup>. Note
that we will make this a function, since it will need to be recomputed
when iterating over train test splits (since it relies on tuning alpha
with respect to the known labels)

    generate_D_alpha_med = function(tfulldistanceMatrix, pfulldistanceMatrix, all_labels, K)
    {
      alphas = seq(0.01, 0.99, length = 20)
      homogeneity = 0
      
      
      # get lower tri distance matricies 
      tdistanceMatrix = tfulldistanceMatrix
      tdistanceMatrix[upper.tri(tdistanceMatrix)] = 0
      
      pdistanceMatrix = pfulldistanceMatrix
      pdistanceMatrix[upper.tri(pdistanceMatrix)] = 0
      
      for (ind.alphas in 1:length(alphas))
      {
        alpha = alphas[ind.alphas]
        fulldistanceMatrix = alpha*tfulldistanceMatrix + (1-alpha)*pfulldistanceMatrix
        distanceMatrix = alpha* tdistanceMatrix + (1-alpha)*pdistanceMatrix
        #### compute homogeneity of groups of classified objects
        homogeneity[ind.alphas] = 0
        for (ind.k in 1:K)
        {
          ## find the medoids of the groups from the classification labels
          med = pam(as.dist(distanceMatrix[which(all_labels == (ind.k-1)),which(all_labels == (ind.k-1))]),k = 1)
          homogeneity[ind.alphas] = homogeneity[ind.alphas] + sum(fulldistanceMatrix[med$medoids,which(all_labels == (ind.k-1))])               
        }
      }
      alpha = alphas[which.min(homogeneity)]
      
      return(alpha)
    }

    alpha = generate_D_alpha_med(tfulldistanceMatrix, pfulldistanceMatrix, all_labels, K = K)

    fulldistanceMatrix_alpha_fixed_classif = alpha*tfulldistanceMatrix + (1-alpha)*pfulldistanceMatrix
    distanceMatrix_alpha_fixed_classif = alpha*tdistanceMatrix + (1-alpha)*pdistanceMatrix

## Classification – corresponds to Section 3.1

### Leave One Out Cross Validation (LOOCV)

    classif_results = matrix(data = 0, nrow = 5, ncol = 2)
    colnames(classif_results) = c("LOOCV", "70-30 train test")
    rownames(classif_results) = c("KNN - tDistance",
                                  "KNN - pDistance",
                                  "KNN - D_alpha_|chi|*",
                                  "KNN - D_alpha_med*",
                                  "KNN - Classic"
    )

    # tDistance 
    classify = KNN_dist_matrix(tfulldistanceMatrix, 1:n, 1:n, all_labels = all_labels, n.neighbors)
    classif_results[1,1] = sum(classify == all_labels)/n

    # pDistance
    classify = KNN_dist_matrix(pfulldistanceMatrix, 1:n, 1:n, all_labels = all_labels, n.neighbors)
    classif_results[2,1] = sum(classify == all_labels)/n

    # D_alpha_|chi|*
    classify = KNN_dist_matrix(fulldistanceMatrix_indv_alpha_classif, 1:n, 1:n, all_labels = all_labels, n.neighbors)
    classif_results[3,1] = sum(classify == all_labels)/n

    # D_alpha_med*
    # need to compute alpha with respect to ONLY training data

    D_alpha_med_alphas = sapply(1:n, FUN = function(test_ind)
    {
      alpha = generate_D_alpha_med(tfulldistanceMatrix[-test_ind, -test_ind],pfulldistanceMatrix[-test_ind, -test_ind], all_labels[-test_ind], K = K)
    })

    classify = numeric(n)

    # get LOOCV classification results, with each alpha for D_alpha_med computed on current training data
    for(i in 1:n)
    {
      fulldistanceMatrix_alpha_fixed_classif_curr = D_alpha_med_alphas[i]*tfulldistanceMatrix + (1-D_alpha_med_alphas[i])*pfulldistanceMatrix
      
      classify[i] = KNN_dist_matrix(fulldistanceMatrix_alpha_fixed_classif_curr,
                                    train_sample_index = (1:n)[-i], 
                                    test_sample_index = i,  
                                    all_labels = all_labels, 
                                    n.neighbors = n.neighbors)
    }

    classif_results[4,1] = sum(classify == all_labels)/n

    # Classic KNN 

    data_set_all = t(sapply(vlmcs, "[[", 3))
    classify = knn(data_set_all, data_set_all, all_labels, n.neighbors)
    classif_results[5,1] = sum(classify == all_labels)/n

### 100 70-30 Train Test splits

    n.random_splits = 100       
    for (index.split in 1:n.random_splits)
    {   
      train_sample_index = sample(1:n,n*0.70)
      test_sample_index = (1:n)[-train_sample_index]
      
      train_vlmcs = vlmcs[train_sample_index]
      test_vlmcs = vlmcs[test_sample_index]
      true_labels_train = do.call("rbind", lapply(train_vlmcs, "[[", 1))
      true_labels_test = do.call("rbind", lapply(test_vlmcs, "[[", 1))
      
      # tDistance 
      classify = KNN_dist_matrix(tfulldistanceMatrix, train_sample_index, test_sample_index, all_labels = all_labels, n.neighbors)
      classif_results[1,2] = classif_results[1,2] + sum(classify == true_labels_test)/length(test_sample_index)
      
      # pDistance
      classify = KNN_dist_matrix(pfulldistanceMatrix, train_sample_index, test_sample_index, all_labels = all_labels, n.neighbors)
      classif_results[2,2] = classif_results[2,2] + sum(classify == true_labels_test)/length(test_sample_index)
      
      # D_alpha_|chi|*
      classify = KNN_dist_matrix(fulldistanceMatrix_indv_alpha_classif, train_sample_index, test_sample_index, all_labels = all_labels, n.neighbors)
      classif_results[3,2] = classif_results[3,2] + sum(classify == true_labels_test)/length(test_sample_index)
      
      # D_alpha_med*
      # need to compute alpha with respect to ONLY training data
      curr_alpha = generate_D_alpha_med(tfulldistanceMatrix[-test_sample_index, -test_sample_index], 
                                        pfulldistanceMatrix[-test_sample_index, -test_sample_index], 
                                        all_labels = all_labels[-test_sample_index], K = K)
      
      fulldistanceMatrix_alpha_fixed_classif_curr = curr_alpha*tfulldistanceMatrix + (1-curr_alpha)*pfulldistanceMatrix
      
      classify = KNN_dist_matrix(fulldistanceMatrix_alpha_fixed_classif_curr, train_sample_index, test_sample_index, all_labels = all_labels, n.neighbors)
      classif_results[4,2] = classif_results[4,2] + sum(classify == true_labels_test)/length(test_sample_index)
      
      # Classic KNN 
      
      # get training state sequences 
      
      data_set_train = t(sapply(train_vlmcs, "[[", 3))
      data_set_test = t(sapply(test_vlmcs, "[[", 3))
      
      classify = knn(data_set_train, data_set_test, true_labels_train, n.neighbors)
      classif_results[5,2] = classif_results[5,2] + sum(classify == true_labels_test)/length(test_sample_index)
    }


    # divide by simulation replicates for 100 70-30 train test splits 
    classif_results[,2] = classif_results[,2] / n.random_splits

    cat(paste0("Classification Results for n.time_series = ", n.time_series, ", n = ",n, ":\n"))

    ## Classification Results for n.time_series = 500, n = 80:

    classif_results

    ##                       LOOCV 70-30 train test
    ## KNN - tDistance      0.9875        0.9870833
    ## KNN - pDistance      0.9875        0.9879167
    ## KNN - D_alpha_|chi|* 1.0000        1.0000000
    ## KNN - D_alpha_med*   1.0000        0.9975000
    ## KNN - Classic        0.7000        0.5216667

## Clustering – corresponds to Section 3.2

    # number of simulation replicates.  Note that multiple simulation replicates are 
    # run to average over Kmedoid and K means cluster performance from random initialization 
    # of centroids 
    n.iter = 100

    # the following are initialized for storing results
    clustering_results = matrix(0, nrow = 5, ncol = 3)
    colnames(clustering_results) = c("rate", "rand.index_rate","mutual_info rate")
    rownames(clustering_results) = c("K-Medoids - tDistance", 
                                     "K-Medoids - pDistance",
                                     "K-Medoids - D_{alpha_{|chi|}}*", 
                                     "K-Medoids - D_{alpha_{WCSS}}*", 
                                     "Classical KNN")



    # to compare with K means clustering, we need to extract the state sequences.  
    # They will be treated as points with n.time_series dimensions in KNN. 
    state_sequences = t(sapply(vlmcs, "[[", 3))


    # Iterate over n.iter iterations -- accounting for randomness in initial centroid selection 
    for (iter in 1:n.iter) 
    {
      # K medoids using tDistance
      k_medoids = pam(x = as.dist(tdistanceMatrix),k = K, medoids = sample(1:n,K))
      clustering_results[1,1] = clustering_results[1,1] + (n - Mismatch(k_medoids$clustering, (as.numeric(all_labels)+1), K = K))/n
      clustering_results[1,2] = clustering_results[1,2] + rand.index(k_medoids$clustering, (as.numeric(all_labels)+1))
      clustering_results[1,3] = clustering_results[1,3] + mutinformation(k_medoids$clustering, (as.numeric(all_labels)+1))
      
      # K medioids using pDistance 
      k_medoids = pam(x = as.dist(pdistanceMatrix),k = K, medoids = sample(1:n,K))
      clustering_results[2,1] = clustering_results[2,1] + (n - Mismatch(k_medoids$clustering, (as.numeric(all_labels)+1), K = K))/n
      clustering_results[2,2] = clustering_results[2,2] + rand.index(k_medoids$clustering, (as.numeric(all_labels)+1))
      clustering_results[2,3] = clustering_results[2,3] + mutinformation(k_medoids$clustering, (as.numeric(all_labels)+1))
      
      # K medoids using D_{\alpha_{|\chi|}}^*
      k_medoids = pam(as.dist(distanceMatrix_indv_alpha),k = K, medoids = sample(1:n,K))
      clustering_results[3,1] = clustering_results[3,1] + (n - Mismatch(k_medoids$clustering, (as.numeric(all_labels)+1), K = K))/n
      clustering_results[3,2] = clustering_results[3,2] + rand.index(k_medoids$clustering, (as.numeric(all_labels)+1))
      clustering_results[3,3] = clustering_results[3,3] + mutinformation(k_medoids$clustering, (as.numeric(all_labels)+1))
      
      # K medoids using D_{\alpha_{WCSS}}^*
      k_medoids = pam(as.dist(distanceMatrix_alpha_fixed),k = K, medoids = sample(1:n,K))
      clustering_results[4,1] = clustering_results[4,1] + (n - Mismatch(k_medoids$clustering, (as.numeric(all_labels)+1), K = K))/n
      clustering_results[4,2] = clustering_results[4,2] + rand.index(k_medoids$clustering, (as.numeric(all_labels)+1))
      clustering_results[4,3] = clustering_results[4,3] + mutinformation(k_medoids$clustering, (as.numeric(all_labels)+1))
      
      # Classical K means 
      
      classic_kmeans = kmeans(state_sequences, centers = K)
      clustering_results[5,1] = clustering_results[5,1] + (n - Mismatch(classic_kmeans$cluster, (as.numeric(all_labels)+1), K = K))/n
      clustering_results[5,2] = clustering_results[5,2] + rand.index(classic_kmeans$cluster, (as.numeric(all_labels)+1))
      clustering_results[5,3] = clustering_results[5,3] + mutinformation(classic_kmeans$cluster, (as.numeric(all_labels)+1))
    }

    # divide by simulation replicates 
    clustering_results = clustering_results / n.iter

    cat(paste0("Clustering Results for n.time_series = ", n.time_series, ", n = ",n, ":\n"))

    ## Clustering Results for n.time_series = 500, n = 80:

    clustering_results 

    ##                                    rate rand.index_rate mutual_info rate
    ## K-Medoids - tDistance          0.990000       0.9807595       0.65963721
    ## K-Medoids - pDistance          1.000000       1.0000000       0.69314718
    ## K-Medoids - D_{alpha_{|chi|}}* 1.000000       1.0000000       0.69314718
    ## K-Medoids - D_{alpha_{WCSS}}*  1.000000       1.0000000       0.69314718
    ## Classical KNN                  0.527625       0.4968576       0.00323722

# Generation of all Tables and scenarios

We will now perform all the simulations for Scenarios 1, 2, and 3 in the
main manuscript. First, we define some general functions for the
clustering and classification simulations.

    run_VLMC_simulation = function(vlmcA, vlmcB,
                                   vlmcC = NULL,
                                   n.time_series = 100, n.neighbors = 7, 
                                   K = 2, n = 80, 
                                   n.random_splits = 100, 
                                   n.iter = 100)
    {
      # n.random_splits -- number of train test splits to assess 
      # n.iter -- number of times to repeat K medoid clustering
      ##############################################################################################
      # generate VLMCs 
      ##############################################################################################
      
      vlmcs = vector(mode = "list", length = n)
      
      # if only two generating VLMCs -- Scenario 1 and 2 
      if(is.null(vlmcC))
      {
        
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
      }else # if three generating VLMCs -- Scenario 3
      {
         for (index in 1:(n/3))
        {
          simulated_data = simulate(vlmcA, n.time_series)
          vlmcs[[index]] = list(0,vlmc(simulated_data, threshold.gen=10, alpha=0.05), as.numeric(simulated_data))
        }
        for (index in (n/3+1):(2*n/3))
        {
          simulated_data = simulate(vlmcB, n.time_series)
          vlmcs[[index]] = list(1,vlmc(simulated_data, threshold.gen=10, alpha=0.05), as.numeric(simulated_data))
        }   
          for (index in (2*n/3 + 1):n)
        {
          simulated_data = simulate(vlmcC, n.time_series)
          vlmcs[[index]] = list(2,vlmc(simulated_data, threshold.gen=10, alpha=0.05), as.numeric(simulated_data))
        }   
        
      }
      
      # shuffle the observations 
      vlmcs = vlmcs[sample(n)]
      
      # true labels of shuffled VLMCs
      all_labels = do.call("rbind", lapply(vlmcs, "[[", 1))
      
      ##############################################################################################
      # compute distances
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
      # pDistance is the distance based on the probabilities of the reduced subgraph
      pdistanceMatrix = matrix(0, n, n)
      for(index in 2:n) #first entry a_11 = 0
      {
        for(j in 1:index)
          pdistanceMatrix[index,j]= pDistance(vlmcs[[index]][[2]], vlmcs[[j]][[2]], vlmcs[[index]][[3]], vlmcs[[j]][[3]])
        
      }
      pdistanceMatrix[is.na(pdistanceMatrix)] = 0
      pfulldistanceMatrix = pdistanceMatrix + t(pdistanceMatrix)
      
      ##############################################################################################
      # generate D_alpha_|chi|* -- Classification
      alpha_divisor = length(c(0,1,2,3,4))^2/4 
      alpha = pmin(tfulldistanceMatrix/alpha_divisor,1)
      
      fulldistanceMatrix_indv_alpha_classif = alpha*tfulldistanceMatrix + (1-alpha)*pfulldistanceMatrix
      distanceMatrix_indv_alpha_classif = alpha*tdistanceMatrix + (1-alpha)*pdistanceMatrix
      
      ##############################################################################################
      # D_alpha_med* -- generated on training loop since hyperparameter is tuned  -- Classification
      ##############################################################################################
      # Clustering D_alpha_|chi|* -- Clustering
      alpha_divisor = length(c(0,1,2,3,4))^2/4 
      
      alpha = pmin(tfulldistanceMatrix/alpha_divisor,1)
      
      fulldistanceMatrix_indv_alpha = alpha*tfulldistanceMatrix + (1-alpha)*pfulldistanceMatrix
      distanceMatrix_indv_alpha = alpha*tdistanceMatrix + (1-alpha)*pdistanceMatrix
      ##############################################################################################
      # D_alpha_WCSS -- Clustering 
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
      # CLASSIFICATION
      ##############################################################################################
      ##############################################################################################
      
      ##############################################################################################
      # LOOCV
      ##############################################################################################
      classif_results = matrix(data = 0, nrow = 5, ncol = 2)
      colnames(classif_results) = c("LOOCV", "70-30 train test")
      rownames(classif_results) = c("KNN - tDistance",
                                    "KNN - pDistance",
                                    "KNN - D_alpha_|chi|*",
                                    "KNN - D_alpha_med*",
                                    "KNN - Classic"
      )
      
      # tDistance 
      classify = KNN_dist_matrix(tfulldistanceMatrix, 1:n, 1:n, all_labels = all_labels, n.neighbors)
      classif_results[1,1] = sum(classify == all_labels)/n
      
      # pDistance
      classify = KNN_dist_matrix(pfulldistanceMatrix, 1:n, 1:n, all_labels = all_labels, n.neighbors)
      classif_results[2,1] = sum(classify == all_labels)/n
      
      
      # D_alpha_|chi|*
      classify = KNN_dist_matrix(fulldistanceMatrix_indv_alpha_classif, 1:n, 1:n, all_labels = all_labels, n.neighbors)
      classif_results[3,1] = sum(classify == all_labels)/n
      
      # D_alpha_med*
      # need to compute alpha with respect to ONLY training data
      D_alpha_med_alphas = sapply(1:n, FUN = function(test_ind)
      {
        alpha = generate_D_alpha_med(tfulldistanceMatrix[-test_ind, -test_ind],pfulldistanceMatrix[-test_ind, -test_ind], all_labels[-test_ind], K = K)
      })
      
      classify = numeric(n)
      
      # get LOOCV classification results, with each alpha for D_alpha_med computed on current training data
      for(i in 1:n)
      {
        fulldistanceMatrix_alpha_fixed_classif_curr = D_alpha_med_alphas[i]*tfulldistanceMatrix + (1-D_alpha_med_alphas[i])*pfulldistanceMatrix
        
        classify[i] = KNN_dist_matrix(fulldistanceMatrix_alpha_fixed_classif_curr,
                                      train_sample_index = (1:n)[-i], 
                                      test_sample_index = i,  
                                      all_labels = all_labels, 
                                      n.neighbors = n.neighbors)
      }
      
      classif_results[4,1] = sum(classify == all_labels)/n
      
      # Classic KNN 
      
      data_set_all = t(sapply(vlmcs, "[[", 3))
      classify = knn(data_set_all, data_set_all, all_labels, n.neighbors)
      classif_results[5,1] = sum(classify == all_labels)/n
      
      ##############################################################################################
      # Train/Test Splits
      ##############################################################################################
      
      for (index.split in 1:n.random_splits)
      { 
        train_sample_index = sample(1:n,n*0.70)
        test_sample_index = (1:n)[-train_sample_index]
        
        train_vlmcs = vlmcs[train_sample_index]
        test_vlmcs = vlmcs[test_sample_index]
        true_labels_train = do.call("rbind", lapply(train_vlmcs, "[[", 1))
        true_labels_test = do.call("rbind", lapply(test_vlmcs, "[[", 1))
        
        # tDistance 
        classify = KNN_dist_matrix(tfulldistanceMatrix, train_sample_index, test_sample_index, all_labels = all_labels, n.neighbors)
        classif_results[1,2] = classif_results[1,2] + sum(classify == true_labels_test)/length(test_sample_index)
        
        # pDistance
        classify = KNN_dist_matrix(pfulldistanceMatrix, train_sample_index, test_sample_index, all_labels = all_labels, n.neighbors)
        classif_results[2,2] = classif_results[2,2] + sum(classify == true_labels_test)/length(test_sample_index)
        
        # D_alpha_|chi|*
        classify = KNN_dist_matrix(fulldistanceMatrix_indv_alpha_classif, train_sample_index, test_sample_index, all_labels = all_labels, n.neighbors)
        classif_results[3,2] = classif_results[3,2] + sum(classify == true_labels_test)/length(test_sample_index)
        
        # D_alpha_med*
        # need to compute alpha with respect to ONLY training data
        curr_alpha = generate_D_alpha_med(tfulldistanceMatrix[-test_sample_index, -test_sample_index], 
                                          pfulldistanceMatrix[-test_sample_index, -test_sample_index], 
                                          all_labels = all_labels[-test_sample_index], K = K)
        
        fulldistanceMatrix_alpha_fixed_classif_curr = curr_alpha*tfulldistanceMatrix + (1-curr_alpha)*pfulldistanceMatrix
        
        classify = KNN_dist_matrix(fulldistanceMatrix_alpha_fixed_classif_curr, train_sample_index, test_sample_index, all_labels = all_labels, n.neighbors)
        classif_results[4,2] = classif_results[4,2] + sum(classify == true_labels_test)/length(test_sample_index)
        
        # Classic KNN 
        
        # get training state sequences 
        
        data_set_train = t(sapply(train_vlmcs, "[[", 3))
        data_set_test = t(sapply(test_vlmcs, "[[", 3))
        
        classify = knn(data_set_train, data_set_test, true_labels_train, n.neighbors)
        classif_results[5,2] = classif_results[5,2] + sum(classify == true_labels_test)/length(test_sample_index)
      }
      
      
      # divide by simulation replicates for 100 70-30 train test splits 
      classif_results[,2] = classif_results[,2] / n.random_splits
      
      # reshape -- add names
      classif_results = c(classif_results[,1], classif_results[,2])
      names(classif_results)[1:5] = paste0("LOOCV- ", names(classif_results)[1:5])
      names(classif_results)[6:10] = paste0("Train/Test- ", names(classif_results)[1:5])
      
      
      ##############################################################################################
      ##############################################################################################
      # CLUSTERING
      ##############################################################################################
      ##############################################################################################
      
      # the following are initialized for storing results
      clustering_results = matrix(0, nrow = 5, ncol = 3)
      colnames(clustering_results) = c("rate", "rand.index_rate","mutual_info rate")
      rownames(clustering_results) = c("K-Medoids - tDistance", 
                                       "K-Medoids - pDistance",
                                       "K-Medoids - D_{alpha_{|chi|}}*", 
                                       "K-Medoids - D_{alpha_{WCSS}}*", 
                                       "Classical KNN")
      
      
      
      # to compare with K means clustering, we need to extract the state sequences.  
      # They will be treated as points with n.time_series dimensions in KNN. 
      state_sequences = t(sapply(vlmcs, "[[", 3))
      
      # Iterate over n.iter iterations -- accounting for randomness in initial centroid selection 
      for (iter in 1:n.iter) 
      {
        # K medoids using tDistance
        k_medoids = pam(x = as.dist(tdistanceMatrix),k = K, medoids = sample(1:n,K))
        clustering_results[1,1] = clustering_results[1,1] + (n - Mismatch(k_medoids$clustering, (as.numeric(all_labels)+1), K = K))/n
        clustering_results[1,2] = clustering_results[1,2] + rand.index(k_medoids$clustering, (as.numeric(all_labels)+1))
        clustering_results[1,3] = clustering_results[1,3] + mutinformation(k_medoids$clustering, (as.numeric(all_labels)+1))
        
        # K medioids using pDistance 
        k_medoids = pam(x = as.dist(pdistanceMatrix),k = K, medoids = sample(1:n,K))
        clustering_results[2,1] = clustering_results[2,1] + (n - Mismatch(k_medoids$clustering, (as.numeric(all_labels)+1), K = K))/n
        clustering_results[2,2] = clustering_results[2,2] + rand.index(k_medoids$clustering, (as.numeric(all_labels)+1))
        clustering_results[2,3] = clustering_results[2,3] + mutinformation(k_medoids$clustering, (as.numeric(all_labels)+1))
        
        # K medoids using D_{\alpha_{|\chi|}}^*
        k_medoids = pam(as.dist(distanceMatrix_indv_alpha),k = K, medoids = sample(1:n,K))
        clustering_results[3,1] = clustering_results[3,1] + (n - Mismatch(k_medoids$clustering, (as.numeric(all_labels)+1), K = K))/n
        clustering_results[3,2] = clustering_results[3,2] + rand.index(k_medoids$clustering, (as.numeric(all_labels)+1))
        clustering_results[3,3] = clustering_results[3,3] + mutinformation(k_medoids$clustering, (as.numeric(all_labels)+1))
        
        # K medoids using D_{\alpha_{WCSS}}^*
        k_medoids = pam(as.dist(distanceMatrix_alpha_fixed),k = K, medoids = sample(1:n,K))
        clustering_results[4,1] = clustering_results[4,1] + (n - Mismatch(k_medoids$clustering, (as.numeric(all_labels)+1), K = K))/n
        clustering_results[4,2] = clustering_results[4,2] + rand.index(k_medoids$clustering, (as.numeric(all_labels)+1))
        clustering_results[4,3] = clustering_results[4,3] + mutinformation(k_medoids$clustering, (as.numeric(all_labels)+1))
        
        # Classical K means 
        
        classic_kmeans = kmeans(state_sequences, centers = K)
        clustering_results[5,1] = clustering_results[5,1] + (n - Mismatch(classic_kmeans$cluster, (as.numeric(all_labels)+1), K = K))/n
        clustering_results[5,2] = clustering_results[5,2] + rand.index(classic_kmeans$cluster, (as.numeric(all_labels)+1))
        clustering_results[5,3] = clustering_results[5,3] + mutinformation(classic_kmeans$cluster, (as.numeric(all_labels)+1))
      }
      
      # divide by simulation replicates 
      clustering_results = clustering_results / n.iter
      
      return_list = list(classif_results = classif_results, 
                         clustering_results = clustering_results[,1])
      
      return(return_list)
    }

## Scenario 1

    dataA = c(1, 3, 3, 3, 0, 0, 2, 1, 3, 3, 3, 2, 1, 0, 0, 3, 0, 4, 3, 2, 1, 2, 1, 0, 0, 4, 2, 1, 0, 3, 3, 0, 3, 3, 3, 3, 0, 0, 2, 1, 3, 0, 0, 3, 3, 3, 2, 1, 0, 0, 3, 0, 3, 3, 0, 3, 0, 3, 3, 0, 3, 0, 3, 0, 3, 4, 3, 2, 1, 0, 0, 3, 2, 1, 2, 1, 0, 3, 3, 3, 3, 2, 1, 3, 3, 0, 2, 1, 4, 3, 0, 3, 3, 3, 3, 4, 3, 0, 3, 3)

    dataB = c(3, 1, 3, 1, 3, 3, 3, 0, 2, 1, 2, 1, 3, 3, 3, 3, 3, 4, 3, 0, 3, 3, 0, 3, 3, 3, 0, 0, 3, 0, 3, 0, 3, 2, 1, 0, 0, 3, 3, 3, 0, 2, 0, 1, 0, 2, 1, 0, 0, 3, 2, 0, 1, 2, 0, 1, 0, 4, 3, 0, 3, 0, 2, 1, 4, 3, 1, 1, 0, 3, 3, 3, 3, 2, 1, 3, 3, 3, 3, 3, 0, 3, 3, 2, 1, 3, 2, 1, 3, 0, 3, 0, 0, 3, 4, 3, 2, 1, 4, 3, 0, 3, 3, 0, 4, 3, 0, 0, 3, 0, 3, 3, 2, 1, 0, 3, 2, 1, 3, 3, 3, 2, 1, 0, 3, 3, 4, 2, 1, 3, 4, 3, 2, 1, 3, 2, 0, 1, 3, 2)

    # Create pop VLMCs
    vlmcA = vlmc(dataA)
    vlmcB = vlmc(dataB)

    draw(vlmcA)
    draw(vlmcB)

    n.time_series = c(50, 100,500,1000, 2000) ## how many observations of the chain we are generating
    scenario_1_classif_results = matrix(data = 0, nrow = length(n.time_series), 
                                        ncol = 10)
    scenario_1_clustering_results = matrix(data = 0, nrow = length(n.time_series), 
                                        ncol = 5)

    rownames(scenario_1_classif_results) = paste0("T = ", n.time_series)
    rownames(scenario_1_clustering_results) = paste0("T = ", n.time_series)


    n.simulations = 500
    n = 80
    n.neighbors = 7
    K = 2

    future_seeds = c(8432, 4901,  219, 6553,  138)


    for(i in 1:length(n.time_series))
    {
      
      plan(multisession)
      current_sim_res = future_lapply(1:n.simulations, FUN = function(x)
      {
        return(run_VLMC_simulation(vlmcA = vlmcA, 
                                   vlmcB = vlmcB, 
                                   n.time_series = n.time_series[i], 
                                   n.neighbors = n.neighbors, 
                                   K = K,
                                   n = n,
                                   n.random_splits = 100, 
                                   n.iter = 100))
      },future.seed = future_seeds[i])
      plan(sequential)
      
      scenario_1_classif_results[i, ] = rowMeans(sapply(current_sim_res, "[[", 1))
      scenario_1_clustering_results[i, ] = rowMeans(sapply(current_sim_res, "[[", 2))
      
      if(i == 1)
      {
        colnames(scenario_1_classif_results) = names(rowMeans(sapply(current_sim_res, "[[", 1)))
        colnames(scenario_1_clustering_results) = names(rowMeans(sapply(current_sim_res, "[[", 2)))
      }
      
      cat(paste("Finished sim", i))
      
    }

    scenario_1_classif_results[,1:5]
    scenario_1_classif_results[,6:10]

    scenario_1_clustering_results


    # SAVE Results 
    Scenario_1_tables = list(scenario_1_classif_results = scenario_1_classif_results, 
                             scenario_1_clustering_results = scenario_1_clustering_results)

    saveRDS(Scenario_1_tables, "Scenario_1_tables.RDS")

# Scenario 2

    read_time_series = function(file)
    {
      suppressWarnings(data <- readLines(file))
      data = paste(data, collapse = "")
      data = unlist(strsplit(data, split = ","))
      data = as.numeric(data)
      return(data)
    }

    dataA = read_time_series("Scenario2_dataA.txt")
    dataB = read_time_series("Scenario2_dataB.txt")


    # Create pop VLMCs
    vlmcA = vlmc(dataA)
    vlmcB = vlmc(dataB)

    draw(vlmcA)
    draw(vlmcB)

    n.time_series = c(50, 100,500,1000, 2000) ## how many observations of the chain we are generating
    scenario_2_classif_results = matrix(data = 0, nrow = length(n.time_series), 
                                        ncol = 10)
    scenario_2_clustering_results = matrix(data = 0, nrow = length(n.time_series), 
                                        ncol = 5)
                                        
    rownames(scenario_2_classif_results) = paste0("T = ", n.time_series)
    rownames(scenario_2_clustering_results) = paste0("T = ", n.time_series)

    n.simulations = 500
    n = 80
    n.neighbors = 7
    K = 2

    future_seeds = c(8432, 4901,  219, 6553,  138)


    for(i in 1:length(n.time_series))
    {
      
      plan(multisession)
      current_sim_res = future_lapply(1:n.simulations, FUN = function(x)
      {
        return(run_VLMC_simulation(vlmcA = vlmcA, 
                                   vlmcB = vlmcB, 
                                   n.time_series = n.time_series[i], 
                                   n.neighbors = n.neighbors, 
                                   K = K,
                                   n = n,
                                   n.random_splits = 100, 
                                   n.iter = 100))
      },future.seed = future_seeds[i])
      plan(sequential)
      
      scenario_2_classif_results[i, ] = rowMeans(sapply(current_sim_res, "[[", 1))
      scenario_2_clustering_results[i, ] = rowMeans(sapply(current_sim_res, "[[", 2))
      
      if(i == 1)
      {
        colnames(scenario_2_classif_results) = names(rowMeans(sapply(current_sim_res, "[[", 1)))
        colnames(scenario_2_clustering_results) = names(rowMeans(sapply(current_sim_res, "[[", 2)))
      }
      
      cat(paste("Finished sim", i))
      
    }

    scenario_2_classif_results[,1:5]
    scenario_2_classif_results[,6:10]

    scenario_2_clustering_results

    # SAVE Results 
    Scenario_2_tables = list(scenario_2_classif_results = scenario_2_classif_results, 
                             scenario_2_clustering_results = scenario_2_clustering_results)

    saveRDS(Scenario_2_tables, "Scenario_2_tables.RDS")

# Scenario 3

    dataA = read_time_series("Scenario3_dataA.txt")
    dataB = read_time_series("Scenario3_dataB.txt")
    dataC = read_time_series("Scenario3_dataC.txt")


    # Create pop VLMCs
    vlmcA = vlmc(dataA)
    vlmcB = vlmc(dataB)
    # changed threshold gen solely to influence construction of population vlmc 
    vlmcC = vlmc(dataC, threshold.gen = 100)

    draw(vlmcA)
    draw(vlmcB)
    draw(vlmcC)


    n.time_series = c(50, 100,500,1000, 2000) ## how many observations of the chain we are generating
    scenario_3_classif_results = matrix(data = 0, nrow = length(n.time_series), 
                                        ncol = 10)
    scenario_3_clustering_results = matrix(data = 0, nrow = length(n.time_series), 
                                        ncol = 5)
    rownames(scenario_3_classif_results) = paste0("T = ", n.time_series)
    rownames(scenario_3_clustering_results) = paste0("T = ", n.time_series)
    n.simulations = 500
    n = 120
    n.neighbors = 7
    K = 3
    future_seeds = c(8432, 4901,  219, 6553,  138)


    for(i in 1:length(n.time_series))
    {
      
      plan(multisession)
      current_sim_res = future_lapply(1:n.simulations, FUN = function(x)
      {
        return(run_VLMC_simulation(vlmcA = vlmcA, 
                                   vlmcB = vlmcB,
                                   vlmcC = vlmcC, 
                                   n.time_series = n.time_series[i], 
                                   n.neighbors = n.neighbors, 
                                   K = K,
                                   n = n,
                                   n.random_splits = 100, 
                                   n.iter = 100))
      },future.seed = future_seeds[i])
      plan(sequential)
      
      scenario_3_classif_results[i, ] = rowMeans(sapply(current_sim_res, "[[", 1))
      scenario_3_clustering_results[i, ] = rowMeans(sapply(current_sim_res, "[[", 2))
      
      if(i == 1)
      {
        colnames(scenario_3_classif_results) = names(rowMeans(sapply(current_sim_res, "[[", 1)))
        colnames(scenario_3_clustering_results) = names(rowMeans(sapply(current_sim_res, "[[", 2)))
      }
      
      cat(paste("Finished sim", i))
      
    }

    scenario_3_classif_results

    scenario_3_clustering_results

    # SAVE Results 
    Scenario_3_tables = list(scenario_3_classif_results = scenario_3_classif_results, 
                             scenario_3_clustering_results = scenario_3_clustering_results)

    saveRDS(Scenario_3_tables, "Scenario_3_tables.RDS")

# Load Results

    Scenario_1_tables <- readRDS("~/Desktop/VLMC_files/Github-Files/Scenario_1_tables.RDS")
    Scenario_2_tables <- readRDS("~/Desktop/VLMC_files/Github-Files/Scenario_2_tables.RDS")
    Scenario_3_tables <- readRDS("~/Desktop/VLMC_files/Github-Files/Scenario_3_tables.RDS")

    Scenario_1_tables

    ## $scenario_1_classif_results
    ##          LOOCV- KNN - tDistance LOOCV- KNN - pDistance
    ## T = 50                 0.538925               0.499425
    ## T = 100                0.798400               0.555900
    ## T = 500                0.987600               0.998400
    ## T = 1000               0.999025               1.000000
    ## T = 2000               0.999100               1.000000
    ##          LOOCV- KNN - D_alpha_|chi|* LOOCV- KNN - D_alpha_med*
    ## T = 50                      0.543100                  0.489775
    ## T = 100                     0.944425                  0.654025
    ## T = 500                     0.999650                  0.997200
    ## T = 1000                    1.000000                  0.999925
    ## T = 2000                    1.000000                  1.000000
    ##          LOOCV- KNN - Classic Train/Test- LOOCV- KNN - tDistance
    ## T = 50               0.666875                          0.4933250
    ## T = 100              0.669725                          0.7626375
    ## T = 500              0.661425                          0.9826475
    ## T = 1000             0.644800                          0.9979333
    ## T = 2000             0.603225                          0.9985942
    ##          Train/Test- LOOCV- KNN - pDistance
    ## T = 50                            0.4849525
    ## T = 100                           0.5554058
    ## T = 500                           0.9963983
    ## T = 1000                          0.9999908
    ## T = 2000                          1.0000000
    ##          Train/Test- LOOCV- KNN - D_alpha_|chi|*
    ## T = 50                                 0.4921633
    ## T = 100                                0.9224750
    ## T = 500                                0.9990025
    ## T = 1000                               1.0000000
    ## T = 2000                               0.9999833
    ##          Train/Test- LOOCV- KNN - D_alpha_med* Train/Test- LOOCV- KNN - Classic
    ## T = 50                               0.4874417                        0.5027650
    ## T = 100                              0.6789992                        0.5059075
    ## T = 500                              0.9964650                        0.5237833
    ## T = 1000                             0.9999483                        0.5312508
    ## T = 2000                             0.9999983                        0.5300342
    ## 
    ## $scenario_1_clustering_results
    ##          K-Medoids - tDistance K-Medoids - pDistance
    ## T = 50               0.5369530             0.5157363
    ## T = 100              0.7815798             0.5170448
    ## T = 500              0.9257895             0.9867865
    ## T = 1000             0.9968767             0.9977690
    ## T = 2000             0.9971320             0.9997500
    ##          K-Medoids - D_{alpha_{|chi|}}* K-Medoids - D_{alpha_{WCSS}}*
    ## T = 50                        0.5227022                     0.5188152
    ## T = 100                       0.8981570                     0.6990467
    ## T = 500                       0.9965750                     0.9897000
    ## T = 1000                      0.9999750                     0.9999750
    ## T = 2000                      1.0000000                     0.9998250
    ##          Classical KNN
    ## T = 50       0.5437675
    ## T = 100      0.5452742
    ## T = 500      0.5456190
    ## T = 1000     0.5482075
    ## T = 2000     0.5521608

    Scenario_2_tables

    ## $scenario_2_classif_results
    ##          LOOCV- KNN - tDistance LOOCV- KNN - pDistance
    ## T = 50                 0.551300                0.49715
    ## T = 100                0.586350                0.71785
    ## T = 500                0.907925                0.74620
    ## T = 1000               0.925725                0.99875
    ## T = 2000               0.979325                1.00000
    ##          LOOCV- KNN - D_alpha_|chi|* LOOCV- KNN - D_alpha_med*
    ## T = 50                      0.556150                  0.486400
    ## T = 100                     0.790325                  0.656325
    ## T = 500                     0.922450                  0.902825
    ## T = 1000                    0.999825                  0.958675
    ## T = 2000                    1.000000                  0.996750
    ##          LOOCV- KNN - Classic Train/Test- LOOCV- KNN - tDistance
    ## T = 50               0.667500                          0.5040525
    ## T = 100              0.667350                          0.5268083
    ## T = 500              0.657025                          0.8765383
    ## T = 1000             0.641950                          0.9078342
    ## T = 2000             0.610150                          0.9461800
    ##          Train/Test- LOOCV- KNN - pDistance
    ## T = 50                            0.4819842
    ## T = 100                           0.6313500
    ## T = 500                           0.7039717
    ## T = 1000                          0.9981617
    ## T = 2000                          1.0000000
    ##          Train/Test- LOOCV- KNN - D_alpha_|chi|*
    ## T = 50                                 0.4954467
    ## T = 100                                0.7041350
    ## T = 500                                0.8819292
    ## T = 1000                               0.9995875
    ## T = 2000                               0.9999975
    ##          Train/Test- LOOCV- KNN - D_alpha_med* Train/Test- LOOCV- KNN - Classic
    ## T = 50                               0.4857958                        0.5057583
    ## T = 100                              0.6594367                        0.5083600
    ## T = 500                              0.9036458                        0.5162383
    ## T = 1000                             0.9623233                        0.5216092
    ## T = 2000                             0.9937283                        0.5235550
    ## 
    ## $scenario_2_clustering_results
    ##          K-Medoids - tDistance K-Medoids - pDistance
    ## T = 50               0.5372260             0.5176978
    ## T = 100              0.5581182             0.5208182
    ## T = 500              0.8427905             0.6000762
    ## T = 1000             0.8228310             0.9783508
    ## T = 2000             0.7702577             0.9998250
    ##          K-Medoids - D_{alpha_{|chi|}}* K-Medoids - D_{alpha_{WCSS}}*
    ## T = 50                        0.5266687                     0.5134725
    ## T = 100                       0.5999775                     0.5212412
    ## T = 500                       0.7655147                     0.8513825
    ## T = 1000                      0.9900043                     0.9549853
    ## T = 2000                      0.9998750                     0.9999000
    ##          Classical KNN
    ## T = 50       0.5456635
    ## T = 100      0.5454192
    ## T = 500      0.5471882
    ## T = 1000     0.5476665
    ## T = 2000     0.5502243

    Scenario_3_tables

    ## $scenario_3_classif_results
    ##          LOOCV- KNN - tDistance LOOCV- KNN - pDistance
    ## T = 50                0.5548833                0.28595
    ## T = 100               0.7932667                0.21125
    ## T = 500               0.7272500                0.97330
    ## T = 1000              0.7560000                1.00000
    ## T = 2000              0.7906667                1.00000
    ##          LOOCV- KNN - D_alpha_|chi|* LOOCV- KNN - D_alpha_med*
    ## T = 50                     0.4744167                 0.3236333
    ## T = 100                    0.8288333                 0.3519333
    ## T = 500                    0.9996167                 0.9632500
    ## T = 1000                   0.9999667                 0.9918000
    ## T = 2000                   0.9977000                 0.9999667
    ##          LOOCV- KNN - Classic Train/Test- LOOCV- KNN - tDistance
    ## T = 50              0.4491667                          0.5143767
    ## T = 100             0.4553333                          0.7453428
    ## T = 500             0.4516667                          0.6579194
    ## T = 1000            0.4245833                          0.6659911
    ## T = 2000            0.3853000                          0.6806056
    ##          Train/Test- LOOCV- KNN - pDistance
    ## T = 50                            0.2760100
    ## T = 100                           0.1964317
    ## T = 500                           0.9719233
    ## T = 1000                          1.0000000
    ## T = 2000                          1.0000000
    ##          Train/Test- LOOCV- KNN - D_alpha_|chi|*
    ## T = 50                                 0.3903528
    ## T = 100                                0.7917022
    ## T = 500                                0.9986733
    ## T = 1000                               0.9999689
    ## T = 2000                               0.9864039
    ##          Train/Test- LOOCV- KNN - D_alpha_med* Train/Test- LOOCV- KNN - Classic
    ## T = 50                               0.3235694                        0.3602233
    ## T = 100                              0.3503217                        0.3691950
    ## T = 500                              0.9600061                        0.3866228
    ## T = 1000                             0.9923717                        0.3793056
    ## T = 2000                             0.9997417                        0.3596383
    ## 
    ## $scenario_3_clustering_results
    ##          K-Medoids - tDistance K-Medoids - pDistance
    ## T = 50               0.5582463             0.3817453
    ## T = 100              0.6334448             0.3807857
    ## T = 500              0.6358165             0.8273597
    ## T = 1000             0.6827180             0.9967833
    ## T = 2000             0.5904812             0.9998667
    ##          K-Medoids - D_{alpha_{|chi|}}* K-Medoids - D_{alpha_{WCSS}}*
    ## T = 50                        0.4041017                     0.3481343
    ## T = 100                       0.6877273                     0.3460660
    ## T = 500                       0.9883325                     0.6962727
    ## T = 1000                      1.0000000                     0.9843167
    ## T = 2000                      0.9939220                     0.9999000
    ##          Classical KNN
    ## T = 50       0.4029032
    ## T = 100      0.4052450
    ## T = 500      0.4276772
    ## T = 1000     0.4542878
    ## T = 2000     0.4938758

# Process LaTeX table

    library(magrittr)

    # function to write code for latex tables 
    to_latex_table = function(scenario_tables)
    {
      # process classification 
      
      latex_classif_code = scenario_tables[[1]] %>% formatC(format = "f", digits = 2) %>% 
        apply(MARGIN = 1, FUN = paste, collapse = " & ")
      latex_classif_code = paste(paste("$",names(latex_classif_code), "$ & "), latex_classif_code) %>% 
        paste(collapse = " \\ ")
      
      # process clustering
      
      latex_cluster_code = scenario_tables[[2]] %>%  formatC(format = "f", digits = 2) %>% 
        apply(MARGIN = 1, FUN = paste, collapse = " & ")
      latex_cluster_code = paste(paste("$",names(latex_cluster_code), "$ & "), latex_cluster_code) %>% 
        paste(collapse = " \\ ")
      
      latex_code = list(latex_classif_code = latex_classif_code, 
                        latex_cluster_code = latex_cluster_code)
      
      return(latex_code)
      
    }

    # Scenario 1
    to_latex_table(Scenario_1_tables)

    # Scenario 2
    to_latex_table(Scenario_2_tables)

    # Scenario 3
    to_latex_table(Scenario_3_tables)
