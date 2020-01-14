### RE-WORK COMPUTECODE FUNCTION TO WORK ON LISTS AND GIVE CONSISTENT OUTPUT ###

#this file is for computation of mapping codes between a list of DNA motifs that are given
#one application is that the DNA motifs are obtained from RNA and DNA and the pairing and the mapping code are estimated


#projectD generates the ideal PWM partner of a given PWM using a given code matrix
projectD=function(R,C){
  return(sapply(1:ncol(R), function(j){sapply(1:nrow(R), function(i){sum(C[,i] * R[,j])})}))
}

#readMotifs parses the combined.meme motif file containing all significant motifs as output from MEME-ChIP into a
#list of motifs which can be used in the Optimize pipeline
#it takes the path of a combined.meme file as input and returns a list of formatted PWMs

readMotifs = function(filename) {
  file = read.delim(file = filename, stringsAsFactors = F)
  file_arr = file[,1]
  file_index = suppressWarnings(which(!(is.na(as.integer(file_arr))), arr.ind = T))
  start = c(1, which(diff(file_index) != 1 & diff(file_index) != 0) + 1)
  end = c(start - 1, length(file_index))
  individual_motifs = split(file_index, cumsum(c(0, diff(file_index) > 1)))
  motif_list = list()
  for (i in 1:length(individual_motifs)) {
    motif = matrix(as.numeric(file[individual_motifs[[i]],1]), nrow = 4)
    dimnames(motif) = list(c('A', 'C', 'G', 'T'), c(1:ncol(motif)))
    motif_list[[i]] = motif
    names(motif_list)[i] = paste('motif_', i, sep = '')
  }
  return(motif_list)
}

#createRandomPWM generates a random PWM of a given length and of a given alphabet size
createRandomPWM = function(alph, length, dp = 7, lowEntropy=TRUE){
  #randomly generate a code matrix with the constraints
  #added ?variable for minimum entropy requirement per position of PWM to produce more discrete motifs
  library(entropy)
  rownames = c('rA', 'rC', 'rG', 'rT')
  code = matrix(nrow = alph, ncol = length, dimnames = list(rownames[1:alph]))
  for (i in 1:length) {
    if(isTRUE(lowEntropy)){
      entropy = 1
      while(entropy >= 0.5) {
        ord_col = XLessThan(alph, 1, roundTo = dp)
        entropy = entropy.empirical(y = ord_col)
      }
    } else {
      ord_col = XLessThan(alph, 1, roundTo = dp)
    }
    rand_ord = sample(1:alph, alph)
    rand_col = c()
    for (j in 1:alph) {
      rand_col = c(rand_col, ord_col[rand_ord[j]])
    }
    code[,i] = rand_col
  }
  return(code)
}

#createRandomPWMList will create two lists of matrices related by a given code, which may either be of
#uniform length or random lengths. It takes the same input as createRandomPWM, along with a vector specifying whether lengths
#of matrices should be uniform or randomly spaced between two values given in a numeric vector.
createRandomPWMList = function(alph, lengthRange, dp, listLengths, code, lowEntropy){
  listR = list()
  listD = list()
  for(i in 1:listLengths) {
    randomR = createRandomPWM(alph = alph, length = round(runif(min = min(lengthRange), n = 1, max = max(lengthRange)), digits = 0),dp = 2, lowEntropy = lowEntropy)
    listR[[i]] = randomR
    randomD = projectD(R = randomR, C = code)
    listD[[i]] = randomD
  }
  return(list(listR = listR, listD = listD))
}

#creates a random position of a motif with high entropy to simulate a non-discrete position in a motif
randomColumn = function(iterations){
  randomColumns = c()
  for (i in 1:iterations) {
    entropy = 0
    while(entropy == 0){
      randomColumn = XLessThan(X = 4, LessThan = 1)
      ent = entropy.empirical(y = randomColumn)
      if(ent > 1.35){
        entropy = 1
      }
    }
    randomColumns = c(randomColumns,randomColumn)
  }

  return(matrix(randomColumns, nrow = 4))
}

#padMotifList adds a random number (drawn from Poisson distribution) of high entropy positions to each end of the motifs
#in a given list
padMotifList = function(motifList,lambda){
  for(i in 1:length(motifList)){
    motifList[[i]] = matrix(c(randomColumn(iterations = rpois(n = 1, lambda = lambda)),
                              motifList[[i]], randomColumn(iterations = rpois(n = 1, lambda = lambda))),nrow = 4)
  }
  return(motifList)
}

#createRandomPairing produces a random pairing of motifs from which a code can be computed as a starting point to the
#optimize function
createRandomPairing = function(listR, listD) {
  condition = length(listR) >= length(listD)
  if(condition){
    pairing = matrix(c(c(sample(1:length(listR), size = length(listD))), c(1:length(listD))), nrow = length(listD))
  } else {
    pairing = matrix(c(c(1:length(listR), c(sample(1:length(listD), size = length(listR))))), nrow = length(listR))
  }
  return(pairing)
}

#createRandomCode2 produces a random code from the random pairing produced by the above function
createRandomCode2 = function(listR, listD){
  code = getCodeFromList(listR = listR, listD = listD, pairing = createRandomPairing(listR, listD), code = createRandomCode(alph = 4, dp = 2))
  return(code)
}



## objective calculates the absolute difference between an ideal and real partner of a motif and nomralises the value to the length of the motif
objective = function(R,D,C){
  objective=0
  objective =  objective + sum(abs(D - projectD(R = R, C = C)) / ncol(R))
  return(objective)
}


#sliding_objective calculates the objective value of a pair of matrices of unequal lengths, by sliding the shorter motif
#along the longer motif and calculating at each position. The function takes the same input as the regular objective function
#but outputs a list containing both the minumum objective value obtained and the reduced profile of the longer matrix
sliding_objective = function(R,D,C){
  shortR = ncol(R) < ncol(D)
  if(ncol(R)==ncol(D)){
    return(objective(R,D,C))
  }
  else if(shortR){
    tract = ncol(R)-1
    ranges = lapply(c(1:(ncol(D)-tract)), function(x) c(x:(x+tract)))
    objective_list = sapply(ranges, function(x) objective(R = R, D = D[,x], C = C))
    int_range = D[,ranges[[which.min(objective_list)]]]
    return(list(Minimum_objective_value=min(objective_list), Motif_section=int_range))
  } else {
    tract = ncol(D)-1
    ranges = lapply(c(1:(ncol(R)-tract)), function(x) c(x:(x+tract)))
    objective_list = sapply(ranges, function(x) objective(R = D, D = R[,x], C = C))
    int_range = R[,ranges[[which.min(objective_list)]]]
    return(list(Minimum_objective_value=min(objective_list), Motif_section=int_range))
  }
}

#getBestPartner calculates the objective value of a given code for a single motif of one list against all possible partners in a second list, and returns
#the best partner motif given the supplied code
getBestPartner = function(R,listD,C,Value=FALSE){
  obj_list = sapply(listD, function(y) sliding_objective(R = R, D = y, C = C)[[1]])
  if(isTRUE(Value)){
    return(list(partner = which.min(obj_list), objective = min(obj_list)))
  } else {
    return(which.min(obj_list))
  }
}

#getBestPair returns for each matrix in R the best fit in D using the current code matrix
#it returns a list of indices which represent the best fit to each of the matrices in listR
#getBestPair returns multpile-matching indices, unlike getBestPair2 which returns only unique matches
getBestPair = function(listR,listD,C){
  if(length(listR)>length(listD)){
    bestPairs = matrix(c(sapply(listD, function(x){getBestPartner(R = x, listD = listR, C = C)}),c(1:length(listD))),nrow = length(listD))
  } else {
    bestPairs = matrix(c(c(1:length(listR)),sapply(listR, function(x){getBestPartner(R = x, listD = listD, C = C)})),nrow = length(listR))
  }
  return(bestPairs)
}

#computeCode given a pairing of matrices in R and D compute the code
#this function needs to convert the matrices in R and D according to their pairing into
#matrices that can be used in the solve.QP function
computeCode = function(R,D){

  #Now the trick is to add zeros in case that, we are looking at the first column of D for the c and d variables
  #and vice versa for the a and b if we are looking at the second column of D:
  a=c(t(R)[,1],rep(0,ncol(R)*3))
  b=c(t(R)[,2],rep(0,ncol(R)*3))
  c=c(t(R)[,3],rep(0,ncol(R)*3))
  d=c(t(R)[,4],rep(0,ncol(R)*3))
  e=c(rep(0,ncol(R)),t(R)[,1], rep(0,ncol(R)*2))
  f=c(rep(0,ncol(R)),t(R)[,2], rep(0,ncol(R)*2))
  g=c(rep(0,ncol(R)),t(R)[,3], rep(0,ncol(R)*2))
  h=c(rep(0,ncol(R)),t(R)[,4], rep(0,ncol(R)*2))
  i=c(rep(0,ncol(R)*2),t(R)[,1],rep(0,ncol(R)))
  j=c(rep(0,ncol(R)*2),t(R)[,2],rep(0,ncol(R)))
  k=c(rep(0,ncol(R)*2),t(R)[,3],rep(0,ncol(R)))
  l=c(rep(0,ncol(R)*2),t(R)[,4],rep(0,ncol(R)))
  m=c(rep(0,ncol(R)*3),t(R)[,1])
  n=c(rep(0,ncol(R)*3),t(R)[,2])
  o=c(rep(0,ncol(R)*3),t(R)[,3])
  p=c(rep(0,ncol(R)*3),t(R)[,4])

  #the response is obtained by appending both columns in D
  y=c(t(D)[,1],t(D)[,2], t(D)[,3], t(D)[,4])

  library('quadprog')
  library('Matrix')
  library('matrixcalc')
  Q=cbind(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p)
  A = matrix(c(rep(c(1,0,0,0),4), rep(c(0,1,0,0),4), rep(c(0,0,1,0),4), rep(c(0,0,0,1),4),rep(c(1,rep(0,16)), 15), 1), nrow = 16)
  b0 = c(1,1,1,1,rep(0,16))
  if(is.positive.definite(t(Q)%*%Q)){
    sol=solve.QP(t(Q) %*% Q, y %*% Q,A,b0,meq=4)
  } else {
    pd = nearPD(t(Q)%*%Q)
    sol=solve.QP(pd$mat, y %*% Q,A,b0,meq=4)
  }


  return(round(matrix(sol$solution, nrow = 4), digits = 2))

}

## getCodeFromList iterates over the current pairing with computeCode and returns the mean code which explains the pairings
getCodeFromList = function(listR, listD, pairing, code) {

  codeList = lapply(X = c(1:nrow(pairing)), FUN = function(i){
    if(ncol(listR[[pairing[i,1]]]) < ncol(listD[[pairing[i,2]]])){
      computeCode(R = listR[[pairing[i,1]]], D = sliding_objective(R=listR[[pairing[i,1]]],D=listD[[pairing[i,2]]],C=code)[[2]])
    } else if(ncol(listR[[pairing[i,1]]]) > ncol(listD[[pairing[i,2]]])){
      computeCode(R = sliding_objective(R=listR[[pairing[i,1]]],D=listD[[pairing[i,2]]],C = code)[[2]], D = listD[[pairing[i,2]]])
    } else if(ncol(listR[[pairing[i,1]]]) == ncol(listD[[pairing[i,2]]])){
      computeCode(R = listR[[pairing[i,1]]], D = listD[[pairing[i,2]]])
    }})
  meanCode = Reduce("+", codeList) / length(codeList)
  return(meanCode)
}


#computeOverallObjective for a given code and pairing
#this function is for evaluation whether the objective value of all matrices is improving during optimization
computeOverallObjective = function(listR,listD,pairing,code){
  objectives = sapply(X = c(1:nrow(pairing)), FUN = function(i){
    if(ncol(listR[[pairing[i,1]]]) == ncol(listD[[pairing[i,2]]])){
      objective(listR[[pairing[i,1]]],listD[[pairing[i,2]]],code)
    } else {
      sliding_objective(listR[[pairing[i,1]]],listD[[pairing[i,2]]],code)[[1]]
    }})
  return(sum(objectives))
}


# XLessThan creates a vector of length X which sums to a value specified by LessThan
# This will be used in the createRandomCode function
XLessThan = function(X, LessThan,roundTo = 7) {
  element_list = c()
  for (i in 1:X) {
    if (i == 1) {
      element_list[i] = round(runif(n = 1, min = 0, max = LessThan),digits = roundTo)
    } else if (i != X) {
      element_list[i] = round(runif(n = 1, min = 0, max = LessThan - sum(element_list)),digits = roundTo)
    } else {
      element_list[i] = LessThan - sum(element_list)
    }
  }
  return(element_list)
}

#createRandomCode this function creates a random code matrix as a starting point
createRandomCode = function(alph, dp = 7){
  #randomly generate a code matrix with the constraints
  rownames = c('rA', 'rC', 'rG', 'rT')
  colnames = c('dA', 'dC', 'dG', 'dT')
  code = matrix(nrow = alph, ncol = alph, dimnames = list(rownames[1:alph], colnames[1:alph]))
  for (i in 1:alph) {
    ord_row = XLessThan(alph, 1, roundTo = dp)
    rand_ord = sample(1:alph, alph)
    rand_row = c()
    for (j in 1:alph) {
      rand_row = c(rand_row, ord_row[rand_ord[j]])
    }
    code[i,] = rand_row
  }
  return(code)
}


## assignPairingsToCodes clusters motif pairings to one of a user-provided list of codes
assignPairingsToCodes = function(pairings, codes, listR, listD){
  objs = sapply(X = c(1:nrow(pairings)), FUN = function(i){
    sapply(X = c(1:length(codes)), FUN = function(j){
      R = listR[[pairings[i,1]]]
      D = listD[[pairings[i,2]]]
      C = codes[[j]]
      if(ncol(R) == ncol(D)){
        objective(R,D,C)
      } else {
        sliding_objective(R,D,C)[[1]]
      }
    })
  })
  return(apply(objs,2,which.min))
}

#optimize(listR,listD) does the complete optimization
optimize =function(listR,listD, unique_pairs = F, dp = 7, plot = FALSE){
  #create random code
  iteration = -1
  iteration_list = c()
  objective_list = c()
  oldEstimate = createRandomCode2(listR, listD)
  if(unique_pairs == T) {
    pairing = getBestPair2(listR,listD,oldEstimate)
  } else {
    pairing = getBestPair(listR,listD,oldEstimate)
  }
  if(plot == TRUE){
    pair_error = c()
    pair_error = c(pair_error, length(which(pairing[,2] != c(1:nrow(pairing)))))
  }
  #while the objective value of the pairing and the code is improving (gets smaller) do
  differenceSmaller = 1
  while(differenceSmaller == 1){
    iteration = iteration + 1
    iteration_list = c(iteration_list, iteration)
    currentObjectiveValue = computeOverallObjective(listR,listD,pairing,oldEstimate)
    objective_list = c(objective_list, currentObjectiveValue)
    newCode = getCodeFromList(listR, listD, pairing, code = oldEstimate)
    if(unique_pairs == T) {
      pairing = getBestPair2(listR,listD,oldEstimate)
    } else {
      pairing = getBestPair(listR,listD,oldEstimate)
    }
    if(plot == TRUE){
      pair_error = c(pair_error, length(which(pairing[,2] != c(1:nrow(pairing)))))
    }
    newObjectiveValue = computeOverallObjective(listR,listD,pairing,newCode)
    if(currentObjectiveValue > newObjectiveValue){
      differenceSmaller =1
    }else{
      differenceSmaller = 0
    }
    oldEstimate=newCode
  }
  if (plot == TRUE) {
    return(data.frame(Iteration = iteration_list, Objective = objective_list, Pair_error = pair_error[2:length(pair_error)]))
  } else {
    if(length(listR) > length(listD)){
      return(list(pairing = matrix(pairing, nrow = length(listD)), objective = currentObjectiveValue/nrow(pairing), code = newCode))#, dimnames = list(c(1:20), c('RNA', 'DNA'))))
    } else {
      return(list(pairing = matrix(pairing, nrow = length(listR)), objective = currentObjectiveValue/nrow(pairing), code = newCode))#, dimnames = list(c(1:20), c('RNA', 'DNA'))))
    }
  }
}

#iterateOptimize iterates the optimize function a set number of times with desired settings and returns several lists which can
#be used for plotting the performance of the function in retreiving the code linking the motif sets
iterateOptimize = function(iterations, listR,listD, unique_pairs = F, dp = 2, plot = T) {
  objective_list = c()
  code_list = list()
  pairing_list = list()
  for(i in 1:iterations) {
    opt = optimize(listR = listR, listD = listD, unique_pairs = unique_pairs, dp = dp, plot = plot)
    objective_list = c(objective_list, opt$objective)
    pairing_list[[i]] = opt$pairing
    code_list[[i]] = opt$code
    print(paste(i, '%', sep = ""))
  }
  iterations = c(1:i)
  return(list(objective = objective_list, code = code_list, pairing = pairing_list, iterations = iterations))
}


###plotOptimization automatically uses ggplot to plot the optimize function as objective value vs iteration
plotOptimization = function(listR, listD, dp = 7, unique_pairs = F, starts) {
  for (i in 1:starts) {
    if (i == 1) {
      current = optimize(listR = listR, listD = listD, dp = dp, unique_pairs = unique_pairs, plot = T)
      current$group = i
    } else {
      new = optimize(listR = listR, listD = listD, dp = dp, unique_pairs = unique_pairs, plot = T)
      new$group = i
      current = rbind(current, new)
    }
  }
  colnames(current) = c('Iteration', 'Objective', 'group')

  library(ggplot2)

  theme_set(theme_bw())

  p <- ggplot(data = current, mapping = aes(x = Iteration, y = Objective, group = group, color = factor(group))) +
    geom_line(size=1.5) +
    labs(x = "Iteration", y = "Objective value", colour = 'Random starting point') +
    coord_cartesian(ylim=c(0, 2)) +
    geom_point(size=3) +
    # ggtitle('Training dataset: Subtle code') +
    theme(plot.title = element_text(size=22))

  return(current)
  #return(p)

}

##get original index of motif in list

originalIndex = function(motif,originalList){
  for(i in 1:length(originalList)){
    if(all.equal(motif,originalList[[i]]) == TRUE){
      originalIndex = i
    }
  }

  return(originalIndex)
}

##get unique codes from list including many (possibly duplicated) codes

findUniqueCodes = function(code_results){
  matrix_list = list()
  for(i in 1:length(code_results)){
    matrix_list[[i]] = code_results[[i]]$code
  }

  unique_code = list()
  m = 1
  while (length(matrix_list)>0) {
    unique_code[[m]] = matrix_list[[1]]
    index = c()
    for (i in 1:length(matrix_list)) {
      if(all.equal(unique_code[[m]],matrix_list[[i]])==TRUE){
        index = c(index,i)
      }
    }
    matrix_list = matrix_list[-(index)]
    m = m+1
  }
  return(unique_code)
}

allMatricesEqual = function(matrixList1,matrixList2){
  all_equal = rep(0,length(matrixList1))
  for(i in 1:length(matrixList1)){
    equal = rep(0,length(matrixList2))
    for(j in 1:length(matrixList2)){
      if(isTRUE(all.equal(matrixList1[[i]],matrixList2[[j]]))){
        equal[j] = 1
      }
    }
    if(sum(equal) > 0){
      all_equal[i] = 1
    }
  }
  if(all(all_equal == 1)){
    return(TRUE)
  } else {
    return(FALSE)
  }
}

whichMatrixEqual = function(target,options){
  test = lapply(options, function(x){isTRUE(all.equal(x,target))})
  res = which(test == TRUE)
  return(res)
}

collateResults = function(listR,listD,group,code){
  result = list()
  result$pairing = group
  result$code = getCodeFromList(listR = listR, listD = listD, pairing = group, code = code)
  result$objective = computeOverallObjective(listR = listR, listD = listD, pairing = group, code = result$code)
  return(result)
}

getBestAssignment = function(R,listD,codes){
  code_pairs = c(1:length(codes))
  code_objs = c(1:length(codes))
  for(n in 1:length(codes)){
    res = getBestPartner(R, listD, C = codes[[n]],Value = T)
    code_pairs[n] = res$partner
    code_objs[n] = res$objective
  }
  ass = which.min(code_objs)
  partner = code_pairs[which.min(code_objs)]
  return(c(partner,ass))
}


## findMultipleCodes is a wrapper of the optimize function which will report a user-defined number of codes. The function produces a user-defined number of random codes which is
## are used to cluster motif pairings and subsequently optimize these codes to produce the results. Each motif's best partner per code is discovered and then used to assign motif
## pairings to codes and form groups of motifs. Each group of motif pairs is then optimized to produce one code/pairing result per code requested by the user.

findMultipleCodes = function(listR,listD,initiations,cores,nCodes){
  library(foreach)
  library(doParallel)
  start = Sys.time()
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  clusterEvalQ(cl,library('computeCode'))
  #results = foreach(q=1:initiations) %dopar% {
    pairsChanged = 1
    iteration = 1
    old_groups = matrix(c(1:length(listD),1:length(listD)),nrow = length(listD))
    while (pairsChanged == 1) {
      if(iteration == 1){
        codes = vector(mode = 'list', length = nCodes)
        for(i in 1:nCodes){
          codes[[i]] = createRandomCode(alph = 4, dp = 2)
        }
      } else {
        codes = vector(mode = 'list', length = nCodes)
        for (i in 1:length(results)) {
          codes[[i]] = results[[i]]$code
        }
      }
      ass = lapply(listR, function(x){getBestAssignment(x,listD,codes)})
      pairings = matrix(c(c(1:length(listR)), unlist(ass)[c(T,F)], unlist(ass)[c(F,T)]),nrow = length(listR))
      group_index = sort(unique(pairings[,3]))
      groups = lapply(group_index, function(x){subset(pairings[,1:2], pairings[,3] == x)})
      if(iteration > 1 & allMatricesEqual(old_groups, groups)){
        pairsChanged = 0
      }
      results = lapply(groups, function(x){collateResults(listR = listR, listD = listD, group = x, code = codes[[whichMatrixEqual(target = x,groups)]])})
      iteration = iteration + 1
      old_groups = groups
    }
    show(iteration)
    show(results)
  #}
  stopCluster(cl)
  end = Sys.time() - start
  show(end)
  objectives = c()
  for (i in 1:length(results)) {
    objectives[i] = results[[i]]$obj
  }
  #best_res = results[[which.min(objectives)]][[1]]
  #return(best_res)
  return(results)
}

## splitMotifListByLength splits a motif list of mixed length motifs into several sublists of motifs with the same lengths for comparison with
## other identical length motifs

splitMotifListByLength = function(motifList){
  lengths = sapply(motifList, function(x) ncol(x))
  names(lengths) = NULL
  splitMotifs = lapply(c(min(lengths):max(lengths)), function(x) which(lengths == x))
  names(splitMotifs) = c(min(lengths):max(lengths))
  return(splitMotifs)
}


## compareMotifsAgainstDatabases compares a list of parsed motifs against a database of motifs (parsed into a text file, e.g. JASPAR, AtTRACT) for the purpose
## of removing possible artefacts arising from TF-binding/RNA-binding proteins. The function takes text files for both the database motifs and user-provided motifs,
## and a minimum length of motif to compare between the sets. The function will output the motifs which are present in both sets of motifs (and hence could be
## removed from the set before further analysis). This function is used during pre-processing of the real data.

compareMotifsAgainstDatabases = function(motifs, database, min) {

  library(computeCode)
  #file = read.table(file = database, fill = T, sep = '\t', stringsAsFactors = F)
  file = read.table(file = database, fill = T, stringsAsFactors = F)
  file = file[,1:4]
  file_arr = file[,1]
  file_index = suppressWarnings(expr = which(!(is.na(as.integer(file_arr))), arr.ind = T))
  start = c(1, which(diff(file_index) != 1 & diff(file_index) != 0) + 1)
  end = c(start - 1, length(file_index))
  individual_motifs = split(file_index, cumsum(c(0, diff(file_index) > 1)))
  motif_list = lapply(individual_motifs, function(x) matrix(as.numeric(t(file[x,])),nrow = 4))
  #user = readMotifs(filename = motifs)
  user = motifs
  db_split = splitMotifListByLength(motifList = motif_list)
  db_split = db_split[sapply(db_split, function(x) length(x) > 0)]
  user_split = splitMotifListByLength(motifList = user)
  user_split = user_split[sapply(user_split, function(x) length(x) > 0)]
  db_split = db_split[which(names(db_split) %in% names(user_split))]
  user_split = user_split[which(names(user_split) %in% names(db_split))]
  IDcode = matrix(c(1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1),nrow = 4)
  strip = c()
  for(i in 1:length(user_split)){
    objs = sapply(user_split[[i]], function(x) getBestPartner(R = user[[x]], listD = motif_list[db_split[[i]]], C = IDcode, Value = T)$objective)
    names(objs) = user_split[[i]]
    strip = c(strip, as.numeric(names(objs[which(objs < min)])))
  }
  return(strip)

}

## stripTerminalHighEntropy removes high-entropy (i.e. non-discrete) positions in motifs if they are present at either the first or last
## position of the motif. This results in the remaining motifs having a higher proportion of discrete positions which are well suited for
## pairing in the code learning algorithm and should avoid the problem of promiscuous motifs pairing with many other motifs due to the
## presence of non-discrete positions within the PWMs. It is used in pre-processing of the real data.

stripTerminalHighEntropy = function(motif, entropy = 0.5) {

  library(entropy)
  disc = findDiscretePositions(matrix = motif, entropyLimit = entropy)
  if(isFALSE(ncol(motif) %in% disc)){
    motif = motif[,1:(ncol(motif)-1)]
  }
  if(isFALSE(1 %in% disc)) {
    motif = motif[,2:ncol(motif)]
  }
  return(motif)
}























