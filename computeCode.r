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



## objective calculates the absolute difference between an ideal and real partner of a motif and normalises the value to the length of the motif
objective = function(R,D,C){
  return(sum(abs(D - projectD(R = R, C = C)) / ncol(R)))
}


### slide finds the window of lowest objective in a motif pair of unequal ncol()
slideValue = function(long,short,C){
  tract = ncol(short)
  len = ncol(long) - tract
  return(min(sapply(lapply(0:len, function(y) 1:tract + y), function(x) objective(R = short, D = long[,x], C = C))))
}

### slide finds the window of lowest objective in a motif pair of unequal ncol()
slideRegion = function(long,short,C){
  tract = ncol(short)
  len = ncol(long) - tract
  regions = lapply(0:len, function(y) 1:tract + y)
  cols = regions[[which.min(sapply(regions, function(x) objective(R = short, D = long[,x], C = C)))]]
  return(long[,cols])
}


#sliding_objective calculates the objective value of a pair of matrices of unequal lengths, by sliding the shorter motif
#along the longer motif and calculating at each position. The function takes the same input as the regular objective function
#but outputs a list containing both the minumum objective value obtained and the reduced profile of the longer matrix
sliding_objective = function(R,D,C){
  equal = ncol(R) == ncol(D)
  shortR = ncol(R) < ncol(D)
  return(ifelse(equal, yes = objective(R,D,C), no = ifelse(shortR, yes = slideValue(long = D, short = R, C = C), no = slideValue(long = R, short = D, C = C))))
}


#getBestPartner calculates the objective value of a given code for a single motif of one list against all possible partners in a second list, and returns
#the best partner motif given the supplied code
getBestPartner = function(R,listD,C){
  which.min(sapply(listD, function(y) sliding_objective(R = R, D = y, C = C)))
  return(which.min(sapply(listD, function(y) sliding_objective(R = R, D = y, C = C))))
}

# #getBestPair returns for each matrix in R the best fit in D using the current code matrix
# #it returns a list of indices which represent the best fit to each of the matrices in listR
# #getBestPair returns multpile-matching indices, unlike getBestPair2 which returns only unique matches
# getBestPair = function(listR,listD,C){
#   start = Sys.time()
#   if(length(listR)>length(listD)){
#     bestPairs = matrix(c(sapply(listD, function(x){getBestPartner(R = x, listD = listR, C = C)}),c(1:length(listD))),nrow = length(listD))
#   } else {
#     bestPairs = matrix(c(c(1:length(listR)),sapply(listR, function(x){getBestPartner(R = x, listD = listD, C = C)})),nrow = length(listR))
#   }
#   Sys.time() - start
#   return(bestPairs)
# }

#getBestPair returns for each matrix in R the best fit in D using the current code matrix
#it returns a list of indices which represent the best fit to each of the matrices in listR
#getBestPair returns multpile-matching indices, unlike getBestPair2 which returns only unique matches
getBestPair = function(listR,listD,C){
  start = Sys.time()
  longR = length(listR) > length(listD)
  ifelse(test = longR, yes = matrix(c(sapply(listD, function(x){getBestPartner(R = x, listD = listR, C = C)}),c(1:length(listD))),nrow = length(listD)),
         no = matrix(c(c(1:length(listR)),sapply(listR, function(x){getBestPartner(R = x, listD = listD, C = C)})),nrow = length(listR)))
  Sys.time() - start
  return(bestPairs)
}

#getBestPair returns for each matrix in R the best fit in D using the current code matrix
#it returns a list of indices which represent the best fit to each of the matrices in listR
#getBestPair returns multpile-matching indices, unlike getBestPair2 which returns only unique matches
getBestPair = function(listR,listD,C,cores){
  start = Sys.time()
  library(foreach)
  library(doParallel)
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  clusterEvalQ(cl,library('computeCode'))
  longR = length(listR) > length(listD)
  if(longR) {
    res = foreach(q=1:length(listD)) %dopar% {
      c(getBestPartner(R = listD[[q]], listD = listR, C = C),q)
    }
  } else {
    res = foreach(q=1:length(listR)) %dopar% {
      return(c(q,getBestPartner(R = listR[[q]], listD = listD, C = C)))
    }
  }
  stopCluster(cl)
  bestPairs = do.call('rbind', res)
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
  sol = ifelse(test = is.positive.definite(t(Q)%*%Q), yes = solve.QP(t(Q) %*% Q, y %*% Q,A,b0,meq=4), no = solve.QP(nearPD(t(Q)%*%Q)$mat, y %*% Q,A,b0,meq=4))
  return(matrix(sol[[1]], nrow = 4))
}

## getCodeFromList iterates over the current pairing with computeCode and returns the mean code which explains the pairings
getCodeFromList = function(listR, listD, pairing, code) {
  codeList = lapply(X = c(1:nrow(pairing)), FUN = function(i){
    if(ncol(listR[[pairing[i,1]]]) < ncol(listD[[pairing[i,2]]])){
      computeCode(R = listR[[pairing[i,1]]], D = slideRegion(short = listR[[pairing[i,1]]],long = listD[[pairing[i,2]]],C = code))
    } else if(ncol(listR[[pairing[i,1]]]) > ncol(listD[[pairing[i,2]]])){
      computeCode(R = slideRegion(long = listR[[pairing[i,1]]],short = listD[[pairing[i,2]]],C = code), D = listD[[pairing[i,2]]])
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
    equalCols = ncol(listR[[pairing[i,1]]]) == ncol(listD[[pairing[i,2]]])
    return(ifelse(test = equalCols, yes = objective(listR[[pairing[i,1]]],listD[[pairing[i,2]]],code), no = sliding_objective(listR[[pairing[i,1]]],listD[[pairing[i,2]]],code)))
  })
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
      equalCols = ncol(R) == ncol(D)
      return(ifelse(test = equalCols, yes = objective(R,D,C), no = sliding_objective(R,D,C)))
    })
  })
  return(apply(objs,2,which.min))
}

#optimize(listR,listD) does the complete optimization
optimize =function(listR,listD, unique_pairs = F, dp = 7, plot = FALSE, cores){
  #create random code
  iteration = -1
  iteration_list = c()
  objective_list = c()
  oldEstimate = createRandomCode2(listR, listD)
  if(unique_pairs == T) {
    pairing = getBestPair2(listR,listD,oldEstimate)
  } else {
    pairing = getBestPair(listR,listD,oldEstimate, cores = cores)
  }
  if(plot == TRUE){
    pair_error = c()
    pair_error = c(pair_error, length(which(pairing[,2] != c(1:nrow(pairing)))))
  }
  #while the objective value of the pairing and the code is improving (gets smaller) do
  differenceSmaller = 1
  while(differenceSmaller == 1){
    message(iteration)
    iteration = iteration + 1
    iteration_list = c(iteration_list, iteration)
    currentObjectiveValue = computeOverallObjective(listR,listD,pairing,oldEstimate)
    objective_list = c(objective_list, currentObjectiveValue)
    newCode = getCodeFromList(listR, listD, pairing, code = oldEstimate)
    if(unique_pairs == T) {
      pairing = getBestPair2(listR,listD,oldEstimate)
    } else {
      pairing = getBestPair(listR,listD,oldEstimate, cores = cores)
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
  return(which(sapply(originalList, function(x) isTRUE(all.equal(motif,x))) == TRUE))
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

getBestAssignment = function(R,long,codes,cores){
  code_pairs = c(1:length(codes))
  code_objs = c(1:length(codes))
  for(n in 1:length(codes)){
    library(foreach)
    library(doParallel)
    cl <- makeCluster(cores)
    registerDoParallel(cl)
    clusterEvalQ(cl,library('computeCode'))
    res = foreach(q=1:length(long)) %dopar% {
      sliding_objective(R = R, D = long[[q]], C = codes[[n]])
    }
    stopCluster(cl)
    code_objs[n] = min(unlist(res))
    code_pairs[n] = which.min(unlist(res))
  }
  ass = which.min(code_objs)
  partner = code_pairs[which.min(code_objs)]
  return(c(partner,ass))
}


### score matrix per code ###
scoreMatrix = function(short, long, codes, cores){
  library(foreach)
  library(doParallel)
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  clusterEvalQ(cl,library('computeCode'))
  mats = vector(mode = 'list', length = length(codes))
  for (z in 1:length(codes)){
    res = foreach(q=1:length(short)) %dopar% {
      lapply(long, function(x) sliding_objective(R = short[[q]], D = x, C = codes[[z]]))
    }
    mats[[z]] = matrix(data = unlist(res), byrow = T, nrow = length(res))
  }
  stopCluster(cl)
  return(mats)
}

### rowMins is a function to return the minimum value per row of a numeric matrix ###
rowMins = function(mat) {
  return(sapply(1:nrow(mat), function(x) min(mat[x,])))
}

### rowWhichMins is a function to return the minimum value per row of a numeric matrix ###
rowWhichMins = function(mat) {
  return(sapply(1:nrow(mat), function(x) which.min(mat[x,])))
}

assignPairingsToCodes = function(scoreMatrix) {
  mins = matrix(unlist(lapply(scoreMatrix,rowMins)), ncol = length(scoreMatrix))
  whichmins = matrix(unlist(lapply(scoreMatrix,rowWhichMins)), ncol = length(scoreMatrix))
  ass = lapply(1:nrow(mins), function(x){
    cd = which.min(mins[x,])
    return(c(whichmins[x,cd], cd))
  })
  ass = cbind(1:length(short),matrix(unlist(ass), nrow = length(ass), byrow = T))
  return(ass)
}

optimizeGroup =function(short,long,group,cores,code){
  #create random code
  subShort = short[group[,1]]
  subLong = long[group[,2]]
  iteration = 1
  iteration_list = c()
  objective_list = c()
  #oldEstimate = getCodeFromList(listR = short, listD = long, pairing = group, code = code)
  oldEstimate = code
  pairing = getBestPair(subShort,subLong,oldEstimate, cores = cores)
  #while the objective value of the pairing and the code is improving (gets smaller) do
  differenceSmaller = 1
  while(differenceSmaller == 1){
    message(paste('Optimization iteration', iteration))
    iteration_list = c(iteration_list, iteration)
    currentObjectiveValue = computeOverallObjective(subShort,subLong,pairing,oldEstimate)
    objective_list = c(objective_list, currentObjectiveValue)
    newCode = getCodeFromList(subShort, subLong, pairing, code = oldEstimate)
    pairing = getBestPair(subShort,subLong,oldEstimate, cores = cores)
    newObjectiveValue = computeOverallObjective(subShort,subLong,pairing,newCode)
    if(currentObjectiveValue > newObjectiveValue){
      differenceSmaller =1
    }else{
      differenceSmaller = 0
    }
    oldEstimate=newCode
    iteration = iteration + 1
  }
  ori1 = sapply(pairing[,1], function(x) originalIndex(motif = subShort[[x]], originalList = short))
  ori2 = sapply(pairing[,2], function(x) originalIndex(motif = subLong[[x]], originalList = long))
  return(list(pairing = matrix(cbind(ori1,ori2), nrow = length(subShort)), objective = currentObjectiveValue/nrow(pairing), code = newCode))
}

collateResults = function(short,long,groups,codes,cores){
  opt = lapply(1:length(groups), function(x) {
    message(paste('Optimizing group', x))
    optimizeGroup(short = short, long = long, group = groups[[x]], cores = cores, code = codes[[x]])
  })
  return(opt)
}

## findMultipleCodes is a wrapper of the optimize function which will report a user-defined number of codes. The function produces a user-defined number of random codes which is
## are used to cluster motif pairings and subsequently optimize these codes to produce the results. Each motif's best partner per code is discovered and then used to assign motif
## pairings to codes and form groups of motifs. Each group of motif pairs is then optimized to produce one code/pairing result per code requested by the user.

findMultipleCodes = function(short,long,initiations,cores,nCodes,reportAll = F,iterationLimit, fileName){
  start = Sys.time()
  finalRes = vector(mode = 'list', length = initiations)
  for(q in 1:initiations) {
    message(paste('INITIATION', q))
    pairsChanged = 1
    iteration = 0
    old_groups = matrix(c(1:length(short),1:length(short)),nrow = length(short))
    res = list()
    while (pairsChanged == 1) {
      iteration = iteration + 1
      message(paste('Iteration', iteration))
      if(iteration == 1){
        codes = vector(mode = 'list', length = nCodes)
        for(i in 1:nCodes){
          #codes[[i]] = createRandomCode(alph = 4, dp = 2)
          codes[[i]] = t(createRandomPWM(alph = 4, dp = 2, length = 4, lowEntropy = T))
        }
      } else {
        codes = vector(mode = 'list', length = nCodes)
        for (i in 1:length(results)) {
          codes[[i]] = results[[i]]$code
        }
      }
      message('Generating objective matrix...')
      scoreMat = scoreMatrix(short = short, long = long, codes = codes, cores = cores)
      ass = assignPairingsToCodes(scoreMatrix = scoreMat)
      group_index = sort(unique(ass[,3]))
      groups = lapply(group_index, function(x){subset(ass[,1:2], ass[,3] == x)})
      results = collateResults(short = short, long = long, groups = groups, codes = codes, cores = cores)
      if(isTRUE(reportAll)){
        res[[iteration]] = list()
        res[[iteration]]$iteration = iteration
        res[[iteration]]$codes = lapply(results, function(x) x$code)
        res[[iteration]]$groups = lapply(results, function(x) x$pairing)
        res[[iteration]]$objectives = sapply(results, function(x) x$objective)
        res[[iteration]]$complexity = sapply(results, function(x) length(unique(x$pairing[,2])) / length(unique(x$pairing[,1]))*100)
        res[[iteration]]$iterations = iteration
      }
      if(iteration > 1 & allMatricesEqual(old_groups, groups)){
        message('Groups unchanged, terminating...')
        pairsChanged = 0
      } else if(iteration > iterationLimit){
        message('Iteration limit reached, terminating...')
        pairsChanged = 0
      } else if(iteration > 1){
        if(sum(sapply(results, function(x) x$objective)) > oldObj){
          message('Objectives increased, terminating...')
          pairsChanged = 0
        } else if(isTRUE(reportAll)){
          if(sum(res[[iteration]]$objectives) > sum(res[[iteration - 1]]$objectives)){
            message('Objectives increased, terminating...')
            pairsChanged = 0
          }
        } else if(isFALSE(allMatricesEqual(old_groups, groups))){
          message('Groups changed, continuing...')
          old_groups = groups
        }
      }
      oldObj = sum(sapply(results, function(x) x$objective))
      old_groups = groups
    }
    if(isTRUE(reportAll)){
      finalRes[[q]] = res
    } else {
      results$Iterations = iteration
      finalRes[[q]] = results
    }
    message('Saving results...')
    saveRDS(object = finalRes, file = fileName)
  }
  return(finalRes)
}

shuffleCode = function(code){
  shuffled = do.call('rbind', lapply(1:nrow(code), function(x) code[x,][sample(1:4)]))[,sample(1:4)]
  dimnames(shuffled) = list(c('A', 'C', 'G', 'T'), c('A', 'C', 'G', 'T'))
  return(shuffled)
}

addN = function(code){
  code = cbind(code,rep(0,4))
  code = rbind(code, rep(0,5))
  code[code < 0] = 0
  dimnames(code) = list(c('A', 'C', 'G', 'T', 'N'), c('A', 'C', 'G', 'T', 'N'))
  return(as.matrix(code))
}



