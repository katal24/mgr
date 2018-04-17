########################### Methods for incomlete matrices #########################

#' @title Saaty index (CI)
#' @description Counts the value of inconsistency for the matrix
#' @param matrix - PC matrix (could be incomplete)
#' @return inconsistency of the matrix
#' @export
saaty <- function(matrix){
  n <- nrow(matrix)
  if( 0 %in% matrix ){
    matrix <- createHarkerMatrix(matrix)
  }
  alpha <- principalEigenValue(matrix)
  chopV((alpha - n)/(n-1))
}


#' @title Geometric consistency index (GCI)
#' @description Counts the value of inconsistency for the matrix
#' @param matrix - PC matrix (could be incomplete)
#' @return inconsistency of the matrix
#' @export
geometric <- function(matrix){
  dim =nrow(matrix)
  w <- geometricRankForIncomplete(matrix);
  e <- matrix(nrow = dim, ncol = dim)
  
  for(i in 1:dim){
    for(j in 1:dim){
      e[i,j] = matrix[i,j]*(w[j]/w[i])
    }
  }
  
  sum = 0
  
  for(i in 1:(dim-1)){
    for(j in (i+1):dim){
      if(e[i,j]!=0){
        sum = sum + log(e[i,j])^2
      }
    }
  }
  
  chopV(2/((dim-1)*(dim-2))*sum)
}


#' @title Koczkodaj matrix inconsistency
#' @description Counts the value of inconsistency for the matrix
#' @param matrix - PC matrix (could be incomplete)
#' @return inconsistency of the matrix
#' @export
koczkodaj <- function(matrix){
  triadsAndIdxs <- countIndexesForTriads("koczkodajForTriad", matrix)
  chopV(max(triadsAndIdxs))
}


#' @title Kazibudzki matrix inconsistency (LTI1)
#' @description Counts the value of inconsistency for the matrix
#' @param matrix - PC matrix (could be incomplete)
#' @return inconsistency of the matrix
#' @export
kazibudzkiLTI1 <- function(matrix){
  triadsAndIdxs <- countIndexesForTriads("grzybowskiLTI1ForTriad", matrix)
  chopV(mean(triadsAndIdxs))
}


#' @title Kazibudzki matrix inconsistency (LTI2)
#' @description Counts the value of inconsistency for the matrix
#' @param matrix - PC matrix (could be incomplete)
#' @return inconsistency of the matrix
#' @export
kazibudzkiLTI2 <- function(matrix){
  triadsAndIdxs <- countIndexesForTriads("grzybowskiLTI2ForTriad", matrix)
  chopV(mean(triadsAndIdxs))
}


#' @title Kazibudzki matrix inconsistency (CMLTI2)
#' @description Counts the value of inconsistency for the matrix
#' @param matrix - PC matrix (could be incomplete)
#' @return inconsistency of the matrix
#' @export
kazibudzkiCMLTI2 <- function(matrix){
  triadsAndIdxs <- countIndexesForTriads("kazibudzkiLTI2ForTriad", matrix)
  chopV(mean(triadsAndIdxs)/(1+max(triadsAndIdxs)))
}


#' @title Index of determinants (Pelaez and Lamata index)
#' @description Counts the value of inconsistency for the matrix
#' @param matrix - PC matrix (could be incomplete)
#' @return inconsistency of the matrix
#' @export
pelaeLamata <- function(matrix){
  triadsAndIdxs <- countIndexesForTriads("pelaeLamataForTriad", matrix)
  chopV(mean(triadsAndIdxs))
}


#' @title Kulakowski and Szybowski consistency index (I1)
#' @description Counts the value of inconsistency for the matrix
#' @param matrix - PC matrix (could be incomplete)
#' @return inconsistency of the matrix
#' @export
kulakowskiSzybowski <- function(matrix){
  triadsAndIdxs <- countIndexesForTriads("koczkodajForTriad", matrix)
  chopV(mean(triadsAndIdxs))
}


#' @title Kulakowski and Szybowski consistency index (I2)
#' @description Counts the value of inconsistency for the matrix
#' @param matrix - PC matrix (could be incomplete)
#' @return inconsistency of the matrix
#' @export
kulakowskiSzybowski2 <- function(matrix){
  triadsAndIdxs <- countIndexesForTriads("koczkodajForTriad", matrix)
  chopV((sqrt(sum(triadsAndIdxs**2)))/length(triadsAndIdxs))
}


#' @title Kulakowski and Szybowski consistency index (Ia)
#' @description Counts the value of inconsistency for the matrix
#' @param matrix - PC matrix (could be incomplete)
#' @return inconsistency of the matrix
#' @export
kulakowskiSzybowskiIa <- function(matrix, alfa, beta=0){
  triadsAndIdxs <- countIndexesForTriads("koczkodajForTriad", matrix)
  chopV(alfa*max(triadsAndIdxs) + (1-alfa)*mean(triadsAndIdxs))
}


#' @title Kulakowski and Szybowski consistency index (Iab)
#' @description Counts the value of inconsistency for the matrix
#' @param matrix - PC matrix (could be incomplete)
#' @return inconsistency of the matrix
#' @export
kulakowskiSzybowskiIab <- function(matrix, alfa, beta){
  triadsAndIdxs <- countIndexesForTriads("koczkodajForTriad", matrix)
  chopV(alfa*max(triadsAndIdxs) + beta*mean(triadsAndIdxs) + (1-alfa-beta)*(sqrt(sum(triadsAndIdxs**2)))/length(triadsAndIdxs))
}


#' @title Harmonic consistency index (HCI)
#' @description Counts the value of inconsistency for the matrix
#' @param matrix - PC matrix (could be incomplete)
#' @return inconsistency of the matrix
#' @export
harmonic <- function(matrix){
  n <- nrow(matrix)
  s <- sumValuesInColumns(matrix)
  sReverse <- 1/s
  hm <- n/(sum(sReverse))
  hm
  hci <- ((hm - n)*(n+1))/(n*(n-1))
  hci
}


#' @title Golden and Wand consistency index (GW)
#' @description Counts the value of inconsistency for the matrix
#' @param matrix - PC matrix (could be incomplete)
#' @return inconsistency of the matrix
#' @export
goldenWang <- function(matrix){
  n <- nrow(matrix)
  matrix[matrix==0] = 1
  g <- normalizeVector(geometricRankForIncomplete(matrix))
  a <- normalizeColumnsInMatrix(matrix)
  sum = 0
  for(i in 1:n){
    for(j in 1:n){
      sum = sum + abs(a[i,j]-g[i])
    }
  }
  gw <- 1/n*sum
  chopV(gw)
}


#' @title Salo and Hamalainen consistency index (SH)
#' @description Counts the value of inconsistency for the matrix
#' @param matrix - PC matrix (could be incomplete)
#' @return inconsistency of the matrix
#' @export
saloHamalainen <- function(matrix){
  n <- nrow(matrix)
  sum <- 0
  
  rMin <- matrix(nrow = n, ncol = n, data = 0)
  rMax <- matrix(nrow = n, ncol = n, data = 0)
  
  for(i in 1:n){
    for(j in 1:n){
      c <- rep(0,n)
      for(k in 1:n){
        c[k] = matrix[i,k]*matrix[k,j]
      }
      rMin[i,j] = min(c)
      rMax[i,j] = max(c)
    }
  }
  
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      sum <- sum + (rMax[i,j]-rMin[i,j])/((1+rMax[i,j])*(1+rMin[i,j]))
    }
  }
  
  cm <- 2/(n*(n-1))*sum
  chopV(cm)
}


#' @title Cavallo D'Apuzzo consistency index (CD)
#' @description Counts the value of inconsistency for the matrix
#' @param matrix - PC matrix (could be incomplete)
#' @return inconsistency of the matrix
#' @export
cavalloDapuzzo <- function(matrix){
  n <- nrow(matrix)
  prod <- 1
  licz <- 0
  for(i in 1:(n-2))
    for(j in (i+1):(n-1))
      for(k in (j+1):n){
        if(matrix[i,k]!=0 && (matrix[i,j]*matrix[j,k])!=0){
          licz <- licz+1
          prod <- prod * max( matrix[i,k]/(matrix[i,j]*matrix[j,k]), (matrix[i,j]*matrix[j,k])/matrix[i,k] )
        }
      }
  
  prod <- prod ^ (1/licz)
  chopV(prod)
}


#' @title Relative error (RE)
#' @description Counts the value of inconsistency for the matrix
#' @param matrix - PC matrix (could be incomplete)
#' @return inconsistency of the matrix
#' @export
relativeError <- function(matrix){
#  matrix[matrix==0] = 1
  n <- nrow(matrix)
 # matrix <- log(matrix)  - ostatnia zmiana
 # matrix[matrix==-Inf]=-100000 
  w <- apply(matrix, 1, function(row){mean(row)})
  C <- matrix(nrow = n, ncol = n, data = 0)
  E <- matrix(nrow = n, ncol = n, data = 0)
  
  for(i in 1:n)
    for(j in 1:n){
      C[i,j] = w[i] - w[j]
      E[i,j] = matrix[i,j] - C[i,j]
    }

  sumE = 0;
  sumA = 0;
  
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      sumE <- sumE + (E[i,j])^2
      sumA <- sumA + (matrix[i,j])^2
    }
  }
  
  chopV(sumE/sumA)
}



########################### Methods for triads #########################

#' @title Koczkodaj innconsistency for one triad 
#' @description Returns the Koczkodaj triad inconsistency
#' @param triad - vector of 3 elements
#' @return the Koczkodaj triad inconsistency
koczkodajForTriad <- function(triad) {
  min(abs(1-(triad[3]/(triad[1]*triad[2]))), abs(1-(triad[1]*triad[2])/triad[3]))
}


#' @title kazibudzki (LTI1) innconsistency for one triad 
#' @description Returns the Kazibudzki (LTI1) triad inconsistency
#' @param triad - vector of 3 elements
#' @return the Kazibudzki (LTI1)  triad inconsistency
kazibudzkiLTI1ForTriad <- function(triad) {
  abs(log((triad[1]*triad[2])/triad[3]))
}


#' @title kazibudzki (LTI2) innconsistency for one triad 
#' @description Returns the Kazibudzki (LTI2) triad inconsistency
#' @param triad - vector of 3 elements
#' @return the Kazibudzki (LTI2)  triad inconsistency
kazibudzkiLTI2ForTriad <- function(triad) {
  log((triad[1]*triad[2])/triad[3])**2
}


#' @title Index of determinants for one triad 
#' @description Returns Index of determinants for one triad 
#' @param triad - vector of 3 elements
#' @return index of determinants
pelaeLamataForTriad <- function(triad) {
  triad[3]/(triad[1]*triad[2]) + (triad[1]*triad[2])/triad[3] - 2 
}



########## Auxiliary functions ##########

#' @title Checks if the number is much greater than 0
#' @description return 0 if the number is close to zero else value
#' @param value - numeric value
#' @return 0 if the number is close to zero else value
chopV <- function(value){
  if(value>0.000001 || value < -0.000001){
    return(value)
  } else{
    return(0)
  }
}


#' @title The largest principal eigenvalue
#' @description Returns the largest principal eigenvalue of matrix
#' @param matrix - PC matrix
#' @return the largest principal eigenvalue of matrix
principalEigenValue <- function(matrix){
  maxEigenVal <- max(Mod(eigen(apply(matrix, 2, as.numeric))$values))
  chopV(maxEigenVal)
}


#' @title Rank list given as geometric means
#' @description Returns rank list given as geometric means of rows of the matrix
#' @param matrix - PC matrix
#' @return The rank list given as geometric means
geometricRank <- function(matrix){
  apply(matrix, 1, function(row){
    prod(row)^(1/length(row))
  })
}


#' @title Create Harker matrix
#' @description Takes incompletePC matrix and sets on diagonal the number of zeros (incremented by 1)
#' @param incompleteMatrix - incomplete PC matrix
#' @return Harker matrix
createHarkerMatrix <- function(incompleteMatrix){
  
  nonZerosElements <- apply(incompleteMatrix, 1, function(row){
    sum(row==0) + 1
  })
  
  diag(incompleteMatrix) <- nonZerosElements
  incompleteMatrix
}


#' @title Recipropal matrix
#' @description Recreates recipropal matrix of the basis of upper-triagle of matrix
#' @param matrix - PC matrix
#' @return recipropal PC matrix
recreatePCMatrix <- function(matrix){
  for(r in 1:nrow(matrix)-1)
    for(c in (r+1):nrow(matrix))
      matrix[c,r] <- 1/matrix[r,c]
    
    diag(matrix) <- 1
    matrix
}


#' @title Generates random priority vector
#' @param numOfElements - size of priority vector
#' @return normalized priority vector
generatePriorityVector <- function(numOfElements){
  randomVector <- runif(numOfElements)
  priorityVector <- normalizeVector(randomVector)
  priorityVector
}


#' @title Generates PC Matrix on the basis of the priority vector
#' @param priorityVector - the priotiry vector
#' @return reciprocal PC matrix
generatePCMatrixFromPV <- function(priorityVector){
  dim <- length(priorityVector)
  matrix <- matrix(0, nrow = dim, ncol = dim)
  for(r in 1:dim-1)
    for(c in (r+1):dim)
      matrix[r,c] <- priorityVector[r]/priorityVector[c]
  
  recreatePCMatrix(matrix)
}


#' @title Generates PC Matrix on the basis of the size
#' @param numOfElements - number of elements of the matrix
#' @return reciprocal PC matrix
generatePCMatrix<- function(numOfElements){
  priorityVector <- generatePriorityVector(numOfElements)
  generatePCMatrixFromPV(priorityVector)
}


#' @title Disturbs the PC matrix
#' @description Disturbs the PC matrix in order to obtain inconsistency
#' @param matrix - PC matrix
#' @param scale - extend of disorders. This parametr is the upper limit of the interval that is used to scale the elements. The lower limit is defined as 1 / scale
#' @return disturbed PC matrix
disturbPCMatrix<- function(matrix, scale){
  dim = nrow(matrix)
  
  for(r in 1:dim-1)
    for(c in (r+1):dim)
      matrix[r,c] <- runif(1, min = 1/scale, max = scale) * matrix[r,c]
    
    diag(matrix) <- 1
    recreatePCMatrix(matrix)
}


#' @title Generates disturbed PC Matrix
#' @param numOfElements - size of PC matrix
#' @param matrix - PC matrix
#' @return disturbed PC matrix
generateDisturbedPCMatrix<- function(numOfElements, scale){
  matrix <- generatePCMatrix(numOfElements)
  disturbPCMatrix(matrix, scale)
}


#' @title Generate broken PC Matrix
#' @description 
#' @param numOfElements -
#' @param grade - level of incompleteness [%]
#' @return 
generateBrokenPCMatrix<- function(numOfElements, scale, grade){
  matrix <- generateDisturbedPCMatrix(numOfElements, scale)
  breakPCMatrix(matrix, grade)
}


#' @title Breaks PC Matrix
#' @description Breaks PC Matrix to obtain reciprocal, incomplete matrix
#' @param  matrix - PC matrix
#' @param grade - percentage of the value to be removed (not applicable to the diagonal)
#' @return incomplete PC matrix
breakPCMatrix<- function(matrix, grade){
  
  dim <- nrow(matrix)
  numOfEmptyElements <- grade*0.01*(dim*(dim-1))
  
  while(numOfEmptyElements > 1){
    rowToClear <- sample(1:dim, 1)
    colToClear <- sample(c(1:dim)[-rowToClear], 1)
    
    if(matrix[rowToClear, colToClear] != 0) {
      matrix[rowToClear, colToClear] <- matrix[colToClear, rowToClear] <- 0
      numOfEmptyElements <- numOfEmptyElements - 2
    }
  }
  
  matrix
  
}

#' @title Rank list for the incomplete matrix given as geometric means
#' @description Returns rank list given as geometric means of rows of the incomplete matrix
#' @param matrix - PC matrix
#' @return The rank list given as geometric means
geometricRankForIncomplete <- function(matrix){
  apply(matrix, 1, function(row){
    prod(row[row!=0])^(1/length(row[row!=0]))
  })
}


#' @title Sum columns
#' @description Sums up values in each matrix column
#' @param  matrix - PC matrix
#' @return vector of sum
sumValuesInColumns <- function(matrix){
  apply(matrix, 2, function(col){
    sum(col)
  })
}


#' @title Normalize vector
#' @description Normalize vector to add up to 1.
#' @param  vector - vector to normalize
#' @return normalized vector
normalizeVector <- function(vector){
  vector/sum(vector)
}


#' @title Normalize columns in matrix
#' @description Normalize columns in matrix to add up to 1.
#' @param  matrix - matrix to normalize
#' @return matrix of normalized columns
normalizeColumnsInMatrix <- function(matrix){
  apply(matrix, 2, function(col){
    normalizeVector(col)
  })
}



#' @title Generates triad from the matrix on the basis of tuple
#' @param matrix - PC matrix
#' @param tuple - vector of three numbers
#' @return triad
makeATriad <- function(matrix, tuple){
  c(
    matrix[tuple[1],tuple[2]],
    matrix[tuple[2],tuple[3]],
    matrix[tuple[1],tuple[3]]
  )
}


#' @title Generates triades from the matrix
#' @param matrix - PC matrix
#' @return triads
generateTriads <- function(matrix) {
  dim <- nrow(matrix)
  tuples <- combn(dim ,3)
  triads <- apply(tuples, 2, function(x) makeATriad(matrix=matrix, tuple=x))
  numOfTriads <- dim(triads)[2]
  numOfTriadsToDelete <- c()
  
  for( c in 1:numOfTriads ){
    if ( triads[1,c]==0 || triads[2,c]==0 || triads[3,c]==0) {
      numOfTriadsToDelete <- c(numOfTriadsToDelete, c)
    }
  }
  
  if(!is.null(numOfTriadsToDelete)){
    triads <- triads[,-numOfTriadsToDelete]
  } 
  
  triads
}


#' @title Generates inconsistency indexes for each triad from the matrix
#' @param methodName - name of the method which computes inconsistency index for triad
#' @param matrix - PC matrix
#' @return vector of inconsistency indexes for triads
countIndexesForTriads <- function(methodName, matrix){
  triads <- generateTriads(matrix)
  
  triadIdxs <- apply(triads, 2, methodName)
  triadIdxs
}



######################### TESTS - Monte Carlo ###########################

#' @title Explore the incomplete PC matrixes for every method and scale <1.1, 1.2, ... , 4.0>
#' @description Examines what is the relative error between the full matrixes and indomplete matrixes.
#' @param numOfElements - dimension of tested matrix
#' @param gradeOfIncomplete -  percentage of the value to be removed (not applicable to the diagonal)
#' @param numOfAttempts - number of tested matrixes
#' @param numOfAttemptsForOneMatrix - number of test cases for each matrix
#' @param alfa - a parameter for kulakowskiSzybowskiIa method
#' @param beta - a parameter for kulakowskiSzybowskiIab method
#' @return average value of the relative error between the full matrix and indomplete matrix
test <- function(numOfElements, gradeOfIncomplete, numOfAttempts, numOfAttemptsForOneMatrix, alfa=0, beta=0){
  
  results <- matrix(nrow=31, ncol=16, data=0)
  counter <- 1
  
  for(i in seq(1.1, 4, 0.1)){
    print(counter)
    results[counter,] <- monteCarloOnTheSameMatrix(numOfElements, i, gradeOfIncomplete, numOfAttempts, numOfAttemptsForOneMatrix, alfa=0, beta=0)
    counter <- counter+1
  }
  
  for(i in 1:16){
    results[31, i] = sum(results[,i])/30
  }
  
  results
}


#' @title Explore the incomplete PC matrix
#' @description Examines what is the relative error between the full matrix and indomplete matrix
#' @param methodName - a name of the method which is tested
#' @param scale - extend of disorders. This parametr is the upper limit of the interval that is used to scale the elements. The lower limit is defined as 1 / scale
#' @param numOfElements - dimension of tested matrix
#' @param gradeOfIncomplete -  percentage of the value to be removed (not applicable to the diagonal)
#' @param numOfAttempts - number of test cases
#' @param alfa - a parameter for kulakowskiSzybowskiIa method
#' @param beta - a parameter for kulakowskiSzybowskiIab method
#' @return average value of the relative error between the full matrix and indomplete matrix
exploreMatrix <- function(methodName, scale, numOfElements, gradeOfIncomplete, numOfAttempts, alfa=0, beta=0) {

  matrix <- generateDisturbedPCMatrix(numOfElements, scale)
  dim <- ncol(matrix)
  if(alfa==0 && beta==0){
    realIdx <- methodName(matrix)
  } else {
    realIdx <- methodName(matrix, alfa, beta)
  }
  
  vectorOfIdsx <- integer(numOfAttempts)
  
  for( i in 1:numOfAttempts ) {
    brokenMatrix <- matrix(nrow = dim, ncol = dim, data = 0)
    
    while( !(0 %in% (matrix %^% (n-1))) ){
      brokenMatrix <- breakPCMatrix(matrix, gradeOfIncomplete)
    }

    if(alfa==0 && beta==0){
      idx <- methodName(brokenMatrix)
    } else {
      idx <- methodName(brokenMatrix, alfa, beta)
    }
    
    vectorOfIdsx[i] <- abs(realIdx - idx)
  }
  
  incompleteIdx <- mean(vectorOfIdsx)
  
  abs(incompleteIdx/realIdx*100)
}


#' @title Explore the incomplete PC matrixes for one method
#' @description Examines what is the relative error between the full matrix and indomplete matrix
#' @param methodName - name of the method which is tested
#' @param scale - extend of disorders. This parametr is the upper limit of the interval that is used to scale the elements. The lower limit is defined as 1 / scale
#' @param numOfElements - dimension of tested matrix
#' @param gradeOfIncomplete -  percentage of the value to be removed (not applicable to the diagonal)
#' @param numOfAttempts - number of tested matrixes
#' @param numOfAttemptsForOneMatrix - number of test cases for each matrix
#' @param alfa - a parameter for kulakowskiSzybowskiIa method
#' @param beta - a parameter for kulakowskiSzybowskiIab method
#' @return average value of the relative error between the full matrix and indomplete matrix for one method
monteCarlo <- function(methodName, scale, numOfElements, gradeOfIncomplete, numOfAttempts, numOfAttemptsForOneMatrix, alfa=0, beta=0) {
  vectorOfIdsx <- integer(numOfAttempts)
  for( i in 1:numOfAttempts ) {
    vectorOfIdsx[i] <- exploreMatrix(methodName, scale, numOfElements, gradeOfIncomplete, numOfAttemptsForOneMatrix, alfa, beta)
  }
  
  incompleteIdx <- mean(vectorOfIdsx)
  incompleteIdx
}


#' @title Explore the incomplete PC matrixes
#' @description Examines what is the relative error between the full matrixes and indomplete matrixes.
#' @param numOfElements - dimension of tested matrix
#' @param scale - extend of disorders. This parametr is the upper limit of the interval that is used to scale the elements. The lower limit is defined as 1 / scale
#' @param gradeOfIncomplete -  percentage of the value to be removed (not applicable to the diagonal)
#' @param numOfAttempts - number of tested matrixes
#' @param numOfAttemptsForOneMatrix - number of test cases for each matrix
#' @param alfa - a parameter for kulakowskiSzybowskiIa method
#' @param beta - a parameter for kulakowskiSzybowskiIab method
#' @return average value of the relative error between the full matrix and indomplete matrix
monteCarloOnTheSameMatrix <- function(numOfElements, scale, gradeOfIncomplete, numOfAttempts, numOfAttemptsForOneMatrix, alfa=0, beta=0) {
  numOfMethods <- 16
  vectorOfIdsx <- integer(numOfAttempts)
  idxs = integer(numOfMethods)
  
  for( i in 1:numOfAttempts ) {
    idxs <- idxs + exploreMatrixOnTheSameMatrix(numOfElements, scale, gradeOfIncomplete, numOfAttemptsForOneMatrix, alfa, beta)
  }
  
  incompleteIdx <- idxs/numOfAttempts
  incompleteIdx
}


#' @title Explore the incomplete PC matrixes
#' @description Examines what is the relative error between the full matrixes and indomplete matrixes for each method
#' @param numOfElements - dimension of tested matrix
#' @param scale - extend of disorders. This parametr is the upper limit of the interval that is used to scale the elements. The lower limit is defined as 1 / scale
#' @param gradeOfIncomplete -  percentage of the value to be removed (not applicable to the diagonal)
#' @param numOfAttempts - number of test cases
#' @param alfa - a parameter for kulakowskiSzybowskiIa method
#' @param beta - a parameter for kulakowskiSzybowskiIab method
#' @return average value of the relative error between the full matrix and indomplete matrix for each method
exploreMatrixOnTheSameMatrix <- function(numOfElements, scale, gradeOfIncomplete, numOfAttempts, alfa=0, beta=0) {
  numOfMethods <- 16;
  
  realIdx <- integer(numOfMethods)
  brokenIdx <- integer(numOfMethods)
  differences <- integer(numOfMethods)
  
  matrix <- generateDisturbedPCMatrix(numOfElements, scale)
  
  # oblicza prawidłowe wartości wsp. niesp.
  for( i in 1:numOfMethods ) {
    realIdx[i] <- runMethod(i, matrix, alfa, beta)
  }

  # dekompletuje macierz i oblicza wszystkie współczynniki niespójności. Powtarza się to numOfAttemps razy
  for( i in 1:numOfAttempts ) {
    brokenMatrix <- breakPCMatrix(matrix, gradeOfIncomplete)
    for(j in 1:numOfMethods) {
      brokenIdx[j] <- brokenIdx[j] + abs(realIdx[j] - runMethod(j, brokenMatrix, alfa, beta))
    }
  }
  
  
  for( k in 1:numOfMethods ) {
    brokenIdx[k] <- brokenIdx[k]/numOfAttempts
    differences[k] <- brokenIdx[k]/realIdx[k]*100
  }
  
  differences
}


#' @title Runs method which counts inconsistency index
#' @param nr - the number of method
#' @param matrix - PC matrix (could be incomplete)
#' @param alfa - a parameter for kulakowskiSzybowskiIa method
#' @param beta - a parameter for kulakowskiSzybowskiIab method
#' @return the inconsistency index
runMethod <- function(nr, matrix, alfa=0, beta=0){
  switch(nr,
         "1"={
           saaty(matrix)
         },
         "2"={
           geometric(matrix)
         },
         "3"={
           koczkodaj(matrix)
         },
         "4"={
           kazibudzkiLTI1(matrix)
         },
         "5"={
           kazibudzkiLTI2(matrix)
         },
         "6"={
           kazibudzkiCMLTI2(matrix)
         },
         "7"={
           pelaeLamata(matrix)
         },
         "8"={
           kulakowskiSzybowski(matrix)
         },
         "9"={
           kulakowskiSzybowski2(matrix)
         },
         "10"={
           kulakowskiSzybowskiIa(matrix, alfa)
         },
         "11"={
           kulakowskiSzybowskiIab(matrix, alfa, beta)
         },
         "12"={
           harmonic(matrix)
         },
         "13"={
           goldenWang(matrix)
         },
         "14"={
           saloHamalainen(matrix)
         },
         "15"={
           cavalloDapuzzo(matrix)
         },
         "16"={
           relativeError(matrix)
         }
  )
}