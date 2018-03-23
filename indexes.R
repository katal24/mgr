# przypisuje 0 jesli wartosc jest bliska 0
#' @export
chopV <- function(value){
  if(value>0.000001){
    return(value)
  } else{
    return(0)
  }
}


#' @title The largest principal eigenvalue
#' @description Returns the largest principal eigenvalue of matrix
#' @param matrix - PC matrix
#' @return the largest principal eigenvalue of matrix
#' @export
principalEigenValue <- function(matrix){
  chopM(matrix);
  maxEigenVal <- max(Mod(eigen(apply(matrix, 2, as.numeric))$values))
  chopV(maxEigenVal)
}


#' @title Rank list given as geometric means
#' @description Returns rank list given as geometric means of rows of the matrix
#' @param matrix - PC matrix
#' @return The rank list given as geometric means
#' @export
geometricRank <- function(matrix){
  apply(matrix, 1, function(row){
    prod(row)^(1/length(row))
  })
}


#' @title Recipropal matrix
#' @description Recreates recipropal matrix of the basis of upper-triagle of matrix
#' @param matrix - matrix
#' @return recipropal PC matrix
#' @export
recreatePCMatrix <- function(matrix){
  for(r in 1:dim(matrix)[1]-1)
    for(c in (r+1):dim(matrix)[1])
      matrix[c,r] <- 1/matrix[r,c]
    diag(matrix) <- 1
    matrix
}


#' @title Koczkodaj triad innconsistency
#' @description Returns the Koczkodaj triad inconsistency
#' @param triad - vector of 3 elements
#' @return the Koczkodaj triad inconsistency
#' @export
#koczkodajTriadIdx <- function(triad) {
#  min(abs(1-(triad[3]/(triad[1]*triad[2]))), abs(1-(triad[1]*triad[2])/triad[3]))
#}

# generuje macierz z trojkami z niepowtarzajacymi sie elementami
getUniqueTuples <- function(n){
  combn(n,3)
}

# generuje jedną trójkę (triad) z wartościami z macierzy
makeATriad <- function(matrix, tuple){
  c(
    matrix[tuple[1],tuple[2]],
    matrix[tuple[2],tuple[3]],
    matrix[tuple[1],tuple[3]]
  )
}

uniqueTriadsTuplesAndIdx <- function(matrix){
  tuples <- getUniqueTuples(dim(matrix)[1])
  triads <- apply(tuples, 2, function(x) makeATriad(matrix=matrix,tuple=x))
  triadIdxs <- apply(triads, 2, koczkodajTriadIdx)
  t(rbind(tuples, triads, c(triadIdxs)))
}



########## Auxiliary functions ##########

#### Anonymous function syntax
#(function(x) x * 10)(10)

############################## MGR ##########################################################
#' @title Generate priority vector
#' @description 
#' @param numOfElements - 
#' @return 
generatePriorityVector <- function(numOfElements){
  randomVector <- runif(numOfElements)
  priorityVector <- randomVector/sum(randomVector)
  priorityVector
}


#' @title Generate PC Matrix
#' @description 
#' @param priorityVector - 
#' @return 
generatePCMatrixFromPV <- function(priorityVector){
  rDim = length(priorityVector)
  cDim = rDim
  matrix <- matrix(0, nrow = rDim, ncol = cDim)
  for(r in 1:rDim-1)
    for(c in (r+1):cDim)
      matrix[r,c] <- priorityVector[r]/priorityVector[c]
    
  diag(matrix) <- 1
  recreatePCMatrix(matrix)
}

generateMatrixE <- function(matrix, priorityVector){
  rDim = dim(matrix)[1]
  cDim = rDim
  for(r in 1:rDim)
    for(c in 1:cDim)
      matrix[r,c] <- matrix[r,c]*(priorityVector[c]/priorityVector[r])
  
  matrix
}

#' @title Generate Priority Vector and PCMatrix 
#' @description 
#' @param numOfElements - 
#' @return 
generatePCMatrix<- function(numOfElements){
  priorityVector <- generatePriorityVector(numOfElements)
  generatePCMatrixFromPV(priorityVector)
}


#' @title Distribute PC Matrix
#' @description 
#' @param  matrix - PC matrix
#' @return 
distributePCMatrix<- function(matrix, scale){
  rDim = dim(matrix)[1]
  cDim = dim(matrix)[2]
 
  for(r in 1:rDim-1)
    for(c in (r+1):cDim)
      matrix[r,c] <- runif(1, min = 1/scale, max = scale) * matrix[r,c]
  
  # for(r in 1:rDim)
  #  for(c in 1:cDim)
  #    matrix[r,c] <- runif(1, min = 0.5, max = 1.5) * matrix[r,c]
  
  diag(matrix) <- 1
  matrix <- recreatePCMatrix(matrix)
  matrix
}


#' @title Generate distributed PC Matrix
#' @description 
#' @param numOfElements - 
#' @return 
generateDistributedPCMatrix<- function(numOfElements, scale){
  matrix <- generatePCMatrix(numOfElements)
  matrix <- distributePCMatrix(matrix, scale)
  matrix
}


#' @title Break PC Matrix
#' @description 
#' @param  matrix - PC matrix
#' @return 
breakPCMatrix<- function(matrix, grade){
  dim <- dim(matrix)[1]
  numOfEmptyElements <- round(grade*0.01*(dim*dim-dim))
  while(numOfEmptyElements > 1){
    rowToClear <- sample(1:dim, 1)
    colToClear <- sample(1:dim, 1)
    if(rowToClear != colToClear && matrix[rowToClear, colToClear] != 0 && matrix[colToClear, rowToClear] != 0) {
      matrix[rowToClear, colToClear] <- 0
      matrix[colToClear, rowToClear] <- 0
      numOfEmptyElements <- numOfEmptyElements - 2
    }
  }
  matrix
}


#' @title Generate broken PC Matrix
#' @description 
#' @param numOfElements -
#' @param grade - level of incompleteness [%]
#' @return 
generateBrokenPCMatrix<- function(numOfElements, grade){
  matrix <- generateDistributedPCMatrix(numOfElements)
  breakPCMatrix(matrix, grade)
}


######################### Auxiliary methods #################################
sumValuesInColumns <- function(matrix){
  apply(matrix, 2, function(col){
    col[col==0]=1
    sum(col)
  })
}

normalizeVector <- function(vector){
  vector/sum(vector)
}

normalizeColumnsInMatrix <- function(matrix){
  apply(matrix, 2, function(col){
    normalizeVector(col)
  })
}

generateTriadsForTuples <- function(matrix, tuples) {
  triads <- apply(tuples, 2, function(x) makeATriad(matrix=matrix,tuple=x))
  numOfTriads <- dim(triads)[2]
  numOfTriadsToDelete <- c()
  for( c in 1:numOfTriads ){
    if ( triads[1,c]==0 || triads[2,c]==0 || triads[3,c]==0) {
      numOfTriadsToDelete <- c(numOfTriadsToDelete, c)
    }
  }
  if(!is.null(numOfTriadsToDelete)){
    triads <- triads[,-numOfTriadsToDelete]
    tuples <- tuples[, -numOfTriadsToDelete]
  } 
  
  triads
}

getRealTriadIdx <- function(triadRow){
  triadRow[4]
}


uniqueTriadsTuplesAndIdxForInComplete <- function(methodName, matrix){
  tuples <- getUniqueTuples(dim(matrix)[1])
  triads <- generateTriadsForTuples(matrix, tuples)
  
  triadIdxs <- apply(triads, 2, methodName)
  t(rbind(triads, c(triadIdxs)))
}

# count average of vector
avg <- function(vector){
  w <- sum(vector)/length(vector)
  w
}

avgPositive <- function(vector){
  vector <- vector[vector!=-Inf]
  sum(vector)/length(vector)
}

########################### Methods for incomlete matrices #########################

#' @title Koczkodaj incomplete matrix inconsistency
#' @description Counts the value of inconsistency for the matrix
#' @param matrix - PC matrix
#' @return inconsistency
#' @export ----
#koczkodaj <- function(matrix){
#  triadsAndIdxs <- uniqueTriadsTuplesAndIdxForInComplete("koczkodajForTriad", matrix)
#  sortMatrix <- triadsAndIdxs[order(triadsAndIdxs[,4],decreasing = TRUE),]
#  if(is.null(dim(sortMatrix))){
#    chopV(getRealTriadIdx(sortMatrix))  
#  } else {
#    chopV(getRealTriadIdx(sortMatrix[1,]))
#  }
#}

koczkodaj <- function(matrix){
  triadsAndIdxs <- uniqueTriadsTuplesAndIdxForInComplete("koczkodajForTriad", matrix)
  chopV(max(triadsAndIdxs[,4]))
}

#' @title Koczkodaj triad innconsistency
#' @description Returns the Koczkodaj triad inconsistency
#' @param triad - vector of 3 elements
#' @return the Koczkodaj triad inconsistency
#' @export
koczkodajForTriad <- function(triad) {
  min(abs(1-(triad[3]/(triad[1]*triad[2]))), abs(1-(triad[1]*triad[2])/triad[3]))
}

#### metohod !!!
grzybowski <- function(matrix){
  triadsAndIdxs <- uniqueTriadsTuplesAndIdxForInComplete("koczkodajForTriad", matrix)
  chopV(avg(triadsAndIdxs[,4]))
}


#### method !!!
kazibudzkiLTI1 <- function(matrix){
  triadsAndIdxs <- uniqueTriadsTuplesAndIdxForInComplete("grzybowskiLTI1ForTriad", matrix)
  v = chopV(avg(triadsAndIdxs[,4]))
  v
}

kazibudzkiLTI1ForTriad <- function(triad) {
  abs(log((triad[1]*triad[2])/triad[3]))
}

#### method !!!
kazibudzkiLTI2 <- function(matrix){
  triadsAndIdxs <- uniqueTriadsTuplesAndIdxForInComplete("grzybowskiLTI2ForTriad", matrix)
  v = chopV(avg(triadsAndIdxs[,4]))
  v
}

kazibudzkiLTI2ForTriad <- function(triad) {
  log((triad[1]*triad[2])/triad[3])**2
}


#### method !!!
kazibudzkiCMLTI2 <- function(matrix){
  triadsAndIdxs <- uniqueTriadsTuplesAndIdxForInComplete("kazibudzkiCMLTI2ForTriad", matrix)
  chopV(avg(triadsAndIdxs[,4])/(1+max(triadsAndIdxs[,4])))
}

kazibudzkiCMLTI2ForTriad <- function(triad) {
  log((triad[1]*triad[2])/triad[3])**2
}


#### method !!!
pelaeLamata <- function(matrix){
  triadsAndIdxs <- uniqueTriadsTuplesAndIdxForInComplete("pelaeLamataForTriad", matrix)
  chopV(avg(triadsAndIdxs[,4]))
}

pelaeLamataForTriad <- function(triad) {
  triad[3]/(triad[1]*triad[2]) + (triad[1]*triad[2])/triad[3] - 2 
}

#### method !!!
pelaeLamata2 <- function(matrix){
  n <- dim(matrix)[1]
  sum <- 0
  for(i in 1:(n-2))
    for(j in (i+1):(n-1))
      for(k in (j+1):n){
       sum <- sum + (matrix[i,k]/(matrix[i,j]*matrix[j,k]) + (matrix[i,j]*matrix[j,k])/matrix[i,k] -2)
      }
  sum <- sum / choose(n,3)
  chopV(sum)
}

#### method !!!
kulakowskiSzybowski <- function(matrix){
  grzybowski(matrix)
}

#### method !!!
kulakowskiSzybowski2 <- function(matrix){
  triadsAndIdxs <- uniqueTriadsTuplesAndIdxForInComplete("koczkodajForTriad", matrix)
  chopV((sqrt(sum(triadsAndIdxs[,4]**2)))/length(triadsAndIdxs[,4]))
}

#### method !!!
kulakowskiSzybowskiIa <- function(matrix, alfa, beta=0){
  triadsAndIdxs <- uniqueTriadsTuplesAndIdxForInComplete("koczkodajForTriad", matrix)
  chopV(alfa*max(triadsAndIdxs[,4]) + (1-alfa)*avg(triadsAndIdxs[,4]))
}


#### method !!!
kulakowskiSzybowskiIab <- function(matrix, alfa, beta=0){
  triadsAndIdxs <- uniqueTriadsTuplesAndIdxForInComplete("koczkodajForTriad", matrix)
  chopV(alfa*max(triadsAndIdxs[,4]) + beta*avg(triadsAndIdxs[,4]) + (1-alfa-beta)*(sqrt(sum(triadsAndIdxs[,4]**2)))/length(triadsAndIdxs[,4]) )
}


#### method !!!
geometric <- function(matrix){
  dim =dim(matrix)[1]
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

geometricRankForIncomplete <- function(matrix){
  apply(matrix, 1, function(row){
    prod(row[row!=0])^(1/length(row[row!=0]))
  })
}

#### method !!!
harmonic <- function(matrix){
  n <- dim(matrix)[1]
  #filled <- length(matrix[matrix!=0])
  #n <- (filled/(n^2))*n
  s <- sumValuesInColumns(matrix)
  sReverse <- 1/s
  hm <- n/(sum(sReverse))
  hm
  hci <- ((hm - n)*(n+1))/(n*(n-1))
  hci
}


#### method !!!
saaty <- function(matrix){
  matrix[matrix==0] = 1
  matrix <- apply(matrix, 2, as.numeric)
  n <- nrow(matrix)
  alpha <- principalEigenValueSym(matrix)
  chopV((alpha - n)/(n-1))
}


#### method !!!
goldenWang <- function(matrix){
  n <- dim(matrix)[1]
  g <- normalizeVector(geometricRankForIncomplete(matrix))
  a <- normalizeColumnsInMatrix(matrix)
  sum = 0
  for(i in 1:n){
    for(j in 1:n){
      sum = sum + abs(a[i,j]-g[i])
    }
  }
  gw <- 1/n*sum
  gw
}


#### method !!!
saloHamalainen <- function(matrix){
  matrix[matrix==0] = 1
  n <- dim(matrix)[1]
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

#### method !!!
cavalloDapuzzo <- function(matrix){
  n <- dim(matrix)[1]
  sum <- 1
  licz <- 0
  for(i in 1:(n-2))
    for(j in (i+1):(n-1))
      for(k in (j+1):n){
        if(matrix[i,k]!=0 && (matrix[i,j]*matrix[j,k])!=0){
          licz = licz+1
          sum <- sum * max(matrix[i,k]/(matrix[i,j]*matrix[j,k]), (matrix[i,j]*matrix[j,k])/matrix[i,k])
        }
      }
  sum <- sum ^ (1/licz)
  chopV(sum)
}


#### method !!!
relativeError <- function(matrix){
  matrix[matrix==0] = 1
  n <- dim(matrix)[1]
  matrix <- log(matrix)
 # matrix[matrix==-Inf]=-100000 
  w <- apply(matrix, 1, function(row){avg(row)})
  C <- matrix(nrow = n, ncol = n, data = 0)
  E <- matrix(nrow = n, ncol = n, data = 0)
  
  n <- dim(matrix)[1]
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




######################### Monte Carlo ###########################


# generuje macierz zaburzoną, a potem przeprowadza x (numOfAttempts) prób zdekompletowania i liczy średnią wartość niespójności,
# potem znajduje różnicę tej średniej od prawdziwej wartości
exploreMatrix <- function(methodName, scale, numOfElements, gradeOfIncomplete, numOfAttempts, alfa=0, beta=0) {
#exploreMatrix <- function(methodName, numOfElements, gradeOfIncomplete, numOfAttempts, alfa=0, beta=0) {
  matrix <- generateDistributedPCMatrix(numOfElements, scale)
  
  if(alfa==0 && beta==0){
    realIdx <- methodName(matrix)
  } else {
    realIdx <- methodName(matrix, alfa, beta)
  }
  
  vectorOfIdsx <- integer(numOfAttempts)
  
  for( i in 1:numOfAttempts ) {
    brokenMatrix <- breakPCMatrix(matrix, gradeOfIncomplete)

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

#tests
test <- function(numOfElements, gradeOfIncomplete, numOfAttempts, numOfAttemptsForOneMatrix, alfa=0, beta=0){
  results <- matrix(nrow=31, ncol=17, data=0)
  counter <- 1
  
  for(i in seq(1.1, 4, 0.1)){
    print(counter)
    results[counter,] <- monteCarloOnTheSameMatrix(numOfElements, i, gradeOfIncomplete, numOfAttempts, numOfAttemptsForOneMatrix, alfa=0, beta=0)
    counter <- counter+1
  }
  
  for(i in 1:17){
    results[31, i] = sum(results[,i])/30
  }
  
  results
}

# Wykonuje exploreMatrixKoczkodaj x (numOfAttempts) razy i z tego bierze średnią
monteCarlo <- function(methodName, scale, numOfElements, gradeOfIncomplete, numOfAttempts, numOfAttemptsForOneMatrix, alfa=0, beta=0) {
  vectorOfIdsx <- integer(numOfAttempts)
  for( i in 1:numOfAttempts ) {
    vectorOfIdsx[i] <- exploreMatrix(methodName, scale, numOfElements, gradeOfIncomplete, numOfAttemptsForOneMatrix, alfa, beta)
  }
  
  incompleteIdx <- mean(vectorOfIdsx)
  incompleteIdx
}


# monte carlo dla wielu metod, ale na tych samych macierzach
monteCarloOnTheSameMatrix <- function(numOfElements, scale, gradeOfIncomplete, numOfAttempts, numOfAttemptsForOneMatrix, alfa=0, beta=0) {
  numOfMethods <- 17
  vectorOfIdsx <- integer(numOfAttempts)
  idxs = integer(numOfMethods)
  
  for( i in 1:numOfAttempts ) {
    idxs <- idxs + exploreMatrixOnTheSameMatrix(numOfElements, scale, gradeOfIncomplete, numOfAttemptsForOneMatrix, alfa, beta)
  }
  
  incompleteIdx <- idxs/numOfAttempts
  incompleteIdx
}

exploreMatrixOnTheSameMatrix <- function(numOfElements, scale, gradeOfIncomplete, numOfAttempts, alfa=0, beta=0) {
  numOfMethods <- 17;
  
  realIdx <- integer(numOfMethods)
  brokenIdx = integer(numOfMethods)
  differences <- integer(numOfMethods)
  
  matrix <- generateDistributedPCMatrix(numOfElements, scale)
  
  # oblicza prawidłowe wartości wsp. niesp.
  for(i in 1:numOfMethods) {
    realIdx[i] <- runMethod(i, matrix, alfa, beta)
  }

  # dekompletuje macierz i oblicza wszystkie współczynniki niespójności. Powtarza się to numOfAttemps razy
  for( i in 1:numOfAttempts ) {
    brokenMatrix <- breakPCMatrix(matrix, gradeOfIncomplete)
    for(j in 1:numOfMethods) {
      brokenIdx[j] <- brokenIdx[j] + abs(realIdx[j] - runMethod(j, brokenMatrix, alfa, beta))
    }
  }
  
  
  for(k in 1:numOfMethods) {
    brokenIdx[k] <- brokenIdx[k]/numOfAttempts
    differences[k] <- brokenIdx[k]/realIdx[k]*100
  }
  
  differences
}

runMethod <- function(i, matrix, alfa, beta){
  switch(i,
         "1"={
           koczkodaj(matrix)
         },
         "2"={
           grzybowski(matrix)
         },
         "3"={
           kazibudzkiLTI1(matrix)
         },
         "4"={
           kazibudzkiLTI2(matrix)
         },
         "5"={
           kazibudzkiCMLTI2(matrix)
         },
         "6"={
           pelaeLamata(matrix)
         },
         "7"={
           kulakowskiSzybowski(matrix)
         },
         "8"={
           kulakowskiSzybowski2(matrix)
         },
         "9"={
           kulakowskiSzybowskiIa(matrix, alfa)
         },
         "10"={
           kulakowskiSzybowskiIab(matrix, alfa, beta)
         },
         "11"={
           geometric(matrix)
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
         },
         "17"={
           saaty(matrix)
         }
  )
}