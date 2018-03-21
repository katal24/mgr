# przypisuje 0 jesli wartosc jest bliska 0
#' @export
chopV <- function(value){
  if(value>0.000001){
    return(value)
  } else{
    return(0)
  }
}

# przypisuje 0 do kazdego elementu macierzy bliskiego 0
chopM <- function(matrix){
  matrix[matrix<0.000001] <- 0
  return(matrix)
}

#################### Eigenvalue Rank Methods ###################
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

#' @title The largest principal eigenvalue (symbolic version)
#' @description A symbolic version of function principalEigenValue
#' @param matrix - PC matrix
#' @return the largest principal eigenvalue of matrix
#' @export
principalEigenValueSym <- function(matrix){
  max(Mod(eigen(matrix)$values))
}

#' @title Eigenvector of matrix
#' @description Returns the eigenvector of matrix corresponding to itd principal eigenvalue
#' @param matrix - PC matrix
#' @return the eigenvector of matrix corresponding to itd principal eigenvaluex
#' @export
principalEigenVector <- function(matrix){
  chopM(matrix);
  maxEigenVector <- Mod(eigen(apply(matrix, 2, as.numeric))$vectors)[,1]
  chopM(maxEigenVector)
}

#' @title Eigenvector of matrix (symbolic version)
#' @description A symbolic version of function principalEigenVector
#' @param matrix - PC matrix
#' @return the eigenvector of matrix corresponding to itd principal eigenvaluex
#' @export
principalEigenVectorSym <- function(matrix){
  Mod(eigen(matrix)$vectors)[,1]
}

#' @title Value of the Saaty Inconsistency Index
#' @description Returns the value of the Saaty Inconsistency Index computed for the matrix
#' @param matrix - PC matrix
#' @return the value of the Saaty Inconsistency Index computed for the matrix
#' @export
saatyIdx <- function(matrix){
  chopM(matrix)
  matrix <- apply(matrix, 2, as.numeric)
  n <- nrow(matrix)
  alpha <- principalEigenValueSym(matrix)
  chopV((alpha - n)/(n-1))
}

#' @title Value of the Saaty Inconsistency Index (symbolic version)
#' @description A symbolic version of function saatyIdxSym
#' @param matrix - PC matrix
#' @return the value of the Saaty Inconsistency Index computed for the matrix
#' @export
saatyIdxSym <- function(matrix){
  n <- nrow(matrix)
  alpha <- principalEigenValueSym(matrix)
  return((alpha - n)/(n-1))
}

#' @title Rescaled principal eigen vector of matrix
#' @description Returns the principal eigen vector of matrix rescaled in a way that the sum of its entries is 1
#' @param matrix - PC matrix
#' @return the rescaled principal eigen vector of matrix
#' @export
eigenValueRank <- function(matrix){
  chopM(matrix)
  matrix <- apply(matrix, 2, as.numeric)
  eigenVector <- chopM(principalEigenVectorSym(matrix));
  chopM(eigenVector/sum(eigenVector))
}

#' @title Rescaled principal eigen vector of matrix (symbolic version)
#' @description   A symbolic version of function eigenValueRank
#' @param matrix - PC matrix
#' @return the principal eigen vector of matrix rescaled in a way that the sum of its entries is 1
#' @export
eigenValueRankSym <- function(matrix){
  eigenVector <- principalEigenVectorSym(matrix);
  eigenVector/sum(eigenVector)
}

#' @title AHP ranking using eigenvalue based method
#' @description
#' Computes milticriteria ranking using eigenvalue based method using criteria matrix and alternatives matrixes.
#' This is basic three levels AHP.
#' @param M - n x n criteria matrix
#' @param ... - matrixes with the comparisons of alternatives with respect to the given criteria
#' @return multicriteria ranking using eigenvalue based method
#' @export ----
ahp <- function(M, ...){
  counter <- 0
  eigenVectorAc <- eigenValueRankSym(M)
  vectors <- lapply(list(...), eigenValueRank)
  verctors2 <- lapply(vectors, function(m){
    counter <<- counter+1
    m*eigenVectorAc[counter]
  })
  chopM(rowSums(data.frame(verctors2)))
}

#java
ahpFromVector <- function(M, matrices){
  counter <- 0
  eigenVectorAc <- eigenValueRankSym(M)
  vectors <- lapply(getListOfMatrices(matrices), eigenValueRank)
  verctors2 <- lapply(vectors, function(m){
    counter <<- counter+1
    m*eigenVectorAc[counter]
  })
  chopM(rowSums(data.frame(verctors2)))
}


################### Geometric mean rankings ###################
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

#' @title Rescaled rank list given as geometric means
#' @description Returns rank list given as geometric means of rows of the matrix rescaled in way
#' that the sum of its entries is 1
#' @param matrix - PC matrix
#' @return Rescaled rank list given as geometric means of rows of the matrix
#' @export
geometricRescaledRank <- function(matrix){
  geometricMeanVector <- apply(matrix, 1, function(row){
    prod(row)^(1/length(row))
  })
  geometricMeanVector/sum(geometricMeanVector)
}



################### HRE - Rankings withe the reference value ####################

#' @title Element of matrix
#' @description Returns [r,c] element of matrix
#' @param matrix - PC matrix
#' @param r - row
#' @param c - column
#' @return [r,c] element of matrix
#' @export
getMatrixEntry <- function(matrix, r, c){
  element <- matrix[r,c]
  element
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

#' @title Delete rows
#' @description Delete rows from matrix
#' @param martix - PC matrix
#' @param listOfRows - indices of rows which need to be deleted
#' @return matrix without deleted rows
#' @export ----
deleteRows <- function(matrix, listOfRows){
  matrix <- matrix[-listOfRows,]
  matrix
}

#' @export
removeRows <- function(matrix, listOfRows){
  matrix <- matrix[-listOfRows,]
  matrix
}

#' @title Delete columns
#' @description Delete columns from matrix
#' @param martix - matrix
#' @param listOfColumns - indices of columns which need to be deleted
#' @return matrix without deleted columns
#' @export ----
deleteColumns <- function(matrix, listOfColumns){
  matrix <- matrix[,-listOfColumns]
  matrix
}

#' @export
removeColumns <- function(matrix, listOfColumns){
  matrix <- matrix[,-listOfColumns]
  matrix
}

#' @title Delete rows and columns
#' @description Delete columns and rows from matrix
#' @param martix - matrix
#' @param listOfRowsAndColumns - indices of rows columns which need to be deleted
#' @return matrix without deleted rows and columns
#' @export ----
deleteRowsAndColumns <- function(matrix, listOfRowsAndColumns){
  matrix <- matrix[-listOfRowsAndColumns,-listOfRowsAndColumns]
  matrix
}

#' @export
zerosIndices <- function(knownVector){
  which(knownVector == 0)
}

#pokaz
#' @export
nonZerosIndices <- function(knownVector){
  which(knownVector != 0)
}

#pokaz
#' @export
nonZerosValues <- function(knownVector){
  knownVector[knownVector != 0]
}

#' @export
zerosValues <- function(knownVector){
  knownVector[knownVector == 0]
}

#' @title Set diagonal of matrix
#' @description Sets diagonal of matrix to specifed value
#' @param matrix - matrix
#' @param value - value which need to be set on the diagonal
#' @return Matrix with the specifed value on the diagonal
#' @export
setDiagonal <- function(matrix, value){
  diag(matrix) <- value
  matrix
}

#' @title Matrix for the HRE method
#' @description HRE matrix together with HRE constant term vector forms the linear equation system Au=b
#' @param matrix - PC matrix
#' @param knownVector - vector of known alternatives and others are marked by value 0
#' @return Matrix A for the HRE method
#' @export
HREmatrix <- function(matrix, knownVector){
  nonZerosIndices <- nonZerosIndices(knownVector)
  A <- (-1/(dim(matrix)[1]-1))*(deleteRowsAndColumns(matrix,nonZerosIndices))
  A <- setDiagonal(A,1)
  A
}

#' @title constant term vector for the HRE method
#' @description HRE constant term vector together with HRE matrix forms the linear equation system Au=b
#' @param matrix - PC matrix
#' @param knownVector - vector of known alternatives and others are marked by value 0
#' @return constant term vector b for the HRE method
#' @export
HREconstantTermVector <- function(matrix, knownVector){
  nonZerosIndices <- nonZerosIndices(knownVector)
  zerosIndices <- zerosIndices(knownVector)
  matrixWithoutRows <- deleteRows(matrix, nonZerosIndices)
  matrixWithoutRowsandColumns <- deleteColumns(matrixWithoutRows, zerosIndices)
  nonZerosValues <- nonZerosValues(knownVector)
  vectorB <- (1/(dim(matrix)[1]-1))* (matrixWithoutRowsandColumns %*% nonZerosValues)
  vectorB
}

#' @title Unknown HRE values
#' @description Counts values for the unknown alternatives
#' @param matrix - PC matrix
#' @param knownVector - vector of known alternatives and others are marked by value 0
#' @return values for the unknown alternatives
#' @export
HREpartialRank <- function(matrix, knownVector){
  solve(HREmatrix(matrix, knownVector), HREconstantTermVector(matrix, knownVector))
}

#' @export
createFullRank <- function(knownVector, partialRank){
  zerosIndices <- zerosIndices(knownVector)
  for(i in zerosIndices){
    knownVector[i] <- head(partialRank,1)
    partialRank <- partialRank[-1]
  }
  knownVector
}

#' @title full rank HRE
#' @description Counts values for the unknown alternatives and add their to known alternatives
#' @param matrix - PC matrix
#' @param knownVector - vector of known alternatives and others are marked by value 0
#' @return values for the both: known and unknown alternatives
#' @export
HREfullRank <- function(matrix, knownVector){
  partialRank <- HREpartialRank(matrix, knownVector);
  fullVector <- createFullRank(knownVector, partialRank)
  fullVector
}


#' @title Rescaled HRE full rank
#' @description Counts the rescaled HRE full rank list. Rescaled in the entries sum up to 1
#' @param matrix - PC matrix
#' @param knownVector - vector of known alternatives and others are marked by value 0
#' @return the rescaled HRE full rank list
#' @export
HRErescaledRank <- function(matrix, knownVector){
  fullRank <- HREfullRank(matrix,knownVector);
  chopM(fullRank/sum(fullRank))
}



################### HRE Geometric - Rankings with the reference values ###################

numberOfUnknownValues <- function(knownVector){
  length(zerosValues(knownVector))
}

#' @title Matrix for the HRE geometric method
#' @description HRE matrix together with HRE constant term vector forms the linear equation system Au=b
#' @param matrix - PC matrix
#' @param knownVector - vector of known alternatives and others are marked by value 0
#' @return Matrix A for the HRE method
#' @export
HREgeomMatrix <- function(matrix, knownVector){
  dimOfMatrix <- dim(matrix)[1]
  numberOfUnknownValues <- numberOfUnknownValues(knownVector)
  A <- matrix(-1,numberOfUnknownValues,numberOfUnknownValues)
  diag(A) <- (dimOfMatrix-1)
  A
}

#' @title constant term vector for the HRE geometric method
#' @description HRE constant term vector together with HRE matrix forms the linear equation system Au=b
#' @param matrix - PC matrix
#' @param knownVector - vector of known alternatives and others are marked by value 0
#' @return constant term vector b for the HRE method
#' @export
HREgeomConstantTermVector <- function(matrix, knownVector){
  nonZerosIndices <- nonZerosIndices(knownVector)
  nonZerosValues <- nonZerosValues(knownVector)
  zerosIndices <- zerosIndices(knownVector)
  matrixWithoutRows <- deleteRows(matrix, nonZerosIndices)
  matrixWithoutRows
  log10(apply(matrixWithoutRows, 1, prod)*prod(nonZerosValues))
}

#' @title The base of unknown HRE values
#' @description Counts values for the unknown alternatives befor they are raised by the power 10
#' @param matrix - PC matrix
#' @param knownVector - vector of known alternatives and others are marked by value 0
#' @return vector of values for the unknown alternatives
#' @export
HREgeomIntermediateRank <- function(matrix, knownVector){
  solve(HREgeomMatrix(matrix, knownVector), HREgeomConstantTermVector(matrix, knownVector))
}

#' @title Unknown HRE values
#' @description Counts values for the unknown alternatives using HRE geometric method
#' @param matrix - PC matrix
#' @param knownVector - vector of known alternatives and others are marked by value 0
#' @return vector of values for the unknown alternatives
#' @export
HREgeomPartialRank <- function(matrix, knownVector){
  sapply(HREgeomIntermediateRank(matrix, knownVector), FUN=function(x) 10^x)
}

#' @title full rank HRE using geometric method
#' @description Counts values for the unknown alternatives and add their to known alternatives
#' @param matrix - PC matrix
#' @param knownVector - vector of known alternatives and others are marked by value 0
#' @return values for the both: known and unknown alternatives
#' @export
HREgeomFullRank <- function(matrix, knownVector){
  partialRank <- HREgeomPartialRank(matrix, knownVector);
  fullVector <- createFullRank(knownVector, partialRank)
  fullVector
}

#' @title Rescaled HRE full rank using geometric method
#' @description Counts the rescaled HRE full rank list. Rescaled in the entries sum up to 1
#' @param matrix - PC matrix
#' @param knownVector - vector of known alternatives and others are marked by value 0
#' @return the rescaled HRE full rank list
#' @export
HREgeomRescaledRank <- function(matrix, knownVector){
  fullRank <- HREgeomFullRank(matrix,knownVector);
  chopM(fullRank/sum(fullRank))
}


################### Koczkodaj inconsistency methods ###################

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


#' @title The most inconsistent triad
#' @description Counts the worst triad in matrix according to Koczkodaj inconsistency criterion
#' @param matrix - PC Matrix
#' @return the worst triad in matrix according to Koczkodaj inconsistency criterion
#' @export
koczkodajTheWorstTriad <- function(matrix){
  sortMatrix <- sortKoczkodajMatrix(matrix)
  if(is.null(dim(sortMatrix))){
    sortMatrix  
  } else {
    sortMatrix[1,]
  }
}

#' @title The most inconsistent triads
#' @description Counts the worst triads in matrix according to Koczkodaj inconsistency criterion
#' @param matrix - PC Matrix
#' @return n worst triads in matrix according to Koczkodaj inconsistency criterion
#' @export  ---
koczkodajTheWorstTriads <- function(matrix, n){
  sortMatrix <- sortKoczkodajMatrix(matrix)
  head(sortMatrix,n)
}

#' @title Koczkodaj inconsistency
#' @description Counts the value of Koczkodaj inconsistency for the matrix
#' @param matrix - PC matrix
#' @return Koczkodaj inconsistency
#' @export ----
koczkodajIdx <- function(matrix){
  sortMatrix <- sortKoczkodajMatrix(matrix)
  if(is.null(dim(sortMatrix))){
    chopV(getTriadIdx(sortMatrix))  
  } else {
    chopV(getTriadIdx(sortMatrix[1,]))
  }
}

getTriad <- function(triadRow){
  c(triadRow[4],triadRow[5],triadRow[6])
}

getTuple<- function(triadRow){
  c(triadRow[1],triadRow[2],triadRow[3])
}

getTriadIdx <- function(triadRow){
  triadRow[7]
}

sortKoczkodajMatrix <- function(M){
  matrixFullKoczkodaj <- uniqueTriadsTuplesAndIdxForInComplete(M)
  matrixFullKoczkodaj[order(matrixFullKoczkodaj[,7],decreasing = TRUE),]
}

#' @title Consistent triad
#' @description Counts the closest consistenr triad
#' @param triad - vector of three values
#' @return the closest consistend triad
#' @export
koczkodajConsistentTriad <- function(triad){
  q <- triad[1]
  r <- triad[2]
  s <- triad[3]
  c(q^(2/3)*r^(-1/3)*s^(1/3),
    q^(-1/3)*r^(2/3)*s^(1/3),
    q^(1/3)*r^(1/3)*s^(2/3))
}

#' @title More consistent matrix
#' @description Counts improved (more consistent) matrix in which the most inconsistent triad is replaced by tje closest consistent one
#' @param matric - PC matrix
#' @return more consistent PC matrix
#' @export
koczkodajImprovedMatrixStep <- function(matrix){
  worstTriad <- koczkodajTheWorstTriad(matrix)
  improvedTriad <- koczkodajConsistentTriad(getTriad(worstTriad))
  mik = getTuple(worstTriad)[-3]
  mkj = getTuple(worstTriad)[-1]
  mij = getTuple(worstTriad)[-2]
  
  mikval = improvedTriad[1]
  mkjval = improvedTriad[2]
  mijval = improvedTriad[3]
  
  matrix[t(mij)] <- mijval
  matrix[mij[2],mij[1]] <- 1/mijval
  matrix[t(mkj)] <- mkjval
  matrix[mkj[2],mkj[1]] <- 1/mkjval
  matrix[t(mik)] <- mikval
  matrix[mik[2],mik[1]] <- 1/mikval
  matrix
}



################### Aggregation of Individual Judgments (AIJ) ###################
#' @title Aggregation matrix/vector (aritmetic means)
#' @description Computes aggregation matrix/vector whose elements are aritmethic means
#' of elements from input matrices/vectors (whose need to be the same dimension)
#' @param ... - matrices orvectors the same dimension
#' @return aggregation matrix/vector
#' @export
AIJadd <- function(...){
  matrixes <- list(...)
  mean <- Reduce('+', matrixes)/length(matrixes)
  mean
}

#' @title Aggregation matrix/vector (geometric means)
#' @description Computes aggregation matrix/vector whose elements are aritmethic means
#' of elements from input matrices/vectors (whose need to be the same dimension)
#' @param ... - matrices or vectors the same dimension
#' @return aggregation matrix/vector
#' @export
AIJgeom <- function(...){
  matrixes <- list(...)
  meanGeom <- Reduce('*', matrixes)^(1/length(matrixes))    # Reduce - wywoluje funkcje po kolei na elementach listy
  meanGeom
}

#java
AIJgeomFromVector<- function(matrixes){
  matrixes1 <- getListOfMatrices(matrixes)
  meanGeom <- Reduce('*', matrixes1)^(1/length(matrixes1))    # Reduce - wywoluje funkcje po kolei na elementach listy
  meanGeom
}
#java
AIJaddFromVector<- function(matrixes){
  matrixes1 <- getListOfMatrices(matrixes)
  mean <- Reduce('+', matrixes1)/length(matrixes1)
  mean
}
#java
AIJvectorsAddFromVector<- function(vectors, len){
  vectors <- getListOfVectors(len, vectors)
  mean <- Reduce('+', vectors)/length(vectors)
  mean
}
#java
AIJvectorsGeomFromVector<- function(vectors, len){
  vectors <- getListOfVectors(len, vectors)
  meanGeom <- Reduce('*', vectors)^(1/length(vectors))    # Reduce - wywoluje funkcje po kolei na elementach listy
  meanGeom
}

getListOfMatrices <- function(matrix){
  matrices <- list()
  dim <- dim(matrix)[2]
  dim1 <- dim(matrix)[1]
  count <- 1;
  for(i in 1:dim1){
    if(i %% dim == 0 ){
      matrices[count] <- list(matrix[(i-dim+1):i,])
      count <- count+1;
    }
  }
  matrices
}

getListOfVectors <- function(len, vs){
  vectors <- list()
  dim <- length(vs)
  
  count <- 1;
  for(i in 1:dim){
    if(i %% len == 0 ){
      vectors[count] <- list(vs[(i-len+1):i])
      count <- count+1;
    }
  }
  vectors
}
################### Incomplete pairwise compairsons matrices ###################

#' @title Number of zeros elements in row
#' @description Counts number of zeros elements in given row in matrix
#' @param matrix - matrix
#' @param row - number of row which will be checked
#' @return number of zeros elements in row
#' @export
harkerMatrixPlaceHolderCount <- function(matrix, row){
  length(which(matrix[row,] != 0))   # which - zwraca elementy spelniajace watrunek
}

#' @title Harher matrix
#' @description Computes harker matrix A for a given incomplete PC matrix.
#' Then matrix A is ready to use in eigenvalue based method.
#' Empty and incorrect values are replaced by 0 and diagonal values are replaced by 1.
#' @param matrix - matrix
#' @return matrix which is ready to use in eigenvalue based method
#' @export
harkerMatrix <- function(matrix){
  matrix[matrix < 0] <- 0
  diag(matrix) <- 1
  matrix
}


#################### Ranking discrepancy ###################

#' @title Ranking discrepancy
#' @description Computes discrepancy for each element of matrix
#' and create matrix containing entries in form e_ij = m_ji*(mju_i/mju_j).
#' When the M matrix is consistent every r_ij equals 1.
#' @param M - matrix
#' @param mju - ranking of matrix
#' @return error matrix E = [e_ij]
#' @export
errorMatrix <- function(matrix, mju){
  E <- matrix
  for(i in 1:dim(matrix)[1])
    for(j in 1:dim(matrix)[2])
      E[i,j] = matrix[j,i]*mju[i]/mju[j]
    # E = outer(1:nrow(mm), 1:ncol(mm), FUN=function(r,c) {print(r); print(c); mm[r,][c]})
    E
}

#' @title Local discrepancy
#' @description Compute matrix with entries d_ij = max{e_ij-1, 1/e_ij-1}
#' @param matrix - matrix
#' @param mju - ranking of matrix
#' @return matrix with locals discrepancy
#' @export
localDiscrepancyMatrix <- function(matrix, mju){
  errorMatrix <- errorMatrix(matrix,mju)
  for(i in 1:dim(errorMatrix)[1])
    for(j in 1:dim(errorMatrix)[2])
      errorMatrix[i,j] = max( errorMatrix[i,j]-1, 1/errorMatrix[i,j] -1)
    chopM(errorMatrix)
}

#' @title Global discrepancy
#' @description Finds maximal entry of local discrepancy matrix
#' @param matrix - matrix
#' @param mju - ranking of matrix
#' @return maximal value of local discrepancy
#' @export
globalDiscrepancy <- function(matrix, mju){
  max(localDiscrepancyMatrix(matrix, mju))
}



########## Condition of order preservation ##########
# tej jedej nie robie - bo to lista
#' @title List of COP1
#' @description counts the list of indices that violate the first Condition of Order Preservation
#' postulate (formulated by Bana e Costa and Vansnick)
#' @param matrix - PC matrix
#' @param resultList - ranking of matrix
#' @return the matrix of indices that violate the first Condition of Order Preservation (every row includes one pair)
#' @export
cop1ViolationList <- function(matrix, resultList){
  errorList <- matrix(c(0,0),1)
  for(i in 1:dim(matrix)[1]){
    for(j in 1:dim(matrix)[2]){
      element <- matrix[i,j]
      if(i!=j && (element > 1 && resultList[i]<resultList[j]) || (element < 1 && resultList[i]>resultList[j])){
        errorList <- rbind(errorList, c(i,j))
      }
    }
  }
  errorList[-1,]
}

# tego jednego nie mam - bo to boolean
#' @title check COP1
#' @description checks if the first Condition of Order Preservation  is fulfilled
#' @param matrix - PC matrix
#' @param resultList - ranking of matrix
#' @return true if the first Condition of Order Preservation  is fulfilled, else false
#' @export
cop1Check <- function(matrix, resultList){
  length(cop1ViolationList(matrix, resultList)) == 0
}

# ## JAVA
#' @title check COP1Details
#' @description checks if the first Condition of Order Preservation  is fulfilled
#' @param matrix - PC matrix
#' @param resultList - ranking of matrix
#' @return true if the first Condition of Order Preservation  is fulfilled, else false
#' @export
cop1CheckDetails <- function(matrix, resultList){
  if(length(cop1ViolationList(matrix, resultList)) == 0){
    return(1)
  } else {
    return(0);
  }
}

#tego jesdengo nie - bo to lista
#' @title List of COP2
#' @description counts the list of indices that violate the second Condition of Order Preservation
#' postulate (formulated by Bana e Costa and Vansnick)
#' @param matrix - PC matrix
#' @param resultList - ranking of matrix
#' @return the matrix of indices that violate the second Condition of Order Preservation (every row includes one pair)
#' @export
cop2ViolationList <- function(matrix, resultList){
  errorList <- matrix(c(0,0,0,0),1)
  # errorList <- list()
  rows <- dim(matrix)[1]
  columns <- dim(matrix)[2]
  for(i in 1:rows)
    for(j in 1:columns)
      for(k in 1:rows)
        for(l in 1:columns){
          if(i!=j && k!=l){
            a <- matrix[i,j]
            b <- matrix[k,l]
            if((a>b && resultList[i]/resultList[j] < resultList[k]/resultList[l]) || (a<b && resultList[i]/resultList[j] > resultList[k]/resultList[l])){
              #    errorList[[length(errorList)+1]] <- list(c(i,j),c(k,l))
              errorList <- rbind(errorList, c(i,j,k,l))
            }
          }
        }
  errorList
}


#' @title check COP2
#' @description checks if the second Condition of Order Preservation  is fulfilled
#' @param matrix - PC matrix
#' @param resultList - ranking of matrix
#' @return true if the second Condition of Order Preservation is fulfilled, else false
#' @export
cop2Check <- function(matrix, resultList){
  length(cop2ViolationList(matrix, resultList)) == 0
}

# ## JAVA
#' @title check COP2
#' @description checks if the second Condition of Order Preservation  is fulfilled
#' @param matrix - PC matrix
#' @param resultList - ranking of matrix
#' @return true if the second Condition of Order Preservation is fulfilled, else false
#' @export
cop2CheckDetails <- function(matrix, resultList){
  if(length(cop2ViolationList(matrix, resultList)) == 0){
    return(1)
  } else {
    return(0);
  }
}




#################### Comparing ranks ####################
kendallTauDistanceForPair <- function(list1, list2, element1, element2){
  posInList1Elem1 <- match(element1, list1)
  posInList1Elem2 <- match(element2, list1)
  posInList2Elem1 <- match(element1, list2)
  posInList2Elem2 <- match(element2, list2)
  if(sign(posInList1Elem1 - posInList1Elem2) == sign(posInList2Elem1 - posInList2Elem2)){
    return(0)
  } else{
    return(1)
  }
}

kendallTauDistanceForPairAsVector <- function(list1, list2, c){
  kendallTauDistanceForPair(list1, list2, c[1], c[2])
}

#' @title Kendall Tau distance for two vectors
#' @description Computes Kendall (bubble sort) distance between two rank vectors
#' @param list1 - first rank to compare
#' @param list2 - second rank to compare
#' @return number of swap which need to be make to order of elements will be identical
#' @export
kendallTauDistance <- function(list1, list2){
  sum(apply(combn(list1,2),2,kendallTauDistanceForPairAsVector,list1=list1, list2=list2))
}

#' @title Normalized Kendall Tau distance for two vectors
#' @description Computes Kendall (bubble sort) distance between two rank vectors and divited it by numbers of all possible pairs
#' @param list1 - first rank to compare
#' @param list2 - second rank to compare
#' @return proportion of Kendall Tau distance to numbers of all possible pairs
#' @export
normalizedKendallTauDistance <- function(list1,list2){
  n <- length(list1)
  kendallTauDistance <- kendallTauDistance(list1,list2)
  kendallTauDistance/(n*(n-1)/2)
}



########## Auxiliary functions ##########

#' @title Consistent matrix from rank
#' @description Creates consistent PC matrix from the rank vector V so that m_ij = v[i]/v[j]
#' @param rankList - ranking of matrix
#' @return  - full PC matrix
#' @export
consistentMatrixFromRank <- function(rankList){
  length <- length(rankList)
  m <- matrix(, nrow=length, ncol=length)
  for(i in 1:length)
    for(j in 1:length)
      m[i,j] <- rankList[i]/rankList[j]
  m
}


#' @title Sort rank
#' @description Sorts elements of rank
#' @param rankList - ranking of matrix
#' @return ranking - sorted rank
#' @export
rankOrder <- function(rankList){
  sort(rankList, decreasing=TRUE)
}

###################
# Anonymous function syntax
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
distributePCMatrix<- function(matrix){
  rDim = dim(matrix)[1]
  cDim = dim(matrix)[2]
 
  for(r in 1:rDim-1)
    for(c in (r+1):cDim)
      matrix[r,c] <- runif(1, min = 0.5, max = 1.5) * matrix[r,c]
  
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
generateDistributedPCMatrix<- function(numOfElements){
  matrix <- generatePCMatrix(numOfElements)
  matrix <- distributePCMatrix(matrix)
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



########################### Methods for incomlete matrices #########################

#' @title Koczkodaj incomplete matrix inconsistency
#' @description Counts the value of inconsistency for the matrix
#' @param matrix - PC matrix
#' @return inconsistency
#' @export ----
koczkodaj <- function(matrix, alfa=0, beta=0){
  triadsAndIdxs <- uniqueTriadsTuplesAndIdxForInComplete("koczkodajForTriad", matrix)
  sortMatrix <- triadsAndIdxs[order(triadsAndIdxs[,4],decreasing = TRUE),]
  if(is.null(dim(sortMatrix))){
    chopV(getRealTriadIdx(sortMatrix))  
  } else {
    chopV(getRealTriadIdx(sortMatrix[1,]))
  }
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
grzybowski <- function(matrix, alfa=0, beta=0){
  triadsAndIdxs <- uniqueTriadsTuplesAndIdxForInComplete("koczkodajForTriad", matrix)
  chopV(avg(triadsAndIdxs[,4]))
}


#### method !!!
kazibudzkiLTI1 <- function(matrix, alfa=0, beta=0){
  triadsAndIdxs <- uniqueTriadsTuplesAndIdxForInComplete("grzybowskiLTI1ForTriad", matrix)
  v = chopV(avg(triadsAndIdxs[,4]))
  v
}

kazibudzkiLTI1ForTriad <- function(triad) {
  abs(log((triad[1]*triad[2])/triad[3]))
}

kazibudzkiLTI2 <- function(matrix, alfa=0, beta=0){
  triadsAndIdxs <- uniqueTriadsTuplesAndIdxForInComplete("grzybowskiLTI2ForTriad", matrix)
  v = chopV(avg(triadsAndIdxs[,4]))
  v
}

kazibudzkiLTI2ForTriad <- function(triad) {
  log((triad[1]*triad[2])/triad[3])**2
}


#### method !!!
kazibudzkiCMLTI2 <- function(matrix, alfa=0, beta=0){
  triadsAndIdxs <- uniqueTriadsTuplesAndIdxForInComplete("kazibudzkiCMLTI2ForTriad", matrix)
  chopV(avg(triadsAndIdxs[,4])/(1+max(triadsAndIdxs[,4])))
}

kazibudzkiCMLTI2ForTriad <- function(triad) {
  log((triad[1]*triad[2])/triad[3])**2
}


#### method !!!
pelaeLamata <- function(matrix, alfa=0, beta=0){
  triadsAndIdxs <- uniqueTriadsTuplesAndIdxForInComplete("pelaeLamataForTriad", matrix)
  chopV(avg(triadsAndIdxs[,4]))
}

pelaeLamataForTriad <- function(triad) {
  triad[3]/(triad[1]*triad[2]) + (triad[1]*triad[2])/triad[3] - 2 
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
kulakowskiSzybowski2 <- function(matrix, alfa=0, beta=0){
  triadsAndIdxs <- uniqueTriadsTuplesAndIdxForInComplete("koczkodajForTriad", matrix)
  chopV((sqrt(sum(triadsAndIdxs[,4]**2)))/length(triadsAndIdxs[,4]))
}

#### method !!!
geometric <- function(matrix, alfa=0, beta=0){
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

sumValuesInColumns <- function(matrix){
  apply(matrix, 2, function(col){
    col[col==0]=1
    sum(col)
  })
}

#### method !!!
saaty <- function(matrix){
  matrix[matrix==0] = 1
  matrix <- apply(matrix, 2, as.numeric)
  n <- nrow(matrix)
  alpha <- principalEigenValueSym(matrix)
  chopV((alpha - n)/(n-1))
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




######################### Monte Carlo ###########################

# count average of vector
avg <- function(vector){
  sum(vector)/length(vector)
}

# generuje macierz zaburzoną, a potem przeprowadza x (numOfAttempts) prób zdekompletowania i liczy średnią wartość niespójności,
# potem znajduje różnicę tej średniej od prawdziwej wartości
exploreMatrix <- function(methodName, numOfElements, gradeOfIncomplete, numOfAttempts, alfa=0, beta=0) {
  matrix <- generateDistributedPCMatrix(numOfElements)
  if(alfa==0 && beta==0){
    realIdx <- methodName(matrix)
  } else {
    realIdx <- methodName(matrix, alfa, beta)
  }
  vectorOfIdsx <- integer(numOfAttempts)
  
  for( i in 1:numOfAttempts ) {
    brokenMatrix <- breakPCMatrix(matrix, gradeOfIncomplete)
    #idx <- methodName(brokenMatrix)
    if(alfa==0 && beta==0){
      idx <- methodName(brokenMatrix)
    } else {
      idx <- methodName(brokenMatrix, alfa, beta)
    }
    vectorOfIdsx[i] <- idx
  }
  
  incompleteIdx <- avg(vectorOfIdsx)
  abs((realIdx - incompleteIdx)/realIdx*100)
}

# Wykonuje exploreMatrixKoczkodaj x (numOfAttempts) razy i z tego bierze średnią
monteCarlo <- function(methodName, numOfElements, gradeOfIncomplete, numOfAttempts, numOfAttemptsForOneMatrix, alfa=0, beta=0) {
  vectorOfIdsx <- integer(numOfAttempts)
  for( i in 1:numOfAttempts ) {
    idx <- exploreMatrix(methodName, numOfElements, gradeOfIncomplete, numOfAttemptsForOneMatrix, alfa, beta)
    vectorOfIdsx[i] <- idx
  }
  
  incompleteIdx <- avg(vectorOfIdsx)
  incompleteIdx
}


# monte carlo dla wielu metod, ale na tych samych macierzach
monteCarloOnTheSameMatrix <- function(numOfElements, gradeOfIncomplete, numOfAttempts, numOfAttemptsForOneMatrix, alfa=0, beta=0) {
  vectorOfIdsx <- integer(numOfAttempts)
  idxs = integer(10)
  for( i in 1:numOfAttempts ) {
    idxs <- idxs + exploreMatrixOnTheSameMatrix(numOfElements, gradeOfIncomplete, numOfAttemptsForOneMatrix, alfa, beta)
  }
  
  incompleteIdx <- idxs/numOfAttempts
  incompleteIdx
}

exploreMatrixOnTheSameMatrix <- function(numOfElements, gradeOfIncomplete, numOfAttempts, alfa=0, beta=0) {
  matrix <- generateDistributedPCMatrix(numOfElements)
  methods <- c(koczkodaj, grzybowskiATI, grzybowskiALTI1, grzybowskiALTI2, kazibudzkiCMLTI2, pelaeLamata, kulakowskiSzybowskiIa, kulakowskiSzybowskiIa, kulakowskiSzybowskiIab, kulakowskiSzybowskiIab)
  realIdx <- integer(10)
  
  for(i in 1:10) {
    realIdx[i] <- runMethod(i, matrix, alfa, beta)
    #realIdx[i] <- methods[i](matrix, alfa, beta)
  }

  brokenIdx = integer(10)
  idxs <- integer(10)
  
  for( i in 1:numOfAttempts ) {
    brokenMatrix <- breakPCMatrix(matrix, gradeOfIncomplete)
    for(i in 1:10) {
      brokenIdx[i] <- brokenIdx[i] + runMethod(i, brokenMatrix, alfa, beta)
    }
  }
  
  for(i in 1:10) {
    brokenIdx[i] <- brokenIdx[i]/numOfAttempts
    idxs[i] <- abs((realIdx[i] - brokenIdx[i])/realIdx[i]*100)
  }
  
  idxs
}

runMethod <- function(i, matrix, alfa, beta){
  switch(i,
         "1"={
           koczkodaj(matrix, alfa, beta)
         },
         "2"={
           grzybowski(matrix, alfa, beta)
         },
         "3"={
           kazibudzkiLTI1(matrix, alfa, beta)
         },
         "4"={
           kazibudzkiLTI2(matrix, alfa, beta)
         },
         "5"={
           kazibudzkiCMLTI2(matrix, alfa, beta)
         },
         "6"={
           pelaeLamata(matrix, alfa, beta)
         },
         "7"={
           kulakowskiSzybowskiIa(matrix, 0.3, 0.7)
         },
         "8"={
           kulakowskiSzybowskiIa(matrix, 0.8, 0.2)
         },
         "9"={
           kulakowskiSzybowskiIab(matrix, 0.3, 0.4)
         },
         "10"={
           kulakowskiSzybowskiIab(matrix, 0.1, 0.3)
         }
  )
}

test <- function() {
  b=15;
  numOfMatrices=100;
  numOfAttempts=100;
  for(n in c(4,5,6,7,8,10,15,20,50)){
    print("=====================================")
    print(n)
    print(monteCarlo(koczkodaj,n,b,numOfMatrices,numOfAttempts))
    print(monteCarlo(grzybowskiATI,n,b,numOfMatrices,numOfAttempts))
    print(monteCarlo(grzybowskiALTI1,n,b,numOfMatrices,numOfAttempts))
    print(monteCarlo(grzybowskiALTI2,n,b,numOfMatrices,numOfAttempts))
    print(monteCarlo(kazibudzkiCMLTI2,n,b,numOfMatrices,numOfAttempts))
    print(monteCarlo(pelaeLamata,n,b,numOfMatrices,numOfAttempts))
  }
}


testOnTheSameMatrix <- function(b, numOfMatrices, numOfAttempts) {
  # for(n in c(4,5,6,7,8,10,15,20,50)){
  for(n in c(20,50)){
    print("=====================================")
    print(n)
    print(monteCarloOnTheSameMatrix(n, b, 100, 100, 0.4, 0.6))
  }
}

test2OnTheSameMatrix <- function(n, numOfMatrices, numOfAttempts) {
  # for(n in c(4,5,6,7,8,10,15,20,50)){
  for(b in c(4,5,8,10,20)){
    print("=====================================")
    print(b)
    print(monteCarloOnTheSameMatrix(n, b, 100, 100, 0.4, 0.6))
  }
}

