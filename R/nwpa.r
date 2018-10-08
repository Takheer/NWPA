#' nwpa package provides pairwise alignment for two nucleotide
#' sequences. It also can help you to parse .fas file into a matrix of strings
#' 
#'
#' @name nwpa
#' @docType package


#' This function is for parsing .fas file into a matrix.
#' It returns a matrix with two columns: label and sequence.
#' 
#' @export
#' @param text a .fas file with unaligned sequences
read.fas <- function(text) {
  raw <- readLines(text)
  result <- 0

  # Delete all spaces in "raw"
  gsub(" ", "", raw)

  # forming the matrix
  i <- 1
  while (i < length(raw)){
    label <- raw[i]
    i <- i + 1
    sequence <- ""
    while (!(substr(raw[i], 1, 1) == ">") & (i <= length(raw))) {
      sequence <- stringr::str_c(sequence, raw[i])
      i <- i + 1
    }

    # In this place labels and sequences must be united into a matrix.
    if (length(result) == 1){
      result <- matrix(c(label, sequence), nrow = 1)
    }
    else{
      t <- matrix(c(label, sequence), nrow = 1)
      result <- rbind(result, t)
    }
  }
  return(result)
}

#' Used for writing aligned sequences into .fas file
#' 
#' @export
#' @param aligned a matrix of sequences in format "label - sequence"
#' @param file output file in .fas format
write.fas <- function(aligned, file){
  result <- ""
  for (i in c(1:2)){
    result <- stringr::str_c(result, aligned[i, 1], "\n")
    for (j in c(1:stringr::str_length(aligned[i,2]))){
      if (j %% 60 == 0){
        result <- stringr::str_c(result, "\n")
      }
      else {
        result <- stringr::str_c(result, substr(aligned[i, 2], j, j))
      }
    }
    result <- stringr::str_c(result, "\n")
  }
  return(writeLines(result, file))
}


#' Calculates similarity score. It supports gapped alignments (The gap 
#' character will just always count as "not a match") and sequences with 
#' ambiguous base count
#' 
#' @param char.a first nucleotide for evaluation
#' @param char.b second nucleotide for evaluation
#' @param mismatch used to evaluate a cell of similarity matrix in case of 
#' matching of two nucleotides
#' @param gap used to evaluate a cell of similarity matrix in case when you
#' should add a gap to your alignment
score <- function(char.a=" ", char.b=" ", match, mismatch) {
  # I believe I can do it in not such a bad way
  vec <- function(c="") {
    if(c == "A"){
      return(c("A"))
    } else if(c == "C"){
      return(c("C"))
    } else if(c == "G"){
      return(c("G"))
    } else if(c == "T" | c == "U"){
      return(c("T"))
    } else if(c == "N"){
      return(c("A", "C", "G", "T"))
    } else if(c == "R"){
      return(c("A", "G"))
    } else if(c == "Y"){
      return(c("C","T"))
    } else if(c == "S"){
      return(c("C", "G"))
    } else if(c == "W"){
      return(c("A", "T"))
    } else if(c == "K"){
      return(c("G", "T"))
    } else if(c == "M"){
      return(c("A", "C"))
    } else if(c == "B"){
      return(c("C", "G", "T"))
    } else if(c == "D"){
      return(c("A", "G", "T"))
    } else if(c == "H"){
      return(c("A", "C", "T"))
    } else if(c == "V"){
      return(c("A", "C", "G"))
    }
  }

  if(char.a == "-" | char.b == "-" | char.a == "." | char.b == ".") {
    result <- mismatch
  } else {
    v1 <- vec(char.a)
    v2 <- vec(char.b)
    result <- ifelse(length(intersect(v1, v2)) == 0, mismatch, match)
  }
  return(result)
}

#' This function is used to calculate matrix for finding best-scoring alignment
#' 
#' @param first a sequence to align
#' @param second another sequence
#' @param match used to evaluate a cell of similarity matrix in case of 
#' matching of two nucleotides
#' @param mismatch used to evaluate a cell of similarity matrix in case of 
#' matching of two nucleotides
#' @param gap used to evaluate a cell of similarity matrix in case when you
#' should add a gap to your alignment
grid <- function(first, second, match=1, mismatch=0, gap=-1) {
  gr <- matrix(0, ncol = stringr::str_length(second) + 1, nrow = stringr::str_length(first) + 1)

  # prior assignment of the grid
  j <- 1
  while (j <= stringr::str_length(second) + 1) {
    gr[1, j] <- ((j - 1) * gap)
    j <- j + 1
  }
  i <- 1
  while (i <= stringr::str_length(first) + 1) {
    gr[i, 1] <- ((i - 1) * gap)
    i <- i + 1
  }

  # this code used for printing current state of calculations
  total <- stringr::str_length(first) * stringr::str_length(second)
  cat(stringr::str_c("building matrix: ", total, " iterations in total\n"))
  p <- total %/% 10 # ten percent of total nucleotides
  count <- 0
  d <- 1 # iterators for current state printing
  # main assignment of the grid
  for (i in c(2:(stringr::str_length(first) + 1))) {
    for (j in c(2:(stringr::str_length(second) + 1))) {
      match0 <- gr[i - 1, j - 1] + score(substr(first, i, i), substr(second, j, j), match, mismatch)
      delete <- gr[i - 1, j] + gap
      insert <- gr[i, j - 1] + gap
      gr[i, j] <- max(match0, delete, insert)
      count <- count + 1
      if (count == p * d){
        cat(stringr::str_c(d * 10, "%\n")) # prints a state for every 10 percent of analysed nucleotides
        d <- d + 1
      }
    }
  }

  return(gr)
}


#' This function provides pairwise alignment with simple similarity matrix
#' Notice that the length of your sequence's name should be one string 
#' 
#' @export
#' @param input a file with .fas extension. These are your DNA sequences. 
#' The function will align only first two sequences
#' @param match used to evaluate a cell of similarity matrix in case of 
#' matching of two nucleotides
#' @param mismatch used to evaluate a cell of similarity matrix in case of 
#' matching of two nucleotides
#' @param gap used to evaluate a cell of similarity matrix in case when you
#' should add a gap to your alignment
align <- function(input, output, match=1, mismatch=0, gap=-1){
  # reading sequences
  sequences <- read.fas(input)
  seq.a <- stringr::str_c("-", sequences[1, 2])
  seq.b <- stringr::str_c("-", sequences[2, 2])
  count <- 0

  # matrix of scores
  grid <- grid(seq.a, seq.b, match, mismatch, gap)

  # the algorithm
  alignment.a <- ""
  alignment.b <- ""
  i <- stringr::str_length(seq.a)
  j <- stringr::str_length(seq.b)
  while ((i > 0) & (j > 0)) {
    score <- grid[i + 1, j + 1]
    score.diag <- grid[i, j]
    score.up <- grid[i + 1, j]
    score.left <- grid[i, j + 1]
    if (score == score.diag + score(substr(seq.a, i + 1, i + 1), substr(seq.b, j + 1, j + 1), match, mismatch)){
      alignment.a <- stringr::str_c(substr(seq.a, i + 1, i + 1), alignment.a)
      alignment.b <- stringr::str_c(substr(seq.b, j + 1, j + 1), alignment.b)
      i <- i - 1
      j <- j - 1
    }
    else if (score == score.left + gap){
      alignment.a <- stringr::str_c(substr(seq.a, i + 1, i + 1), alignment.a)
      alignment.b <- stringr::str_c("-", alignment.b)
      i <- i - 1
    }
    else if (score == score.up + gap){
      alignment.a <- stringr::str_c("-", alignment.a)
      alignment.b <- stringr::str_c(substr(seq.b, j, j), alignment.b)
      j <- j - 1
    }
  }

  while (i > 0){
    alignment.a <- stringr::str_c(substr(seq.a, i + 1, i + 1), alignment.a)
    alignment.b <- stringr::str_c("-", alignment.b)
    i <- i - 1
  }
  while (j > 0){
    alignment.a <- stringr::str_c("-", alignment.a)
    alignment.b <- stringr::str_c(substr(seq.b, j + 1, j + 1), alignment.b)
    j <- j - 1
  }

  cat(stringr::str_c("\n\nMaximum score: ", grid[stringr::str_length(seq.a), stringr::str_length(seq.b)], "\n\n\n"))
  r1 <- matrix(c(sequences[1,1], alignment.a), nrow=1, ncol=2)
  r2 <- matrix(c(sequences[2,1], alignment.b), nrow=1, ncol=2)
  result <- rbind(r1, r2)
  write.fas(result, output)
  if(output != 0){
    cat(stringr::str_c("Aligned sequences have been written to ", getwd(), "/", output))
  }
  system(paste("open ", output))
}
NULL
