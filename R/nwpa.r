#' needleman.wunsch package provides pairwise alignment for two nucleotide
#' sequences. 
#'
#'
#' @name needleman.wunsch
#' @docType package


#' This is an internal function for parsing .fas file into a matrix 
#' with label and sequence
#' 
#' @param text a .fas file with unaligned sequences
read.fas <- function(text) {
  raw <- readLines(text)
  result <- 0

  # Delete all the spaces in "raw"
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
    if (result == 0){
      result <- matrix(c(label, sequence), nrow = 1)
    }
    else{
      t <- matrix(c(label, sequence), nrow = 1)
      result <- rbind(result, t)
    }
  }
  return(result)
}


#' You don't need to use it
#' 
#' @param char.a first nucleotide for evaluation
#' @param char.b second nucleotide for evaluation
#' 
score <- function(char.a=" ", char.b=" ", match, mismatch) {
  result <- ifelse(char.a == char.b, match, mismatch)
  return(result)
}

#' This function is used to calculate matrix finding best-scoring alignment
#' 
#' @param match used to evaluate a cell of similarity matrix in case of 
#' matching of two nucleotides
#' @param mismatch used to evaluate a cell of similarity matrix in case of 
#' matching of two nucleotides
#' @param gap used to evaluate a cell of similarity matrix in case when you
#' should add a gap to your alignment
grid <- function(first, second, match=1, mismatch=0, gap=-1) {
  # номер эл-та в матрице == номер эл-та в посл-ти + 1
  gr <- matrix(0, ncol = str_length(first) + 1, nrow = str_length(second) + 1)

  # prior assignment of the grid
  j <- 1
  while (j <= str_length(first) + 1) {
    gr[1, j] <- ((j - 1) * mismatch)
    j <- j + 1
  }
  i <- 1
  while (i <= str_length(second) + 1) {
    gr[i, 1] <- ((i - 1) * mismatch)
    i <- i + 1
  }

  # main assignment of the grid
  for (i in c(2:(str_length(first) + 1))) {
    for (j in c(2:(str_length(second) + 1))) {
      match0 <- gr[i - 1, j - 1] + score(substr(first, i, i), substr(second, j, j), match, mismatch)
      delete <- gr[i - 1, j] + gap
      insert <- gr[i, j - 1] + gap
      gr[i, j] <- max(match0, delete, insert)
    }
  }

  return(gr)
}


#' This function provides pairwise alignment with simple similarity matrix
#' Notice that the length of your sequence's name should be one string 
#' 
#' @export
#' @param file a file with .fas extension. These are your DNA sequences. The function will align only first two sequences
#' @param match used to evaluate a cell of similarity matrix in case of 
#' matching of two nucleotides
#' @param mismatch used to evaluate a cell of similarity matrix in case of 
#' matching of two nucleotides
#' @param gap used to evaluate a cell of similarity matrix in case when you
#' should add a gap to your alignment
align <- function(file, match=1, mismatch=0, gap=-1){
  
  # reading sequences
  sequences <- read.fas(file)
  seq.a <- sequences[1, 2]
  seq.b <- sequences[2, 2]

  grid <- grid(seq.a, seq.b, match, mismatch, gap)

  alignment.a <- ""
  alignment.b <- ""
  i <- str_length(seq.a) + 1
  j <- str_length(seq.b) + 1
  while ((i > 1) & (j > 1)) {
    score <- grid[i, j]
    score.diag <- grid[i - 1, j - 1]
    score.up <- grid[i, j - 1]
    score.left <- grid[i - 1, j]
    if (score == score.diag + score(substr(seq.a, i, i), substr(seq.b, j, j), match, mismatch))
    {
      alignment.a <- str_c(substr(seq.a, i, i), alignment.a)
      alignment.b <- str_c(substr(seq.b, j, j), alignment.b)
      i <- i - 1
      j <- j - 1
    }
    else if (score == score.left + gap)
    {
      alignment.a <- str_c(substr(seq.a, i, i), alignment.a)
      alignment.b <- str_c("-", alignment.b)
      i <- i - 1
    }
    else if (score == score.up + gap)
    {
      alignment.a <- str_c("-", alignment.a)
      alignment.b <- str_c(substr(seq.b, j, j), alignment.b)
      j <- j - 1
    }
  }

  while (i > 0)
  {
    alignment.a <- str_c(substr(seq.a, i, i), alignment.a)
    alignment.b <- str_c("-", alignment.b)
    i <- i - 1
  }
  while (j > 0)
  {
    alignment.a <- str_c("-", alignment.a)
    alignment.b <- str_c(substr(seq.b, j, j), alignment.b)
    j <- j - 1
  }

  print(alignment.a)
  print(alignment.b)
}
NULL