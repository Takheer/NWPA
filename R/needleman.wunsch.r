#' needleman.wunsch package provides pairwise alignment for two nucleotide
#' sequences. 
#'
#'
#' @name needleman.wunsch
#' @docType package


#' This is an internal function for parsing .fas file into a dataset 
#' with label and sequence
#' @export # should be removed before release
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
    print("Label: ")
    print(label)
    sequence <- ""
    while (!(substr(raw[i], 1, 1) == ">") & (i <= length(raw))) {
      sequence <- stringr::str_c(sequence, raw[i])
      i <- i + 1
    }
    print("sequence: ")
    print(sequence)

    # In this place labels and sequences must be united into a matrix. It should work well
    if (result == 0){
      result <- matrix(c(label, sequence), nrow = 1)
    }
    else{
      t <- matrix(c(label, sequence), nrow = 1)
      result <- rbind(result, t)
    }
  }
  print(result)
  return(result)
}


#' This function provides alignment with simple similarity matrix
#' 
#' @export
#' @param sequence a file with .fas extension. These are your DNA sequences.
#' @param match used to evaluate a cell of similarity matrix in case of 
#' matching of two nucleotides
#' @param mismatch used to evaluate a cell of similarity matrix in case of 
#' matching of two nucleotides
#' @param gap used to evaluate a cell of similarity matrix in case when you
#' should add a gap to your alignment
align <- function(sequence, match=1, mismatch=-1, gap=-1){
  # reading sequences


  # Forming a similarity matrix
  r1 <- c(match, mismatch, mismatch, mismatch)
  r2 <- c(mismatch, match, mismatch, mismatch)
  r3 <- c(mismatch, mismatch, match, mismatch)
  r4 <- c(mismatch, mismatch, mismatch, match)
  s <- rbind(r1, r2, r3, r4)
}
NULL