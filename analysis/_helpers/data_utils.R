# analysis/_helpers/data_utils.R
# Shared data loading and string manipulation utilities
#
# Usage: source("analysis/_helpers/data_utils.R") from the repository root.

#' Split a string and extract field(s)
#'
#' Vectorized string splitting helper. Extracts one or more fields from each
#' element of a character vector after splitting on a separator.
#'
#' @param myStr Character vector to split
#' @param mySep Separator string (passed to strsplit)
#' @param myField Integer vector of field indices to extract
#' @return Character vector of extracted fields (collapsed if length(myField) > 1)
#'
#' @examples
#' chop("ACGT-B73_v5", "-", 1)        # "ACGT"
#' chop("ACGT-B73_v5", "-", 2)        # "B73_v5"
#' chop("A-B-C-D", "-", c(2, 3))      # "B-C"
#'
#' Shared by the Fig S1 script (analysis/supplementary/figS1.R) and the Socrates QC engine
#' (workflows/05_qc_and_embedding/common/1_1_QC_scifiATAC_data.R).
chop <- function(myStr, mySep, myField) {
  choppedString <- sapply(strsplit(myStr, mySep), "[", myField)
  if (length(myField) > 1) {
    choppedString <- apply(choppedString, 2, function(x) {
      paste0(x[!is.na(x)], collapse = mySep)
    })
  }
  return(choppedString)
}
