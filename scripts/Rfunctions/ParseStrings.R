# Jake Yeung
# ParseStrings.R
# Functions to parse strings and convert them to different data types  
# 2018-12-22

StrToVector <- function(s, delim = "-"){
  # WT-KO-Kidney -> c("WT", "KO", "Kidney")
  return(strsplit(s, delim)[[1]])
}

StrToNumeric <- function(s){
  # string to numeric, if fails then keep string as string
  num <- as.numeric(s)
  if (is.na(num)){
    if (s == "NA"){
      return(NA)
    } else {
      return(s)
    }
  } else {
    return(num)
  }
}   


StrToBool <- function(s, default){
  if (s == "TRUE"){
    return(TRUE)
  } else if (s == "FALSE"){
    return(FALSE)
  } else if (s == "NA"){
    return(NA)  
  } else {
    warning(paste("String didnt match TRUE or FALSE, setting as default", default))
    return(default)
  }
}

StrToNA <- function(s){
  if (s == "NA"){
    return(NA)
  } else {
    return(s)
  }
}
