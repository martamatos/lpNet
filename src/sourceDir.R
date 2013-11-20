#
# sources all the R files in the given directory
#
sourceDir <- function(path, trace = TRUE) {
  for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
    if(trace) cat(nm,":")
#    print(file.path(path, nm))          
    source(file.path(path, nm))
    
    if(trace) cat("\n")
  }
}
