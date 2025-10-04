# See https://stackoverflow.com/questions/5577221/can-i-load-a-saved-r-object-into-a-new-object-namehttps://stackoverflow.com/questions/5577221/can-i-load-a-saved-r-object-into-a-new-object-name
loadRData <- function(fileName) {
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}


