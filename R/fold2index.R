#' Fold to Index
#' @description Creates index and indexOut for caret::trainControl from a fold vector
#'
#' @param folds vector with fold labels
#'
#' @return list, with index and indexOut
#'
#' @import purrr
#' @import dplyr
#'
#' @export
#'
#' @author Marvin Ludwig
#'
#'
#'

fold=kndm$clusters
fold2index = function(fold){
    
    fold = data.frame(fold = fold)

    indOut = fold %>% dplyr::group_by(fold) %>%
        attr('groups') %>% dplyr::pull(.rows)

    ind = purrr::map(seq(length(indOut)), function(x){
        s = seq(nrow(fold))
        s = s[!s %in% indOut[[x]]]
        return(s)
    })
    return(
        list(
            index = ind,
            indexOut = indOut
        )
    )

}

x=data.frame(clust=kndm$clusters)
spacevar = "clust"
k = 10
timevar=NA
class=NA
seed=sample(1:1000, 1)

CreateSpacetimeFolds <- function(x,spacevar=NA,timevar=NA,
                                 k=10,class=NA,seed=sample(1:1000, 1)){
  x <- data.frame(x)
  
  
  if(!is.na(spacevar)){
    if(k>length(unique(x[,spacevar]))){
      k <- length(unique(x[,spacevar]))
      print(paste0("warning: k is higher than number of unique locations. k is set to ",k))
    }
  }
 
  #split space into k folds
  if(!is.na(spacevar)){
    set.seed(seed)
    spacefolds <- lapply(caret::createFolds(1:length(unique(x[,spacevar])),k),function(y){
      unique(x[,spacevar])[y]})
  }
  
  # combine space and time folds
  cvindices_train <- list()
  cvindices_test <- list()
  for (i in 1:k){
    if(is.na(timevar)&!is.na(spacevar)){
      cvindices_test[[i]]<- which(x[,spacevar]%in%spacefolds[[i]])
      cvindices_train[[i]]<- which(!x[,spacevar]%in%spacefolds[[i]])
    }
  }
  
  return(list("index"=cvindices_train,"indexOut"=cvindices_test))
}


spf <- list("index"=cvindices_train,"indexOut"=cvindices_test)
ind <- fold2index(kndm$clusters)

# unterschiede !!!
tail(spf$indexOut[[10]])
tail(ind$indexOut[[10]])
spf$indexOut[[10]]==ind$indexOut[[10]]
spf$index[[10]]==ind$index[[10]]
