#' Evaluate 'global' cross-validation
#' @description Calculate validation metric using all held back predictions at once
#' @param model an object of class \code{\link{train}} or \code{\link{ResampleResult}}
#' @return regression (\code{\link{postResample}}) or classification  (\code{\link{confusionMatrix}}) statistics
#' @details Relevant when folds are not representative for the entire area of interest.
#' In this case, metrics like R2 are not meaningful since it doesn't reflect the general ability of
#' the model to explain the entire gradient of the response.
#' Comparable to LOOCV, predictions from all held back folds are used here together to calculate validation statistics.
#' @author Hanna Meyer, Jan Linnenbrink
#' @seealso \code{\link{CreateSpacetimeFolds}}
#' @examples
#' dat <- get(load(system.file("extdata","Cookfarm.RData",package="CAST")))
#' dat <- dat[sample(1:nrow(dat),500),]
#' indices <- CreateSpacetimeFolds(dat,"SOURCEID","Date")
#' ctrl <- caret::trainControl(method="cv",index = indices$index,savePredictions="final")
#' model <- caret::train(dat[,c("DEM","TWI","BLD")],dat$VW, method="rf", trControl=ctrl, ntree=10)
#' global_validation(model)
#' @export global_validation
#' @aliases global_validation

global_validation <- function(model){

  if (inherits(model, "data.frame")) {
    obs <- model$observed
    pred <- model$predicted
  } else if (inherits(model, "ResampleResult")) {
    predictions <- model$predictions()
    predictions <- lapply(predictions, mlr3::as.data.table)
    predictions <- do.call("rbind", predictions)
    obs <- predictions$truth
    pred <- predictions$response

  } else {
    predictions <- model$pred
    if(is.null(predictions)){stop("Global performance could not be estimated because predictions were not saved.
                                Train model with savePredictions='final'")}

    ### only use predictions of best tune:
    for (i in 1:length(model$bestTune)){
      predictions <- predictions[predictions[,names(model$bestTune)[i]]==model$bestTune[,i],]
    }

    obs <- predictions$obs
    pred <- predictions$pred
  }

  if(model$modelType=="Regression" || grepl("regr", model$learner$man)){
    out <- caret::postResample(pred = pred, obs = obs)
  }else{
    out <- caret::confusionMatrix(pred, obs)$overall
  }
  return(out)
}

