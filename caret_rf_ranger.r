loadNamespace("caret")
## Required interface functions to use ranger::ranger Random Forest implementation with caret::rfe feature selection

## our constant, not part of the caret interface

make.rangerCaretRfeFuncs <- function(rerank=F,
                                     dependent.varID = ".ranger.dependent.varID",
                                     selectSize=caret::pickSizeBest,
                                     fit.args.first=list(),
                                     fit.args.other=fit.args.first,
                                     fit.args.last=fit.args.first) {
  list(
  

summary = caret::defaultSummary
,

fit = 
function (x, y, first, last, ...) 
{
  loadNamespace("ranger")
  x = as.data.frame(x)
  x[[dependent.varID]]= y
  ##need to write.forest for predict() to work
  ##need to always compute with importance in rerank mode of rfeControl,
  ##otherwise get "out of bound" error. This happens with buil-in rfFuncs
  ##object too.
  importance = if(rerank || first || last) "impurity" else "none"
  #print(paste("first=",first,"last=",last,"dim=",ncol(x),"importance=",importance))
  
  if(first) fit.args = fit.args.first
  else if(last) fit.args = fit.args.last
  else fit.args = fit.args.other
  
  do.call(ranger::ranger,
          c(
          list(data = x, dependent.variable.name = dependent.varID, 
                  importance = importance, write.forest = T, 
                 ...),
          fit.args
          )
  )

}
,
  
pred =
function (object, x) 
{
  x = as.data.frame(x)
  x[[dependent.varID]] = 0
  #print(paste("pred dim=",ncol(x)))
  ##ranger::predict looks at dependent.varID inside its tree object, so we need to add such column,
  ##otherwise it cuts away one of the predictors
  pred=predict(object,x,verbose=F)$predictions
}
,
  
rank =
function (object, x, y) 
{
  vimp <- ranger::importance(object)
  #print(paste("rank dim=",length(vimp)))
  vimp <- vimp[order(vimp, decreasing = TRUE)]
  data.frame(var=names(vimp),Overall=vimp)
}
,

selectSize = selectSize
,
  
selectVar = caret::pickVars

)
}