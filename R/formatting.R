reshape_to_docs <- function(doc,vocab){

  doc <- doc[doc != 0]
  idx <- match(names(doc),vocab)

  return(rbind(idx=idx,count=doc))

}
