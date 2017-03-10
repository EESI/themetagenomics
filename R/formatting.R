format_to_docs <- function(doc,vocab){

  doc <- doc[doc != 0]
  idx <- match(names(doc),vocab)

  return(rbind(idx=idx,count=doc))

}

format_gene_metadata <- function(gene_metadata){

  gene_metadata_list <- vector(mode='list')
  for (i in seq_along(gene_metadata)){

    line <- gene_metadata[i]
    values <- unlist(strsplit(gene_metadata[i],'\t'))
    type <- values[1]
    values <- values[-1]

    if (grepl('\\;',line)){
      if (grepl('\\|',line)){
        gene_metadata_list[[type]] <- lapply(values, function(x) lapply(strsplit(trimws(unlist(strsplit(x,'\\|'))),'\\;'),trimws))
      }else{
        gene_metadata_list[[type]] <- lapply(values,function(x) trimws(unlist(strsplit(x,'\\;'))))
      }
    }else{
      gene_metadata_list[[list]] <- values
    }
  }

  return(gene_metadata_list)

}
