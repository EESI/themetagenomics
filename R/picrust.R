picrust <- function(x,rows_are_taxa,reference_path,drop=TRUE){

  if (rows_are_taxa == TRUE){

    otu_table <- t(otu_table)

  }

  out <- picrust_otu(reference_path,colnames(x))
  fxn_mapping <- out$genome_table_out
  rownames(fxn_mapping) <- out$matches
  colnames(fxn_mapping) <- out$gene_ids

  overlap <- intersect(rownames(fxn_mapping),colnames(x))
  x <- x[,overlap]
  fxn_mapping <- fxn_mapping[overlap,]

  fxn_table <- round(x %*% fxn_mapping)
  fxn_meta <- format_gene_metadata(out)
  pi_meta <- out$pimeta_table_out
  rownames(pi_meta) <- fit$vocab
  colnames(pi_meta) <- out$pimeta_ids

  if (drop){
    fxn_table <- fxn_table[,colSums(fxn_table)>0]
    fxn_meta <- lapply(fxn_meta,function(x) x[colnames(fxn_table)])
  }

  predictions <- list(fxn_table=fxn_table,fxn_meta=fxn_meta,pi_meta=pi_meta)

  return(predictions)

}
