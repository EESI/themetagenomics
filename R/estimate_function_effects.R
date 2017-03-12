estimate_function_effects <- function(functions,topics_subset,level=2,method='ML',verbose=FALSE){

  fxn_table <- functions$fxn_table
  fxn_meta <- functions$fxn_meta

  fxn_table <- fxn_table[topics_subset,]
  fxn_table <- fxn_table[,colSums(fxn_table)>0]

  fxn_meta <- lapply(fxn_meta,function(x) x[colnames(fxn_table)])

  functions$fxn_table <- fxn_table
  functions$fxn_meta <- fxn_meta

  gene_table <- format_gene_table(functions,level=level)

  if (method == 'ML'){
    mm <- lme4::glmer.nb(count ~ (1|pw) + (1|topic) + (1|pw:topic),
                         data=gene_table,
                         verbose=verbose)
  }
  if (method == 'HMC'){

    mm <- fit_stan_model(gene_table,iters=1000,chains=4)

  }

}
