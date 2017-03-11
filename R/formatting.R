format_to_docs <- function(doc,vocab){

  doc <- doc[doc != 0]
  idx <- match(names(doc),vocab)

  return(rbind(idx=idx,count=doc))

}

format_gene_metadata <- function(picrust_out){

  gene_metadata <- picrust_out$genemeta
  gene_ids <- picrust_out$gene_ids

  gene_metadata_list <- vector(mode='list')
  for (i in seq_along(gene_metadata)){

    line <- gene_metadata[i]
    values <- unlist(strsplit(gene_metadata[i],'\t'))
    type <- values[1]
    values <- values[-1]

    if (grepl('\\;',line)){
      if (grepl('\\|',line)){
        gene_metadata_list[[type]] <- lapply(values, function(x) strsplit(unlist(strsplit(x,'\\|')),'\\;'))
      }else{
        gene_metadata_list[[type]] <- lapply(values,function(x) unlist(strsplit(x,'\\;')))
      }
    }else{
      gene_metadata_list[[type]] <- values
    }

    names(gene_metadata_list[[i]]) <- gene_ids
    gene_metadata_list[[i]] <- lapply(gene_metadata_list[[i]],
                                      function(x) lapply(x, function(y) trimws(tolower(y))))

  }

  return(gene_metadata_list)

}

format_gene_table <- function(functions,level,pw_targets,keep){

  if (missing(pw_targets)){
    pw_targets <- list(
      c('metabolism',
        'environmental information processing',
        'organismal systems',
        'cellular processes'),
      c('carbohydrate metabolism',
        'amino acid metabolism',
        'lipid metabolism',
        'cell motility',
        'glycan biosynthesis and metabolism',
        'xenobiotics biodegradation and metabolism',
        'metabolism of cofactors and vitamins',
        'nucleotide metabolism',
        'metabolism of terpenoids and polyketides',
        'biosynthesis of other secondary metabolites',
        'energy metabolism',
        'metabolism',
        'membrane transport',
        'cellular processes and signaling',
        'genetic information processing'),
      c('none'))
    keep <- c(TRUE,TRUE,FALSE)
  }

  gene_counts <- functions$fxn_table
  pathways <- functions$fxn_meta$KEGG_Pathways

  K <- nrow(gene_counts)

  pathways_level <- vector(mode='list',length=level)
  for (i in seq_len(level)){
    pathways_level[[i]] <- lapply(pathways,function(x) sapply(x,function(y) y[i]))
  }

  for (i in seq_len(level)){
    if (keep[i] == TRUE){
      next
    }else{
      pw_tmp <- unique(unlist(pathways_level[[i]]))
      pw_targets[[i]] <- pw_tmp[!(pw_tmp %in% pw_targets[[i]])]
    }
  }


  gene_filter <- vector(mode='list',length=level-1)
  for (i in seq_len(level-1)){
    gene_filter[[i]] <- unlist(lapply(pw_targets[[i]],function(pw) names(which(sapply(pathways_level[[i]], function(x) any(x %in% pw))))))
  }
  gene_filter <- Reduce(intersect, gene_filter)

  if (!is.null(gene_filter)){
    pathways_level_tmp <- names(pathways_level[[level]])
    pathways_level[[level]] <- pathways_level[[level]][pathways_level_tmp %in% gene_filter]
  }

  gene_targets <- lapply(pw_targets[[level]],function(pw) names(which(sapply(pathways_level[[level]], function(x) any(x %in% pw)))))
  names(gene_targets) <- pw_targets[[level]]
  gene_targets <- gene_targets[sapply(gene_targets,length)>0]

  j1 <- 1
  j2 <- 0
  gene_table <- as.data.frame(matrix(0,sum(sapply(gene_targets,length)) * K,4,
                                     dimnames=list(NULL,c('topic','pw','ko','count'))))
  for (i in seq_along(gene_targets)){
    genes <- gene_targets[[i]]
    for (k in seq_len(K)){
      j2 <- j2 + length(genes)
      gene_table[j1:j2,1] <- k
      gene_table[j1:j2,2] <- pw_targets[[level]][i]
      gene_table[j1:j2,3] <- genes
      gene_table[j1:j2,4] <- gene_counts[k,genes]
      j1 <- j2+1
    }
  }

  return(gene_table)

}

