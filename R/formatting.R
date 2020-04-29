# sample vector of counts to document format
format_to_docs <- function(sample,vocab){

  sample <- sample[sample != 0]
  idx <- match(names(sample),vocab)

  doc <- matrix(0L,2,length(sample),
                dimnames=list(c('idx','count'),names(sample)))
  doc[1,] <- idx
  doc[2,] <- sample

  return(doc)

}

# functional metadata strings to lists of lists
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

# format functional output to gene table for inference
format_gene_table <- function(functions,level,pw_targets,keep=NULL){

  description_idx <- grepl('Description',names(functions$fxn_meta))
  prefix <- gsub('^(.*)\\_Description$','\\1',names(functions$fxn_meta)[description_idx])
  names(functions$fxn_meta)[description_idx] <- 'Description'
  names(functions$fxn_meta)[!description_idx] <- 'Category'

  if (missing(pw_targets)){

    if (prefix == 'COG'){
      pw_targets <- list(
        c('poorly characterized'),
        c('[r] general function prediction only',
          '[s] function unknown')
      )
      keep <- c(FALSE,FALSE)
    }
    if (prefix == 'KEGG'){
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
  }

  gene_counts <- functions$fxn_table
  category <- functions$fxn_meta[['Category']]
  description <- sapply(functions$fxn_meta[['Description']],paste0,collapse='; ')

  K <- nrow(gene_counts)

  category_level <- vector(mode='list',length=level)
  for (i in seq_len(level)){
    category_level[[i]] <- lapply(category,function(x) sapply(x,function(y) y[i]))
  }

  for (i in seq_len(level)){
    if (keep[i] == TRUE){
      next
    }else{
      pw_tmp <- unique(unlist(category_level[[i]]))
      pw_targets[[i]] <- pw_tmp[!(pw_tmp %in% pw_targets[[i]])]
    }
  }


  gene_filter <- vector(mode='list',length=level-1)
  for (i in seq_len(level-1)){
    gene_filter[[i]] <- unlist(lapply(pw_targets[[i]],function(pw) names(which(sapply(category_level[[i]], function(x) any(x %in% pw))))))
  }
  gene_filter <- Reduce(intersect,gene_filter)

  if (!is.null(gene_filter)){
    category_level_tmp <- names(category_level[[level]])
    category_level[[level]] <- category_level[[level]][category_level_tmp %in% gene_filter]
  }

  gene_targets <- lapply(pw_targets[[level]],function(pw) names(which(sapply(category_level[[level]], function(x) any(x %in% pw)))))
  names(gene_targets) <- pw_targets[[level]]
  gene_targets <- gene_targets[sapply(gene_targets,length)>0]

  j1 <- 1
  j2 <- 0
  gene_table <- as.data.frame(matrix(0,sum(sapply(gene_targets,length)) * K,5,
                                     dimnames=list(NULL,c('topic','pw','ko','description','count'))))
  for (i in seq_along(gene_targets)){
    genes <- gene_targets[[i]]
    for (k in seq_len(K)){
      j2 <- j2 + length(genes)
      gene_table[j1:j2,1] <- k
      gene_table[j1:j2,2] <- names(gene_targets)[i]
      gene_table[j1:j2,3] <- genes
      gene_table[j1:j2,4] <- description[genes]
      gene_table[j1:j2,5] <- gene_counts[k,genes]
      j1 <- j2+1
    }
  }

  return(gene_table)

}

# generate pretty taxa names
pretty_taxa_names <- function(x){

  taxonomy <- c('Species','Genus','Family','Order','Class','Phylum','Kingdom')

  taxon <- gsub('[a-z]__','',x[,taxonomy[1]])
  taxon <- ifelse(taxon != '',paste(gsub('[a-z]__','',x[,taxonomy[2]]),taxon),'')
  bad_names <- taxon == ''

  for (tax in taxonomy[-1]){
    clean <- gsub('[a-z]__','',x[bad_names,tax])
    taxon[bad_names] <- ifelse(clean != '',paste(clean,tolower(tax)),'')
    bad_names <- taxon == ''
  }

  taxon

}


# group_by summarize function for taxa without dplyr
sum_taxa_by_group <- function(otu_ids,taxa,otu_table,metadata,cov_list,group=c('Phylum','Class','Order','Family','Genus'),sample_norm=FALSE){

  if (sample_norm) z <- cov_list$coverage else {z <- sapply(cov_list$ids,length); z <- z/stats::median(z)}

  group_list <- vector(mode='list',length=length(group))
  names(group_list) <- group
  for (g in names(group_list)){

    group_df <- rbind(data.frame(taxon=taxa[otu_ids,g],abundance=colSums(otu_table[cov_list$ids[[1]],otu_ids])/z[1],cov=names(cov_list$ids)[1],group=g,
                                 stringsAsFactors=TRUE),
                      data.frame(taxon=taxa[otu_ids,g],abundance=colSums(otu_table[cov_list$ids[[2]],otu_ids])/z[2],cov=names(cov_list$ids)[2],group=g,
                                 stringsAsFactors=TRUE))

    group_list[[g]] <- with(group_df,aggregate(abundance,by=list(taxon=taxon,cov=cov,group=group),FUN=sum))

  }

  group_list <- do.call('rbind',group_list)
  colnames(group_list) <- gsub('^x$','abundance',colnames(group_list))

  group_list <- group_list[order(group_list$abundance,decreasing=TRUE),]
  group_list$taxon <- with(group_list,factor(taxon,levels=c(as.character(unique(taxon)[unique(taxon) != 'Other']),'Other'),ordered=TRUE))
  group_list$group <- with(group_list,factor(group,levels=colnames(taxa),ordered=TRUE))

  return(group_list)

}

# rename taxa below rank to other
rename_taxa_to_other <- function(x,taxa,top_n=7,group=c('Phylum','Class','Order','Family','Genus'),type=c('otu_table','docs'),as_factor=FALSE){

  type <- match.arg(type)

  taxa <- as.data.frame(taxa)

  if (type == 'otu_table') taxa <- taxa[colnames(x),]

  for (g in group){

    taxa[,g] <- gsub('^[a-z]__','',taxa[,g])
    taxa_temp <- taxa[taxa[,g] != '',g]
    taxa_temp_ids <- rownames(taxa)[taxa[,g] != '']

    if (type == 'otu_table'){

      group_df <- data.frame(count=colSums(x[,taxa_temp_ids]),taxon=taxa_temp,stringsAsFactors=TRUE)

    }else{

      taxa_temp_ids_total <- numeric(length(taxa_temp_ids))
      names(taxa_temp_ids_total) <- taxa_temp_ids
      for (d in seq_along(x)){

        taxa_ids_in_doc <- taxa_temp_ids[taxa_temp_ids %in% colnames(x[[d]])]
        taxa_temp_ids_total[taxa_ids_in_doc] <- taxa_temp_ids_total[taxa_ids_in_doc] + x[[d]][2,taxa_ids_in_doc]

      }

      group_df <- data.frame(count=taxa_temp_ids_total,taxon=taxa_temp,stringsAsFactors=TRUE)
      group_df <- group_df[group_df$count > 0,]

    }
    group_df <- with(group_df,aggregate(count,by=list(taxon=taxon),FUN=sum))
    taxa_top <- as.character(group_df[order(group_df$x,decreasing=TRUE),'taxon'][1:top_n])

    taxa[,g] <- ifelse(taxa[,g] %in% taxa_top,taxa[,g],'Other')

    if (as_factor) taxa[,g] <- factor(taxa[,g],levels=unique(taxa[,g]),ordered=TRUE)

  }

  return(taxa)

}
