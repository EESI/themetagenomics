#' @import ggplot2 shiny
#' @importFrom plotly plotlyOutput renderPlotly ggplotly
NULL

#' Generate interactive binary covariate effect figure
#'
#' TBD
#'
#' @param topics Output of \code{\link{find_topics}} that contains the STM object.
#' @param topic_effects Output of \code{\link{estimate_topic_effects}} that contains regression weights.
#' @param function_effects Output of \code{\link{estimate_function_effects}} that contains the results from either HMC or ML.
#' @param taxa Dataframe or matrix containing the taxonomy information.
#' @param beta_min (optional) Minimum probability in topics over taxa distribution to set to 0.
#' @param gene_min (optional) Mininum count for gene set table.
#'
#' @export

vis_function_effects <- function(topics,topic_effects,otu_table,taxa,covariate,metadata){


  cov1_name <- unique(metadata[,covariate])[1]
  cov2_name <- unique(metadata[,covariate])[2]

  cov1 <- rownames(otu_table)[metadata[,covariate] == cov1_name]
  cov2 <- rownames(otu_table)[metadata[,covariate] == cov2_name]


  pretty_names <- pretty_taxa_names(taxa)

  taxa_top <- as.data.frame(taxa)
  taxa_top[,2] <- gsub('^[a-z]__','',taxa_top[,2])
  taxa_top[,3] <- gsub('^[a-z]__','',taxa_top[,3])
  taxa_top[,4] <- gsub('^[a-z]__','',taxa_top[,4])

  phyla_top <- names(sort(table(taxa_top[,2]),decreasing=TRUE)[1:7])
  taxa_top[,2] <- ifelse(taxa_top[,2] %in% phyla_top,taxa_top[,2],'Other')

  class_top <- names(sort(table(taxa_top[,3]),decreasing=TRUE)[1:7])
  taxa_top[,3] <- ifelse(taxa_top[,3] %in% class_top,taxa_top[,3],'Other')

  order_top <- names(sort(table(taxa_top[,4]),decreasing=TRUE)[1:7])
  taxa_top[,4] <- ifelse(taxa_top[,4] %in% order_top,taxa_top[,4],'Other')

  fit <- topics$fit
  K <- fit$settings$dim$K
  beta <- t(exp(fit$beta$logbeta[[1]]))
  rownames(beta) <- fit$vocab
  colnames(beta) <- paste0('T',1:K)

  est_mat <- do.call('rbind',topic_effects$est)

  df <- data.frame(topic=paste0('T',1:K),
                   est=est_mat[,1],
                   lower=est_mat[,2],
                   upper=est_mat[,3],
                   sig=ifelse(1:K %in% topic_effects$sig,'1','0'))[order(topic_effects$rank),]
  df$sig <- as.character(sign(df$est) * as.numeric(as.character.factor(df$sig)))
  df$topic <- factor(df$topic,levels=df$topic,ordered=TRUE)
  p_est <- ggplot(df,aes(topic,y=est,ymin=lower,ymax=upper,color=sig)) +
    geom_hline(yintercept=0,linetype=3) +
    geom_pointrange(size=2) +
    theme_void() +
    labs(x='',y='Estimate') +
    scale_color_manual(values=c('gray','darkred','darkblue')) +
    scale_fill_brewer(type='qual',palette=6,direction=-1) +
    theme(legend.position='none')

  df0 <- data.frame(count=c(1,max(otu_table)),
                    otu=c('x','y'),
                    p=c(0,1),
                    covariate=unique(metadata[,covariate]),
                    taxon=c('x','y'))


  shinyApp(


    ui <- fluidPage(

      titlePanel('Thematic Structure of Predicted Functional Effects'),

      fixedRow(
        column(12,
               plotlyOutput('est',height='200px'),
               fixedRow(
                 column(2,
                        sliderInput('z', label='exp(x): ',
                                    min=-6,max=-1,value=-2,step=.5)
                 ),
                 column(10,
                        plotlyOutput('dis')
                 )
               )
        )
      ),


      fixedRow(
        column(12,plotOutput('highlight'))
      )

    ),


    server = function(input,output){






      output$est <- renderPlotly({
        ggplotly(p_est,source='reg_est')
      })






      output$dis <- renderPlotly({

        z <- reactive(as.integer(input$z))

        s <- event_data('plotly_click',source='reg_est')

        if (length(s)){

          k <- levels(df$topic)[s[['x']]]

          beta_subset <- beta[,k]
          beta_subset <- beta_subset[beta_subset > 10^(z())]

          otu_subset <- t(otu_table[,names(beta_subset)])
          suppressWarnings(
          df1 <- data.frame(count=matrix(otu_subset,ncol=1)+1,
                            otu=rownames(otu_subset),
                            p=beta_subset,
                            sample=rep(colnames(otu_subset),each=length(nrow(otu_subset))),
                            covariate=metadata[colnames(otu_subset),covariate])
          )
          df1$taxon <- paste0(pretty_names[df1$otu],' (',df1$otu,')')


          df2 <- df1[df1$count > 1,]

          suppressWarnings(
          p_dis <- ggplot(data=df0,aes(covariate,count,color=p)) +
            geom_blank() +
            geom_violin(data=df2,aes(covariate,count,color=p),fill=NA) +
            geom_jitter(data=df1,aes(covariate,count,color=p,labels=taxon),alpha=.5,size=2) +
            scale_y_log10() +
            viridis::scale_color_viridis('Probability') +
            theme_classic() +
            labs(x='',y='Count')
          )

          p_dis <- ggplotly(p_dis,source='dis_cov',tooltip=c('taxon','otu'))
          layout(p_dis,dragmode='select')

        }else{

          p_dis <- ggplot(data=df0,aes(covariate,count,color=p)) +
            geom_blank() +
            scale_y_log10() +
            viridis::scale_color_viridis('Probability') +
            theme_classic() +
            labs(x='',y='Count')

          p_dis

        }

      })




      output$highlight <- renderPlot({

        h <- event_data('plotly_selected',source='dis_cov')

        if (length(h)){

          h_otu <- gsub('^taxon:.*\\(([0-9]+)\\)$','\\1',p_dis$x$data[[3]]$text)[h$pointNumber]

          h_otu <- sample(colnames(otu_table),25)

          df_tax1 <- data.frame(taxon=taxa_top[h_otu,2],count=colSums(otu_table[cov1,h_otu]),cov=cov1_name)
          df_tax2 <- data.frame(taxon=taxa_top[h_otu,2],count=colSums(otu_table[cov2,h_otu]),cov=cov2_name)
          df_tax_p <- rbind(aggregate(df_tax1$count,by=list(taxon=df_tax1$taxon,cov=df_tax1$cov),FUN=sum),
                          aggregate(df_tax2$count,by=list(taxon=df_tax2$taxon,cov=df_tax2$cov),FUN=sum))

          df_tax1 <- data.frame(taxon=taxa_top[h_otu,3],count=colSums(otu_table[cov1,h_otu]),cov=cov1_name)
          df_tax2 <- data.frame(taxon=taxa_top[h_otu,3],count=colSums(otu_table[cov2,h_otu]),cov=cov2_name)
          df_tax_c <- rbind(aggregate(df_tax1$count,by=list(taxon=df_tax1$taxon,cov=df_tax1$cov),FUN=sum),
                          aggregate(df_tax2$count,by=list(taxon=df_tax2$taxon,cov=df_tax2$cov),FUN=sum))

          df_tax1 <- data.frame(taxon=taxa_top[h_otu,4],count=colSums(otu_table[cov1,h_otu]),cov=cov1_name)
          df_tax2 <- data.frame(taxon=taxa_top[h_otu,4],count=colSums(otu_table[cov2,h_otu]),cov=cov2_name)
          df_tax_o <- rbind(aggregate(df_tax1$count,by=list(taxon=df_tax1$taxon,cov=df_tax1$cov),FUN=sum),
                          aggregate(df_tax2$count,by=list(taxon=df_tax2$taxon,cov=df_tax2$cov),FUN=sum))

          df_tax <- data.frame(rbind(df_tax_p,df_tax_c,df_tax_o),
                               group=rep(c('phylum','class','order'),c(nrow(df_tax_p),nrow(df_tax_c),nrow(df_tax_o))))
          colnames(df_tax)[3] <- 'count'

          df_tax <- df_tax[order(df_tax$count,decreasing=TRUE),]
          df_tax$taxon <- factor(df_tax$taxon,levels=c(as.character(unique(df_tax$taxon)[unique(df_tax$taxon) != 'Other']),'Other'),ordered=TRUE)
          df_tax$group <- factor(df_tax$group,levels=c('phylum','class','order'),ordered=TRUE)

          p_tax <- ggplot(df_tax,aes(taxon,count,fill=cov)) +
            geom_bar(color='black',stat='identity',position='dodge') +
            facet_grid(.~group,scales='free_x') +
            theme_classic() +
            theme(axis.text.x=element_text(angle=45,hjust=1)) +
            labs(x='',y='Count')

          p_tax

        }else{

          ggplot() + geom_blank() + theme_classic() + labs(x='',y='')

        }

      })





    }



  )

}
