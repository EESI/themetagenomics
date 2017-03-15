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

vis_covariate_effects_binary <- function(topics,topic_effects,otu_table,taxa,covariate,metadata){

  cov_names <- c(unique(metadata[,covariate])[1],unique(metadata[,covariate])[2])
  cov_ids <- list(rownames(otu_table)[metadata[,covariate] == cov_names[1]],rownames(otu_table)[metadata[,covariate] == cov_names[2]])
  names(cov_ids) <- cov_names
  cov_coverage <- c(sum(otu_table[cov_ids[[1]],]),sum(otu_table[cov_ids[[2]],]))
  names(cov_coverage) <- cov_names
  cov_list <- list(ids=cov_ids,coverage=cov_coverage)

  pretty_names <- pretty_taxa_names(taxa)
  taxa_other <- rename_taxa_to_other(otu_table,taxa,top_n=7)

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
                            sample=rep(colnames(otu_subset),each=nrow(otu_subset)),
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

          df_tax <- sum_taxa_by_group(h_otu,taxa_other,otu_table,metadata,cov_list)

          p_tax <- ggplot(df_tax,aes(taxon,abundance,fill=cov)) +
            geom_bar(color='black',stat='identity',position='dodge') +
            facet_grid(.~group,scales='free_x') +
            scale_fill_manual(values=c('darkred','darkblue')) +
            theme_classic() +
            theme(axis.text.x=element_text(angle=45,hjust=1,size=16),
                  strip.text=element_text(face='bold',size=16)) +
            labs(x='',y='Count')

          p_tax

        }else{

          ggplot() + geom_blank() + theme_classic() + labs(x='',y='')

        }

      })





    }



  )

}
