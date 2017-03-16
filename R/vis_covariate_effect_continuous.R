#' @import ggplot2 shiny
#' @importFrom plotly plotlyOutput renderPlotly ggplotly
NULL

#' Generate interactive continuous covariate effect figure
#'
#' TBD
#'
#' @param topics Output of \code{\link{find_topics}} that contains the STM object.
#' @param topic_effects Output of \code{\link{estimate_topic_effects}} that contains regression weights.
#' @param function_effects Output of \code{\link{estimate_function_effects}} that contains the results from either HMC or ML.
#' @param taxa Dataframe or matrix containing the taxonomy information.
#' @param beta_min (optional) Minimum probability in topics over taxa distribution to set to 0.
#' @param gene_min (optional) Mininum count for gene set table.

vis_covariate_effects_continuous <- function(topics,topic_effects,otu_table,taxa,covariate,metadata,taxa_n=12){

  otu_table <- otu_table + 1

  pretty_names <- pretty_taxa_names(taxa)

  fit <- topics$fit
  K <- fit$settings$dim$K

  beta <- t(exp(fit$beta$logbeta[[1]]))
  rownames(beta) <- fit$vocab
  colnames(beta) <- paste0('T',1:K)

  theta <- fit$theta
  rownames(theta) <- rownames(metadata)
  colnames(theta) <- paste0('T',1:K)

  df0 <- data.frame(p=c(0,max(theta)),
                    covariate=range(metadata[,covariate]))

  beta_rank <- apply(beta,2,function(x) rownames(beta)[order(x,decreasing=TRUE)])

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
    scale_color_manual(values=c('gray','indianred3','dodgerblue3')) +
    scale_fill_brewer(type='qual',palette=6,direction=-1) +
    theme(legend.position='none')



  shinyApp(

    ui <- fluidPage(

      titlePanel('Topic-Covariate Effects'),

    fixedRow(
      column(12,plotlyOutput('est',height='200px'),
             fixedRow(
               column(2,checkboxGroupInput('z',
                                           label='',
                                           c('Ignore zeros'=1,'Normalize'=2),
                                           selected=c(1,2))),
               fixedRow(
                 column(8,htmlOutput('text'))
               ),
               column(12,plotlyOutput('ss')
               )
             )
      )
    ),

    fixedRow(
      column(12,plotlyOutput('full'))
    )

    ),



    server = function(input,output){

      output$est <- renderPlotly({

        ggplotly(p_est,source='reg_est')

      })

      TAB <- reactive({

        if (any(input$z %in% '2')) x <- otu_table/rowSums(otu_table) else x <- otu_table

        df_temp <- df0
        df_temp$abundance <- c(min(x),max(x))

        list(table=x,df=df_temp)

      })

      P <- reactive({

        suppressWarnings({

        s <- event_data('plotly_click',source='reg_est')

        validate(need(!is.null(s),'Please choose a topic by clicking a point.'))

        k <- levels(df$topic)[s[['x']]]

        beta_subset <- beta_rank[1:taxa_n,k]

        otu_subset <- t(TAB()$table[,beta_subset])

        df_temp <- data.frame(abundance=matrix(otu_subset,ncol=1),
                          otu=rownames(otu_subset),
                          sample=rep(colnames(otu_subset),each=nrow(otu_subset)),
                          covariate=metadata[colnames(otu_subset),covariate])

        df_temp$p <- theta[df_temp$sample,k]
        df_temp$taxon <- paste0(pretty_names[df_temp$otu],' (',df_temp$otu,')')

        p_full <- ggplot(data=TAB()$df,aes(covariate,abundance,color=p)) +
          geom_blank() +
          geom_point(data=df_temp,aes(covariate,abundance,color=p),alpha=.5,size=1.2) +
          scale_y_log10(labels=scales::comma) +
          viridis::scale_color_viridis('Beta Prob.') +
          theme_classic() +
          theme(aspect.ratio=.5) +
          labs(x=covariate,y='Abundance')

        p_ss <- p_full +
          facet_wrap(~taxon,ncol=4)

        list(df=df_temp,min=min(df_temp$abundance),p_ss=p_ss,p_full=p_full)

        })

      })


      observeEvent(event_data('plotly_click',source='reg_est'),
                   {
                     output$text <- renderUI({
                       HTML(sprintf("The check box on the <b>left</b> sets whether to ignore zeros for the best fit lines and whether raw counts or
                                     relative abundances will be shown. The scatter plots <b>below</b> show the top %s taxa in terms of their
                                     probability in the topics over taxa distributions versus %s. Each point represents the abundance of that taxa
                                     in that sample. Each point is colored based on the probability of that sample occuring in the chosen topic. The
                                     large scatter plot shows all %s taxa combined.",
                                    taxa_n,covariate,taxa_n))
                     })
                   }
      )

      output$ss <- renderPlotly({

        suppressWarnings({

        if (any(input$z %in% '1')){

          df2 <- P()$df[P()$df$abundance > round(P()$min,5),]

          p_ss <- P()$p_ss +
            stat_smooth(data=df2,aes(covariate,abundance),color='black',size=1.1,method='loess',se=FALSE) +
            stat_smooth(data=df2,aes(covariate,abundance),color='darkred',size=1.1,method='lm',linetype=2,se=FALSE)

        }else{

          p_ss <- P()$p_ss +
            stat_smooth(data=P()$df,aes(covariate,abundance),color='black',size=1.1,method='loess',se=FALSE) +
            stat_smooth(data=P()$df,aes(covariate,abundance),color='darkred',size=1.1,method='lm',linetype=2,se=FALSE)

        }

        ggplotly(p_ss)

        })

      })

      output$full <- renderPlotly({

        suppressWarnings({

          if (any(input$z %in% '1')){

            df2 <- P()$df[P()$df$abundance > round(P()$min,5),]

            p_full <- P()$p_full +
              stat_smooth(data=df2,aes(covariate,abundance),color='black',size=1.1,method='loess',se=FALSE) +
              stat_smooth(data=df2,aes(covariate,abundance),color='darkred',size=1.1,method='lm',linetype=2,se=FALSE)

          }else{

            p_full <- P()$p_full +
              stat_smooth(data=P()$df,aes(covariate,abundance),color='black',size=1.1,method='loess',se=FALSE) +
              stat_smooth(data=P()$df,aes(covariate,abundance),color='darkred',size=1.1,method='lm',linetype=2,se=FALSE)
          }

          ggplotly(p_full)

        })

      })


    }


  )

}
