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

vis_covariate_effects_continuous <- function(topics,topic_effects,otu_table,taxa,covariate,metadata,sample_norm=FALSE,taxa_n=10){

  pretty_names <- pretty_taxa_names(taxa)

  fit <- topics$fit
  K <- fit$settings$dim$K

  beta <- t(exp(fit$beta$logbeta[[1]]))
  rownames(beta) <- fit$vocab
  colnames(beta) <- paste0('T',1:K)

  theta <- fit$theta
  rownames(theta) <- rownames(metadata)
  colnames(theta) <- paste0('T',1:K)

  beta_rank <- apply(beta,2,function(x) rownames(beta)[order(x,decreasing=TRUE)])

  ra_table <- otu_table/rowSums(otu_table)

  est_mat <- do.call('rbind',topic_effects$est)

  if (sample_norm) df0_range <- c(1e-10,max(ra_table)) else df0_range <- c(1,max(otu_table))
  df0 <- data.frame(abundance=df0_range,
                    p=c(0,max(theta)),
                    covariate=range(metadata[,covariate]))


  shinyApp(


    ui <- fluidPage(

      titlePanel('Thematic Structure of Predicted Functional Effects'),

      fixedRow(
        column(12,plotlyOutput('est',height='200px'),
               fixedRow(
                 column(12,plotlyOutput('dis')
                 )
               )
        )
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

          beta_subset <- beta_rank[1:taxa_n,k]

          if (sample_norm) x <- ra_table + 1e-10 else x <- otu_table + 1

          otu_subset <- t(x[,beta_subset])
          suppressWarnings(
            df1 <- data.frame(abundance=matrix(otu_subset,ncol=1),
                              otu=rownames(otu_subset),
                              sample=rep(colnames(otu_subset),each=nrow(otu_subset)),
                              covariate=metadata[colnames(otu_subset),covariate])
          )
          df1$p <- theta[df1$sample,k]
          df1$taxon <- paste0(pretty_names[df1$otu],' (',df1$otu,')')

          if (sample_norm) df2 <- df1[df1$abundance > 1e-10,] else df2 <- df1[df1$abundance > 1,]

          suppressWarnings({
            p_reg <- ggplot(data=df0,aes(covariate,abundance,color=p)) +
              geom_blank() +
              geom_point(data=df1,aes(covariate,abundance,color=p),alpha=.3,size=2) +
              stat_smooth(data=df2,aes(covariate,abundance),color='black',size=1.1,method='loess',se=FALSE) +
              stat_smooth(data=df2,aes(covariate,abundance),color='darkred',size=1.1,method='lm',linetype=2,se=FALSE) +
              facet_wrap(~taxon,ncol=5) +
              scale_y_log10() +
              viridis::scale_color_viridis('Probability') +
              theme_classic() +
              theme(aspect.ratio=.5) +
              labs(x=covariate,y='Abundance')

          ggplotly(p_reg)
          })

        }else{

          p_reg <- ggplot(data=df0,aes(covariate,abundance,color=p)) +
            geom_blank() +
            scale_y_log10() +
            viridis::scale_color_viridis('Probability') +
            theme_classic() +
            theme(aspect.ratio=.5) +
            labs(x='',y='')

          p_reg

        }

      })






    }



  )

}
