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
#' @param gene_min (optional) Mininum abundance for gene set table.

vis_covariate_effects_binary <- function(topics,topic_effects,otu_table,taxa,covariate,metadata,taxa_n=7){

  otu_table <- otu_table + 1
  ra_table <- otu_table/rowSums(otu_table)

  cov_names <- c(unique(metadata[,covariate])[1],unique(metadata[,covariate])[2])
  cov_ids <- list(rownames(otu_table)[metadata[,covariate] == cov_names[1]],rownames(otu_table)[metadata[,covariate] == cov_names[2]])
  names(cov_ids) <- cov_names
  cov_coverage <- c(sum(otu_table[cov_ids[[1]],]),sum(otu_table[cov_ids[[2]],]))
  names(cov_coverage) <- cov_names
  cov_list <- list(ids=cov_ids,coverage=cov_coverage)

  pretty_names <- pretty_taxa_names(taxa)
  taxa_other <- rename_taxa_to_other(otu_table,taxa,top_n=taxa_n)

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
    scale_color_manual(values=c('gray','indianred3','dodgerblue3')) +
    scale_fill_brewer(type='qual',palette=6,direction=-1) +
    theme(legend.position='none')

  df0 <- data.frame(otu=c('x','y'),
                    p=c(0,1),
                    covariate=unique(metadata[,covariate]),
                    taxon=c('x','y'))

  s <- NULL

  shinyApp(


    ui <- fluidPage(

      titlePanel('Topic-Covariate Effects'),

      fixedRow(
        column(12,plotlyOutput('est',height='200px'),
               fixedRow(
                 column(2,sliderInput('z', label='Min. beta prob.: ',min=-6,max=-1,value=-2,step=.25,
                                      pre='10^(',post=')')),

                        fixedRow(
                          column(8,htmlOutput('text'))
                        ),
                        fixedRow(
                          column(2,checkboxInput('norm',label='Normalize',value=TRUE))
                        ),
                 column(12,plotlyOutput('dis')
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


      DF1 <- reactive({

        suppressWarnings({
        s <- event_data('plotly_click',source='reg_est')

        validate(need(!is.null(s),'Please choose a topic by clicking a point.'))

        k <- levels(df$topic)[s[['x']]]

        beta_subset <- beta[,k]
        beta_subset <- beta_subset[beta_subset > 10^(as.integer(input$z))]

        if (input$norm) x <- ra_table else x <- otu_table
        df0$abundance <- c(0,max(x))

        otu_subset <- t(x[,names(beta_subset)])

        df1 <- data.frame(abundance=matrix(otu_subset,ncol=1),
                          otu=rownames(otu_subset),
                          p=beta_subset,
                          sample=rep(colnames(otu_subset),each=nrow(otu_subset)),
                          covariate=metadata[colnames(otu_subset),covariate])

        df1$covariate <- factor(df1$covariate,levels=cov_names,ordered=TRUE)
        df1$taxon <- paste0(pretty_names[df1$otu],' (',df1$otu,')')

        df2 <- df1[df1$abundance > min(df1$abundance),]

        p_dis <- ggplot(data=df0,aes(covariate,abundance,color=p)) +
          geom_blank() +
          geom_violin(data=df2,
                      aes(covariate,abundance,color=p),fill=NA) +
          geom_jitter(data=df1,
                      aes(covariate,abundance,color=p,labels=taxon),alpha=.5,size=2) +
          scale_y_log10(labels=scales::comma) +
          viridis::scale_color_viridis('Beta Prob.') +
          theme_classic() +
          labs(x='',y='Abundance')

        class(p_dis)

        p_dis <- ggplotly(p_dis,source='dis_cov',tooltip=c('taxon','abundance','p'))
        p_dis <- layout(p_dis,dragmode='select')

        list(p_dis=p_dis,df1=df1,df0=df0)
        })

      })


      observeEvent(event_data('plotly_click',source='reg_est'),
                   {
                     output$text <- renderUI({
                       HTML('The slider on the <b>left</b> adjusts the value in which taxa are filtered based on their probability in the topics over
                             taxa distribution beta. The violin plots <b>below</b> are interactive. Select a subset of points to visualize their
                             distribution in the raw data as a function on taxonomy.')
                       })
                     }
                   )


      output$dis <- renderPlotly({

        if (length(DF1()$p_dis)){

          DF1()$p_dis

        }else{

          ggplot(data=DF1()$df0,aes(covariate,abundance,color=p)) +
            geom_blank() +
            scale_y_log10() +
            viridis::scale_color_viridis('Probability') +
            theme_classic() +
            labs(x='',y='Abundance')

        }

      })

      output$highlight <- renderPlot({

        h <- event_data('plotly_selected',source='dis_cov')

        if (length(h)){

          h_otu <- gsub('^taxon:.*\\(([0-9]+)\\)$','\\1',DF1()$df1$otu)[h$pointNumber + 1] # indexed from 0

          df_tax <- try({

            df_tax <- sum_taxa_by_group(h_otu,taxa_other,otu_table,metadata,cov_list,sample_norm=input$norm)
            df_tax$cov <- factor(df_tax$cov,levels=cov_names,ordered=TRUE)

            df_tax

          },silent=TRUE)

          validate(need(!(class(df_tax) == 'try-error'),'Please select more points.'))

          p_tax <- ggplot(df_tax,aes(taxon,abundance,fill=cov)) +
            geom_bar(color='black',stat='identity',position='dodge') +
            facet_grid(.~group,scales='free_x') +
            scale_fill_manual(values=c('indianred3','dodgerblue3')) +
            theme_classic() +
            theme(axis.text.x=element_text(angle=45,hjust=1,size=16),
                  strip.text=element_text(face='bold',size=16)) +
            labs(x='',y='abundance',fill='')

          p_tax


        }else{

          ggplot() + geom_blank() + theme_classic() + labs(x='',y='')

        }

      })


    }


  )

}
