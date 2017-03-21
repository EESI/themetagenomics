#' @import ggplot2 shiny plotly
NULL

#' Generate interactive functional effect heatmap
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

vis_function_effects <- function(topics,topic_effects,function_effects,taxa,beta_min=1e-5,ui_level=.8,gene_min=10,pw_min=20){

  if (!(ui_level %in% c(.99,.95,.9,.5))){
    warning(sprintf('ui_level must be .5, .8, .9, .95, or .99 -- defaulting to .8.'))
    ui_level <- .8
  }

  ui_interval <- paste0(round(100*c((1-ui_level)/2,1-(1-ui_level)/2),1),'%')

  topic_effects <- topic_effects$topic_effects

  covariates <- lapply(names(topic_effects),identity)
  names(covariates) <- tolower(names(topic_effects))
  est_range <- range(c(sapply(topic_effects, function(x) unlist(x$est))))

  pretty_names <- pretty_taxa_names(taxa)

  fit <- topics$fit
  K <- fit$settings$dim$K
  N_pw <- length(unique(function_effects$model$summary$b_pwxtopic$pw))
  pws <- unique(function_effects$model$summary$b_pwxtopic$pw)

  beta <- t(exp(fit$beta$logbeta[[1]]))
  rownames(beta) <- fit$vocab
  colnames(beta) <- paste0('T',1:K)
  beta[beta < beta_min] <- beta_min
  logbeta_global <- log(beta)

  uncertainty <- any(grepl('%',colnames(function_effects$model$summary$mu)))

  int_mat <- matrix(0.0,N_pw,K,dimnames=list(pws,paste0('T',1:K)))

  for (i in seq_len(nrow(function_effects$model$summary$b_pwxtopic))){
    pw <- as.character(function_effects$model$summary$b_pwxtopic$pw[i])
    k <- paste0('T',as.character(function_effects$model$summary$b_pwxtopic$topic[i]))
    int_mat[pw,k] <- function_effects$model$summary$b_pwxtopic$mean[i]
  }

  dd_row <- as.dendrogram(hclust(dist(int_mat,method='euclidean'),method='ward.D2'))
  pw_ord <- data.frame(pw=rownames(int_mat)[order.dendrogram(dd_row)],
                       pw_rank=1:nrow(int_mat))

  df1 <- data.frame(weight=matrix(int_mat,ncol=1),
                    pw=pws,
                    topic=rep(colnames(int_mat),each=N_pw))
  df1$pw <- factor(df1$pw,levels=unique(df1$pw),ordered=TRUE)


  if (uncertainty){

    sig_mat <- matrix(0,N_pw,K,dimnames=list(pws,paste0('T',1:K)))

    sig <- rowSums(sign(function_effects$model$summary$b_pwxtopic[,ui_interval]))

    for (i in seq_len(nrow(function_effects$model$summary$b_pwxtopic))){
      pw <- as.character(function_effects$model$summary$b_pwxtopic$pw[i])
      k <- paste0('T',as.character(function_effects$model$summary$b_pwxtopic$topic[i]))
      sig_mat[pw,k] <- sig[i]
    }

    df2 <- data.frame(sig=matrix(sig_mat,ncol=1),
                      weight=matrix(int_mat,ncol=1),
                      pw=pws,
                      topic=rep(colnames(int_mat),each=N_pw))
    df2$pw <- factor(df2$pw,levels=unique(df2$pw),ordered=TRUE)

    if (length(unique(df2$pw)) > pw_min){

      fxn_subset <- df2$pw %in% unique(df2[df2$sig != 0,'pw'])

      df1 <- df1[fxn_subset,]
      df2 <- df2[fxn_subset,]

    }

  }else{

    if (length(unique(df2$pw)) > pw_min){

      upper <- quantile(df1$weight,ui_level + (1-ui_level)/2)
      lower <- quantile(df1$weight,(1-ui_level)/2)

      df1 <- df1[df1$weight > upper | df1$weight < lower,]

    }

  }







  shinyApp(

    ui <- fluidPage(

      titlePanel('Topic-Function Effects'),

      fixedRow(
        column(2,selectInput('choose', label='Covariate',
                             choices=covariates,selected=covariates[[1]])),
        column(10,plotlyOutput('est',height='200px'))
        ),

      fixedRow(
        column(1,''),
        column(10,htmlOutput('text')),
        column(1,'')
        ),
      fixedRow(
        column(2,radioButtons('dist_tax',label=strong('Distance'),
                              choices=list('Bray Curtis'='bray','Jaccard'='jaccard','Euclidean'='euclidean',
                                           'Hellinger'='hellinger','Chi Squared'='chi2','Jensen Shannon'='jsd'),
                              selected='bray')),
        column(10,plotlyOutput('tax'))
      ),

      plotlyOutput('fxn'),

      div(tableOutput('tbl'),style='font-size:80%;line-height:50%')

    ),





    server = function(input,output){

      EST <- reactive({
        suppressWarnings({

          covariate <- input$choose

          topic_order <- names(sort(topic_effects[[covariate]]$rank))

          est_mat <- topic_effects[[covariate]][['est']]

          df <- data.frame(topic=paste0('T',1:K),
                           est=est_mat[,1],
                           lower=est_mat[,2],
                           upper=est_mat[,3],
                           sig=ifelse(1:K %in% topic_effects[[covariate]][['sig']],'1','0'))[order(topic_effects[[covariate]][['rank']]),]
          df$sig <- as.character(sign(df$est) * as.numeric(as.character.factor(df$sig)))
          df$topic <- factor(df$topic,levels=df$topic,ordered=TRUE)

          bars <- c(which(max(df[df$sig == -1,'topic']) == levels(df$topic)) + .5,
                    which(min(df[df$sig == 1,'topic']) == levels(df$topic)) - .5)

          p_est <- ggplot(df,aes(topic,y=est,ymin=lower,ymax=upper,color=sig)) +
            geom_hline(yintercept=0,linetype=3) +
            geom_pointrange(size=2) +
            theme_minimal() +
            ylim(est_range[1],est_range[2]) +
            labs(x='',y='Estimate') +
            scale_color_manual(values=c('gray','indianred3','dodgerblue3')) +
            scale_fill_brewer(type='qual',palette=6,direction=-1) +
            theme(legend.position='none',
                  axis.text.y=element_blank(),
                  axis.ticks.y=element_blank())


          list(p_est=p_est,covariate=covariate,bars=bars,topic_ordered=topic_order)

        })
      })

      output$est <- renderPlotly({

        ggplotly(EST()$p_est)

      })

      output$tax <- renderPlotly({

        if (input$dist_tax == 'hellinger'){

          d <- vegan::vegdist(vegan::decostand(beta,'norm'),method='euclidean')

        }else if (input$dist_tax == 'chi2'){

          d <- vegan::vegdist(vegan::decostand(beta,'chi.square'),method='euclidean')

        }else if (input$dist_tax == 'jsd'){

          if (min(beta) == 0){
            min_beta <- min(beta[beta > 0])
            pbeta <- beta + min_beta
            pbeta <- pbeta/colSums(pbeta)
          }

          d <- proxy::dist(pbeta,jsd)

        }else{

          d <- vegan::vegdist(beta,method=input$dist_tax)

        }

        dd_row <- as.dendrogram(hclust(d,method='ward.D2'))

        logbeta <- logbeta_global[order.dendrogram(dd_row),]
        df <- data.frame(probability=matrix(logbeta,ncol=1),
                         otu=rownames(logbeta),
                         topic=rep(colnames(logbeta),each=nrow(logbeta)))
        df$taxon <- pretty_taxa_names(taxa[df$otu,])
        df$taxon <- paste0(df$taxon,' (',df$otu,')')

        df$topic <- factor(df$topic,levels=EST()$topic_order,ordered=TRUE)
        df$otu <- factor(df$otu,levels=unique(df$otu),ordered=TRUE)
        df$taxon <- factor(df$taxon,levels=unique(df$taxon),ordered=TRUE)

        p_tax <- ggplot(df,aes(topic,taxon)) +
          geom_raster(aes(fill=probability)) +
          viridis::scale_fill_viridis(name='Probability Rank') +
          labs(x='',y='',fill='Probability Rank') +
          geom_vline(xintercept=EST()$bars,color='white',size=2) +
          theme(legend.position='none',
                axis.title.y=element_text(face='bold',size=20),
                axis.ticks.y=element_blank(),
                axis.text.y=element_blank(),
                axis.text.x=element_text(angle=-90,hjust=0,vjust=.5))

        ggplotly(p_tax)
      })

      output$fxn <- renderPlotly({

        df1$topic <- factor(df1$topic,levels=EST()$topic_order,ordered=TRUE)
        df2$topic <- factor(df2$topic,levels=EST()$topic_order,ordered=TRUE)

        p_fxn <- ggplot(df1,aes(topic,pw,fill=weight)) +
          geom_raster() +
          geom_point(data=df2[df2$sig > 0,],
                     aes(topic,pw),color='red',shape=3)

        if (uncertainty)  p_fxn <- p_fxn + geom_point(data=df2[df2$sig < 0,],aes(topic,pw),color='cyan',shape=3)

        p_fxn <- p_fxn + viridis::scale_fill_viridis(name='Coefficient') +
          labs(x='',y='',fill='Coefficient') +
          theme(axis.text.x=element_text(angle=-90,hjust=0,vjust=.5),
                legend.position='left') +
          theme(legend.position='none') +
          geom_vline(xintercept=EST()$bars,color='white',size=2)

        ggplotly(p_fxn,source='hm_fxn')

      })


      output$text <- renderUI({

        HTML(sprintf('The choices on the <b>left</b> sets the type of distance metric to use for hierarchical clustering in the taxa heatmap <b>below</b>.
                      The functional heatmap is interactive; by clicking a given cell, a table will appear at the bottom of the page reporting the gene
                      set that makes up that topic-pathway combination. Crosses on the functional heatmap cells indicates weights that did not enclose 0
                      at the %s%% uncertaintly level.',round(100*ui_level,0)))

      })


      show_table <- reactiveValues(gene_table=NULL)

      observeEvent(event_data('plotly_click',source='hm_fxn'),{

        s <- event_data('plotly_click',source='hm_fxn')

        k <- levels(df1$topic)[s[['x']]]
        pw <- levels(df1$pw)[s[['y']]]

        gene_table <- function_effects$gene_table[paste0('T',function_effects$gene_table$topic) == k & function_effects$gene_table$pw == pw,c('count','ko','description')]
        gene_table <- gene_table[gene_table$count >= 10,]
        gene_table <- gene_table[order(gene_table$count,decreasing=TRUE),]
        colnames(gene_table) <- c('Count','ID','Description')

        show_table$gene_table <- gene_table

      })

      observeEvent(input$choose,{

        show_table$gene_table <- NULL

      })


      output$tbl <- renderTable({

        validate(need(!is.null(show_table$gene_table),'Please explore a gene set by clicking a cell in the heatmap.'))

        show_table$gene_table


      },digits=0,width='900',align='l',rownames=FALSE,striped=TRUE,hover=TRUE,bordered=FALSE,spacing='s')

    }



  )

}
