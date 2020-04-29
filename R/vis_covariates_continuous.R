#' @rdname vis
#'
#' @param taxa_reg_n Number of most relevant taxa within topic to regress. Defaults to 8.
#'
#' @export
vis.continuous <- function(object,lambda_step=.1,taxa_reg_n=8,...){

  topics <- object$topics
  topic_effects <- object$topic_effects
  otu_table <- object$topics$otu_table
  tax_table <- object$topics$tax_table
  metadata <- object$modelframe

  covariates <- colnames(metadata)
  cov_f <- sapply(metadata,class) == 'factor'

  cov_cont <- covariates[!cov_f]

  if (length(cov_cont) == 0) stop('For continuous, formula provided must have contained a continuous, numeric, or integer covariate.')

  cov_fact <- covariates[cov_f]
  names(cov_cont) <- tolower(unlist(cov_cont))

  mod_fact <- names(topic_effects[[cov_cont[1]]]$fitted_switch)

  pretty_names <- pretty_taxa_names(tax_table)
  otu_table <- otu_table + 1

  fit <- topics$fit
  vocab <- fit$vocab
  K <- fit$settings$dim$K

  doc_lengths <- sapply(topics$docs,function(x) sum(x[2,]))

  theta <- fit$theta
  rownames(theta) <- names(topics$docs)
  colnames(theta) <- paste0('T',1:K)
  if (min(theta) == 0){
    min_theta <- min(theta[theta > 0])
    theta <- theta + min_theta
    theta <- theta/rowSums(theta)
  }

  beta <- exp(fit$beta$logbeta[[1]])
  colnames(beta) <- vocab
  rownames(beta) <- paste0('T',1:K)
  if (min(beta) == 0){
    min_beta <- min(beta[beta > 0])
    beta <- beta + min_beta
    beta <- beta/rowSums(beta)
  }

  topic_freq <- colSums(theta*doc_lengths)
  topic_marg <- topic_freq/sum(topic_freq)
  term_topic_freq <- beta*topic_freq

  # compute term frequencies as column sums of term.topic.frequency
  # we actually won't use the user-supplied term.frequency vector.
  # the term frequencies won't match the user-supplied frequencies exactly
  # this is a work-around to solve the bug described in Issue #32 on github:
  # https://github.com/cpsievert/LDAvis/issues/32
  term_freq <- colSums(term_topic_freq) # topic_freq <- colSums(theta*doc_lengths)

  # marginal distribution over terms (width of blue bars)
  term_marg <- term_freq/sum(term_freq)

  beta <- t(beta)

  # compute the distinctiveness and saliency of the terms:
  # this determines the R terms that are displayed when no topic is selected
  topic_g_term <- beta/rowSums(beta)  # (W x K)
  kernel <- topic_g_term*log(t(t(topic_g_term)/topic_marg))
  distinctiveness <- rowSums(kernel)
  saliency <- term_marg*distinctiveness

  # Order the terms for the "default" view by decreasing saliency:
  default_terms <- vocab[order(saliency,decreasing=TRUE)][1:taxa_reg_n]
  counts <- as.integer(term_freq[match(default_terms,vocab)])
  taxa_n_rev <- rev(seq_len(taxa_reg_n))

  default <- data.frame(Term=default_terms,
                        logprob=taxa_n_rev,
                        loglift=taxa_n_rev,
                        Freq=counts,
                        Total=counts,
                        Category='Default',
                        stringsAsFactors=FALSE)
  default$Term <- factor(default$Term,levels=rev(default$Term),ordered=TRUE)

  topic_seq <- rep(seq_len(K),each=taxa_reg_n)
  category <- paste0('T',topic_seq)
  lift <- beta/term_marg

  # Collect taxa_reg_n most relevant terms for each topic/lambda combination
  # Note that relevance is re-computed in the browser, so we only need
  # to send each possible term/topic combination to the browser
  find_relevance <- function(i){
    relevance <- i*log(beta) + (1 - i)*log(lift)
    idx <- apply(relevance,2,function(x) order(x,decreasing=TRUE)[seq_len(taxa_reg_n)])
    indices <- cbind(c(idx),topic_seq)
    data.frame(Term=vocab[idx],
               Category=category,
               logprob=round(log(beta[indices]),4),
               loglift=round(log(lift[indices]),4),
               stringsAsFactors=FALSE)
  }

  lambda_seq <- seq(0,1,by=lambda_step) # 0.01
  tinfo <- lapply(as.list(lambda_seq),find_relevance)
  tinfo <- unique(do.call('rbind', tinfo))
  tinfo$Total <- term_freq[match(tinfo$Term,vocab)]
  rownames(term_topic_freq) <- paste0('T',seq_len(K))
  colnames(term_topic_freq) <- vocab
  tinfo$Freq <- term_topic_freq[as.matrix(tinfo[c('Category','Term')])]
  tinfo <- rbind(default[,-7],tinfo)

  new_order <- default


  shinyApp(

    ui <- fluidPage(

      tags$head(tags$style(type='text/css','.side{font-size: 10px;} .side{color: gray;} .side{font-weight: bold;}')),
      tags$head(tags$style(type='text/css','.below{font-size: 10px;} .below{color: gray;} .below{font-weight: bold;}')),
      tags$head(tags$style(type='text/css','.capt{font-size: 9px;} .capt{color: gray;} .capt{font-weight: bold;} .capt{margin-top: -10px;}')),

      titlePanel('Topic-Covariate Effects'),

      fixedRow(
        column(2,selectInput('choose_cov', label='Covariate',
                             choices=cov_cont,selected=cov_cont[[1]]),
               fixedRow(column(1,''),
                        column(11,tags$div(paste0('Choosing a covariate determines which weight estimates will shown',
                                                  ' The order of the topics will be adjusted accordingly. By clicking',
                                                  ' an estimate, all figures below will rerender.'),class='side')))),
        column(10,plotlyOutput('est',height='200px'))
      ),

      fixedRow(
        column(2,''),
        column(8,htmlOutput('text')),
        column(2,'')
      ),

      plotlyOutput('theta'),

      fixedRow(
        column(3,
               conditionalPanel(condition='output.show_mods',
                               selectInput('choose_mod',label='Group',
                                           choices=mod_fact,selected=mod_fact[[1]]))),
        column(6,
               conditionalPanel(condition='output.show',
                                sliderInput('lambda',label=strong('Lambda'),min=0,max=1,value=1,step=lambda_step))),
        column(3,
               conditionalPanel(condition='output.show',
                                checkboxGroupInput('z',
                                                   label='',
                                                   c('Ignore zeros'=1,'Normalize'=2,'Split'=3),
                                                   selected=c(1,2))))
      ),

      fixedRow(
        column(3,conditionalPanel(condition='output.show_mods',
                                  tags$div(paste0('Choice of factor to split scatter plots.'),class='capt'))),
        column(6,conditionalPanel(condition='output.show',
                                  tags$div(paste0('Relative weighting that influences taxa shown in the scatterplots',
                                 ' If equal to 1, p(taxa|topic)l if 0, p(taxa|topic)/p(taxa).'),class='capt'))),
        column(3,conditionalPanel(condition='output.show',
                                  tags$div(paste0('Adjust scatter plots by ignoring zeros for smoothing, normalizing to relative abundances,',
                                 ' or splitting the plots as a function of the selected factor (if present).'),class='capt')))
      ),

      plotlyOutput('ss'),

      plotlyOutput('full')

    ),


    server = function(input,output,session){


      current_mods <- reactiveValues(mods=mod_fact)

      observeEvent(input$choose_cov,{

        covariate <- input$choose_cov
        mods <- names(topic_effects[[covariate]]$fitted_switch)
        current_mods$mods <- mods
        updateNumericInput(session,'k_in',value=mods)

      })


      EST <- reactive({
        suppressWarnings({

          covariate <- input$choose_cov

          df0 <- data.frame(p=c(0,max(theta)),
                            covariate=range(metadata[,covariate]),
                            stringsAsFactors=TRUE)

          est_mat <- topic_effects[[covariate]]$est

          df <- data.frame(topic=paste0('T',1:K),
                           est=est_mat[,1],
                           lower=est_mat[,2],
                           upper=est_mat[,3],
                           sig=ifelse(1:K %in% topic_effects[[covariate]]$sig,'1','0'),
                           stringsAsFactors=TRUE)[order(topic_effects[[covariate]]$rank),]
          df$sig <- factor(as.character(sign(df$est) * as.numeric(as.character.factor(df$sig))),levels=c('0','1','-1'),ordered=TRUE)
          df$topic <- factor(df$topic,levels=df$topic,ordered=TRUE)

          p_est <- ggplot(df,aes_(~topic,y=~est,ymin=~lower,ymax=~upper,color=~sig)) +
            geom_hline(yintercept=0,linetype=3)
          p_est <- p_est +
            geom_pointrange(size=2) +
            theme_minimal() +
            labs(x='',y='Estimate') +
            scale_color_manual(values=c('gray','indianred3','dodgerblue3'),drop=FALSE) +
            scale_fill_brewer(type='qual',palette=6,direction=-1) +
            theme(legend.position='none',
                  axis.text.x=element_text(angle=-90,hjust=0,vjust=.5))

          list(p_est=p_est,k_levels=levels(df$topic),df0=df0,covariate=covariate)

        })
      })


      output$est <- renderPlotly({

        ggplotly(EST()$p_est,source='reg_est')

      })

      TAB <- reactive({

        if (any(input$z %in% '2')) x <- otu_table/rowSums(otu_table) else x <- otu_table

        df_temp <- EST()$df0
        df_temp$abundance <- c(min(x),max(x))

        list(table=x,df=df_temp)

      })


      P <- reactive({

        suppressWarnings({

          s <- event_data('plotly_click',source='reg_est')

          validate(need(!is.null(s),''))

          current_k <- EST()$k_levels[s[['x']]]

          l <- input$lambda

          tinfo_k <- tinfo[tinfo$Category == current_k,] #subset(tinfo,Category == current_k)
          rel_k <- l*tinfo_k$logprob + (1-l)*tinfo_k$loglift
          new_order <- tinfo_k[order(rel_k,decreasing=TRUE)[1:taxa_reg_n],]
          new_order$Term <- as.character.factor(new_order$Term)

          otu_subset <- t(TAB()$table[,new_order$Term])

          df_temp <- data.frame(abundance=matrix(otu_subset,ncol=1),
                                otu=rownames(otu_subset),
                                sample=rep(colnames(otu_subset),each=nrow(otu_subset)),
                                stringsAsFactors=FALSE)
          df_temp$covariate <- metadata[df_temp$sample,EST()$covariate]
          df_temp$switch <- metadata[df_temp$sample,input$choose_mod]
          df_temp$p <- theta[df_temp$sample,current_k]
          df_temp$taxon <- paste0(pretty_names[df_temp$otu],' (',df_temp$otu,')')
          df_temp$taxon <- factor(df_temp$taxon,levels=unique(df_temp$taxon),ordered=TRUE)

          p_full <- ggplot(data=TAB()$df,aes_(~covariate,~abundance,color=~p)) +
            geom_blank() +
            geom_point(data=df_temp,aes_(~covariate,~abundance,color=~p),alpha=.5,size=1.2) +
            scale_y_log10(labels=scales::comma) +
            viridis::scale_color_viridis('Theta Prob.') +
            theme_classic() +
            theme(aspect.ratio=.5) +
            labs(x=EST()$covariate,y='Abundance')

          p_ss <- p_full +
            facet_wrap(~taxon,ncol=4) +
            theme(aspect.ratio=.5,strip.text.x=element_text(size=7)) +
            xlab('')

          list(df=df_temp,min=min(df_temp$abundance),p_ss=p_ss,p_full=p_full)

        })

      })


      observeEvent(event_data('plotly_click',source='reg_est'),
                   {
                     output$text <- renderUI({
                       HTML(sprintf("The scatterplot below shows the sample over topic distribution as a function of %s.
                                    If additional factors were present in the model formula,
                                    then their posterior predictive estimates will be shown as a function of the given
                                    factor being set to 1 and 0 with all other covariates held fixed (i.e., averaged across samples).
                                    These factors can be cycled through via the selection box. Each point represents a sample's probability of
                                    containing the selected topic.
                                    The next set of scatterplots show the top %s highest probability taxa in the selected topic -- that is,
                                    p(taxa|k). Each point is the taxa's abundance in the raw data for a given sample, shaded based on the probability of that
                                    sample occuring in the chosen topic. By adjusting lambda, the top %s taxa in terms of p(taxa|k)/p(taxa) will be shown.
                                    The large scatter plot in the bottom row shows all %s taxa combined.",
                                    EST()$covariate,
                                    taxa_reg_n,
                                    taxa_reg_n,
                                    taxa_reg_n))
                     })
                   }
                       )

      output$ss <- renderPlotly({

        suppressWarnings({

          if (any(input$z %in% '1')){

            df2 <- P()$df[P()$df$abundance > round(P()$min,5),]

            p_ss <- P()$p_ss +
              stat_smooth(data=df2,aes_(~covariate,~abundance),color='black',size=1.1,method='loess',se=FALSE) +
              stat_smooth(data=df2,aes_(~covariate,~abundance),color='darkred',size=1.1,method='lm',linetype=2,se=FALSE)

          }else{

            p_ss <- P()$p_ss +
              stat_smooth(data=P()$df,aes_(~covariate,~abundance),color='black',size=1.1,method='loess',se=FALSE) +
              stat_smooth(data=P()$df,aes_(~covariate,~abundance),color='darkred',size=1.1,method='lm',linetype=2,se=FALSE)

          }

          if (any(input$z %in% '3') & input$choose_mod != '')
            p_ss <- p_ss + facet_wrap(switch~taxon,ncol=6) + theme(strip.text.x=element_text(size=7,margin=margin(t=15,b=15)))

          ggplotly(p_ss)

        })

      })

      output$theta <- renderPlotly({

        suppressWarnings({

          s <- event_data('plotly_click',source='reg_est')

          validate(need(!is.null(s),'Please choose a topic by clicking a point.'))

          k <- EST()$k_levels[s[['x']]]

          if (input$choose_mod != ''){
            topic_fitted <- topic_effects[[EST()$covariate]]$fitted[[input$choose_mod]][[k]]
            topic_switch <- topic_effects[[EST()$covariate]]$fitted_switch[[input$choose_mod]][[k]]
            df_th <- data.frame(Value=theta[rownames(metadata),k],Covariate=metadata[,EST()$covariate],stringsAsFactors=TRUE)
            df_fit <- data.frame(topic_fitted,stringsAsFactors=TRUE)
            df_switch <- data.frame(topic_switch,stringsAsFactors=TRUE)
            colnames(df_fit) <- c('Value','Lower','Upper','Covariate')
            colnames(df_switch) <- c('Value','Lower','Upper','Covariate')
            R <- range(c(df_th$Covariate,df_fit$covariate,df_switch$covariate))
          }else{
            topic_fitted <- topic_effects[[EST()$covariate]]$fitted[[1]][[k]]
            df_th <- data.frame(Value=theta[rownames(metadata),k],Covariate=metadata[,EST()$covariate],stringsAsFactors=TRUE)
            df_fit <- data.frame(topic_fitted,stringsAsFactors=TRUE)
            colnames(df_fit) <- c('Value','Lower','Upper','Covariate')
            R <- range(c(df_th$Covariate,df_fit$covariate))
          }


          p1 <- plot_ly(df_th,x=df_th$Covariate)
          p1 <- add_trace(p1,
                          x=~df_th$Covariate,y=~df_th$Value,
                          type='scatter',mode='markers',
                          marker=list(symbol='circle',opacity=.5,sizemode='diameter',size=10,color=I('gray')),
                          showlegend=FALSE)
          p1 <- add_lines(p1,
                          x=~df_fit$Covariate,y=~df_fit$Lower,
                          type='scatter',mode='lines',
                          line=list(color='transparent'),
                          showlegend=FALSE,
                          name=input$choose_mod)
          p1 <- add_lines(p1,
                          x=~df_fit$Covariate,y=~df_fit$Upper,
                          type='scatter',mode='lines',
                          fill='tonexty',
                          fillcolor='rgba(205,85,85,0.2)',
                          line=list(color='transparent'),
                          showlegend=FALSE,
                          name=input$choose_mod)
          p1 <- add_lines(p1,
                          x=~df_fit$Covariate,y=~df_fit$Value,
                          type='scatter',mode='lines',
                          line = list(width=5,color = 'rgba(205,85,85,0.9)'),
                          showlegend=TRUE,
                          name=input$choose_mod)

          if (input$choose_mod != ''){
            p1 <- add_lines(p1,
                            x=~df_switch$Covariate,y=~df_switch$Lower,
                            type='scatter',mode='lines',
                            line=list(color='transparent'),
                            showlegend=FALSE,
                            name='Reference/Control')
            p1 <- add_lines(p1,
                            x=~df_switch$Covariate,y=~df_switch$Upper,
                            type='scatter',mode='lines',
                            fill='tonexty',
                            fillcolor='rgba(24,116,205,.2)',
                            line=list(color='transparent'),
                            showlegend=FALSE,
                            name='Reference/Control')
            p1 <- add_lines(p1,
                            x=~df_switch$Covariate,y=~df_switch$Value,
                            type='scatter',mode='lines',
                            line=list(width=5,color='rgba(24,116,205,.9)'),
                            showlegend=TRUE,
                            name='Reference/Control')
            p1 <- layout(p1,title=k,
                         legend=list(orientation='h'),
                         paper_bgcolor='rgb(255,255,255)',
                         xaxis=list(title=sprintf('%s',EST()$covariate),
                                      gridcolor='rgb(255,255,255)',
                                      showgrid=TRUE,
                                      showline=FALSE,
                                      showticklabels=TRUE,
                                      tickcolor='rgb(127,127,127)',
                                      ticks='outside',
                                      zeroline=FALSE),
                         yaxis=list(title=sprintf('p(topic %s|sample)',k),
                                      gridcolor='rgb(255,255,255)',
                                      showgrid=TRUE,
                                      showline=FALSE,
                                      showticklabels=TRUE,
                                      tickcolor='rgb(127,127,127)',
                                      ticks='outside',
                                      zeroline=FALSE))
          }else{
            p1 <- layout(p1,title=k,
                         showlegend=FALSE,
                         paper_bgcolor='rgb(255,255,255)',
                         xaxis=list(title=sprintf('%s',EST()$covariate),
                                    gridcolor='rgb(255,255,255)',
                                    showgrid=TRUE,
                                    showline=FALSE,
                                    showticklabels=TRUE,
                                    tickcolor='rgb(127,127,127)',
                                    ticks='outside',
                                    zeroline=FALSE),
                         yaxis=list(title=sprintf('p(topic %s|sample)',k),
                                    gridcolor='rgb(255,255,255)',
                                    showgrid=TRUE,
                                    showline=FALSE,
                                    showticklabels=TRUE,
                                    tickcolor='rgb(127,127,127)',
                                    ticks='outside',
                                    zeroline=FALSE))
          }
          p1

        })

      })

      output$show_mods <- reactive({
        s <- event_data('plotly_click',source='reg_est')
        !is.null(s) & !is.null(mod_fact)
      })
      outputOptions(output,'show_mods',suspendWhenHidden=FALSE)

      output$show <- reactive({
        s <- event_data('plotly_click',source='reg_est')
        !is.null(s)
      })
      outputOptions(output,'show',suspendWhenHidden=FALSE)

      output$full <- renderPlotly({

        suppressWarnings({

          if (any(input$z %in% '1')){

            df2 <- P()$df[P()$df$abundance > round(P()$min,5),]

            p_full <- P()$p_full +
              stat_smooth(data=df2,aes_(~covariate,~abundance),color='black',size=1.1,method='loess',se=FALSE) +
              stat_smooth(data=df2,aes_(~covariate,~abundance),color='darkred',size=1.1,method='lm',linetype=2,se=FALSE)

          }else{

            p_full <- P()$p_full +
              stat_smooth(data=P()$df,aes_(~covariate,~abundance),color='black',size=1.1,method='loess',se=FALSE) +
              stat_smooth(data=P()$df,aes_(~covariate,~abundance),color='darkred',size=1.1,method='lm',linetype=2,se=FALSE)
          }

          if (any(input$z %in% '3') & input$choose_mod != '')
            p_full <- p_full + facet_wrap(~switch,ncol=1)

          ggplotly(p_full)

        })

      })


   }


 )

}
