#' @import ggplot2 shiny plotly
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

vis_covariate_effects_continuous <- function(topics,topic_effects,otu_table,taxa,moderator,lambda_step=.01,taxa_n=12){

  metadata <- topic_effects$modelframe
  topic_effects <- topic_effects$topic_effects

  covariates <- colnames(metadata)
  cov_f <- sapply(metadata,class) == 'factor'

  # cov_f <- sapply(topic_effects,function(x) all(is.na(x$fitted)))
  # covariates <- lapply(names(topic_effects),identity)
  cov_fact <- moderator
  cov_cont <- covariates[!cov_f]
  names(cov_fact) <- tolower(cov_fact)
  names(cov_cont) <- tolower(unlist(cov_cont))

  pretty_names <- pretty_taxa_names(taxa)
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
  default_terms <- vocab[order(saliency,decreasing=TRUE)][1:taxa_n]
  counts <- as.integer(term_freq[match(default_terms,vocab)])
  taxa_n_rev <- rev(seq_len(taxa_n))

  default <- data.frame(Term=default_terms,
                        logprob=taxa_n_rev,
                        loglift=taxa_n_rev,
                        Freq=counts,
                        Total=counts,
                        Category='Default',
                        stringsAsFactors=FALSE)
  default$Term <- factor(default$Term,levels=rev(default$Term),ordered=TRUE)

  topic_seq <- rep(seq_len(K),each=taxa_n)
  category <- paste0('T',topic_seq)
  lift <- beta/term_marg

  # Collect taxa_n most relevant terms for each topic/lambda combination
  # Note that relevance is re-computed in the browser, so we only need
  # to send each possible term/topic combination to the browser
  find_relevance <- function(i){
    relevance <- i*log(beta) + (1 - i)*log(lift)
    idx <- apply(relevance,2,function(x) order(x,decreasing=TRUE)[seq_len(taxa_n)])
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

      titlePanel('Topic-Covariate Effects'),

      fixedRow(
        column(2,selectInput('choose', label='Covariate',
                             choices=cov_cont,selected=cov_cont[[1]])),
        column(10,plotlyOutput('est',height='200px'))
      ),

      fixedRow(
        column(2,''),
        column(8,htmlOutput('text')),
        column(2,'')
      ),

      plotlyOutput('theta'),

      fixedRow(
        column(4,
               conditionalPanel(condition='output.show',
                                checkboxGroupInput('z',
                                                   label='',
                                                   c('Ignore zeros'=1,'Normalize'=2,'Split'=3),
                                                   selected=c(1,2)))),
        column(4,
               conditionalPanel(condition='output.show',
                                sliderInput('lambda',label=strong('Lambda'),min=0,max=1,value=1,step=lambda_step))),
        column(4,'')
      ),

      plotlyOutput('ss'),

      plotlyOutput('full')

    ),


    server = function(input,output){

      EST <- reactive({
        suppressWarnings({

          covariate <- input$choose

          df0 <- data.frame(p=c(0,max(theta)),
                            covariate=range(metadata[,covariate]))

          est_mat <- topic_effects[[covariate]]$est

          df <- data.frame(topic=paste0('T',1:K),
                           est=est_mat[,1],
                           lower=est_mat[,2],
                           upper=est_mat[,3],
                           sig=ifelse(1:K %in% topic_effects[[covariate]]$sig,'1','0'))[order(topic_effects[[covariate]]$rank),]
          df$sig <- factor(as.character(sign(df$est) * as.numeric(as.character.factor(df$sig))),levels=c('0','1','-1'),ordered=TRUE)
          df$topic <- factor(df$topic,levels=df$topic,ordered=TRUE)

          p_est <- ggplot(df,aes(topic,y=est,ymin=lower,ymax=upper,color=sig)) +
            geom_hline(yintercept=0,linetype=3)
          p_est <- p_est +
            geom_pointrange(size=2) +
            theme_minimal() +
            labs(x='',y='Estimate') +
            scale_color_manual(values=c('gray','indianred3','dodgerblue3')) +
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

          tinfo_k <- subset(tinfo,Category == current_k)
          rel_k <- l*tinfo_k$logprob + (1-l)*tinfo_k$loglift
          new_order <- tinfo_k[order(rel_k,decreasing=TRUE)[1:taxa_n],]
          new_order$Term <- as.character.factor(new_order$Term)

          otu_subset <- t(TAB()$table[,new_order$Term])

          df_temp <- data.frame(abundance=matrix(otu_subset,ncol=1),
                                otu=rownames(otu_subset),
                                sample=rep(colnames(otu_subset),each=nrow(otu_subset)),
                                covariate=metadata[colnames(otu_subset),EST()$covariate],
                                switch=metadata[colnames(otu_subset),cov_fact])
          df_temp$p <- theta[df_temp$sample,current_k]
          df_temp$taxon <- paste0(pretty_names[df_temp$otu],' (',df_temp$otu,')')
          df_temp$taxon <- factor(df_temp$taxon,levels=unique(df_temp$taxon),ordered=TRUE)

          p_full <- ggplot(data=TAB()$df,aes(covariate,abundance,color=p)) +
            geom_blank() +
            geom_point(data=df_temp,aes(covariate,abundance,color=p),alpha=.5,size=1.2) +
            scale_y_log10(labels=scales::comma) +
            viridis::scale_color_viridis('Beta Prob.') +
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
                       HTML(sprintf("The check box on the <b>left</b> sets whether to ignore zeros for the best fit lines and whether raw counts or
                                    relative abundances will be shown. The scatter plots <b>below</b> show the top %s taxa in terms of their
                                    probability in the topics over taxa distributions versus %s. Each point represents the abundance of that taxa
                                    in that sample. Each point is colored based on the probability of that sample occuring in the chosen topic. The
                                    large scatter plot shows all %s taxa combined.",
                                    taxa_n,EST()$covariate,taxa_n))
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

          if (any(input$z %in% '3'))
            p_ss <- p_ss + facet_wrap(switch~taxon,ncol=6) + theme(strip.text.x=element_text(size=7,margin=margin(t=15,b=15)))

          ggplotly(p_ss)

        })

      })

      output$theta <- renderPlotly({

        suppressWarnings({

          s <- event_data('plotly_click',source='reg_est')

          validate(need(!is.null(s),'Please choose a topic by clicking a point.'))

          k <- EST()$k_levels[s[['x']]]
          topic_fitted <- topic_effects[[EST()$covariate]]$fitted[[k]]
          topic_switch <- topic_effects[[EST()$covariate]]$fitted_switch[[k]]

          df_th <- data.frame(Value=theta[rownames(metadata),k],Covariate=metadata[,EST()$covariate])
          df_fit <- data.frame(topic_fitted)
          df_switch <- data.frame(topic_switch)
          colnames(df_fit) <- c('Value','Lower','Upper','Covariate')
          colnames(df_switch) <- c('Value','Lower','Upper','Covariate')
          R <- range(c(df_th$Covariate,df_fit$covariate,df_switch$covariate))

          p_theta <- ggplot() +
            geom_point(data=df_th,aes(Covariate,Value),color='black',alpha=.7) +
            geom_ribbon(data=df_switch,aes(x=Covariate,ymin=Lower,ymax=Upper),alpha=0.2,color='dodgerblue3') +
            geom_line(data=df_switch,aes(Covariate,Value),color='dodgerblue3',size=1.45) +
            geom_ribbon(data=df_fit,aes(x=Covariate,ymin=Lower,ymax=Upper),alpha=0.2,color='indianred3') +
            geom_line(data=df_fit,aes(Covariate,Value),color='indianred3',size=1.45) +
            xlim(R[1],R[2]) +
            ylim(0,max(theta)) +
            theme_classic() +
            labs(x='',y='Theta') +
            ggtitle(k)

          ggplotly(p_theta)

        })

      })

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
              stat_smooth(data=df2,aes(covariate,abundance),color='black',size=1.1,method='loess',se=FALSE) +
              stat_smooth(data=df2,aes(covariate,abundance),color='darkred',size=1.1,method='lm',linetype=2,se=FALSE)

          }else{

            p_full <- P()$p_full +
              stat_smooth(data=P()$df,aes(covariate,abundance),color='black',size=1.1,method='loess',se=FALSE) +
              stat_smooth(data=P()$df,aes(covariate,abundance),color='darkred',size=1.1,method='lm',linetype=2,se=FALSE)
          }

          if (any(input$z %in% '3'))
            p_full <- p_full + facet_wrap(~switch,ncol=1)

          ggplotly(p_full)

        })

      })


                     }


                       )

                   }
