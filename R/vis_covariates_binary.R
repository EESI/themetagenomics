#' @rdname vis
#'
#' @param taxa_grp_n Number of taxa group names to display (remaining are renamed to other). Defaults to 7.
#'
#' @export
vis.binary <- function(object,taxa_grp_n=7,...){

  topics <- object$topics
  topic_effects <- object$topic_effects
  otu_table <- object$topics$otu_table
  tax_table <- object$topics$tax_table
  metadata <- object$modelframe

  fit <- topics$fit
  K <- fit$settings$dim$K
  beta <- t(exp(fit$beta$logbeta[[1]]))
  rownames(beta) <- fit$vocab
  colnames(beta) <- paste0('T',1:K)

  otu_table <- otu_table + 1
  ra_table <- otu_table/rowSums(otu_table)

  pretty_names <- pretty_taxa_names(tax_table)
  taxa_other <- rename_taxa_to_other(otu_table,tax_table,top_n=taxa_grp_n)

  covariates <- colnames(metadata)
  names(covariates) <- tolower(names(topic_effects))
  cov_f <- sapply(metadata,class) == 'factor'
  covariates <- covariates[cov_f]

  if (length(covariates) == 0) stop('For binary, formula provided must have contained a binary, categorical, or factor covariate.')

  est_range <- range(c(sapply(topic_effects, function(x) unlist(x$est))))

  s <- NULL

  shinyApp(

    ui <- fluidPage(

      tags$head(tags$style(type='text/css','.side{font-size: 10px;} .side{color: gray;} .side{font-weight: bold;}')),
      tags$head(tags$style(type='text/css','.below{font-size: 10px;} .below{color: gray;} .below{font-weight: bold;}')),
      tags$head(tags$style(type='text/css','.capt{font-size: 9px;} .capt{color: gray;} .capt{font-weight: bold;} .capt{margin-top: -20px;}')),

      titlePanel('Topic-Covariate Effects'),

      fixedRow(
        column(2,selectInput('choose', label='Covariate',
                             choices=covariates,selected=covariates[[1]]),
               fixedRow(column(1,''),
                        column(11,tags$div(paste0('Choosing a covariate determines which weight estimates will shown',
                                                  ' The order of the topics will be adjusted accordingly. By clicking',
                                                  ' an estimate, all figures below will rerender.'),class='side')))),
        column(10,plotlyOutput('est',height='200px'))
      ),

      fixedRow(
        column(2,
               conditionalPanel(condition='output.show',sliderInput('z', label='Min. beta prob.',min=-6,max=-1,value=-2,step=.25,pre='10^(',post=')'))
        ),
        column(8,htmlOutput('text')),
        column(2,
               conditionalPanel(condition='output.show',checkboxInput('norm',label=strong('Normalize'),value=TRUE))
        )
      ),

      plotlyOutput('dis'),

      plotOutput('highlight')

    ),



    server = function(input,output){

      EST <- reactive({
        suppressWarnings({

          covariate <- input$choose

          cov_names <- levels(unique(metadata[,covariate]))
          cov_ids <- list(rownames(otu_table)[metadata[,covariate] == cov_names[1]],rownames(otu_table)[metadata[,covariate] == cov_names[2]])
          names(cov_ids) <- cov_names
          cov_coverage <- c(sum(otu_table[cov_ids[[1]],]),sum(otu_table[cov_ids[[2]],]))
          names(cov_coverage) <- cov_names
          cov_list <- list(ids=cov_ids,coverage=cov_coverage)

          est_mat <- topic_effects[[covariate]][['est']]

          df <- data.frame(topic=paste0('T',1:K),
                           est=est_mat[,1],
                           lower=est_mat[,2],
                           upper=est_mat[,3],
                           sig=ifelse(1:K %in% topic_effects[[covariate]][['sig']],'1','0'),
                           stringsAsFactors=TRUE)[order(topic_effects[[covariate]][['rank']]),]
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

          df0 <- data.frame(otu=c('x','y'),
                            p=c(0,1),
                            covariate=unique(metadata[,covariate]),
                            taxon=c('x','y'),
                            stringsAsFactors=TRUE)


          list(p_est=p_est,k_levels=levels(df$topic),df0=df0,cov_list=cov_list,cov_names=cov_names,covariate=covariate)

        })
      })

      output$est <- renderPlotly({

        ggplotly(EST()$p_est,source='reg_est')

      })

      show_slider <- reactiveValues(k=FALSE)

      observeEvent(event_data('plotly_click',source='reg_est'),{

        show_slider$k <- TRUE

      })

      DF1 <- reactive({
        suppressWarnings({

          s <- event_data('plotly_click',source='reg_est')

          validate(need(!is.null(s),'Please choose a topic by clicking a point.'))

          df0 <- EST()$df0

          k <- EST()$k_levels[s[['x']]]

          beta_subset <- beta[,k]
          beta_subset <- beta_subset[beta_subset > 10^(as.integer(input$z))]

          if (input$norm) x <- ra_table else x <- otu_table
          df0$abundance <- c(0,max(x))

          otu_subset <- t(x[,names(beta_subset)])

          df1 <- try(data.frame(abundance=matrix(otu_subset,ncol=1),
                                otu=rownames(otu_subset),
                                p=beta_subset,
                                sample=rep(colnames(otu_subset),each=nrow(otu_subset)),
                                stringsAsFactors=TRUE),
                     silent=TRUE)

          validate(need(!(class(df1)[1] == 'try-error'),'Too many points filtered. Lower the minimum beta probability.'))

          df1$covariate <- metadata[df1$sample,EST()$covariate]
          df1$covariate <- factor(df1$covariate,levels=EST()$cov_names,ordered=TRUE)
          df1$taxon <- paste0(pretty_names[df1$otu],' (',df1$otu,')')

          df2 <- df1[df1$abundance > min(df1$abundance),]

          p_dis <- ggplot(data=df0,aes_(~covariate,~abundance,color=~p)) +
            geom_blank() +
            geom_violin(data=df2,
                        aes_(~covariate,~abundance,color=~p),fill=NA) +
            geom_jitter(data=df1,
                        aes_(~covariate,~abundance,color=~p,labels=~taxon),alpha=.5,size=2) +
            scale_y_log10(labels=scales::comma) +
            viridis::scale_color_viridis('Beta Prob.') +
            theme_classic() +
            labs(x='',y='Abundance') +
            ggtitle(k)

          p_dis <- ggplotly(p_dis,source='dis_cov',tooltip=c('taxon','abundance','p'))
          p_dis <- layout(p_dis,dragmode='select')

          list(p_dis=p_dis,df1=df1,df0=df0)

        })
      })


      output$show <- reactive({
        show_slider$k
      })
      outputOptions(output,'show',suspendWhenHidden=FALSE)


      observeEvent(event_data('plotly_click',source='reg_est'),
                   {
                     output$text <- renderUI({
                       HTML('The slider on the <b>left</b> adjusts the value in which taxa are filtered based on their probability in the topics over
                            taxa distribution beta. The check box sets whether to shown counts or relative abundances. The violin plots <b>below</b>
                            are interactive. Select a subset of points to visualize their distribution in the raw data as a function on taxonomy.')
                     })
                     }
                     )


      output$dis <- renderPlotly({

        validate(need(length(DF1()$p_dis),''))

        DF1()$p_dis

      })

      output$highlight <- renderPlot({

        h <- event_data('plotly_selected',source='dis_cov')

        validate(need(length(h),''))

        h_otu <- gsub('^taxon:.*\\(([0-9]+)\\)$','\\1',DF1()$df1$otu)[h$pointNumber + 1] # indexed from 0

        df_tax <- try({

          df_tax <- sum_taxa_by_group(h_otu,taxa_other,otu_table,metadata,EST()$cov_list,sample_norm=input$norm)
          df_tax$cov <- factor(df_tax$cov,levels=EST()$cov_names,ordered=TRUE)

          df_tax

        },silent=TRUE)

        validate(need(!(class(df_tax)[1] == 'try-error'),'Too few or no points selected. Drag select points in the ribbon plots.'))

        p_tax <- ggplot(df_tax,aes_(~taxon,~abundance,fill=~cov)) +
          geom_bar(color='black',stat='identity',position='dodge') +
          facet_grid(.~group,scales='free_x') +
          scale_fill_manual(values=c('indianred3','dodgerblue3')) +
          theme_classic() +
          theme(axis.text.x=element_text(angle=45,hjust=1,size=16),
                strip.text=element_text(face='bold',size=16)) +
          labs(x='',y='abundance',fill='')

        p_tax


      })


                   }


      )

}
