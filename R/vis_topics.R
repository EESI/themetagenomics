#' @rdname vis
#'
#' @export
vis.topics <- function(object,taxa_bar_n=30,top_n=7,method=c('huge','simple'),corr_thresh=.01,lambda_step=.01,...){

  tax_table <- object$tax_table

  method <- match.arg(method)

  fit <- object$fit
  K <- fit$settings$dim$K
  vocab <- fit$vocab
  taxa <- tax_table[vocab,]
  taxon <- paste0(pretty_taxa_names(taxa),' (',vocab,')')
  taxa <- rename_taxa_to_other(object$docs,taxa,top_n=top_n,type='docs',as_factor=TRUE)
  rownames(taxa) <- taxon

  corr <- stm::topicCorr(fit,method=method,cutoff=corr_thresh)

  colors <- c('gray35','gray60','gray85')

  doc_lengths <- sapply(object$docs,function(x) sum(x[2,]))

  theta <- fit$theta
  rownames(theta) <- names(object$docs)
  colnames(theta) <- paste0('T',1:K)
  if (min(theta) == 0){
    min_theta <- min(theta[theta > 0])
    theta <- theta + min_theta
    theta <- theta/rowSums(theta)
  }

  beta <- exp(fit$beta$logbeta[[1]])
  colnames(beta) <- taxon
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
  default_terms <- taxon[order(saliency,decreasing=TRUE)][1:taxa_bar_n]
  counts <- as.integer(term_freq[match(default_terms,taxon)])
  taxa_n_rev <- rev(seq_len(taxa_bar_n))

  default <- data.frame(Term=default_terms,
                        logprob=taxa_n_rev,
                        loglift=taxa_n_rev,
                        Freq=counts,
                        Total=counts,
                        Category='Default',
                        stringsAsFactors=FALSE)
  default$Taxon <- taxa[default$Term,'Phylum']
  default$Term <- factor(default$Term,levels=rev(default$Term),ordered=TRUE)

  topic_seq <- rep(seq_len(K),each=taxa_bar_n)
  category <- paste0('T',topic_seq)
  lift <- beta/term_marg

  # Collect taxa_bar_n most relevant terms for each topic/lambda combination
  # Note that relevance is re-computed in the browser, so we only need
  # to send each possible term/topic combination to the browser
  find_relevance <- function(i){
    relevance <- i*log(beta) + (1 - i)*log(lift)
    idx <- apply(relevance,2,function(x) order(x,decreasing=TRUE)[seq_len(taxa_bar_n)])
    indices <- cbind(c(idx),topic_seq)
    data.frame(Term=taxon[idx],
               Category=category,
               logprob=round(log(beta[indices]),4),
               loglift=round(log(lift[indices]),4),
               stringsAsFactors=FALSE)
  }

  lambda_seq <- seq(0,1,by=lambda_step) # 0.01
  tinfo <- lapply(as.list(lambda_seq),find_relevance)
  tinfo <- unique(do.call('rbind', tinfo))
  tinfo$Total <- term_freq[match(tinfo$Term,taxon)]
  rownames(term_topic_freq) <- paste0('T',seq_len(K))
  colnames(term_topic_freq) <- taxon
  tinfo$Freq <- term_topic_freq[as.matrix(tinfo[c('Category','Term')])]
  tinfo <- rbind(default[,-7],tinfo)

  new_order <- default


  shinyApp(

    ui <- fluidPage(

      tags$head(tags$style(type='text/css','.side{font-size: 10px;} .side{color: gray;} .side{font-weight: bold;}')),
      tags$head(tags$style(type='text/css','.below{font-size: 10px;} .below{color: gray;} .below{font-weight: bold;}')),
      tags$head(tags$style(type='text/css','.capt{font-size: 9px;} .capt{color: gray;} .capt{font-weight: bold;} .capt{margin-top: -20px;}')),

      titlePanel('Taxa Features'),

      fixedRow(
        column(1,''),
        column(10,htmlOutput('text1')),
        column(1,'')
      ),

      br(),

      fixedRow(
        column(1,radioButtons('dim',label=strong('Dim'),
                              choices=list('2D'='2d','3D'='3d'),
                              selected='2d')),
        column(3,selectInput('dist',label=strong('Method'),
                             choices=list('Bray Curtis'='bray','Jaccard'='jaccard','Euclidean'='euclidean',
                                          'Hellinger'='hellinger','Chi Squared'='chi2','Jensen Shannon'='jsd',
                                          't-SNE'='tsne'),
                             selected='jsd')),
        column(1,style='padding: 25px 0px;',actionButton('reset','Reset')),
        column(2,numericInput('k_in',label=strong('Topic Number'),value=0,min=0,max=K,step=1)),
        column(3,sliderInput('lambda',label=strong('Lambda'),min=0,max=1,value=1,step=lambda_step)),
        column(2,selectInput('taxon',label=strong('Taxon'),
                             choices=list('Phylum'='Phylum','Class'='Class','Order'='Order',
                                          'Family'='Family','Genus'='Genus')))
      ),

      fixedRow(
        column(1,tags$div('Number of components to plot.',class='capt')),
        column(3,tags$div('Type of distance and method for ordination.',class='capt')),
        column(1,tags$div('Reset topic selection.',class='capt')),
        column(2,tags$div('Current selected topic.',class='capt')),
        column(3,tags$div(paste0('Relative weighting that influences taxa shown in barplot.',
                                 ' If equal to 1, p(taxa|topic)l if 0, p(taxa|topic)/p(taxa).'),class='capt')),
        column(2,tags$div('Taxonomic group to dictate bar plot shading',class='capt'))
      ),

      fixedRow(
        column(6,offset=0,height='600px',plotlyOutput('ord')),
        column(6,offset=0,height='600px',plotOutput('bar'))
      ),

      fixedRow(
        column(6,tags$div(paste0('Ordination of the samples over topics distribution theta, colored according to',
                                 ' the weights shown in the scatter plot above. The radius of a given point',
                                 ' represents the marginal topic frequency. The amount of variation explained is',
                                 ' annotated on each axis.'),class='below')),
        column(6,tags$div(paste0('Bar plot representing the taxa frequencies. When no topic is selected, the overall',
                                 ' taxa frequencies are shown, colored based on the selected taxonomy and ordered in',
                                 ' in terms of saliency. When a topic is chosen, the red bars show the margina taxa',
                                 ' frequency within the selected topic, ordered in terms of relevency, which in turn',
                                 ' can be reweighted by adjusting the lambda slider.'),class='below'))
      ),

      br(),

      fixedRow(
        column(1,''),
        column(10,htmlOutput('text2')),
        column(1,'')
      ),

      networkD3::forceNetworkOutput('corr')

    ),


    server = function(input,output,session){

      REL <- reactive({

        if (show_topic$k != 0){

          current_k <- paste0('T',show_topic$k)

          l <- input$lambda

          tinfo_k <- tinfo[tinfo$Category == current_k,]  #subset(tinfo,Category == current_k)
          rel_k <- l*tinfo_k$logprob + (1-l)*tinfo_k$loglift
          new_order <- tinfo_k[order(rel_k,decreasing=TRUE)[1:taxa_bar_n],]
          new_order$Term <- as.character.factor(new_order$Term)
          new_order$Taxon <- taxa[new_order$Term,input$taxon]
          new_order$Term <- factor(new_order$Term,levels=rev(new_order$Term),ordered=TRUE)

        }else{

          new_order <- default
          new_order$Taxon <- taxa[as.character.factor(new_order$Term),input$taxon]

        }

        new_order

      })


      output$text1 <- renderUI({

        HTML(sprintf("Below are the results of a %s topic STM. The ordination of the topics over taxa distribution (left) and the frequencies of
                    the top %s taxa (in terms of saliency) across all topics. By selecting a topic, the relative
                    frequencies of the taxa within that topic are shown in red. The ordination figure can be shown in
                    either 2D or 3D and the ordination method can be adjusted. Lambda adjusts the relevance calculation.
                    Choosing the taxon adjusts the group coloring for the bar plot. Clicking Reset resets the topic selection.",
                     K,taxa_bar_n))
      })

      output$text2 <- renderUI({
        HTML(paste0('Below shows topic-to-topic correlations from the samples over topics distribution. The edges represent positive',
                    ' correlation between two topics, with the size of the edge reflecting to the magnitude of the correlation.',
                    ' The size of the nodes are consistent with the ordination figure, reflecting the marginal topic frequencies.'))
      })

      output$ord <- renderPlotly({

        beta <- t(beta)

        if (input$dist == 'hellinger'){

          d <- cmdscale(vegan::vegdist(vegan::decostand(beta,'norm'),method='euclidean'),3,eig=TRUE)

        }else if (input$dist == 'chi2'){

          d <- cmdscale(vegan::vegdist(vegan::decostand(beta,'chi.square'),method='euclidean'),3,eig=TRUE)

        }else if (input$dist == 'jsd'){

          d <- cmdscale(proxy::dist(beta,jsd),3,eig=TRUE)

        }else if (input$dist == 'tsne'){

          p <- 30
          d <- try(Rtsne::Rtsne(beta,3,theta=.5,perplexity=p),silent=TRUE)
          while(class(d)[1] == 'try-error'){
            p <- p-1
            d <- try(Rtsne::Rtsne(beta,3,theta=.5,perplexity=p),silent=TRUE)
          }
          if (p < 30) cat(sprintf('Performed t-SNE with perplexity = %s.\n',p))
          d$points <- d$Y
          d$eig <- NULL

        }else{

          d <- cmdscale(vegan::vegdist(beta,method=input$dist),3,eig=TRUE)

        }

        eig <- d$eig[1:3]/sum(d$eig)
        colnames(d$points) <- c('Axis1','Axis2','Axis3')
        df <- data.frame(d$points,stringsAsFactors=TRUE)
        df$topic <- paste0('T',1:K)
        df$marg <- topic_marg

        if (input$dim == '2d'){

          p1 <- plot_ly(df,source='ord_click')
          p1 <- add_trace(p1,
                          x=~Axis1,y=~Axis2,size=~marg,
                          type='scatter',mode='markers',sizes=c(5,125),
                          color=I('darkred'),opacity=.5,
                          marker=list(symbol='circle',sizemode='diameter',line=list(width=3,color='#FFFFFF')),
                          text=~paste('<br>Topic:',topic),hoverinfo='text')
          p1 <- layout(p1,
                       showlegend=FALSE,
                       xaxis=list(title=sprintf('Axis 1 [%.02f%%]',eig[1]*100),
                                  showgrid=FALSE),
                       yaxis=list(title=sprintf('Axis 2 [%.02f%%]',eig[2]*100),
                                  showgrid=FALSE),
                       paper_bgcolor='rgb(243, 243, 243)',
                       plot_bgcolor='rgb(243, 243, 243)')
          p1 <- add_annotations(p1,x=df$Axis1,y=df$Axis2,text=df$topic,showarrow=FALSE,
                                font=list(size=10))

        }

        if (input$dim == '3d'){

          p1 <- plot_ly(df,source='ord_click',
                        x=~Axis1,y=~Axis2,z=~Axis3,size=~marg,
                        type='scatter3d',mode='markers',sizes=c(5,125),
                        color=I(df$colors),opacity=.5,
                        marker=list(symbol='circle',sizemode='diameter'),
                        text=~paste('<br>Topic:',topic),hoverinfo='text')

          p1 <- layout(p1,
                       showlegend=FALSE,
                       scene=list(
                         xaxis=list(title=sprintf('Axis 1 [%.02f%%]',eig[1]*100),
                                    showgrid=FALSE),
                         yaxis=list(title=sprintf('Axis 2 [%.02f%%]',eig[2]*100),
                                    showgrid=FALSE),
                         zaxis=list(title=sprintf('Axis 3 [%.02f%%]',eig[3]*100),
                                    showgrid=FALSE)),
                       paper_bgcolor='rgb(243, 243, 243)',
                       plot_bgcolor='rgb(243, 243, 243)')

        }

        p1

      })

      show_topic <- reactiveValues(k=0)

      observeEvent(event_data('plotly_click',source='ord_click'),{

        s <- event_data('plotly_click',source='ord_click')

        if (is.null(s)){

          show_topic$k <- 0
          updateNumericInput(session,'k_in',value=0)

        }else{

          t_idx <- s$pointNumber + 1
          updateNumericInput(session,'k_in',value=t_idx)
          show_topic$k <- t_idx

        }

      })

      observeEvent(input$k_in,{

        show_topic$k <- input$k_in

      })

      observeEvent(input$reset,{

        show_topic$k <- 0
        updateNumericInput(session,'k_in',value=0)

      })

      output$bar <- renderPlot({

        if (show_topic$k != 0){

          p_bar <- ggplot(data=REL()) +
            geom_bar(aes_(~Term,~Total,fill=~Taxon),stat='identity',color='white',alpha=.6) +
            geom_bar(aes_(~Term,~Freq),stat='identity',fill='darkred',color='white')

        }else{

          p_bar <- ggplot(data=REL()) +
            geom_bar(aes_(~Term,~Total,fill=~Taxon),stat='identity',color='white',alpha=1)

        }

        p_bar +
          coord_flip() +
          labs(x='',y='Frequency',fill='') +
          theme(axis.text.x=element_text(angle=-90,hjust=0,vjust=.5),
                legend.position='bottom') +
          viridis::scale_fill_viridis(discrete=TRUE,drop=FALSE) +
          guides(fill=guide_legend(nrow=2))

      })

      output$corr <- networkD3::renderForceNetwork({

        K <- nrow(corr$posadj)

        suppressWarnings({suppressMessages({
          g <- igraph::graph.adjacency(corr$posadj,mode='undirected',
                                       weighted=TRUE,diag=FALSE)

          wc <- igraph::cluster_walktrap(g)
          members <- igraph::membership(wc)

          g_d3 <- networkD3::igraph_to_networkD3(g,group=members)

          g_d3$links$edge_width <- 10*(.1+sapply(seq_len(nrow(g_d3$links)),function(r) corr$poscor[g_d3$links$source[r]+1,g_d3$links$target[r]+1]))
          g_d3$nodes$color <- 1
          g_d3$nodes$node_size <- 25*norm10(c(0,topic_marg))[-1]
          g_d3$nodes$name <- paste0('T',g_d3$nodes$name)

          networkD3::forceNetwork(Links=g_d3$links,Nodes=g_d3$nodes,
                                  Source='source',Target='target',
                                  charge=-25,
                                  opacity=.7,
                                  fontSize=12,
                                  zoom=TRUE,
                                  bounded=TRUE,
                                  NodeID='name',
                                  fontFamily='sans-serif',
                                  opacityNoHover=.7,
                                  Group='color',
                                  Value='edge_width',
                                  Nodesize='node_size',
                                  linkColour='#000000',
                                  linkWidth=networkD3::JS('function(d) {return d.value;}'),
                                  radiusCalculation=networkD3::JS('d.nodesize'),
                                  colourScale=networkD3::JS("color=d3.scaleLinear()\n.domain([1,1])\n.range(['red','red']);"))

        })})

      })

    }

      )

    }
