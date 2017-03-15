#' @import ggplot2 shiny
#' @importFrom plotly plotlyOutput renderPlotly ggplotly
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

vis_function_effects <- function(topics,topic_effects,function_effects,taxa,beta_min=1e-5,gene_min=10){

  shinyApp(


    ui <- fluidPage(

      titlePanel('Thematic Structure of Predicted Functional Effects'),

      plotOutput('est',height='100px'),
      plotlyOutput('tax'),
      plotlyOutput('fxn'),
      tableOutput('tbl')

    ),


    server = function(input,output){

      fit <- topics$fit
      K <- fit$settings$dim$K
      beta <- t(exp(fit$beta$logbeta[[1]]))
      rownames(beta) <- fit$vocab
      colnames(beta) <- paste0('T',1:K)
      beta[beta < beta_min] <- beta_min
      beta <- log(beta)

      dd_row <- as.dendrogram(hclust(dist(beta,method='euclidean'),method='ward.D2'))
      otu_ord <- data.frame(otu=rownames(beta)[order.dendrogram(dd_row)],
                            otu_rank=1:nrow(beta))

      df <- data.frame(topic=paste0('T',1:K),
                       est=do.call('rbind',topic_effects$est)[,1],
                       sig=ifelse(1:K %in% topic_effects$sig,'1','0'))[order(topic_effects$rank),]
      df$topic <- factor(df$topic,levels=df$topic,ordered=TRUE)
      p_est <- ggplot(df,aes(topic,est,color=sig)) +
        geom_hline(yintercept=0,linetype=3) +
        geom_point(size=4) +
        theme_void() +
        labs(x='',y='Estimate') +
        scale_color_brewer(type='qual',palette=6,direction=-1) +
        scale_fill_brewer(type='qual',palette=6,direction=-1) +
        theme(legend.position='none')

      df <- data.frame(probability=matrix(beta,ncol=1),
                       otu=rownames(beta),
                       otu_rank=otu_ord$otu_rank,
                       topic=rep(colnames(beta),each=length(fit$vocab)),
                       topic_rank=rep(topic_effects$rank,each=length(fit$vocab)))
      df$taxon <- pretty_taxa_names(taxa[df$otu,])
      df$taxon <- paste0(df$taxon,' (',df$otu,')')
      df <- df[order(-df$topic_rank,df$otu_rank),]
      df$topic <- factor(df$topic,levels=unique(df$topic),ordered=TRUE)
      df$otu <- factor(df$otu,levels=unique(df$otu),ordered=TRUE)
      df$taxon <- factor(df$taxon,levels=unique(df$taxon),ordered=TRUE)

      p_tax <- ggplot(df,aes(topic,taxon)) +
        geom_raster(aes(fill=probability)) +
        viridis::scale_fill_viridis(name='Probability Rank') +
        labs(x='',y='',fill='Probability Rank') +
        scale_y_discrete(labels=taxon) +
        theme(legend.position='none',
              axis.title.y=element_text(face='bold',size=20),
              axis.ticks.y=element_blank(),
              axis.text.y=element_blank(),
              axis.text.x=element_text(angle=-90,hjust=0,vjust=.5))

      int_mat <- matrix(0.0,
                        length(unique(function_effects$model$summary$b_pwxtopic$pw)),
                        length(unique(function_effects$model$summary$b_pwxtopic$topic)),
                        dimnames=list(unique(function_effects$model$summary$b_pwxtopic$pw),
                                      paste0('T',unique(function_effects$model$summary$b_pwxtopic$topic))))
      sig_mat <- matrix(0,
                        length(unique(function_effects$model$summary$b_pwxtopic$pw)),
                        length(unique(function_effects$model$summary$b_pwxtopic$topic)),
                        dimnames=list(unique(function_effects$model$summary$b_pwxtopic$pw),
                                      paste0('T',unique(function_effects$model$summary$b_pwxtopic$topic))))

      sig <- rowSums(sign(function_effects$model$summary$b_pwxtopic[,c('2.5%','97.5%')]))

      for (i in seq_len(nrow(function_effects$model$summary$b_pwxtopic))){
          pw <- as.character(function_effects$model$summary$b_pwxtopic$pw[i])
          k <- paste0('T',as.character(function_effects$model$summary$b_pwxtopic$topic[i]))
          int_mat[pw,k] <- function_effects$model$summary$b_pwxtopic$mean[i]
          sig_mat[pw,k] <- sig[i]
      }

      dd_row <- as.dendrogram(hclust(dist(int_mat,method='euclidean'),method='ward.D2'))
      pw_ord <- data.frame(pw=rownames(int_mat)[order.dendrogram(dd_row)],
                           pw_rank=1:nrow(int_mat))

      df1 <- data.frame(weight=matrix(int_mat,ncol=1),
                       pw=rownames(int_mat),
                       pw_rank=pw_ord$pw_rank,
                       topic=rep(colnames(int_mat),each=nrow(int_mat)),
                       topic_rank=rep(topic_effects$rank,each=nrow(int_mat)))
      df1 <- df1[order(-df1$topic_rank,df1$pw_rank),]
      df1$topic <- factor(df1$topic,levels=unique(df1$topic),ordered=TRUE)
      df1$pw <- factor(df1$pw,levels=unique(df1$pw),ordered=TRUE)

      df2 <- data.frame(sig=matrix(sig_mat,ncol=1),
                        weight=matrix(int_mat,ncol=1),
                        pw=rownames(int_mat),
                        pw_rank=pw_ord$pw_rank,
                        topic=rep(colnames(int_mat),each=nrow(int_mat)),
                        topic_rank=rep(topic_effects$rank,each=nrow(int_mat)))
      df2 <- df2[order(-df2$topic_rank,df2$pw_rank),]
      df2$topic <- factor(df2$topic,levels=unique(df2$topic),ordered=TRUE)
      df2$pw <- factor(df2$pw,levels=unique(df2$pw),ordered=TRUE)

      p_fxn <- ggplot(df1,aes(topic,pw,fill=weight)) +
        geom_raster() +
        geom_point(data=df2[df2$sig > 0,],
                   aes(topic,pw),color='red',shape=3) +
        geom_point(data=df2[df2$sig < 0,],
                   aes(topic,pw),color='cyan',shape=3) +
        viridis::scale_fill_viridis(name='Coefficient') +
        labs(x='',y='',fill='Coefficient') +
        theme(axis.text.x=element_text(angle=-90,hjust=0,vjust=.5),
              legend.position='left') +
        geom_vline(xintercept=3.5,color='white',size=2) +
        geom_vline(xintercept=22.5,color='white',size=2) +
        theme(legend.position='none')


      output$est <- renderPlot({
        p_est
      })

      output$tax <- renderPlotly({
        ggplotly(p_tax)
      })

      output$fxn <- renderPlotly({
        ggplotly(p_fxn,source='hm_fxn')
      })

      output$tbl <- renderTable({

        s <- event_data('plotly_click',source='hm_fxn')

        if (length(s)){

          k <- levels(df1$topic)[s[['x']]]
          pw <- levels(df1$pw)[s[['y']]]

          gene_tbl <- function_effects$gene_table[paste0('T',function_effects$gene_table$topic) == k & function_effects$gene_table$pw == pw,c('count','ko','description')]
          gene_tbl <- gene_tbl[gene_tbl$count >= 10,]
          gene_tbl <- gene_tbl[order(gene_tbl$count,decreasing=TRUE),]

          gene_tbl

        }else{

          NULL

        }

      },digits=0,width='900',align='l',rownames=FALSE,striped=TRUE,hover=TRUE,bordered=FALSE,spacing='s')




    }



  )

}
