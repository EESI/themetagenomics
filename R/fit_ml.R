#' @rdname est.functions
#'
#' @export

est.ml <- function(object,iters=1000,verbose=FALSE,seed=sample.int(.Machine$integer.max,1),...){

  set.seed(check_seed(seed))

  gene_table <- object$gene_table

  if (verbose) cat('Fitting model via ML.\n')

  mm <- glmer.nb(count ~ (1|pw) + (1|topic) + (1|pw:topic),
                 data=gene_table,
                 verbose=verbose,
                 control=glmerControl(calc.derivs=TRUE,
                                      optCtrl=list(maxfun=iters)),
                 ...)


  if (verbose) cat('Extracting summary.\n')

  summary_pars <- c('mu','phi','b_pw_sigma','b_topic_sigma','b_pwxtopic_sigma','b_pw','b_topic','b_pwxtopic','yhat')
  extract_summary <- vector(mode='list',length=length(summary_pars))
  names(extract_summary) <- summary_pars

  extract_summary[['mu']] <- matrix(fixef(mm),ncol=1,dimnames=list('mu','mean'))
  extract_summary[['phi']] <- matrix(getME(mm,'glmer.nb.theta'),ncol=1,dimnames=list('phi','mean'))
  extract_summary[['b_pw_sigma']] <- matrix(c(VarCorr(mm)$`pw`),ncol=1,dimnames=list('b_pw_sigma','mean'))
  extract_summary[['b_topic_sigma']] <- matrix(c(VarCorr(mm)$`topic`),ncol=1,dimnames=list('b_topic_sigma','mean'))
  extract_summary[['b_pwxtopic_sigma']] <- matrix(c(VarCorr(mm)$`pw:topic`),ncol=1,dimnames=list('b_pwxtopic_sigma','mean'))
  extract_summary[['b_pw']] <- data.frame(pw=rownames(ranef(mm)$pw),
                                          mean=ranef(mm)$pw[,1],
                                          stringsAsFactors=TRUE)
  rownames(extract_summary[['b_pw']]) <- extract_summary[['b_pw']]$pw
  extract_summary[['b_topic']] <- data.frame(topic=1:nrow(ranef(mm)$topic),
                                             mean=ranef(mm)$topic[,1],
                                             stringsAsFactors=TRUE)
  rownames(extract_summary[['b_topic']]) <- extract_summary[['b_topic']]$topic
  extract_summary[['b_pwxtopic']] <- data.frame(pw=gsub('^(.*)\\:([0-9]+)$','\\1',rownames(ranef(mm)$`pw:topic`)),
                                                topic=gsub('^(.*)\\:([0-9]+)$','\\2',rownames(ranef(mm)$`pw:topic`)),
                                                mean=ranef(mm)$`pw:topic`[,1],
                                                stringsAsFactors=TRUE)
  rownames(extract_summary[['b_pwxtopic']]) <- rownames(ranef(mm)$`pw:topic`)
  extract_summary[['yhat']] <- matrix(predict(mm),ncol=1)
  dimnames(extract_summary[['yhat']]) <- list(sprintf('yhat[%s]',1:nrow(extract_summary[['yhat']])),'mean')

  out <- list(summary=extract_summary)

  out[['fit']] <- mm

  return(out)

}
