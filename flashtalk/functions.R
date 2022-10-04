plot_gsea <- function(df, es_col="E", pval_col="adj.P.Val") {

  df <- df %>% arrange(-.data[[pval_col]])
  df$log10_pval <- -log10(df[[pval_col]])

	df$Title <- .cleanup_ids(df$Title)
	df$Title <- factor(df$Title, levels=df$Title)
  
  ggplot(df, 
         aes_string(x="Title", y="log10_pval", fill=es_col)) + geom_bar(stat="identity") +
         xlab("") +
         ylab("-log10(FDR)") +
         coord_flip()


}

.cleanup_ids <- function(ids) {
  ids <- gsub("_", " ", ids)

  #return(ids)

  ids <- strsplit(ids, " ")
  min.l <- min(sapply(ids, length))

  max_prefix <- 0
  for(i in 1:min.l) {

    if(length(unique(sapply(ids, function(x) x[i]))) == 1L) {
      max_prefix <- i
    }
  }

  if(max_prefix > 0) {
    ids <- lapply(ids, function(x) x[ -1:-max_prefix ])
  }

  ids <- unlist(lapply(ids, paste, collapse=" "))

  if(!any(grepl("[a-z]", ids))) {
    ids <- tolower(ids)
    substr(ids, 1, 1) <- toupper(substr(ids, 1, 1))
  }

  ids
}


ggplot_gene <- function(x, group, 
                        groups=list(1:2, 3:4),
                        annotation=NULL, ...) {

  df <- data.frame(Expression=x, Group=group)

  p1 <- ggplot(df, aes(x=Group, y=Expression)) + 
      geom_jitter(width = .2, size=1, alpha=.5) +
      geom_boxplot(fill=NA, outlier.shape = NA) 
      
  if(!is.null(annotation)) {
    p1 <- add_pval(p1, groups, annotation = annotation, ...)
  }

  p1

}

ggplot_disco <- function(df,
                         lfc1="log2FoldChange.g1",
                         lfc2="log2FoldChange.g2",
                         pval1="padj.g1",
                         pval2="padj.g2") {

  df$disco <- 1





}


## return the number of "specific" enrichments in the given file
gsea_specific_num <- function(file_name, gsea_auc_thr=0.5, gsea_pval_thr=0.05) {

  res <- readRDS(file_name)

  map_int(res, ~ {
    n1 <- .x %>% filter(adj.P.Val.g1 < gsea_pval_thr) %>% nrow()
    n2 <- .x %>% filter(adj.P.Val.g2 < gsea_pval_thr) %>% nrow()
    n1 + n2
    })

}

## returns a numeric vector. Each number represents the number of
## "specific" genes (in both groups) for one replicate.
get_specific_num <- function(merged_res, lfc_thr=0, padj_thr=0.05) {

  n_spec <- map_int(merged_res, ~ {
    res <- .x %>%
      mutate(DEG.g1 = !is.na(padj.g1) & abs(log2FoldChange.g1) > lfc_thr & padj.g1 < padj_thr) %>%
      mutate(DEG.g2 = !is.na(padj.g2) & abs(log2FoldChange.g2) > lfc_thr & padj.g2 < padj_thr) 
    sum((res$DEG.g1 & !res$DEG.g2) | (!res$DEG.g1 & res$DEG.g2))
  })

  n_spec
}

## returns a numeric vector. Each number represents the number of
## DEGs (in both groups) for one replicate.
get_total_num <- function(merged_res, lfc_thr=0, padj_thr=0.05) {

  n_spec <- map_int(merged_res, ~ {
    res <- .x %>%
      mutate(DEG.g1 = !is.na(padj.g1) & abs(log2FoldChange.g1) > lfc_thr & padj.g1 < padj_thr) %>%
      mutate(DEG.g2 = !is.na(padj.g2) & abs(log2FoldChange.g2) > lfc_thr & padj.g2 < padj_thr) 
    sum((res$DEG.g1 | res$DEG.g2))
  })

  n_spec
}


