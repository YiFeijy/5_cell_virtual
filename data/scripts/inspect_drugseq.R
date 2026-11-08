args <- commandArgs(trailingOnly=TRUE)
root <- if (length(args)) args[[1]] else 'data'
out <- file.path(root, 'analysis_20261106', 'DRUGseq_R_objects.txt')
con <- file(out, 'wt')
inspect <- function(x, label, depth=0L) {
  cat(label, 'class=', paste(class(x), collapse=','), 'dim=', paste(dim(x), collapse='x'),
      'length=', length(x), 'bytes=', as.numeric(object.size(x)), '\n', file=con)
  if (isS4(x)) {
    cat('slots=', paste(slotNames(x), collapse=','), '\n', file=con)
  } else if (is.data.frame(x)) {
    cat('fields=', paste(names(x), collapse=' | '), '\n', file=con)
    cat('row_names=', paste(head(rownames(x), 3), collapse=' | '), '\n', file=con)
  } else if (is.matrix(x)) {
    cat('row_names=', paste(head(rownames(x), 3), collapse=' | '), '\n', file=con)
    cat('column_names=', paste(head(colnames(x), 3), collapse=' | '), '\n', file=con)
    cat('top_left=', paste(x[seq_len(min(2,nrow(x))), seq_len(min(2,ncol(x)))], collapse=' | '), '\n', file=con)
  } else if (is.list(x) && depth < 4L) {
    for (n in head(names(x), 2)) inspect(x[[n]], paste0(label, '$', n), depth+1L)
  }
}
for (name in c('comparisons_metadata.RData', 'Exp_gzip.RData')) {
  cat('FILE', name, '\n', file=con); flush(con)
  env <- new.env()
  loaded <- load(file.path(root, 'raw', 'Novartis_DRUGseq', name), envir=env)
  cat('objects=', paste(loaded, collapse=','), '\n', file=con)
  for (n in loaded) inspect(env[[n]], n)
  if (name == 'comparisons_metadata.RData') {
    d <- env$comparison_group_member_info_public
    cat('distinct_cmpd_sample_ids=', length(unique(d$cmpd_sample_id)), '\n', file=con)
    cat('time_h_values=', paste(sort(unique(d$hours_post_treatment)), collapse=','), '\n', file=con)
    cat('concentration_units=', paste(unique(d$unit), collapse=','), '\n', file=con)
    write.csv(unique(d[, intersect(c('external_biosample_id','plate_barcode','well_id','cmpd_sample_id','concentration','unit','hours_post_treatment','batch_id'), names(d))]), file.path(root, 'analysis_20261106', 'DRUGseq_sample_design.csv'), row.names=FALSE)
  } else {
    inventory <- list()
    walk <- function(x, path) {
      if (is.matrix(x) || is.data.frame(x)) {
        inventory[[length(inventory)+1L]] <<- data.frame(path=path, class=paste(class(x),collapse=','), rows=nrow(x), columns=ncol(x))
      } else if (is.list(x)) {
        for (n in names(x)) walk(x[[n]], paste0(path,'$',n))
      }
    }
    walk(env$Exp, 'Exp')
    write.csv(do.call(rbind, inventory), file.path(root, 'analysis_20261106', 'DRUGseq_expression_object_dimensions.csv'), row.names=FALSE)
  }
  flush(con)
  rm(env); gc()
}
close(con)
