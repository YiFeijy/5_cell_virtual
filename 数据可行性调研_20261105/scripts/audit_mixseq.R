args <- commandArgs(trailingOnly=TRUE)
root <- if (length(args)) args[[1]] else "数据可行性调研_20261105"
x <- readRDS(file.path(root, "metadata", "mixseq_all_CL_features.rds"))
out <- file.path(root, "audit")
dir.create(out, showWarnings=FALSE)
tables <- list()
summary <- list()
for (n in setdiff(names(x), "metadata")) {
  d <- x[[n]]
  if (!is.data.frame(d)) next
  cols <- intersect(c("DEPMAP_ID", "CCLE_ID", "sens", "in_pool"), names(d))
  a <- d[, cols, drop=FALSE]
  a$comparison <- n
  tables[[n]] <- a
  summary[[n]] <- data.frame(comparison=n, rows=nrow(d),
                            in_pool=sum(d$in_pool %in% TRUE),
                            in_pool_with_sens=sum(d$in_pool %in% TRUE & !is.na(d$sens)))
}
write.csv(do.call(rbind, tables), file.path(out, "mixseq_feature_labels.csv"), row.names=FALSE)
write.csv(do.call(rbind, summary), file.path(out, "mixseq_label_coverage.csv"), row.names=FALSE)
print(do.call(rbind, summary))
