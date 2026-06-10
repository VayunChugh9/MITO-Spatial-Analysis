if (!exists("lisa")) {
  stop("Object `lisa` not found. Run the upstream localmoran calculation first.")
}

str(lisa)
colnames(as.data.frame(lisa))
