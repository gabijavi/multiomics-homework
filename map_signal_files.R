library(data.table)
library(GenomicRanges)
library(rtracklayer)

peak_files <- c(
  h3k27ac = "DATA/ENCFF646OKL.bed.gz",
  h3k27me3 = "DATA/ENCFF179UDS.bed.gz"
)
bw_files <- c(
  ENCFF075TZQ = "DATA/ENCFF075TZQ.bigWig",
  ENCFF937OKW = "DATA/ENCFF937OKW.bigWig"
)

set.seed(1)
for (pf in names(peak_files)) {
  p <- fread(peak_files[pf], nrows = 10000, showProgress = FALSE)
  gr <- GRanges(p$V1, IRanges(p$V2 + 1, p$V3))
  gr <- sample(gr, min(1000, length(gr)))
  cat("\nPeaks:", pf, "\n")
  for (bw in names(bw_files)) {
    vals <- import(bw_files[bw], which = gr, as = "NumericList")
    m <- mean(vapply(vals, function(x) mean(x, na.rm = TRUE), numeric(1)), na.rm = TRUE)
    cat(bw, sprintf("%.4f", m), "\n")
  }
}
