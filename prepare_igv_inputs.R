library(data.table)

dir.create("results/igv", recursive = TRUE, showWarnings = FALSE)

dmr <- fread("DATA/dmr.bed", showProgress = FALSE)
setnames(dmr, c("chr", "start", "end", "L", "value"))
dmr[, name := fifelse(value > 0, "hyper_DMR", "hypo_DMR")]
dmr[, score := 0]
dmr[, strand := "."]
dmr[, thickStart := start]
dmr[, thickEnd := end]
dmr[, itemRgb := fifelse(value > 0, "215,48,39", "69,117,180")]

dmr_bed9 <- dmr[, .(chr, start, end, name, score, strand, thickStart, thickEnd, itemRgb)]
fwrite(dmr_bed9, "results/igv/dmr_colored_hyper_red_hypo_blue.bed",
       sep = "\t", col.names = FALSE)

loci <- fread("results/igv/igv_loci.tsv", header = TRUE)

abs_path <- function(path) normalizePath(path, winslash = "/", mustWork = TRUE)
snapshot_dir <- abs_path("results/igv")

batch <- c(
  "new",
  "genome hg19",
  "maxPanelHeight 2000",
  paste("snapshotDirectory", snapshot_dir),
  paste("load", abs_path("DATA/ENCFF937OKW.bigWig")),
  paste("load", abs_path("DATA/ENCFF075TZQ.bigWig")),
  paste("load", abs_path("results/igv/dmr_colored_hyper_red_hypo_blue.bed")),
  "squish"
)

for (i in seq_len(nrow(loci))) {
  batch <- c(
    batch,
    paste("goto", loci$igv_region[i]),
    paste("snapshot", loci$snapshot_file[i])
  )
}
batch <- c(batch, "exit")

writeLines(batch, "results/igv/hw3_igv_batch.txt")

session_xml <- sprintf(
'<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<Session genome="hg19" hasGeneTrack="true" hasSequenceTrack="true" locus="%s" version="8">
  <Resources>
    <Resource path="%s"/>
    <Resource path="%s"/>
    <Resource path="%s"/>
  </Resources>
  <Panel name="DataPanel">
    <Track attributeKey="ENCFF937OKW.bigWig" id="%s" name="H3K27me3 normal signal"/>
    <Track attributeKey="ENCFF075TZQ.bigWig" id="%s" name="H3K27ac normal signal"/>
    <Track attributeKey="dmr_colored_hyper_red_hypo_blue.bed" id="%s" name="DMRs red hyper blue hypo"/>
  </Panel>
</Session>',
  loci$igv_region[1],
  abs_path("DATA/ENCFF937OKW.bigWig"),
  abs_path("DATA/ENCFF075TZQ.bigWig"),
  abs_path("results/igv/dmr_colored_hyper_red_hypo_blue.bed"),
  abs_path("DATA/ENCFF937OKW.bigWig"),
  abs_path("DATA/ENCFF075TZQ.bigWig"),
  abs_path("results/igv/dmr_colored_hyper_red_hypo_blue.bed")
)
writeLines(session_xml, "results/igv/hw3_igv_session.xml")

cat("Wrote results/igv/dmr_colored_hyper_red_hypo_blue.bed\n")
cat("Wrote results/igv/hw3_igv_batch.txt\n")
cat("Wrote results/igv/hw3_igv_session.xml\n")
