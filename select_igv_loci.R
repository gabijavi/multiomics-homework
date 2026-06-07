library(data.table)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(AnnotationDbi)

dir.create("results/igv", recursive = TRUE, showWarnings = FALSE)

dmr <- fread("DATA/dmr.bed", showProgress = FALSE)
setnames(dmr, c("chr", "start", "end", "L", "value"))
deg <- fread("DATA/DEG_All_Genes.csv", showProgress = FALSE)
deg[, ENTREZID := as.character(ENTREZID)]

ac_dt <- fread("DATA/ENCFF646OKL.bed.gz", showProgress = FALSE)
me3_dt <- fread("DATA/ENCFF179UDS.bed.gz", showProgress = FALSE)

dmr_gr <- GRanges(dmr$chr, IRanges(dmr$start + 1, dmr$end), value = dmr$value)
ac_peaks <- GRanges(ac_dt$V1, IRanges(ac_dt$V2 + 1, ac_dt$V3))
me3_peaks <- GRanges(me3_dt$V1, IRanges(me3_dt$V2 + 1, me3_dt$V3))

gene_gr <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
gene_gr <- gene_gr[names(gene_gr) %in% deg$ENTREZID]

gene_annot <- data.table(ENTREZID = names(gene_gr))
gene_annot[, SYMBOL := mapIds(org.Hs.eg.db, ENTREZID, "SYMBOL", "ENTREZID", multiVals = "first")]
deg_anno <- merge(deg, gene_annot, by = "ENTREZID", all.x = TRUE)

sig_down <- deg_anno[!is.na(padj) & padj < 0.05 & log2FoldChange < 0]
down_gr <- gene_gr[names(gene_gr) %in% sig_down$ENTREZID]
hyper_dmr <- dmr_gr[mcols(dmr_gr)$value > 0]

hyper_ac <- hyper_dmr[unique(queryHits(findOverlaps(hyper_dmr, ac_peaks)))]
nearest_down <- nearest(hyper_ac, down_gr)
locus1_dt <- data.table(
  dmr_i = seq_along(hyper_ac),
  gene_id = names(down_gr)[nearest_down],
  distance = distance(hyper_ac, down_gr[nearest_down]),
  dmr_value = mcols(hyper_ac)$value
)
locus1_dt <- merge(locus1_dt, sig_down[, .(gene_id = ENTREZID, SYMBOL, log2FoldChange, padj)], by = "gene_id")
locus1_dt <- locus1_dt[distance < 50000][order(distance, padj)]
locus1 <- hyper_ac[locus1_dt$dmr_i[1]]

hyper_me3 <- hyper_dmr[unique(queryHits(findOverlaps(hyper_dmr, me3_peaks)))]
nearest_any <- nearest(hyper_me3, gene_gr)
locus2_dt <- data.table(
  dmr_i = seq_along(hyper_me3),
  gene_id = names(gene_gr)[nearest_any],
  distance = distance(hyper_me3, gene_gr[nearest_any]),
  dmr_value = mcols(hyper_me3)$value
)
locus2_dt <- merge(locus2_dt, deg_anno[, .(gene_id = ENTREZID, SYMBOL, log2FoldChange, padj)], by = "gene_id")
locus2_dt <- locus2_dt[distance < 50000][order(distance, -dmr_value)]
locus2 <- hyper_me3[locus2_dt$dmr_i[1]]

make_region <- function(gr) {
  region <- resize(gr, width = max(width(gr) + 40000, 60000), fix = "center")
  paste0(as.character(seqnames(region)), ":", start(region), "-", end(region))
}

loci <- data.table(
  locus = c("active_h3k27ac_hyper_downregulated", "repressed_h3k27me3_hyper"),
  igv_region = c(make_region(locus1), make_region(locus2)),
  snapshot_file = c("igv_locus1_active_hyper_down.png", "igv_locus2_repressed_hyper.png"),
  nearest_gene = c(locus1_dt$SYMBOL[1], locus2_dt$SYMBOL[1]),
  entrez_id = c(locus1_dt$gene_id[1], locus2_dt$gene_id[1]),
  log2FoldChange = c(locus1_dt$log2FoldChange[1], locus2_dt$log2FoldChange[1]),
  padj = c(locus1_dt$padj[1], locus2_dt$padj[1]),
  dmr_value = c(mcols(locus1)$value, mcols(locus2)$value)
)

fwrite(loci, "results/igv/igv_loci.tsv", sep = "\t")
print(loci)
