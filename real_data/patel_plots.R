library(cowplot)
library(RColorBrewer)
library(magrittr)
library(ggplot2)

load("patel_covariates.rda")

counts <- read.table("patel/glioblastoma_raw_rnaseq_SCandbulk_counts_withannots.txt", header=TRUE, stringsAsFactors = FALSE)
info <- as.matrix(counts)[1,-(1:3)]

# Drops "X" rownames column and sample info row
counts <- counts[-1,-1]

gene_symbols <- counts[,1]
ensembl_ids <- counts[,2]

sample_names <- colnames(counts)[-(1:2)]

# Drops gene symbol and ID columns
all.counts <- counts[,-(1:2)]
rownames(all.counts) <- ensembl_ids

metadata <- read.table("patel/SraRunTable.txt", sep=',', stringsAsFactors = FALSE, header=TRUE)
rownames(metadata)<-metadata$Run
metadata <- metadata[sample_names,]

batch <- stringr::str_split_fixed(info, "_", 2)[,2]
batch <- stringr::str_split_fixed(batch, "_", 2)[,1]

keep <- which(grepl("^Single cell", metadata$source_name) &
                grepl("MGH", metadata$patient_id) &
		metadata$subtype != "")
		#metadata$subtype != "" & metadata$subtype != "None")
		
metadata <- metadata[keep,]

batch <- as.factor(batch[keep])
all.counts <- all.counts[,keep]
all.counts <- as.matrix(all.counts)
class(all.counts) <- "numeric"

stopifnot(all(rownames(metadata)==colnames(all.counts)))

col1 <- brewer.pal(9, "Set1")
col2 <- c(brewer.pal(8, "Set2"), brewer.pal(8, "Set3"), brewer.pal(8, "Set1"))

detection_rate <- colSums(all.counts>0)
coverage <- colSums(all.counts)

qc <- cbind(detection_rate, coverage)

level1 <- as.factor(metadata$patient_id)
level2 <- as.factor(metadata$subtype)

library(scater)
sceset <- SingleCellExperiment(list(counts=all.counts))

# Only keeping features with non-zero counts in at least one cell
keep_feature <- rowSums(assay(sceset) > 0) > 0
sceset <- sceset[keep_feature,]

# Deprecated, but we only need revisit this if it doesn't work, but perCellQCMetrics doesn't give us everything we need here
sceset <- calculateQCMetrics(sceset)
# pct_counts_top_100_features
# pct_counts_top_200_features
# pct_counts_top_50_features
# pct_counts_top_500_features
# pct_dropout: have this, but it is in rowData and by feature
# total_counts
# total_features

# Maybe pData is just being used to access the QC metrics? 
# Missing dropout here, let's run it and see what happens
# Ok, everything seems fine. Let's try and calculate pct_dropout by hand, taking 0 as dropout
#qc <- pData(sceset)[,c(1, 4, 7, 11:14)]

pct_dropout=100.0*apply(counts(sceset) == 0, 2, sum)/nrow(counts(sceset))
qc <- data.frame(colData(sceset)[,c("total_features_by_counts", "pct_counts_in_top_50_features", "pct_counts_in_top_100_features", "pct_counts_in_top_200_features", "pct_counts_in_top_500_features", "total_counts")], pct_dropout=pct_dropout)

filter <- rowSums(all.counts>10)>=10
raw <- all.counts[filter,]

colMerged <- col1[level1]
colBatch <- col2[level2]

data.frame(Dim1=pc_tc[,1], Dim2=pc_tc[,2]) %>%
  ggplot(aes(Dim1, Dim2, colour=level1)) + geom_point() +
  scale_color_brewer(palette="Set1") + ggtitle("PCA")-> panel1_pca

data.frame(Dim1=zifa_tc[,1], Dim2=zifa_tc[,2]) %>%
  ggplot(aes(Dim1, Dim2, colour=level1)) + geom_point() +
  scale_color_brewer(palette="Set1")+ ggtitle("ZIFA") -> panel1_zifa

data.frame(Dim1=zinb@W[,1], Dim2=zinb@W[,2]) %>%
  ggplot(aes(Dim1, Dim2, colour=level1)) + geom_point()  +
  scale_color_brewer(palette="Set1")+ ggtitle("ZINB")  -> panel1_zinb

p1 <- plot_grid(panel1_pca + theme(legend.position = "none"),
          panel1_zifa + theme(legend.position = "none"),
          panel1_zinb + theme(legend.position = "none"),
          labels=c("a", "c", "e"), align = "h", ncol=3)

legend <- get_legend(panel1_pca)
upper <- plot_grid(p1, legend, rel_widths = c(3, .6))

data.frame(Dim1=pc_tc[,1], Dim2=pc_tc[,2]) %>%
  ggplot(aes(Dim1, Dim2)) + geom_point(aes(color=detection_rate)) + scale_colour_gradient(low="blue", high="yellow")+ggtitle("PCA") -> panel2_pca
data.frame(Dim1=zifa_tc[,1], Dim2=zifa_tc[,2]) %>%
  ggplot(aes(Dim1, Dim2)) + geom_point(aes(color=detection_rate)) + scale_colour_gradient(low="blue", high="yellow") +ggtitle("ZIFA")-> panel2_zifa
data.frame(Dim1=zinb@W[,1], Dim2=zinb@W[,2]) %>%
  ggplot(aes(Dim1, Dim2)) + geom_point(aes(color=detection_rate)) + scale_colour_gradient(low="blue", high="yellow") +ggtitle("ZINB")-> panel2_zinb

p2 <- plot_grid(panel2_pca + theme(legend.position = "none"),
                panel2_zifa + theme(legend.position = "none"),
                panel2_zinb + theme(legend.position = "none"),
                labels=c("b", "d", "f"), align = "h", ncol=3)

legend2 <- get_legend(panel2_pca)
lower <- plot_grid(p2, legend2, rel_widths = c(3, .6))

fig1 <- plot_grid(upper, lower, ncol=1, nrow=2)

save_plot("patel_plots/patel_fig1.pdf", fig1,
          ncol = 3,
          nrow = 3,
          base_aspect_ratio = 1.3
)

cors <- lapply(1:2, function(i) abs(cor(pc_tc[,i], as.matrix(qc), method="spearman")))
cors <- unlist(cors)
bars <- data.frame(AbsoluteCorrelation=cors,
                   QC=rep(stringr::str_to_lower(colnames(qc)), 2),
                   Dimension=as.factor(rep(1:2, each=ncol(qc))))

bars %>%
  ggplot(aes(Dimension, AbsoluteCorrelation, group=QC, fill=QC)) +
  geom_bar(stat="identity", position='dodge') +
  scale_fill_manual(values=col2) + ylim(0, .5)+ggtitle("PCA") -> panel2_pca

cors <- lapply(1:2, function(i) abs(cor(zifa_tc[,i], as.matrix(qc))))
cors <- unlist(cors)
bars <- data.frame(AbsoluteCorrelation=cors,
                   QC=rep(stringr::str_to_lower(colnames(qc)), 2),
                   Dimension=as.factor(rep(1:2, each=ncol(qc))))

bars %>%
  ggplot(aes(Dimension, AbsoluteCorrelation, group=QC, fill=QC)) +
  geom_bar(stat="identity", position='dodge') +
  scale_fill_manual(values=col2) + ylim(0, .5) +ggtitle("ZIFA")-> panel2_zifa

cors <- lapply(1:2, function(i) abs(cor(zinb@W[,i], as.matrix(qc))))
cors <- unlist(cors)
bars <- data.frame(AbsoluteCorrelation=cors,
                   QC=rep(stringr::str_to_lower(colnames(qc)), 2),
                   Dimension=as.factor(rep(1:2, each=ncol(qc))))

bars %>%
  ggplot(aes(Dimension, AbsoluteCorrelation, group=QC, fill=QC)) +
  geom_bar(stat="identity", position='dodge') +
  scale_fill_manual(values=col2) + ylim(0, .5)+ ggtitle("ZINB")-> panel2_zinb

p2 <- plot_grid(panel2_pca + theme(legend.position = "none"),
                panel2_zifa + theme(legend.position = "none"),
                panel2_zinb + theme(legend.position = "none"),
                labels=c("b", "d", "f"), align = "h", ncol=3)

legend2 <- get_legend(panel2_pca)
lower <- plot_grid(p2, legend2, rel_widths = c(3, 1))

fig1bis <- plot_grid(upper, lower, ncol=1, nrow=2)

save_plot("patel_plots/patel_fig1bis.pdf", fig1bis,
          ncol = 3,
          nrow = 3,
          base_aspect_ratio = 1.3
)

library(cluster)

methods <- list(pc_raw[,1:2], pc_tc[,1:2], pc_tmm[,1:2], pc_fq[,1:2],
                zifa_raw, zifa_tc, zifa_tmm, zifa_fq,
                zinb@W)
names(methods) <- c(paste0("PCA_", c("RAW", "TC", "TMM", "FQ")),
                    paste0("ZIFA_", c("RAW", "TC", "TMM", "FQ")),
                    "ZINB")
met_type <- as.factor(c(rep(c("PCA", "ZIFA"), each=4), "ZINB"))

sil_cl <- sapply(seq_along(methods), function(i) {
  d <- dist(methods[[i]])
  ss <- silhouette(as.numeric(level1), d)
  mean(ss[,3])
})

bars <- data.frame(AverageSilhouette=sil_cl, Method=names(methods), Type=met_type)

bars %>%
  ggplot(aes(Method, AverageSilhouette, group=Type, fill=Type)) +
  geom_bar(stat="identity", position='dodge') +
  scale_fill_manual(values=col1) + coord_flip() +
  theme(legend.position = "none") -> sil

sil2 <- plot_grid(sil, NULL, NULL, ncol=3, nrow=1, labels="G")
fig1_tris <- plot_grid(upper, lower, sil2, ncol=1, nrow=3)
save_plot("patel_plots/patel_fig1tris.pdf", fig1_tris,
          ncol = 3,
          nrow = 3,
          base_aspect_ratio = 1.3
)

methods_sub <- methods[c(2, 6, 9)]
sil_pc <- lapply(seq_along(methods_sub), function(i) {
  d <- dist(methods_sub[[i]])
  ss <- silhouette(as.numeric(level1), d)
  sss <-  summary(ss)
  sss$clus.avg.widths
})

bars <- data.frame(AverageSilhouette=unlist(sil_pc),
                   Method=rep(c("PCA", "ZIFA", "ZINB"), each=nlevels(level1)),
                   Cluster=rep(levels(level1), length(methods_sub)))

library(dplyr)
bars %>%
  dplyr::mutate(ClusterByMethod = paste0(Cluster, " ", Method)) %>%
  ggplot(aes(ClusterByMethod, AverageSilhouette, fill=Cluster)) +
  geom_bar(stat="identity", position='dodge') +
  scale_fill_manual(values=col1) + coord_flip() +
  theme(legend.position = "none", axis.text = element_text(size=8)) -> sil

p1 <- plot_grid(panel1_pca + theme(legend.position = "none"),
                panel1_zifa + theme(legend.position = "none"),
                panel1_zinb + theme(legend.position = "none"),
                labels=c("a", "c", "e"), align = "h", ncol=3)

upper <- plot_grid(p1, sil, labels=c("", "g"), rel_widths = c(3, 1))

fig1_4 <- plot_grid(upper, lower, ncol=1, nrow=2)
fig1_4

save_plot("patel_plots/patel_fig1_v4.pdf", fig1_4,
          ncol = 3,
          nrow = 3,
          base_aspect_ratio = 1.3
)

save_plot("patel_plots/patel_supp_sil.pdf", sil)

data.frame(Dim1=pc_raw[,1], Dim2=pc_raw[,2]) %>%
  ggplot(aes(Dim1, Dim2, colour=level1)) + geom_point() +
  scale_color_brewer(palette="Set1") -> pca_raw
data.frame(Dim1=pc_tc[,1], Dim2=pc_tc[,2]) %>%
  ggplot(aes(Dim1, Dim2, colour=level1)) + geom_point() +
  scale_color_brewer(palette="Set1") -> pca_tc
data.frame(Dim1=pc_tmm[,1], Dim2=pc_tmm[,2]) %>%
  ggplot(aes(Dim1, Dim2, colour=level1)) + geom_point() +
  scale_color_brewer(palette="Set1") -> pca_tmm
data.frame(Dim1=pc_fq[,1], Dim2=pc_fq[,2]) %>%
  ggplot(aes(Dim1, Dim2, colour=level1)) + geom_point() +
  scale_color_brewer(palette="Set1") -> pca_fq

fig_pca <- plot_grid(pca_raw, pca_tc, pca_tmm, pca_fq, labels=c("a", "b", "c", "d"))
save_plot("patel_plots/patel_supp_pca.pdf", fig_pca,
          ncol = 2,
          nrow = 2,
          base_aspect_ratio = 1.3)

data.frame(Dim1=zinb@W[,1], Dim2=zinb@W[,2]) %>%
  ggplot(aes(Dim1, Dim2, colour=batch)) + geom_point()  +
  scale_color_brewer(palette="Set1")+ ggtitle("ZINB")  -> zinb_plot

data.frame(Dim1=zinb_det@W[,1], Dim2=zinb_det@W[,2]) %>%
  ggplot(aes(Dim1, Dim2, colour=batch)) + geom_point()  +
  scale_color_brewer(palette="Set1")+ ggtitle("ZINB with det. rate")  -> zinb_det_plot

detPlots<-plot_grid(zinb_plot+theme(legend.position = "none"), 
		zinb_det_plot+theme(legend.position = "none"), 
		ncol=2)

save_plot("patel_plots/patel_det_cov.pdf", detPlots)
