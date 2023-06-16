library(data.table)
library(dplyr)
library(ggplot2)
library(Matrix)
library(Seurat)

### compare library complexity between native 3' scRNA-seq and amplicon libraries

# load count matrices
seurat.data = as.matrix(Read10X('./scRNAseq/sMDS/'))
seurat.data.targeted = as.matrix(Read10X('./scRNAseq/U2AF1_targeted/'))
seurat.data.additionalPCR = as.matrix(Read10X('./scRNAseq/U2AF1_additionalPCR/'))
# aggregate data
seurat.data.counts = colSums(seurat.data[which(rownames(seurat.data) != 'U2AF1'),])
seurat.data.targeted.counts = colSums(seurat.data.targeted[which(rownames(seurat.data.targeted) != 'U2AF1'),])
seurat.data.additionalPCR.counts = colSums(seurat.data.additionalPCR[which(rownames(seurat.data.additionalPCR) != 'U2AF1'),])
seurat.data.genes = rowSums(seurat.data[which(rownames(seurat.data) != 'U2AF1'),])
seurat.data.targeted.genes = rowSums(seurat.data.targeted[which(rownames(seurat.data.targeted) != 'U2AF1'),])
seurat.data.additionalPCR.genes = rowSums(seurat.data.additionalPCR[which(rownames(seurat.data.additionalPCR) != 'U2AF1'),])

# plot comparison of count distribution
ggplot() +
  geom_histogram(data=as.data.frame(seurat.data.counts), aes(x=seurat.data.counts), fill='grey', binwidth = 50) +
  geom_histogram(data=as.data.frame(seurat.data.targeted.counts), aes(x=seurat.data.targeted.counts), fill='orange', binwidth = 50) +
  geom_histogram(data=as.data.frame(seurat.data.additionalPCR.counts), aes(x=seurat.data.additionalPCR.counts), fill='blue', binwidth = 50) +
  scale_x_continuous('counts per cell',limits = c(1,20000)) +
  scale_y_sqrt('cells') +
  theme_classic() +
  theme(axis.text = element_text('Arial', size=10, color='black'),
        axis.title = element_text('Arial', size=10, color='black'))
ggsave('./figures/QC/counts_per_cell.svg', width = 2.5, height = 2)

# plot comparison of counts per gene distribution
ggplot() +
  geom_histogram(data=as.data.frame(seurat.data.genes), aes(x=seurat.data.genes), fill='grey', binwidth = 50) +
  geom_histogram(data=as.data.frame(seurat.data.targeted.genes), aes(x=seurat.data.targeted.genes), fill='orange', binwidth = 50) +
  geom_histogram(data=as.data.frame(seurat.data.additionalPCR.genes), aes(x=seurat.data.additionalPCR.genes), fill='blue', binwidth = 50) +
  scale_x_continuous('counts per gene',limits = c(1,10000)) +
  scale_y_sqrt('genes', limits = c(0,5000)) +
  theme_classic() +
  theme(axis.text = element_text('Arial', size=10, color='black'),
        axis.title = element_text('Arial', size=10, color='black'))
ggsave('./figures/QC/counts_per_gene.svg', width = 2.5, height = 2)

### read out of PCR handles
#
# import of import_kite_counts function modified from
# https://github.com/caleblareau/asap_reproducibility/blob/master/pbmc_stim_multiome/code/11_setup.R
source ('./analysis/import_kite_counts.R')

# load pseudo-aligned counts of PCR handles

# 3' GEX:
# processed 481,211,466 reads, 9,569 reads pseudoaligned
# targeted:
# processed 4,579,212 reads, 3,975,119 reads pseudoaligned
# additionalPCR:
# processed 9,332,345 reads, 4,137,196 reads pseudoaligned

seq.counts.native = as.matrix(import_kite_counts2('./data/PCR_handles/sMDS/featurecounts/'))
seq.counts.native = seq.counts.native[,which(paste0(colnames(seq.counts.native),'-1') %in% colnames(seurat.data))]
seq.counts.native = as.data.frame(t(seq.counts.native)) %>% tidyr::pivot_longer(cols = c('Seq_1', 'Seq_2', 'Seq_3'))

seq.counts.targeted = as.matrix(import_kite_counts2('./data/PCR_handles/sMDS_targeted/featurecounts/'))
seq.counts.targeted = seq.counts.targeted[,which(paste0(colnames(seq.counts.targeted),'-1') %in% colnames(seurat.data))]
seq.counts.targeted = as.data.frame(t(seq.counts.targeted)) %>% tidyr::pivot_longer(cols = c('Seq_1', 'Seq_2', 'Seq_3'))

seq.counts.additionalPCR = as.matrix(import_kite_counts2('./data/PCR_handles/sMDS_additionalPCR/featurecounts/'))
seq.counts.additionalPCR = seq.counts.additionalPCR[,which(paste0(colnames(seq.counts.additionalPCR),'-1') %in% colnames(seurat.data))]
seq.counts.additionalPCR = as.data.frame(t(seq.counts.additionalPCR)) %>% tidyr::pivot_longer(cols = c('Seq_2', 'Seq_3'))

p=ggplot(seq.counts.native, aes(x=name, y=value)) +
  ggrastr::rasterize(geom_jitter(color='grey', size=0.5, alpha=0.5), dpi=600) +
  geom_violin(scale = 'width', fill=NA, draw_quantiles = c(0.25, .5, .75, 1), color='black') +
  scale_y_sqrt('pseudocounts') +
  theme_classic() +
  theme(axis.text = element_text('Arial', size=10, color='black'),
        axis.title = element_text('Arial', size=10, color='black'),
        axis.title.x = element_blank())
ggsave('./figures/QC/pseudo_counts_native_library_seq123.svg', width = 2, height = 2, plot = p)

p=ggplot(seq.counts.targeted, aes(x=name, y=value)) +
  ggrastr::rasterize(geom_jitter(color='orange', size=0.5, alpha=0.5), dpi=600) +
  geom_violin(scale = 'width', fill=NA, draw_quantiles = c(0.25, .5, .75, 1), color='black') +
  scale_y_sqrt('pseudocounts') +
  theme_classic() +
  theme(axis.text = element_text('Arial', size=10, color='black'),
        axis.title = element_text('Arial', size=10, color='black'),
        axis.title.x = element_blank())
ggsave('./figures/QC/pseudo_counts_targeted_library_seq123.svg', width = 2, height = 2, plot = p)

p=ggplot(seq.counts.additionalPCR, aes(x=name, y=value)) +
  ggrastr::rasterize(geom_jitter(color='blue', size=0.5, alpha=0.5), dpi=600) +
  geom_violin(scale = 'width', fill=NA, draw_quantiles = c(0.25, .5, .75, 1), color='black') +
  scale_y_sqrt('pseudocounts') +
  theme_classic() +
  theme(axis.text = element_text('Arial', size=10, color='black'),
        axis.title = element_text('Arial', size=10, color='black'),
        axis.title.x = element_blank())
ggsave('./figures/QC/pseudo_counts_additionalPCR_library_seq123.svg', width = 1.7, height = 2, plot=p)
