library(dplyr)
library(ggplot2)

sMDS = readRDS(file = './scRNAseq/sMDS_combined.rds')

# load pileup data 
U2AF1.SN = as.data.frame(data.table::fread('./data/pileup/MDS.pileup.U2AF1.targeted.csv.gz'))
U2AF1.SN$library = 'SN'
U2AF1.additionalPCR = as.data.frame(data.table::fread('./data/pileup/MDS.pileup.U2AF1.additionalPCR.csv.gz'))
U2AF1.additionalPCR$library = 'additionalPCR'
U2AF1.original = as.data.frame(data.table::fread('./data/pileup/MDS.pileup.U2AF1.original.csv.gz'))
U2AF1.original$library = 'original'
U2AF1.combined = rbind(U2AF1.SN, U2AF1.additionalPCR)
U2AF1.combined = rbind(U2AF1.combined, U2AF1.original)

# plot knee plot of U2AF1 locus with native library or amplicon libraries
boo = U2AF1.combined %>% group_by(library, bc) %>% summarize(reads = length(umi))
boo = boo %>% group_by(library) %>% mutate(rank = rank(-reads, ties.method = 'random'))

p = ggplot(boo, aes(x=rank, y=reads, color=library)) + geom_line() +
  scale_x_log10('cell rank', limits = c(1,50000)) +
  scale_y_log10('reads', limits = c(1,50000)) +
  scale_color_manual(values = c('SN' = 'orange', 'additionalPCR' = 'blue', 'original' = 'grey')) +
  geom_hline(yintercept = 50, color='red') +
  theme_classic() +
  theme(legend.position = 'none',
        axis.title = element_text('Arial', size=10, color='black'),
        axis.text = element_text('Arial', size=10, color='black'),
        axis.ticks = element_line(color='black'))
ggsave('./figures/comparison/knee_plot_native_targeted_additionalPCR.svg', width = 2, height = 2)

# compare reads obtained per cell barcode using targeted or additionalPCR amplicon library
U2AF1.reads = merge(boo[which(boo$library == 'SN',), c('bc', 'reads')],
                    boo[which(boo$library == 'additionalPCR',), c('bc', 'reads')], all=T, by = 'bc')
colnames(U2AF1.reads) = c('bc', 'SN', 'additionalPCR')
U2AF1.reads$SN.plot = U2AF1.reads$SN
U2AF1.reads$additionalPCR.plot = U2AF1.reads$additionalPCR
U2AF1.reads$SN.plot[which(is.na(U2AF1.reads$SN.plot))] = jitter(rep(0.1, length(which(is.na(U2AF1.reads$SN.plot)))), factor = 10)
U2AF1.reads$additionalPCR.plot[which(is.na(U2AF1.reads$additionalPCR.plot))] = jitter(rep(0.1, length(which(is.na(U2AF1.reads$additionalPCR.plot)))), factor = 10)

p=ggplot(U2AF1.reads, aes(x=SN.plot, y=additionalPCR.plot)) +
  geom_hline(yintercept = 50, color='red') +
  geom_vline(xintercept = 50, color='red') +
  ggrastr::rasterize(ggpointdensity::geom_pointdensity(size=0.5), dpi = 600) +
  scale_x_log10('targeted PCR', breaks = c(0.1, 1, 10, 100, 1000, 10000), labels = c('ND', '1', '10', '100', '1000', '10000'), limits = c(0.05, 20000)) +
  scale_y_log10('additionalPCR PCR', breaks = c(0.1, 1, 10, 100, 1000, 10000), labels = c('ND', '1', '10', '100', '1000', '10000'), limits = c(0.05, 20000)) +
  scale_color_gradientn('density',colours = BuenColors::jdb_palette(name = 'solar_basic')) +
  theme_classic() +
  theme(legend.position = 'none', # change to bottom to get legend
        axis.title = element_text('Arial', size=8, color='black'),
        axis.text = element_text('Arial', size=8, color='black'),
        axis.ticks = element_line(color='black'))
ggsave('./figures/comparison/targeted_additionalPCR_read_comparison.svg', width = 2, height = 2, plot = p)

### compare genotyped cells using targeted and additionalPCR library

source('./analysis/aggregated.celltype.R')

# load data from targeted library
U2AF1.SN = as.data.frame(read.csv2('./data/pileup/U2AF1.targeted.csv', row.names = 1))
U2AF1.SN$bc = paste0('sMDS_', U2AF1.SN$bc)
rownames(U2AF1.SN) = U2AF1.SN$bc
U2AF1.SN = U2AF1.SN[intersect(U2AF1.SN$bc, colnames(sMDS)),]
U2AF1.SN$predicted.celltype = sMDS$predicted.celltype[U2AF1.SN$bc]
U2AF1.SN$aggregated.celltype = sapply(U2AF1.SN$predicted.celltype, aggregated_celltype)
U2AF1.SN.stats = U2AF1.SN %>% group_by(aggregated.celltype) %>% summarize(mutated.cells = length(which(mutated == 'mutated')),
                                                                         wildtype.cells = length(which(mutated == 'wildtype')),
                                                                         total.cells = length(predicted.celltype))
U2AF1.SN.stats$library = 'SN'

# load data from additionalPCR library
U2AF1.additionalPCR = as.data.frame(read.csv2('./data/pileup/U2AF1.additionalPCR.csv', row.names = 1))
U2AF1.additionalPCR$bc = paste0('sMDS_', U2AF1.additionalPCR$bc)
rownames(U2AF1.additionalPCR) = U2AF1.additionalPCR$bc
U2AF1.additionalPCR = U2AF1.additionalPCR[intersect(U2AF1.additionalPCR$bc, colnames(sMDS)),]
U2AF1.additionalPCR$predicted.celltype = sMDS$predicted.celltype[U2AF1.additionalPCR$bc]
U2AF1.additionalPCR$aggregated.celltype = sapply(U2AF1.additionalPCR$predicted.celltype, aggregated_celltype)
U2AF1.additionalPCR.stats = U2AF1.additionalPCR %>% group_by(aggregated.celltype) %>% summarize(mutated.cells = length(which(mutated == 'mutated')),
                                                                                 wildtype.cells = length(which(mutated == 'wildtype')),
                                                                                 total.cells = length(predicted.celltype))
U2AF1.additionalPCR.stats$library = 'additionalPCR'
df = rbind(U2AF1.SN.stats, U2AF1.additionalPCR.stats)
df$aggregated.celltype = factor(df$aggregated.celltype, levels = c('GMP', 'Mono','megakaryopoiesis', 'erythropoiesis',
                                                                   'CD4+ T cell', 'CD8+ T cell', 'NK', 'innate T cell', 'progenitor B cell', 'B cell'))

ggplot(data=df, aes(x=aggregated.celltype, y=total.cells, fill=library)) +
  geom_col(position = 'dodge') +
  scale_fill_manual(values = c('SN' = 'orange', 'additionalPCR' = 'blue')) +
  scale_y_continuous('genotyped cells') +
  theme_classic() +
  theme(legend.position = 'none',
        axis.text = element_text('Arial', size=10, color='black'),
        axis.title = element_text('Arial', size=10, color='black'),
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5),
        axis.title.x = element_blank())
ggsave('./figures/comparison/targeted_additionalPCR_genotyped_cells_comparison.svg', width = 2.5, height = 2.5)

### compare UMIs per cell between targeted and additionalPCR library

U2AF1.targeted = as.data.frame(read.csv2('./data/pileup/U2AF1.targeted.csv', row.names = 1))
U2AF1.targeted$count = U2AF1.targeted$alt + U2AF1.targeted$ref
U2AF1.additionalPCR = as.data.frame(read.csv2('./data/pileup/U2AF1.additionalPCR.csv', row.names = 1))
U2AF1.additionalPCR$count = U2AF1.additionalPCR$alt + U2AF1.additionalPCR$ref

ggplot() +
  geom_histogram(data=U2AF1.targeted, aes(x=count), binwidth=1, fill='orange') +
  geom_histogram(data=U2AF1.additionalPCR, aes(x=count), binwidth=1, fill='blue') +
  scale_x_continuous('UMIs per cell') +
  scale_y_continuous('cells') +
  theme_classic() +
  theme(axis.text = element_text('Arial', size=10, color='black'),
        axis.title = element_text('Arial', size=10, color='black'))
ggsave('./figures/comparison/UMIs_per_cell.svg', width = 2, height = 2)
