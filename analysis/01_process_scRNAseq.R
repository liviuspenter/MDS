library(ComplexHeatmap)
library(dplyr)
library(ggplot2)
library(Seurat)

# functionality for extracting UMIs and generating kneeplots
# can be downloaded like this:
# devtools::install_github('liviuspenter/nanoranger.R')
library(nanoranger.R)

### input libraries:
# MDS: 3' scRNA-seq library at timepoint 31 months
# MDS2.3p: 3' scRNA-seq library at timepoint 33 months

### read in souporcell output from 3 scRNA-seq libraries
sMDS.genotype = as.data.frame(data.table::fread('./data/souporcell/clusters_MDS.tsv'))
rownames(sMDS.genotype) = paste0('sMDS_',sMDS.genotype$barcode)
sMDS.genotype$barcode = paste0('sMDS_',sMDS.genotype$barcode)
sMDS.genotype2.3p = as.data.frame(data.table::fread('./data/souporcell/clusters_MDS2.3p.tsv'))
rownames(sMDS.genotype2.3p) = paste0('sMDS2.3p_',sMDS.genotype2.3p$barcode)
sMDS.genotype2.3p$barcode = paste0('sMDS2.3p_',sMDS.genotype2.3p$barcode)
sMDS.genotype = as.data.frame(dplyr::bind_rows(sMDS.genotype, sMDS.genotype2.3p))

### load count matrices - need to be downloaded from NCBI GEO
# MDS sample
seurat.data = Read10X('./scRNAseq/sMDS/')
sMDS = CreateSeuratObject(seurat.data, project = 'sMDS', min.cells = 3, min.features = 200)
sMDS[['percent.mt']] <- PercentageFeatureSet(sMDS, pattern = '^MT-')
sMDS = RenameCells(sMDS, add.cell.id = 'sMDS')

# MDS progression sample
seurat.data = Read10X('./scRNAseq/sMDS2.3p/')
sMDS2.3p = CreateSeuratObject(seurat.data, project = 'sMDS2.3p', min.cells = 3, min.features = 200)
sMDS2.3p[['percent.mt']] <- PercentageFeatureSet(sMDS2.3p, pattern = '^MT-')
sMDS2.3p = RenameCells(sMDS2.3p, add.cell.id = 'sMDS2.3p')

sMDS = merge(sMDS, list(sMDS2.3p))
VlnPlot(sMDS, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
sMDS <- subset(sMDS, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 20)
sMDS = subset(sMDS, cells = sMDS.genotype$barcode[which(sMDS.genotype$status == 'singlet')])

sMDS = NormalizeData(sMDS)
sMDS = FindVariableFeatures(sMDS)
sMDS = ScaleData(sMDS)
sMDS = RunPCA(sMDS)
sMDS = RunUMAP(sMDS, dims = 1:30)
sMDS = FindNeighbors(sMDS)
sMDS = FindClusters(sMDS, resolution = 0.3)

# bone marrow reference from Seurat
# see https://satijalab.org/seurat/articles/multimodal_reference_mapping.html

bm = readRDS(file = './scRNAseq/reference/bm.reference.rds')
bm[["spca.annoy.neighbors"]] <- LoadAnnoyIndex(object = bm[["spca.annoy.neighbors"]], file = "./scRNAseq/reference/reftmp.idx")

anchors <- FindTransferAnchors(
  reference = bm, query = sMDS, k.filter = NA,
  reference.reduction = "spca", reference.neighbors = "spca.annoy.neighbors", dims = 1:50)

sMDS <- MapQuery(
  anchorset = anchors, query = sMDS,
  reference = bm,
  refdata = list(
    celltype = "celltype.l2",
    predicted_ADT = "ADT"),
  reference.reduction = "spca",
  reduction.model = "wnn.umap"
)
#saveRDS(file='./scRNAseq/sMDS_combined.rds', sMDS)
sMDS = readRDS('./scRNAseq/sMDS_combined.rds')

### plot donor versus recipient
# simple helper function for aggregating meta celltypes as some celltypes have too few events to be meaningful
source('./analysis/aggregated.celltype.R')

sMDS.genotype$predicted.celltype = sMDS$predicted.celltype[sMDS.genotype$barcode]
sMDS.genotype$sample = stringr::str_split_fixed(sMDS.genotype$barcode, pattern = '_', n = 2)[,1]
sMDS.genotype = sMDS.genotype[which(sMDS.genotype$assignment %in% c(0,1)),]
sMDS.genotype$chimerism = 'donor'
# annotation of donor and recipient is based on known T cell chimerism and distribution of souporcell clusters across T cells
sMDS.genotype$chimerism[which(sMDS.genotype$sample == 'sMDS' & sMDS.genotype$assignment == 0)] = 'recipient'
sMDS.genotype$chimerism[which(sMDS.genotype$sample == 'sMDS2.3p' & sMDS.genotype$assignment == 1)] = 'recipient'
sMDS.genotype$aggregated.celltype = sapply(sMDS.genotype$predicted.celltype, aggregated_celltype)

# calculate donor and recipient frequencies
df = as.data.frame(sMDS.genotype %>% filter(!is.na(predicted.celltype)) %>%
                     group_by(sample, aggregated.celltype, chimerism) %>%
                     summarize(count = n()) %>%
                     tidyr::pivot_wider(values_from = 'count', names_from = 'chimerism'))
df[is.na(df)] = 0

df = df %>% group_by(sample, aggregated.celltype) %>%
  summarize(recipient = sum(recipient),
            donor = sum(donor))
df$recipient.freq = df$recipient / (df$recipient + df$donor)
df$donor.freq = df$donor / (df$recipient + df$donor)
df = df %>% tidyr::pivot_longer(cols = c('donor.freq', 'recipient.freq'),
                                names_to = 'chimerism')
df$chimerism = factor(df$chimerism, levels = c('recipient.freq', 'donor.freq'))
df$aggregated.celltype = factor(df$aggregated.celltype, levels = c('HSC','GMP', 'Mono','megakaryopoiesis', 'erythropoiesis',
                                                                   'CD4+ T cell', 'CD8+ T cell', 'NK', 'innate T cell', 'progenitor B cell', 'B cell'))

# plot donor and recipient chimerism across aggregated celltypes for MDS sample timepoint 31 months
ggplot(df[which(df$sample == 'sMDS'),], aes(x=aggregated.celltype, y=100*value, fill=chimerism)) +
  geom_col() +
  scale_y_continuous('% cells') +
  scale_fill_manual(values = c('recipient.freq' = 'firebrick', 'donor.freq' = 'blue')) +
  theme_classic() +
  theme(legend.position = 'none',
        axis.text = element_text('Arial', size=10, color='black'),
        axis.title = element_text('Arial', size=10, color='black'),
        axis.text.x = element_text(angle=90,hjust=1,vjust=0.5),
        axis.title.x = element_blank())
ggsave('./figures/MDS_donor_recipient_freq.svg', width = 2, height = 2.5)

# plot donor and recipient chimerism across aggregated celltypes for MDS sample 2 3' timepoint 33 months
ggplot(df[which(df$sample == 'sMDS2.3p'),], aes(x=aggregated.celltype, y=100*value, fill=chimerism)) +
  geom_col() +
  scale_y_continuous('% cells') +
  scale_fill_manual(values = c('recipient.freq' = 'firebrick', 'donor.freq' = 'blue')) +
  theme_classic() +
  theme(legend.position = 'none',
        axis.text = element_text('Arial', size=10, color='black'),
        axis.title = element_text('Arial', size=10, color='black'),
        axis.text.x = element_text(angle=90,hjust=1,vjust=0.5),
        axis.title.x = element_blank())
ggsave('./figures/MDS2.3p_donor_recipient_freq.svg', width = 2, height = 2.5)

# plot chimerism kinetics for monocytes and megakaryocytes at timepoint 31 and 33 months using both 3' libraries
ggplot() +
  geom_point(data=df[which(df$aggregated.celltype == 'Mono' &
                             df$chimerism == 'recipient.freq' &
                             df$sample %in% c('sMDS', 'sMDS2.3p')),],
             aes(x=sample, y=100*value), size=0.5) +
  geom_line(data=df[which(df$aggregated.celltype == 'Mono' &
                            df$chimerism == 'recipient.freq' &
                            df$sample %in% c('sMDS', 'sMDS2.3p')),],
            aes(x=sample, y=100*value, group=1), color=AML.combined.colors[['CD14 Mono']]) +
  geom_point(data=df[which(df$aggregated.celltype == 'megakaryopoiesis' &
                             df$chimerism == 'recipient.freq' &
                             df$sample %in% c('sMDS', 'sMDS2.3p')),],
             aes(x=sample, y=100*value), size=0.5) +
  geom_line(data=df[which(df$aggregated.celltype == 'megakaryopoiesis' &
                            df$chimerism == 'recipient.freq' &
                            df$sample %in% c('sMDS', 'sMDS2.3p')),],
            aes(x=sample, y=100*value, group=1), color=AML.combined.colors[['Prog_Mk']]) +
  scale_x_discrete(labels = c('diagnosis', 'progression')) +
  scale_y_continuous('% recipient chimerism',limits = c(0,100)) +
  theme_classic() +
  theme(legend.position = 'none',
        axis.text = element_text('Arial', size=10, color='black'),
        axis.title = element_text('Arial', size=10, color='black'),
        axis.text.x = element_text(angle=90,hjust=1,vjust=0.5),
        axis.title.x = element_blank())
ggsave('./figures/chimerism_kinetics.svg', width = 1.1, height = 2.5)

### read out of U2AF1mut cells

# process U2AF1 genotyping information
U2AF1 = extract_mutation(BC.data.file = './data/pileup/MDS.pileup.U2AF1.combined.csv.gz', ALT = 'T', REF = 'G', FILTER = 50)
U2AF1 = as.data.frame(U2AF1)
U2AF1$bc = paste0('sMDS_', U2AF1$bc)
rownames(U2AF1) = U2AF1$bc
# save output for later as processing is time-consuming
# write.csv2(U2AF1, file = './data/pileup/U2AF1.combined.csv', quote = F)

U2AF1 = as.data.frame(read.csv2('./data/pileup/U2AF1.combined.csv', row.names = 1))
U2AF1 = U2AF1[intersect(colnames(sMDS), rownames(U2AF1)),]
U2AF1$chimerism = sMDS.genotype[U2AF1$bc, 'assignment']
U2AF1$predicted.celltype = sMDS$predicted.celltype[rownames(U2AF1)]
U2AF1 = U2AF1[order(U2AF1$chimerism),]

cells.mono = U2AF1$bc[which(U2AF1$predicted.celltype == 'CD14 Mono')]
cells.mk = U2AF1$bc[which(U2AF1$predicted.celltype == 'Prog_Mk')]

# plot sex-specific genes of monocytes and megakaryocytes
df = as.matrix(GetAssayData(sMDS)[ c('RPS4Y1', 'DDX3Y', 'EIF1AY','KDM5D', 'XIST'), c(cells.mono, cells.mk)])
df = df[,c(U2AF1$bc[which(U2AF1$predicted.celltype == 'CD14 Mono' & U2AF1$mutated == 'mutated')],
           U2AF1$bc[which(U2AF1$predicted.celltype == 'CD14 Mono' & U2AF1$mutated == 'wildtype')],
           U2AF1$bc[which(U2AF1$predicted.celltype == 'Prog_Mk' & U2AF1$mutated == 'mutated')],
           U2AF1$bc[which(U2AF1$predicted.celltype == 'Prog_Mk' & U2AF1$mutated == 'wildtype')])]

ha = columnAnnotation(chimerism = as.character(sMDS.genotype[colnames(df), 'assignment']),
                      U2AF1S34Y = as.character(U2AF1[colnames(df), 'mutated']),
                      col = list('celltype' = c('Mono' = as.character(AML.combined.colors['CD14 Mono']),
                                                'Prog_MK' = as.character(AML.combined.colors['Prog_Mk'])),
                                 'chimerism' = c('0' = 'firebrick', '1' = 'blue'),
                                 'U2AF1S34Y' = c('wildtype' = 'grey90', 'mutated' = 'darkgreen')),
                      border=T, simple_anno_size = unit(10, 'pt'),
                      annotation_name_side = 'left', annotation_name_gp = gpar(fontsize=10))

col_fun = circlize::colorRamp2(breaks = seq(0,3,3/8), colors = BuenColors::jdb_palette(name = 'solar_rojos'))
svglite::svglite('./figures/heatmap.svg', width = 4, height = 1.5)
Heatmap(df, top_annotation = ha, cluster_rows = F, cluster_columns = F, show_column_names = F, col = col_fun,
        border=T, column_split = c(rep('Mono', length(cells.mono)), rep('Prog_MK', length(cells.mk))),
        row_names_side = 'left', row_names_gp = gpar(fontsize=10))
dev.off()

### plot UMAPs
p=DimPlot(subset(sMDS, orig.ident %in% c('sMDS', 'sMDS2.3p')), reduction = 'ref.umap', group.by = 'predicted.celltype', cols = AML.combined.colors) +
  NoLegend() + NoAxes() + theme(plot.title = element_blank())
ggsave('./figures/UMAP/predicted_celltype.png', width = 3, height = 3, dpi = 600, plot = p)
q=DimPlot(subset(sMDS, cells = U2AF1$bc), cells.highlight = U2AF1$bc[which(U2AF1$mutated == 'mutated')], reduction = 'ref.umap', cols = 'grey90', cols.highlight = 'darkgreen') +
  NoLegend() + NoAxes() + theme(plot.title = element_blank())
ggsave('./figures/UMAP/U2AF1_mutation.png', width = 3, height = 3, dpi = 600, plot = q)
r=DimPlot(subset(sMDS, orig.ident %in% c('sMDS', 'sMDS2.3p')), cells.highlight = sMDS.genotype$barcode[which(sMDS.genotype$chimerism == 'recipient' & sMDS.genotype$barcode %in% colnames(sMDS))],
          cols = 'grey90', cols.highlight = 'firebrick', reduction = 'ref.umap', sizes.highlight = 0.5) + NoLegend() + NoAxes() + theme(plot.title = element_blank())
ggsave('./figures/UMAP/recipient.png', width = 3, height = 3, dpi = 600, plot = r)
s=DimPlot(subset(sMDS, orig.ident %in% c('sMDS', 'sMDS2.3p')), cells.highlight = sMDS.genotype$barcode[which(sMDS.genotype$chimerism == 'donor')],
          cols = 'grey90', cols.highlight = 'blue', reduction = 'ref.umap', sizes.highlight = 0.5) + NoLegend() + NoAxes() + theme(plot.title = element_blank())
ggsave('./figures/UMAP/donor.png', width = 3, height = 3, dpi = 600, plot = s)
t=FeaturePlot(sMDS, features = 'XIST', order = T,
              cols = c('grey90', BuenColors::jdb_palette(name = 'solar_rojos')), reduction = 'ref.umap') +
  NoLegend() + NoAxes() + theme(plot.title = element_blank())
ggsave('./figures/UMAP/XIST.png', width = 3, height = 3, dpi = 600, plot = t)
u=FeaturePlot(sMDS, features = 'RPS4Y1', order = T,
              cols = c('grey90', BuenColors::jdb_palette(name = 'solar_rojos')), reduction = 'ref.umap') +
  NoLegend() + NoAxes() + theme(plot.title = element_blank())
ggsave('./figures/UMAP/RPS4Y1.png', width = 3, height = 3, dpi = 600, plot = u)
