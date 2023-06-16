# generate coverage plots of chrM, U2AF1 and GAPDH

library(dplyr)
library(ggplot2)

# coverage generated with bedtools coverage -d -a targets.bed -b possorted_genome_bam.bam
MDS.coverage = data.table::fread('./data/coverage/MDS_coverage.bed')
MDS.additionalPCR.coverage = data.table::fread('./data/coverage/MDS_additionalPCR_coverage.bed')
MDS.targeted.coverage = data.table::fread('./data/coverage/MDS_targeted_coverage.bed')

# calculate position on mRNA
MDS.coverage = MDS.coverage %>% mutate(abs.pos = V2+V4-1) %>% group_by(V1) %>% mutate(pos = rank(abs.pos))
MDS.additionalPCR.coverage = MDS.additionalPCR.coverage %>% mutate(abs.pos = V2+V4-1) %>% group_by(V1) %>% mutate(pos = rank(abs.pos))
MDS.targeted.coverage = MDS.targeted.coverage %>% mutate(abs.pos = V2+V4-1) %>% group_by(V1) %>% mutate(pos = rank(abs.pos))

p=ggplot() + 
  ggrastr::rasterize(geom_col(data=MDS.coverage[which(MDS.coverage$V1 == 'chrM'),], aes(x=V4, y=V5/1000), fill='grey30'), dpi=600) +
  scale_x_continuous('chrM') + 
  scale_y_continuous('coverage') +
  theme_classic() + 
  theme(axis.text = element_text('Arial', size=8, color='black'),
        axis.title = element_text('Arial', size=8, color='black'))
q=ggplot() + 
  ggrastr::rasterize(geom_col(data=MDS.targeted.coverage[which(MDS.targeted.coverage$V1 == 'chrM'),], aes(x=V4, y=V5/1000), fill='orange'), dpi=600) +
  scale_x_continuous('chrM') + 
  scale_y_continuous('coverage') +
  theme_classic() + 
  theme(axis.text = element_text('Arial', size=8, color='black'),
        axis.title = element_text('Arial', size=8, color='black'),
        axis.title.y = element_blank())
r=ggplot() + 
  ggrastr::rasterize(geom_col(data=MDS.additionalPCR.coverage[which(MDS.additionalPCR.coverage$V1 == 'chrM'),], aes(x=V4, y=V5/1000), fill='blue'), dpi=600) +
  scale_x_continuous('chrM') + 
  scale_y_continuous('coverage') +
  theme_classic() + 
  theme(axis.text = element_text('Arial', size=8, color='black'),
        axis.title = element_text('Arial', size=8, color='black'),
        axis.title.y = element_blank())

cowplot::plot_grid(plotlist = list(p,q,r), ncol = 3)
ggsave('./figures/coverage/coverage_chrM.svg', width = 5.5, height = 1)

p=ggplot() + 
  ggrastr::rasterize(geom_col(data=MDS.coverage[which(MDS.coverage$V1 == 'chr21'),], aes(x=-pos, y=V5/1000), fill='grey30'), dpi=600) +
  scale_x_continuous('U2AF1', labels = seq(0,927,250), breaks = seq(-927,0,250), limits = c(-927,0)) + 
  scale_y_continuous('coverage') +
  theme_classic() + 
  theme(axis.text = element_text('Arial', size=8, color='black'),
        axis.title = element_text('Arial', size=8, color='black'))
q=ggplot() + 
  ggrastr::rasterize(geom_col(data=MDS.targeted.coverage[which(MDS.targeted.coverage$V1 == 'chr21'),], aes(x=-pos, y=V5/1000), fill='orange'), dpi=600) +
  scale_x_continuous('U2AF1', labels = seq(0,927,250), breaks = seq(-927,0,250), limits = c(-927,0)) + 
  scale_y_continuous('coverage') +
  theme_classic() + 
  theme(axis.text = element_text('Arial', size=8, color='black'),
        axis.title = element_text('Arial', size=8, color='black'),
        axis.title.y = element_blank())
r=ggplot() + 
  ggrastr::rasterize(geom_col(data=MDS.additionalPCR.coverage[which(MDS.additionalPCR.coverage$V1 == 'chr21'),], aes(x=-pos, y=V5/1000), fill='blue'), dpi=600) +
  scale_x_continuous('U2AF1', labels = seq(0,927,250), breaks = seq(-927,0,250), limits = c(-927,0)) + 
  scale_y_continuous('coverage') +
  theme_classic() + 
  theme(axis.text = element_text('Arial', size=8, color='black'),
        axis.title = element_text('Arial', size=8, color='black'),
        axis.title.y = element_blank())

cowplot::plot_grid(plotlist = list(p,q,r), ncol = 3, align = 'v')
ggsave('./figures/coverage/coverage_U2AF1.svg', width = 5.5, height = 1)

p=ggplot() + 
  geom_col(data=MDS.coverage[which(MDS.coverage$V1 == 'chr12'),], aes(x=pos, y=V5/1000), fill='grey30') +
  scale_x_continuous('GAPDH') + 
  scale_y_continuous('coverage') +
  theme_classic() + 
  theme(axis.text = element_text('Arial', size=8, color='black'),
        axis.title = element_text('Arial', size=8, color='black'))
q=ggplot() + 
  geom_col(data=MDS.targeted.coverage[which(MDS.targeted.coverage$V1 == 'chr12'),], aes(x=pos, y=V5/1000), fill='orange') +
  scale_x_continuous('GAPDH') + 
  theme_classic() + 
  theme(axis.text = element_text('Arial', size=8, color='black'),
        axis.title = element_text('Arial', size=8, color='black'),
        axis.title.y = element_blank())
r=ggplot() + 
  geom_col(data=MDS.additionalPCR.coverage[which(MDS.additionalPCR.coverage$V1 == 'chr12'),], aes(x=pos, y=V5/1000), fill='blue') +
  scale_x_continuous('GAPDH') + 
  theme_classic() + 
  theme(axis.text = element_text('Arial', size=8, color='black'),
        axis.title = element_text('Arial', size=8, color='black'),
        axis.title.y = element_blank())

cowplot::plot_grid(plotlist = list(p,q,r), ncol = 3, align = 'v')
ggsave('./figures/coverage/coverage_GAPDH.svg', width = 5.5, height = 1)