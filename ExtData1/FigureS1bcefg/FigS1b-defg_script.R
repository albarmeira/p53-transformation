#CNV profiles of TP53-sAML patients

library(dplyr)
library(stringr)
library(data.table)
library(glue)
library(ggforce)
library(patchwork)
library(survival)
library(survminer)
library(cowplot)
library(ggpubr)
library(ggrepel)
library(readr)
library(dplyr)
library(ggplot2)
library(stringr)
library(readxl)

#Import filtered data
TP53_calls_filtered_full <- read_csv("data/TP53_noCP_IF9.csv")

cytoBand <- read_delim("data/cytoBand.txt",
                       "\t", escape_double = FALSE, col_names = FALSE,
                       trim_ws = TRUE)

chrom_lens <- read_delim("data/chrom_lens.tsv",
                         "\t", escape_double = FALSE, trim_ws = TRUE)

acen <- read_delim("data/acen.tsv",
                   "\t", escape_double = FALSE, col_names = FALSE,
                   trim_ws = TRUE)

chrom_levels = c(as.character(1:22), 'X')

###### CNV plots for all chromosomes in TP53-sAML patients

cytoband = cytoBand %>% setNames(c('chrom', 'start', 'end', 'name', 'stain'))
cytoband = cytoband %>% 
  mutate(chrom = str_remove(chrom, 'chr')) %>%
  mutate(chrom = ifelse(chrom == 'X', 23, chrom)) %>%
  filter(chrom != 'Y') %>%
  mutate(chrom = as.integer(chrom))

colnames(acen) = c('chrom', 'start', 'end', 'arm', 'type')

acen = acen %>% 
  mutate(chrom = str_remove(chrom, 'chr')) %>%
  filter(chrom != 'Y') %>%
  mutate(
    chrom = ifelse(chrom == 'X', 23, chrom),
    chrom = as.integer(chrom)
  ) %>%
  arrange(chrom) %>%
  group_by(chrom) %>%
  summarise(cen_start = min(start), cen_end = max(end)) %>%
  ungroup()

#acen_chrom = acen %>% filter(chrom == gene_loc$chrom)

cyto_colors = c(
  'gpos100'= rgb(0/255.0,0/255.0,0/255.0),
  'gpos'   = rgb(0/255.0,0/255.0,0/255.0),
  'gpos75' = rgb(130/255.0,130/255.0,130/255.0),
  'gpos66' = rgb(160/255.0,160/255.0,160/255.0),
  'gpos50' = rgb(200/255.0,200/255.0,200/255.0),
  'gpos33' = rgb(210/255.0,210/255.0,210/255.0),
  'gpos25' = rgb(200/255.0,200/255.0,200/255.0),
  'gvar'   = rgb(220/255.0,220/255.0,220/255.0),
  'gneg'  = rgb(255/255.0,255/255.0,255/255.0),
  'acen'  = rgb(217/255.0,47/255.0,39/255.0),
  'stalk' = rgb(100/255.0,127/255.0,164/255.0)
)
cnv_pal2 = c('CN-LOH' = 'darkgreen', 'Gain' = 'red', 'Loss' = 'blue', 'Undetermined' = 'darkorange')

all_colors = c(
  'Gain' = 'darkred', 'CN-LOH' = 'darkgreen', 'Loss' = 'blue',
  'Regulatory' = 'darkturquoise', 'Missense' = 'salmon', 'Truncating' = 'purple') 


###CNVs plotted across chromosomes
mocha_calls_filtered2<- TP53_calls_filtered_full
mocha_calls_filtered2 %>% pull(sample_id) %>% unique %>% length %>% paste('unique individuals with CNV')
mocha_calls_filtered2 %>% count(type) %>%mutate(prop = signif(n/sum(n),2))
mocha_calls_filtered2 %>% count(sample_id) %>% count(n) %>% mutate(prop = signif(nn/sum(nn), 2))


TP53_cnvs = mocha_calls_filtered2 %>%
  mutate(chrom = factor(ifelse(chrom == 23, 'X', as.character(chrom)), chrom_levels)) %>%
  group_by(sample_id, chrom) %>%
  mutate(total_length = sum(length)) %>%
  mutate(type = factor(type, c('CN-LOH', 'Gain', 'Loss', 'Undetermined'))) %>%
  arrange(type, round(beg_GRCh37/1e7), -total_length) %>%
  group_by(chrom) %>%
  mutate(index = as.integer(factor(sample_id, unique(sample_id)))) %>% 
  ungroup() %>%
  ggplot(
    aes(x = beg_GRCh37, xend = end_GRCh37, y = index, yend = index, color = type)
  ) +
  geom_rect(
    inherit.aes = F,
    data = chrom_lens %>% mutate(chrom = factor(ifelse(chrom == 23, 'X', as.character(chrom)), chrom_levels)),
    aes(xmin = 0, xmax = length, ymin = -2, ymax = -0.5),
    color = 'black',
    fill = 'white',
    size = 1
  ) +
  geom_rect(
    data = cytoband %>% mutate(chrom = factor(ifelse(chrom == 23, 'X', as.character(chrom)), chrom_levels)),
    inherit.aes = F,
    aes(xmin = start, xmax = end, ymin = -2, ymax = -0.5, fill = stain),
    size = 1,
    show.legend = FALSE
  ) +
  geom_segment(
    lineend = 'round',
    size = 1
  ) +
  scale_fill_manual(values = cyto_colors) +
  theme_void() +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = 'right',
    panel.spacing.y = unit(0, "lines"),
    plot.margin = unit(c(0,0,0,0), "lines"),
    strip.text.y = element_text(size = 10, hjust = 0.5, angle = 0),
    panel.border = element_blank(),
    strip.background = element_blank()
  ) +
  scale_y_discrete(expand = expansion(add = 1)) +
  scale_x_discrete(expand = c(0.02,0)) +
  facet_grid(.~chrom, space = 'free', scales = 'free', switch="both") +
  xlab('') +
  ylab('') +
  scale_color_manual(values = cnv_pal2) +
  labs(color = '')
TP53_cnvs #FigS1b.png

#### TP53 gene CNV plots for TP53-sAML patients

TP53_cna <- read_delim("data/TP53_CNA_tp53pt_noCP_IF9.csv",
                       ",", escape_double = FALSE, trim_ws = TRUE)

TP53_M <- read_excel("data/TP53_M_tp53pt_noCP_IF9.xlsx")

acen_17 <- read_excel("data/acen_17.xlsx")
chr17_lens <- read_csv("data/chr17_lens.csv")
TP53_gene_structure <- read_excel("data/tp53_gene_structure.xlsx")
TP53_loc <- read_excel("data/TP53_loc.xlsx")

TP53_cna_format = TP53_cna %>%
  filter(chrom == chrom) %>%
  mutate(chrom = factor(chrom)) %>%
  group_by(sample_id, chrom) %>%
  mutate(total_length = sum(size)) %>%
  mutate(type = factor(type, c('Gain', 'Loss', 'Undetermined', 'CN-LOH'))) %>%
  arrange(desc(type), round(start/1e7), -size) %>%
  group_by(chrom) %>%
  mutate(index = as.integer(index)) %>%
  ungroup()

TP53_cna_mut = ggplot( #this top bit determines the colours in the top plot
  TP53_cna_format,
  aes(x = start, xend = end, y = index, yend = index, color = type)
) + #this sets the centromere bit
  geom_rect(
    inherit.aes = F,
    data = acen_17,
    aes(xmin = cen_start, xmax = cen_end, ymin = -Inf, ymax = Inf),
    fill = 'gray',
    size = 0,
    alpha = 0.3
  ) + #this determines the size and rounded end to the ins/del bit
  geom_segment(
    size = 1,
    alpha = 0.8,
    lineend = 'round'
  ) + #this sets the bottom bit of the plot
  geom_segment(
    inherit.aes = F,
    data = chr17_lens,
    aes(x = 0, xend = length, y = 0, yend = 0),
    color = 'black',
    size = 0
  ) +
  scale_alpha(
    range = c(0, 1)
  ) +
  theme_classic() +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = 'none',
    panel.spacing.y = unit(1, "mm"),
    plot.margin = unit(c(0.5,0.2,0.2,0.2), "lines"),
    strip.text.y = element_blank(),
    panel.border = element_rect(size = 1, fill = NA, color = 'gray'),
    strip.background = element_rect(size = 0.5, fill = 'seashell3', color = 'seashell3'),
    plot.title = element_text(size = 12)
  ) +
  scale_y_discrete(expand = expansion(add = 1)) +
  xlab('') +
  ylab('') +
  ggtitle("TP53 mutation/CNAs")+
  scale_color_manual(values = all_colors) +
  scale_fill_manual(values = all_colors) +
  geom_point(
    inherit.aes = F,
    mapping = aes(x = start, y = index, fill = VariantClass2),
    data = TP53_M,
    size = 1.5,
    color = 'black',
    pch = 21
  ) +
  geom_point(
    inherit.aes = F,
    mapping = aes(x = start3, y = index, fill = VariantClass3),
    data = TP53_M,
    size = 1.5,
    color = 'black',
    pch = 21
  ) +
  geom_point(
    inherit.aes = F,
    mapping = aes(x = start4, y = index, fill = VariantClass4),
    data = TP53_M,
    size = 1.5,
    color = 'black',
    pch = 21
  ) +
  # exon structures
  geom_segment(
    inherit.aes = F,
    data = TP53_gene_structure %>% mutate(zoom = TRUE),
    aes(x = min(start), xend = max(end), y = -1, yend = -1),
    size = 1,
    color = 'black'
  ) +
  geom_segment(
    inherit.aes = F,
    data = TP53_gene_structure %>% mutate(zoom = TRUE),
    aes(x = start, xend = end, y = -1, yend = -1),
    size = 8,
    color = 'black'
  ) +
  facet_zoom(
    zoom.size = 1,
    zoom.data = zoom,
    xlim = c(TP53_loc$start, TP53_loc$end)
  )
TP53_cna_mut #FigS1d.png

#### JAK2 gene CNV plots for TP53-sAML patients

JAK2_cna <- read_delim("data/JAK2_cna_tp53.txt",
                       "\t", escape_double = FALSE, trim_ws = TRUE)

JAK_M <- read_excel("data/JAK_M_tp53.xlsx")

acen_9 <- read_excel("data/acen_9.xlsx")
chr9_lens <- read_csv("data/chr9_lens.csv")
JAK2_gene_structure <- read_excel("data/JAK2_gene_structure.xlsx")
JAK2_loc <- read_excel("data/JAK2_loc.xlsx")

JAK_M2<- JAK_M
JAK_M2$index<- JAK_M$index2

JAK2_cna_format = JAK2_cna %>%
  filter(chrom == chrom) %>%
  mutate(chrom = factor(chrom)) %>%
  group_by(sample_id, chrom) %>%
  mutate(total_length = sum(size)) %>%
  mutate(type = factor(type, c('Gain', 'Loss', 'Undetermined', 'CN-LOH'))) %>%
  arrange(desc(type), round(start/1e7), -size) %>%
  group_by(chrom) %>%
  mutate(index = as.integer(index)) %>%
  ungroup()

JAK2_cna_mut = ggplot( #this top bit determines the colours in the top plot
  JAK2_cna_format,
  aes(x = start, xend = end, y = index, yend = index, color = type)
) + #this sets the centromere bit
  geom_rect(
    inherit.aes = F,
    data = acen_9,
    aes(xmin = cen_start, xmax = cen_end, ymin = -Inf, ymax = Inf),
    fill = 'gray',
    size = 0,
    alpha = 0.3
  ) + #this determines the size and rounded end to the ins/del bit
  geom_segment(
    size = 1,
    alpha = 0.8,
    lineend = 'round'
  ) + #this sets the bottom bit of the plot
  geom_segment(
    inherit.aes = F,
    data = chr9_lens,
    aes(x = 0, xend = length, y = 0, yend = 0),
    color = 'black',
    size = 0
  ) +
  scale_alpha(
  range = c(0, 1)
) +
  theme_classic() +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = 'none',
    panel.spacing.y = unit(1, "mm"),
    plot.margin = unit(c(0.5,0.2,0.2,0.2), "lines"),
    strip.text.y = element_blank(),
    panel.border = element_rect(size = 1, fill = NA, color = 'gray'),
    strip.background = element_rect(size = 0.5, fill = 'seashell3', color = 'seashell3'),
    plot.title = element_text(size = 12)
  ) +
  scale_y_discrete(expand = expansion(add = 1)) +
  xlab('') +
  ylab('') +
  ggtitle("JAK2 mutation/CNAs")+
  scale_color_manual(values = all_colors) +
  scale_fill_manual(values = all_colors) +
  geom_point(
    inherit.aes = F,
    mapping = aes(x = start, y = index, fill = VariantClass2),
    data = JAK_M2,
    size = 1.5,
    color = 'black',
    pch = 21
  ) +
  # exon structures
  geom_segment(
    inherit.aes = F,
    data = JAK2_gene_structure %>% mutate(zoom = TRUE),
    aes(x = min(start), xend = max(end), y = -1, yend = -1),
    size = 1,
    color = 'black'
  ) +
  geom_segment(
    inherit.aes = F,
    data = JAK2_gene_structure %>% mutate(zoom = TRUE),
    aes(x = start, xend = end, y = -1, yend = -1),
    size = 8,
    color = 'black'
  ) +
  facet_zoom(
    zoom.size = 1,
    zoom.data = zoom,
    xlim = c(JAK2_loc$start, JAK2_loc$end)
  )
JAK2_cna_mut #FigS1e.png

####EZH2 gene CNV plots for TP53-sAML patients

EZH2_cna <- read_csv("data/EZH2_CNA_tp53pt.csv")
EZH2_M <- read_excel("data/EZH2_M_tp53.xlsx")
acen_7 <- read_excel("data/acen_7.xlsx")
chr7_lens <- read_csv("data/chr7_lens.csv")
EZH2_gene_structure <- read_excel("data/EZH2_gene_structure.xlsx")
EZH2_loc <- read_excel("data/EZH2_loc.xlsx")


EZH2_cna_format = EZH2_cna %>%
  filter(chrom == chrom) %>%
  mutate(chrom = factor(chrom)) %>%
  group_by(sample_id, chrom) %>%
  mutate(total_length = sum(size)) %>%
  mutate(type = factor(type, c('Gain', 'Loss', 'Undetermined', 'CN-LOH'))) %>%
  arrange(desc(type), round(start/1e7), -size) %>%
  group_by(chrom) %>%
  mutate(index = as.integer(index)) %>%
  ungroup()



EZH2_cna %>%arrange(desc(type),round(start/1e7), -total_length)


EZH2_cna_mut = ggplot( #this top bit determines the colours in the top plot
  EZH2_cna_format,
  aes(x = start, xend = end, y = index, yend = index, color = type)
) + #this sets the centromere bit
  geom_rect(
    inherit.aes = F,
    data = acen_7,
    aes(xmin = cen_start, xmax = cen_end, ymin = -Inf, ymax = Inf),
    fill = 'gray',
    size = 0,
    alpha = 0.3
  ) + #this determines the size and rounded end to the ins/del bit
  geom_segment(
    size = 1,
    alpha = 0.8,
    lineend = 'round'
  ) + #this sets the bottom bit of the plot
  geom_segment(
    inherit.aes = F,
    data = chr7_lens,
    aes(x = 0, xend = length, y = 0, yend = 0),
    color = 'black',
    size = 0
  ) +
  scale_alpha(
    range = c(0, 1)
  ) +
  theme_classic() +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = 'none',
    panel.spacing.y = unit(1, "mm"),
    plot.margin = unit(c(0.5,0.2,0.2,0.2), "lines"),
    strip.text.y = element_blank(),
    panel.border = element_rect(size = 1, fill = NA, color = 'gray'),
    strip.background = element_rect(size = 0.5, fill = 'seashell3', color = 'seashell3'),
    plot.title = element_text(size = 12)
  ) +
  scale_y_discrete(expand = expansion(add = 1)) +
  xlab('') +
  ylab('') +
  ggtitle("EZH2 mutation/CNAs")+
  scale_color_manual(values = all_colors) +
  scale_fill_manual(values = all_colors) +
  geom_point(
    inherit.aes = F,
    mapping = aes(x = start, y = index, fill = VariantClass2),
    data = EZH2_M,
    size = 1.5,
    color = 'black',
    pch = 21
  ) +
  # exon structures
  geom_segment(
    inherit.aes = F,
    data = EZH2_gene_structure %>% mutate(zoom = TRUE),
    aes(x = min(start), xend = max(end), y = -1, yend = -1),
    size = 1,
    color = 'black'
  ) +
  geom_segment(
    inherit.aes = F,
    data = EZH2_gene_structure %>% mutate(zoom = TRUE),
    aes(x = start, xend = end, y = -1, yend = -1),
    size = 8,
    color = 'black'
  ) +
  facet_zoom(
    zoom.size = 1,
    zoom.data = zoom,
    xlim = c(EZH2_loc$start, EZH2_loc$end)
  )

EZH2_cna_mut #FigS1f.png

####TET2 mutation/CNV plot

TET2_cna <- read_excel("data/TET2_CNA_tp53pt_noCP_IF9.xls")
TET2_M <- read_excel("data/TET2_M_tp53_noCP_IF9.xlsx")

acen_4 <- read_excel("data/acen_4.xlsx")
chr4_lens <- read_csv("data/chr4_lens.csv")
TET2_gene_structure <- read_excel("data/TET2_gene_structure.xlsx")
TET2_loc <- read_excel("data/TET2_loc.xlsx")

TET2_cna_format = TET2_cna %>%
  filter(chrom == chrom) %>%
  mutate(chrom = factor(chrom)) %>%
  group_by(sample_id, chrom) %>%
  mutate(total_length = sum(size)) %>%
  mutate(type = factor(type, c('Gain', 'Loss', 'Undetermined', 'CN-LOH'))) %>%
 arrange(desc(type), round(start/1e7), -size) %>%
  group_by(chrom) %>%
  mutate(index = as.integer(index)) %>%
  ungroup()

TET2_cna_mut = ggplot( #this top bit determines the colours in the top plot
  TET2_cna_format,
  aes(x = start, xend = end, y = index, yend = index, color = type)
) + #this sets the centromere bit
  geom_rect(
    inherit.aes = F,
    data = acen_4,
    aes(xmin = cen_start, xmax = cen_end, ymin = -Inf, ymax = Inf),
    fill = 'gray',
    size = 0,
    alpha = 0.3
  ) + #this determines the size and rounded end to the ins/del bit
  geom_segment(
    size = 1,
    alpha = 0.8,
    lineend = 'round'
  ) + #this sets the bottom bit of the plot
  geom_segment(
    inherit.aes = F,
    data = chr4_lens,
    aes(x = 0, xend = length, y = 0, yend = 0),
    color = 'black',
    size = 0
  ) +
scale_alpha(
  range = c(0, 1)
) +
  theme_classic() +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = 'none',
    panel.spacing.y = unit(1, "mm"),
    plot.margin = unit(c(0.5,0.2,0.2,0.2), "lines"),
    strip.text.y = element_blank(),
    panel.border = element_rect(size = 1, fill = NA, color = 'gray'),
    strip.background = element_rect(size = 0.5, fill = 'seashell3', color = 'seashell3'),
    plot.title = element_text(size = 12)
  ) +
  scale_y_discrete(expand = expansion(add = 1)) +
  xlab('') +
  ylab('') +
  ggtitle("TET2 mutation/CNAs")+
  scale_color_manual(values = all_colors) +
  scale_fill_manual(values = all_colors) +
  geom_point(
    inherit.aes = F,
    mapping = aes(x = start, y = index, fill = VariantClass2),
    data = TET2_M,
    size = 1.5,
    color = 'black',
    pch = 21
  ) +
  # exon structures
  geom_segment(
    inherit.aes = F,
    data = TET2_gene_structure %>% mutate(zoom = TRUE),
    aes(x = min(start), xend = max(end), y = -1, yend = -1),
    size = 1,
    color = 'black'
  ) +
  geom_segment(
    inherit.aes = F,
    data = TET2_gene_structure %>% mutate(zoom = TRUE),
    aes(x = start, xend = end, y = -1, yend = -1),
    size = 8,
    color = 'black'
  ) +
  facet_zoom(
    zoom.size = 1,
    zoom.data = zoom,
    xlim = c(TET2_loc$start, TET2_loc$end)
  )

TET2_cna_mut #FigS1g.png


