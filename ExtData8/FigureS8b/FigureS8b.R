#TP53 CNAs in CP-TP53-MPN patients (no transformation)

library(readxl)
library(readr)
library(dplyr)
library(ggplot2)
library(ggforce)

TP53_cna <- read_delim("data/TP53_CNA_tp53pt_IF9.csv",
                       ",", escape_double = FALSE, trim_ws = TRUE)

TP53_M <- read_excel("data/TP53_M_tp53pt_IF9.xlsx")

View(TP53_cna)

acen_17 <- read_excel("data/acen_17.xlsx")
chr17_lens <- read_csv("data/chr17_lens.csv")
TP53_gene_structure <- read_excel("data/tp53_gene_structure.xlsx")
TP53_loc <- read_excel("data/TP53_loc.xlsx")

D2 = TP53_cna %>%
  filter(chrom == chrom) %>%
  mutate(chrom = factor(chrom)) %>%
  group_by(sample_id, chrom) %>%
  mutate(total_length = sum(size)) %>%
  mutate(type = factor(type, c('Gain', 'Loss', 'Undetermined', 'CN-LOH'))) %>%
  arrange(desc(type), round(start/1e7), -size) %>%
  group_by(chrom) %>%
  mutate(index = as.integer(index)) %>%
  ungroup()


all_colors = c(
  'Gain' = 'darkred', 'CN-LOH' = 'darkgreen', 'Loss' = 'blue',
  'Regulatory' = 'darkturquoise', 'Missense' = 'salmon', 'Truncating' = 'purple') 

p_segs4 = ggplot( #this top bit determines the colours in the top plot
  D2,
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
  # label the mutations
  #         geom_point(
  #             inherit.aes = F,
  #             mapping = aes(x = start, y = index, color = VariantClass2),
  #             data = M_overlap,
  #             size = 6.5,
  #             pch = '*',
  #             color = 'white'
  #         ) +
  #label the mutations by type 
  #         geom_text_repel(
#             inherit.aes = F,
#             mapping = aes(x = start, y = index, color = VariantClass2, label = AAchange),
#             data = M_overlap %>% mutate(zoom = TRUE),
#             size = 1.5,
#             segment.size = 0.2,
#             color = 'black',
#             ylim = c(0,NA)
#         ) +
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

ggsave(filename = "FigS8b.pdf", plot = p_segs4, width = 4, height = 4)
