library(data.table)
library(tidyverse)
library(ggrepel)
library(cowplot)

eig <- fread("~/Downloads/GIAB_1KGP.eigenvec") %>%
  setnames(., "V1", "sample_id")

igsr <- fread("~/Downloads/igsr_samples.tsv")[, c(1, 4, 6)] %>%
  setnames(., c("sample_id", "pop", "superpop"))

eig <- merge(eig, igsr, by = "sample_id", all.x = TRUE)

eig[sample_id %in% c("normal", "tumor"), superpop := "GIAB"]

eig$superpop <- factor(eig$superpop, levels = c("AFR", "AMR", "EAS", "EUR", "SAS", "GIAB"))

a <- ggplot(data = eig, aes(x = V3, y = V4, color = superpop)) + 
  geom_point() +
  geom_text_repel(data = eig[superpop == "GIAB"], aes(label = sample_id), size = 3, max.overlaps = 1000, color = "black") +
  theme_bw() +
  scale_color_brewer(palette = "Set2", name = "") +
  theme(panel.grid = element_blank(), legend.position = "none") +
  xlab("PC1") +
  ylab("PC2")

b <- ggplot(data = eig, aes(x = V3, y = V5, color = superpop)) + 
  geom_point() +
  geom_text_repel(data = eig[superpop == "GIAB"], aes(label = sample_id), size = 3, max.overlaps = 1000, color = "black") +
  theme_bw() +
  scale_color_brewer(palette = "Set2", name = "") +
  theme(panel.grid = element_blank()) +
  xlab("PC1") +
  ylab("PC3")

cplot <- plot_grid(a, b, labels = c("A.", "B."), rel_widths = c(0.44, 0.53))

ggsave("~/Downloads/giab_ancestry.pdf", cplot, device = "pdf", width = 8, height = 3.3)


