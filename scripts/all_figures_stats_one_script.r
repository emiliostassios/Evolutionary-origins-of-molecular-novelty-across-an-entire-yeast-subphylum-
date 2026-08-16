library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(readxl)
library(ggtree)
library(ggpubr)
library(RColorBrewer)
library(ggnewscale)
library(ppcor)
library(ggridges)
library(svglite) 

##Read the large dataset which contains the gene/protein data
df <- read.csv('data/matrices/all_matrix_5_8_26.tsv', sep ='\t')

yeast_genera <- c('Saccharomyces', 'Candida', 'Yarrowia','Wickerhamomyces','Pichia', 'Ogataea', 'Lachancea', 'Lipomyces', 'Ambrosiozyma', 
                  'Kazachstania', 'Barnettozyma', 'Blastobotrys', 'Cyberlindnera', 'Debaryomyces', 'Hanseniaspora', 'Metschnikowia','Saturnispora',
                  'Spathaspora', 'Tetrapisispora', 'Torulaspora')

#Supp Figure 5
plot_data <- df %>%
  dplyr::filter(updated_final_type == "denovo") %>% 
  group_by(species) %>%
  summarise(
    n_denovo = n(),
    outgroup_age = dplyr::first(Max_Outgroup_Distance) 
  ) %>%
  dplyr::filter(!is.na(outgroup_age), !is.na(n_denovo))

cor_res <- cor.test(plot_data$outgroup_age, plot_data$n_denovo, method = "spearman")
rho_val <- round(cor_res$estimate, 2)
p_val <- format(cor_res$p.value, scientific = TRUE, digits = 7)
anno_text <- paste0("rho = ", rho_val, ", p = ", p_val)

ggplot(plot_data, aes(x = outgroup_age, y = n_denovo)) +
  geom_point(color = "#5c5ca8", alpha = 0.8, size = 1.5) + 
  geom_smooth(method = "lm", color = "darkred", fill = "grey75") +
  annotate(
    "text", 
    x = -Inf, y = Inf, 
    label = anno_text, 
    hjust = -0.1, vjust = 1.5, 
    fontface = "bold", 
    size = 6
  ) +
  labs(
    x = "Age of Closest Outgroup (Million Years)",
    y = "Number of De Novo Genes"
  ) +
  theme_bw(base_size = 16) +
  theme(
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "grey90", fill = NA),
    axis.text = element_text(color = "grey30")
  )

###load the candidate denovo matrix and the proteomic similarity natrix
candidate_denovos <- read.csv('data/matrices/denovo_candidates.tsv', sep ='\t')
proteomic_similarity <- read.csv('data/matrices/proteomic_similarity.tsv', sep ='\t')

###all vs alla pairs
my_comparisons <- list(
  c("De novo", "False Positive (change <20Mya)"), 
  c("De novo", "HGT"), 
  c("De novo", "False Positive (change >20Mya)"),
  c("False Positive (change >20Mya)", "False Positive (change <20Mya)"),
  c("False Positive (change >20Mya)", "HGT"),
  c("False Positive (change <20Mya)", "HGT")  
)

####assigning the prot count numbers to the proteins
prot_counts <- proteomic_similarity %>%
  group_by(code_name) %>%
  summarise(prot_hits = n(), .groups = "drop")

prot_distribution <- candidate_denovos %>%
  left_join(prot_counts, by = "code_name") %>%
  mutate(prot_hits = replace_na(prot_hits, 0))


##Sup figure 6
ggplot(candidate_denovos, aes(x = status_col, fill = status_col)) +
  geom_bar(color = "black", linewidth = 0.5, alpha = 0.8) + 
  geom_text(stat = "count", aes(label = after_stat(count)), 
            vjust = -0.5, size = 6, fontface = "bold") +
  scale_fill_manual(values = c("#FAD510FF", "royalblue1", "royalblue3", "forestgreen")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) + 
  labs(
    x = NULL,
    y = "Count"
  ) +
  theme_classic(base_size = 18) + theme_minimal() +
  theme(
    axis.text.x = element_blank(),  
    axis.ticks.x = element_blank(), 
    axis.text.y = element_text(size = 16, color = "black"),
    axis.title.y = element_text(size = 18, face = "bold", margin = margin(r = 15)),
    axis.line = element_line(color = "black", linewidth = 0.8),
    legend.position = "right", 
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 14),
    legend.key.size = unit(1, "cm"),
    legend.key = element_rect(color = "black", linewidth = 0.5) 
  ) + 
  guides(fill = guide_legend(title = 'Gene Type'))


ggplot(prot_distribution, aes(x = status_col, y = prot_hits, fill = status_col)) +
  geom_boxplot(alpha = 0.8, color = "black", linewidth = 0.6, outlier.shape = 21, outlier.alpha = 0.6) + 
  scale_y_log10(expand = expansion(mult = c(0.05, 0.15))) +
  scale_fill_manual(values = c("#FAD510FF", "royalblue1", "royalblue4", "forestgreen")) +
  labs(x = NULL, y = "Proteomic Hits per Gene (log)") +
  theme_classic(base_size = 18) + theme_minimal() +
  theme(
    axis.text.x = element_blank(),  
    axis.ticks.x = element_blank(), 
    axis.text.y = element_text(size = 16, color = "black"),
    axis.title.y = element_text(size = 18, face = "bold", margin = margin(r = 15)),
    axis.line = element_line(color = "black", linewidth = 0.8),
    legend.position = "right", 
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 14),
    legend.key.size = unit(1, "cm"),
    legend.key = element_rect(color = "black", linewidth = 0.5) 
  ) + 
  guides(fill = guide_legend(title = 'Gene Type')) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", size = 5, tip.length = 0.02)


ggplot(candidate_denovos, aes(x = status_col, y = length_x, fill = status_col)) +
  geom_boxplot(alpha = 0.8, color = "black", linewidth = 0.6, outlier.shape = 21, outlier.alpha = 0.6) + 
  scale_y_log10(expand = expansion(mult = c(0.05, 0.15))) +
  scale_fill_manual(values = c("#FAD510FF", "royalblue1", "royalblue4", "forestgreen")) +
  labs(x = NULL, y = "Sequence Length") +
  theme_classic(base_size = 18) + theme_minimal()+
  theme(
    axis.text.x = element_blank(),  
    axis.ticks.x = element_blank(), 
    axis.text.y = element_text(size = 16, color = "black"),
    axis.title.y = element_text(size = 18, face = "bold", margin = margin(r = 15)),
    axis.line = element_line(color = "black", linewidth = 0.8),
    legend.position = "right", 
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 14),
    legend.key.size = unit(1, "cm"),
    legend.key = element_rect(color = "black", linewidth = 0.5) 
  ) + 
  guides(fill = guide_legend(title = 'Gene Type')) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", size = 5, tip.length = 0.02)


###Figure 3 tree

tree <- read.tree('data/trees/332_noded.nwk')
genes_per_species <- read.csv('data/matrices/genes_per_species.txt', sep ='\t', header =F)
clade_metadata <- read.csv('data/metadata/343taxa_speicies-name_clade-name_color-code.txt',sep ='\t')

denovo <- df[df$updated_final_type=='denovo',]
denovo_families <- data.frame(table(denovo$merged_families))
genes_per_species <- genes_per_species[match(tree$tip.label, genes_per_species$V2), ]

tree_data <- ggtree(tree)$data

tip_metadata_simple <- df %>%
  group_by(species, updated_final_type) %>%
  summarise(count = n(), .groups = 'drop') %>%
  pivot_wider(names_from = updated_final_type, values_from = count, values_fill = 0) %>%
  dplyr::rename(label = species)

inner_metadata <- df %>%
  group_by(recalculated_origin_forty, final_age_, updated_final_type) %>%
  summarise(count = n(), .groups = 'drop') %>%
  pivot_wider(names_from = updated_final_type, values_from = count, values_fill = 0) %>%
  dplyr::rename(label = recalculated_origin_forty)

tree_data <- tree_data %>%
  left_join(tip_metadata_simple, by = "label") %>%
  left_join(inner_metadata, by = "label", suffix = c("_tip", "_internal")) %>%
  mutate(denovo_value = coalesce(denovo_tip, denovo_internal))

tree_data$denovo_all_genes <- tree_data$denovo_value / genes_per_species$V1

clade_metadata <- clade_metadata[,c("old_speceis_names", "Major.clade")]
ordered_clades <- unique(clade_metadata$Major.clade)
num_clades <- length(ordered_clades)
clade_colors_vector <- colorRampPalette(brewer.pal(12, "Set3"))(num_clades)
clade_colors <- setNames(clade_colors_vector, ordered_clades)

tip_metadata_bar <- df %>%
  group_by(species, genus, updated_final_type) %>%
  summarise(count = n(), .groups = 'drop') %>%
  pivot_wider(names_from = updated_final_type, values_from = count, values_fill = 0)

tip_metadata_bar <- tip_metadata_bar %>%
  left_join(clade_metadata, by = c("species" = "old_speceis_names")) 

tip_metadata_bar <- tip_metadata_bar[match(tree$tip.label, tip_metadata_bar$species), ]
tip_metadata_bar$denovo_all_genes <- tip_metadata_bar$denovo / genes_per_species$V1
tip_metadata_bar <- tip_metadata_bar %>%
  mutate(Major.clade = factor(Major.clade, levels = ordered_clades))


p <- ggtree(tree, layout = "circular")

p <- p +
  geom_fruit(
    data = tip_metadata_bar,
    geom = geom_bar,
    mapping = aes(y = species, x = denovo_all_genes, fill = Major.clade), 
    stat = "identity",
    orientation = "y",
    width = 0.8, 
    offset = 0.05, 
    size = 0.4,
    grid.params = list(
      vline = TRUE, 
      color = "grey85",  
      size = 0.15,       
      linetype = "dashed"
    ),
    axis.params = list(
      axis = "x", 
      text.size = 2,     
      text.angle = 0, 
      hjust = 0.5,
      color = "grey50"   
    )
  ) +
  scale_fill_manual(
    values = clade_colors, 
    name = "Saccharomycotina Order" 
  )

p <- p %<+% tree_data +
  new_scale_color() + 
  geom_tree(aes(color = denovo_internal), size = 1) + 
  scale_color_gradient(
    low = "grey", 
    high = "black", 
    name = "de novo \n origination events",
    guide = guide_colorbar(position = "left") 
  ) +
  theme(
    legend.position = "right", 
    legend.justification = "top" 
  )

ggsave(
  "figures/figure_3c.svg", 
  plot = p,
  width = 8, 
  height = 6,
  device = svglite)


###figure 3 length vs conserved
df_ridge <- df %>%
  filter(genus %in% yeast_genera) %>%
  filter(updated_final_type %in% c("conserved", "denovo")) %>%
  mutate(length_x = as.numeric(length_x)) %>%
  filter(!is.na(length_x) & is.finite(length_x)) %>%
  mutate(updated_final_type = factor(updated_final_type, levels = c("conserved", "denovo")))

ggplot(df_ridge, aes(x = length_x, y = genus, fill = genus, linetype = updated_final_type)) +
  geom_density_ridges(alpha = 0.7, scale = 1.2, size = 0.6) +
  scale_x_log10(
    breaks = c(10, 100, 1000, 10000),
    labels = c("10", "100", "1000", "10000")
  ) +
  scale_linetype_manual(
    name = "Gene Type",
    values = c("conserved" = "solid", "denovo" = "dashed")
  ) +
  guides(
    fill = guide_legend(title = NULL, order = 1),
    linetype = guide_legend(order = 2, override.aes = list(fill = "gray80"))
  ) +
  labs(
    x = "CDS length (nt)",
    y = NULL
  ) +
  theme_minimal(base_size = 18) +
  theme(
    axis.text.y = element_text(color = "black", size = 16),
    axis.text.x = element_text(color = "black", size = 16),
    axis.title.x = element_text(color = "black", size = 18),
    panel.grid.major.x = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "grey90"),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 18),
    legend.text = element_text(size = 16)
  )


###Figure 4
seq_div <- df[df$updated_final_type =='seq_div',]
size_counts <- seq_div %>%
  group_by(merged_families) %>%
  summarise(family_size = n(), .groups = "drop") %>%
  group_by(family_size) %>%
  summarise(count_of_families = n(), .groups = "drop")

ggplot(size_counts, aes(x = family_size, y = count_of_families)) +
  geom_col(fill = "#b30000", color = "black", width = 1) +
  geom_text(aes(label = count_of_families), vjust = -0.8, fontface = "bold", size = 5, check_overlap = TRUE) +
  scale_y_log10(
    labels = scales::comma,
    expand = expansion(mult = c(0, 0.2)) 
  ) +
  scale_x_continuous(breaks = seq(min(size_counts$family_size), max(size_counts$family_size), by = 1)) +
  labs(
    x = "Number of Divergent Genes per Family",
    y = "Count of Families (Log10 Scale)"
  ) +
  theme_classic(base_size = 18) + 
  theme(
    axis.title.x = element_text(face = "bold", size = 22, margin = margin(t = 15)),
    axis.title.y = element_text(face = "bold", size = 22, margin = margin(r = 15)),
    axis.text = element_text(size = 16, color = "black", face = "bold"),
    axis.line = element_line(linewidth = 1.2, color = "black"),
    axis.ticks = element_line(linewidth = 1.2, color = "black"),
    panel.grid.major.y = element_line(color = "grey80", linetype = "dashed", linewidth = 0.8),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )

###figure 4 c
seq_div_data <- df[df$updated_final_type == 'seq_div', ]
seq_div_families <- data.frame(table(seq_div_data$merged_families)) 

genes_per_species <- genes_per_species[match(tree$tip.label, genes_per_species$V2), ]
tree_data <- ggtree(tree)$data

tip_metadata_simple <- df %>%
  group_by(species, updated_final_type) %>%
  summarise(count = n(), .groups = 'drop') %>%
  pivot_wider(names_from = updated_final_type, values_from = count, values_fill = 0) %>%
  dplyr::rename(label = species)

inner_metadata <- df %>%
  group_by(recalculated_origin_forty, final_age_, updated_final_type) %>%
  summarise(count = n(), .groups = 'drop') %>%
  pivot_wider(names_from = updated_final_type, values_from = count, values_fill = 0) %>%
  dplyr::rename(label = recalculated_origin_forty)

tree_data <- tree_data %>%
  left_join(tip_metadata_simple, by = "label") %>%
  left_join(inner_metadata, by = "label", suffix = c("_tip", "_internal")) %>%
  mutate(seq_div_value = coalesce(seq_div_tip, seq_div_internal))

tree_data$seq_div_all_genes <- tree_data$seq_div_value / genes_per_species$V1

clade_metadata <- read.csv('data/metadata/343taxa_speicies-name_clade-name_color-code.txt',sep ='\t')
clade_metadata <- clade_metadata[,c("old_speceis_names", "Major.clade")]
ordered_clades <- unique(clade_metadata$Major.clade)
num_clades <- length(ordered_clades)
clade_colors_vector <- colorRampPalette(brewer.pal(12, "Set3"))(num_clades)
clade_colors <- setNames(clade_colors_vector, ordered_clades)

tip_metadata_bar <- df %>%
  group_by(species, genus, updated_final_type) %>%
  summarise(count = n(), .groups = 'drop') %>%
  pivot_wider(names_from = updated_final_type, values_from = count, values_fill = 0)

tip_metadata_bar <- tip_metadata_bar %>%
  left_join(clade_metadata, by = c("species" = "old_speceis_names")) 

tip_metadata_bar <- tip_metadata_bar[match(tree$tip.label, tip_metadata_bar$species), ]
tip_metadata_bar$seq_div_all_genes <- tip_metadata_bar$seq_div / genes_per_species$V1
tip_metadata_bar <- tip_metadata_bar %>%
  mutate(Major.clade = factor(Major.clade, levels = ordered_clades))

p <- ggtree(tree, layout = "circular")

p <- p +
  geom_fruit(
    data = tip_metadata_bar,
    geom = geom_bar,
    mapping = aes(y = species, x = seq_div_all_genes, fill = Major.clade), 
    stat = "identity",
    orientation = "y",
    width = 0.8, 
    offset = 0.05, 
    size = 0.4,
    grid.params = list(
      vline = TRUE, 
      color = "grey85",  
      size = 0.15,       
      linetype = "dashed"
    ),
    axis.params = list(
      axis = "x", 
      text.size = 2,     
      text.angle = 0, 
      hjust = 0.5,
      color = "grey50"   
    )
  ) +
  scale_fill_manual(
    values = clade_colors, 
    name = "Saccharomycotina Order" 
  )


p <- p %<+% tree_data +
  new_scale_color() + 
  geom_tree(aes(color = seq_div_internal), size = 1) + 
  scale_color_gradient(
    low = "grey", 
    high = "black", 
    name = "Divergent origination events",
    guide = guide_colorbar(position = "left") 
  ) +
  theme(
    legend.position = "right", 
    legend.justification = "top" 
  )

ggsave(
  "figures/figure_4c.svg", 
  plot = p,
  width = 8, 
  height = 6,
  device = svglite)


###property comparison
properties_to_plot <- c("length_x", "gc", "isoelectric", "cost", "propensity", "diso_pct" ,"hydropathicity", "nuc_coverage")
length_col <- rlang::sym("length_x")
cost_col   <- rlang::sym("cost")
tm_col     <- rlang::sym("propensity")
properties_to_pivot <- c("length_x", "gc", "isoelectric", "propensity", "diso_pct", "hydropathicity", "nuc_coverage", "Normalized Cost")
my_colors <- c("denovo" = "#FAD510FF", "seq_div" = "firebrick")

property_labels <- c(
  "Length (aa) (log scale)", "GC Content (%)", "Isoelectric Point (pH)", 
  "TM propensity", "Disorder Pct (%)", 
  "Hydropathicity (GRAVY)", "Nucleosomal coverage", "Normalized Cost"
)
names(property_labels) <- properties_to_pivot


df_long <- df %>%
  filter(updated_final_type %in% c("denovo", "seq_div")) %>%
  mutate(across(all_of(properties_to_plot), as.numeric)) %>%
  mutate(
    `Normalized Cost` = !!cost_col / !!length_col
  ) %>%
  pivot_longer(
    cols = all_of(properties_to_pivot), 
    names_to = "Raw_Property",
    values_to = "Value"
  ) %>%
  mutate(Property = property_labels[Raw_Property]) %>%
  mutate(Property = factor(Property, levels = property_labels)) %>%
  filter(!is.na(Value)) %>%
  filter(is.finite(Value)) 

ggplot(df_long, aes(x = final_type, y = Value, fill = final_type)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  stat_compare_means(
    method = "wilcox.test",          
    label = "p.format",              
    label.x.npc = "center",          
    vjust = -0.5                     
  ) +
  scale_fill_manual(values = my_colors) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  facet_wrap(
    ~ Property, 
    scales = "free_y", 
    strip.position = "left"
  ) +
  labs(x = NULL, y = NULL, fill = "Gene Type") +
  theme_minimal() +
  theme(
    strip.placement = "outside", 
    strip.text.y.left = element_text(size = 12, color = "black", face = "bold", angle = 90),
    strip.background = element_blank(),
    axis.title = element_text(size = 18, color = "gray30"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.text.x = element_blank(),  
    axis.ticks.x = element_blank(),
    legend.position = "bottom",  
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 18)
  )

###just  for formatting
df_length <- df %>%
  filter(updated_final_type %in% c("denovo", "seq_div")) %>%
  mutate(length_x = as.numeric(length_x)) %>%
  filter(!is.na(length_x) & is.finite(length_x))

ggplot(df_length, aes(x = updated_final_type, y = length_x, fill = updated_final_type)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  stat_compare_means(
    method = "wilcox.test",          
    label = "p.format",              
    label.x.npc = "center",          
    vjust = -0.5                     
  ) +
  scale_fill_manual(values = c("denovo" = "#FAD510FF", "seq_div" = "firebrick")) +
  scale_y_log10(expand = expansion(mult = c(0.05, 0.15))) + 
  labs(
    x = NULL, 
    y = "Length (aa) (log scale)", 
    fill = "Gene Type"
  ) +
  theme_minimal() +
  theme(
    axis.title.y = element_text(size = 18, color = "black", face = 'bold'),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.text.x = element_blank(),  
    axis.ticks.x = element_blank(),
    legend.position = "bottom",  
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 18)
  )


tree <- read.tree('data/trees/332_noded.nwk')
branch_data <- data.frame(
  Ancestor_Node = tree$edge[, 1],
  Descendant_Node = tree$edge[, 2],
  Branch_Length = tree$edge.length
)

all_node_names <- c(tree$tip.label, tree$node.label)
name_mapping <- data.frame(
  Descendant_Node = 1:length(all_node_names),
  Node_Name = all_node_names
)

final_branch_data <- branch_data %>%
  left_join(name_mapping, by = "Descendant_Node")

df_with_lengths <- df %>%
  left_join(final_branch_data, by = c("recalculated_origin_forty" = "Node_Name")) %>%
  filter(updated_final_type %in% c("denovo", "seq_div"))

node_counts <- df_with_lengths %>%
  group_by(recalculated_origin_forty, Branch_Length) %>%
  count(updated_final_type) %>%
  pivot_wider(names_from = updated_final_type, values_from = n, values_fill = 0) %>%
  ungroup()


head(node_counts)
partial_cor <- pcor.test(
  x = node_counts$denovo,
  y = node_counts$seq_div,
  z = node_counts$Branch_Length,
  method = "spearman" 
)

print(partial_cor)
cor_estimate <- round(partial_cor$estimate, 2)
cor_pvalue <- signif(partial_cor$p.value, 3)
cor_label <- paste0("Partial Spearman's rho: ", cor_estimate, "\np-value: ", cor_pvalue)

ggplot(node_counts, aes(x = seq_div, y = denovo)) +
  geom_point(aes(size = Branch_Length), alpha = 0.6, color = "#1f77b4") +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +
  scale_size_continuous(range = c(2, 10)) +
  annotate(
    "text", 
    x = Inf, 
    y = Inf, 
    label = cor_label, 
    hjust = 1.1,
    vjust = 1.5,
    size = 6, 
    color = "black", 
    fontface = "bold"
  ) +
  labs(
    x = "Number of Divergent Events",
    y = "Number of De Novo Events",
    size = "Branch Length"
  ) +
  theme_minimal(base_size = 18) + 
  theme(
    axis.title.x = element_text(face = "bold", size = 22, margin = margin(t = 15)),
    axis.title.y = element_text(face = "bold", size = 22, margin = margin(r = 15)),
    axis.text = element_text(size = 16, color = "black"),
    legend.title = element_text(face = "bold", size = 18),
    legend.text = element_text(size = 16),
    legend.position = "right", 
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey85")
  )

node_family_counts <- df_with_lengths %>%
  group_by(recalculated_origin_forty, Branch_Length, updated_final_type) %>%
  summarise(n = n_distinct(merged_families), .groups = "drop") %>%
  pivot_wider(names_from = updated_final_type, values_from = n, values_fill = 0)

ggplot(node_family_counts, aes(x = seq_div, y = denovo, size = Branch_Length)) +
  geom_point(alpha = 0.5, color = "black", fill = "#1f77b4", shape = 21, stroke = 1) +
  scale_size_continuous(range = c(2, 10)) +
  labs(
    x = "Number of Divergent Gene Families",
    y = "Number of De Novo Gene Families",
    size = "Branch Length"
  ) +
  theme_minimal(base_size = 18) + 
  theme(
    axis.title.x = element_text(face = "bold", size = 22, margin = margin(t = 15)),
    axis.title.y = element_text(face = "bold", size = 22, margin = margin(r = 15)),
    axis.text = element_text(size = 16, color = "black"),
    legend.title = element_text(face = "bold", size = 18),
    legend.text = element_text(size = 16),
    legend.position = "right", 
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey85")
  )


standard_cor <- cor.test(
  x = node_family_counts$denovo,
  y = node_family_counts$seq_div,
  method = "spearman",
  exact = FALSE
)

print(standard_cor)
cor_label <- sprintf("r = %.3f, p = %.2e", standard_cor$estimate, standard_cor$p.value)

ggplot(node_family_counts, aes(x = seq_div, y = denovo)) +
  geom_point(alpha = 0.5, color = "black", fill = "#1f77b4", shape = 21, size = 4, stroke = 1) +
  geom_smooth(method = "lm", color = "red", linetype = "dashed", fill = "grey80") +
  annotate("text", x = Inf, y = Inf, label = cor_label, hjust = 2, vjust = 1.5, size = 7, fontface = "bold", color = "black") +
  labs(
    x = "Number of Divergent Gene Families",
    y = "Number of De Novo Gene Families"
  ) +
  theme_minimal(base_size = 18) + 
  theme(
    axis.title.x = element_text(face = "bold", size = 22, margin = margin(t = 15)),
    axis.title.y = element_text(face = "bold", size = 22, margin = margin(r = 15)),
    axis.text = element_text(size = 16, color = "black"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey85")
  )

###motifs
genome_lengths <- read.csv('data/assemblies/genome_lengths.tsv', header =F, sep ='\t')
genome_lengths <- subset(genome_lengths, select = -c(V1))
colnames(genome_lengths) <- c('species', 'cds_length')

non_coding_overlapping_1st <- read.csv('data/assemblies/1st_motif/all_species_summary_positive_strand.txt', header = F,sep ='\t')
non_coding_overlapping_2nd <- read.csv('data/assemblies/2nd_motif/all_species_summary_positive_strand.txt', header = F,sep ='\t')
non_coding_overlapping_3rd <- read.csv('data/assemblies/3rd_motif/all_species_summary_positive_strand.txt', header = F,sep ='\t')
colnames(non_coding_overlapping_1st) <- c('species', 'non_coding_space_1st')
colnames(non_coding_overlapping_2nd) <- c('species', 'non_coding_space_2nd')
colnames(non_coding_overlapping_3rd) <- c('species', 'non_coding_space_3rd')

genome_sizes <- read.csv('data/assemblies/genome_sizes.txt', header = F, sep ='\t')
colnames(genome_sizes) <- c('species','genome_size')

df_list <- list(
  genome_lengths, 
  non_coding_overlapping_1st, 
  non_coding_overlapping_2nd, 
  non_coding_overlapping_3rd, 
  genome_sizes
)

non_coding_overlapping <- Reduce(function(x, y) merge(x, y, by = "species", all = TRUE), df_list)
non_coding_overlapping$total_non_coding <- non_coding_overlapping$genome_size - non_coding_overlapping$cds_length
non_coding_overlapping$norm_space_1st <- non_coding_overlapping$non_coding_space_1st / non_coding_overlapping$total_non_coding
non_coding_overlapping$norm_space_2nd <- non_coding_overlapping$non_coding_space_2nd / non_coding_overlapping$total_non_coding
non_coding_overlapping$norm_space_3rd <- non_coding_overlapping$non_coding_space_3rd / non_coding_overlapping$total_non_coding

space_sizes <- genome_sizes %>%
  mutate(species = gsub(" ", "_", species)) %>% 
  left_join(
    genome_lengths %>% mutate(species = gsub(" ", "_", species)), 
    by = "species"
  ) %>%
  mutate(
    size_coding = cds_length,
    size_non_coding = genome_size - cds_length
  ) %>%
  dplyr::select(species, size_coding, size_non_coding)

plot_df <- non_coding_overlapping %>%
  pivot_longer(
    cols = starts_with("norm_space"), 
    names_to = "Motif", 
    values_to = "Normalized_Coverage"
  ) %>%
  mutate(Motif = case_when(
    Motif == "norm_space_1st" ~ "1st Motif",
    Motif == "norm_space_2nd" ~ "2nd Motif",
    Motif == "norm_space_3rd" ~ "3rd Motif"
  ))

coding_overlapping_1st <- read.csv('data/assemblies/1st_motif/all_species_coding_summary_positive_strand.txt', header = F,sep ='\t')
coding_overlapping_2nd <- read.csv('data/assemblies/2nd_motif/all_species_coding_summary_positive_strand.txt', header = F,sep ='\t')
coding_overlapping_3rd <- read.csv('data/assemblies/3rd_motif/all_species_coding_summary_positive_strand.txt', header = F,sep ='\t')
colnames(coding_overlapping_1st) <- c('species', 'coding_space_1st')
colnames(coding_overlapping_2nd) <- c('species', 'coding_space_2nd')
colnames(coding_overlapping_3rd) <- c('species', 'coding_space_3rd')

df_list <- list(
  non_coding_overlapping, 
  coding_overlapping_1st, 
  coding_overlapping_2nd, 
  coding_overlapping_3rd)

motif_overlapping <- Reduce(function(x, y) merge(x, y, by = "species", all = TRUE), df_list)
motif_overlapping$norm_coding_space_1st <- motif_overlapping$coding_space_1st / motif_overlapping$cds_length
motif_overlapping$norm_coding_space_2nd <- motif_overlapping$coding_space_2nd / motif_overlapping$cds_length
motif_overlapping$norm_coding_space_3rd <- motif_overlapping$coding_space_3rd / motif_overlapping$cds_length

motif_overlapping <- motif_overlapping %>%
  mutate(genus = word(species, 1, sep = "[ _]"))

plot_comparison_df <- motif_overlapping %>%
  dplyr::select(species, starts_with("norm_space"), starts_with("norm_coding_space")) %>%
  pivot_longer(
    cols = -species,
    names_to = "Original_Column",
    values_to = "Normalized_Coverage"
  ) %>%
  mutate(
    Region = ifelse(grepl("norm_coding", Original_Column), "Coding Space", "Non-Coding Space"),
    Motif = case_when(
      grepl("1st", Original_Column) ~ "1st Motif",
      grepl("2nd", Original_Column) ~ "2nd Motif",
      grepl("3rd", Original_Column) ~ "3rd Motif"
    )
  ) %>%
  mutate(Normalized_Coverage = replace_na(Normalized_Coverage, 0))

normalized_diff_df <- plot_comparison_df %>%
  filter(!species %in% c("all_species_summary.txt", "non_coding_regions.gff")) %>%
  dplyr::select(species, Region, Motif, Normalized_Coverage) %>%
  pivot_wider(names_from = Region, values_from = Normalized_Coverage) %>%
  left_join(genome_sizes, by = "species") %>%
  mutate(
    Coverage_Diff = `Non-Coding Space` - `Coding Space`,
    Size_Normalized_Diff = Coverage_Diff / genome_size
  ) %>%
  drop_na(Size_Normalized_Diff) %>%
  mutate(
    Motif = case_when(
      Motif == "1st Motif" ~ "polyT",
      Motif == "2nd Motif" ~ "TATA",
      Motif == "3rd Motif" ~ "polyA",
      TRUE ~ Motif 
    ),
    Motif = factor(Motif, levels = c("polyT", "TATA", "polyA"))
  )

ggplot(normalized_diff_df, aes(x = Motif, y = Coverage_Diff, fill = Motif)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA, color = "black", linewidth = 0.8) +
  geom_jitter(width = 0.2, height = 0, alpha = 0.4, size = 2, color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 1.2) +
  scale_fill_manual(
    values = c(
      "polyT" = "cadetblue",
      "polyA" = "rosybrown3", 
      "TATA"  = "olivedrab"  
    )
  ) +
  theme_bw(base_size = 18) +
  labs(
    x = NULL, 
    y = "Relative Motif Density Difference\n(Non-Coding - Coding)"
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(face = "bold", size = 18, color = "black"),
    axis.text.y = element_text(size = 16, color = "black"),
    axis.title.y = element_text(face = "bold", size = 20, margin = margin(r = 15)),
    panel.grid.major = element_line(color = "grey85", linewidth = 0.5),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 1.2),
    axis.ticks = element_line(color = "black", linewidth = 1),
    axis.ticks.length = unit(0.25, "cm")
  )

motifs <- read.csv('data/matrices/fixed_genus_only.tsv', sep ='\t')
original_cols <- colnames(motifs)

motifs <- motifs %>%
  left_join(dplyr::select(df, code_name, updated_final_type), by = c("sequence" = "code_name")) %>%
  mutate(type = updated_final_type) %>%
  dplyr::select(-updated_final_type) %>%
  distinct(sequence, .keep_all = TRUE) %>%
  mutate(type = coalesce(type, "conserved"))

long_motifs <- motifs %>%
  pivot_longer(
    cols = c(polyA_div, polyT_div, TATA_div),
    names_to = "motif",
    values_to = "normalized_occurrence"
  ) %>%
  mutate(motif = factor(motif, 
                        levels = c("polyA_div", "polyT_div", "TATA_div"), 
                        labels = c("polyA", "polyT", "TATA")))

my_order <- c("conserved", "seq_div", "trg", "denovo")

long_motifs <- long_motifs %>%
  left_join(dplyr::select(df, code_name,  cleaned_species), by = c("sequence" = "code_name")) %>%
  filter(type %in% my_order) %>% 
  mutate(raw_motif_count = normalized_occurrence * length_x) %>%
  group_by(type, motif) %>%
  summarise(
    total_motifs = sum(raw_motif_count, na.rm = TRUE),
    total_length = sum(length_x, na.rm = TRUE),
    p = total_motifs / total_length,
    se = sqrt((p * (1 - p)) / total_length),
    pooled_percentage = p * 100,
    se_percentage = se * 100,
    .groups = 'drop'
  ) %>%
  mutate(type = factor(type, levels = my_order))

my_colors <- c(
  "conserved" = "springgreen4",
  "seq_div"   = "firebrick",
  "trg"       = "cyan3",
  "denovo"    = "#FAD510FF"
)

ggplot(long_motifs, aes(x = type, y = pooled_percentage, fill = type)) +
  geom_col(color = "black", alpha = 0.8) +
  geom_errorbar(
    aes(ymin = pooled_percentage - se_percentage,
        ymax = pooled_percentage + se_percentage),
    width = 0.2,
    linewidth = 0.8
  ) +
  scale_fill_manual(
    values = my_colors, 
    name = "Gene Type",
    labels = c(
      "conserved" = "Conserved", 
      "seq_div" = "Divergent", 
      "trg" = "TRG", 
      "denovo" = "De novo"
    )
  ) +
  facet_wrap(~ motif, scales = "free_y") +
  theme_bw(base_size = 16) +
  labs(
    title = NULL,
    x = NULL, 
    y = "Motifs per 100 Base Pairs (%)"
  ) +
  theme(
    legend.position = c(0.88, 0.75), 
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 11),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.background = element_rect(fill = "gray90", color = "black"),
    strip.text = element_text(face = "bold", size = 14)
  )

motifs_submatrix <- motifs %>%
  inner_join(df %>% dplyr::select(code_name, updated_final_type), by = c("sequence" = "code_name")) %>%
  mutate(type = updated_final_type) %>%
  dplyr::select(all_of(original_cols))

filtered_motifs <- motifs_submatrix %>%
  filter(genus %in% yeast_genera) %>% 
  filter(type %in% my_order) %>% 
  mutate(type = factor(type, levels = my_order))

ggplot(filtered_motifs, aes(x = genus, y = polyA_div, fill = type)) + ##change the y to polyT_div and polyA_div
  stat_summary(fun = "mean", geom = "bar", position = position_dodge(width = 0.8), color = "black") +
  stat_summary(fun.data = "mean_se", geom = "errorbar", position = position_dodge(width = 0.8), width = 0.2) +
  scale_fill_manual(
    values = my_colors, 
    name = "Gene Type",
    labels = c(
      "conserved" = "Conserved", 
      "seq_div"   = "Divergent", 
      "trg"       = "TRG", 
      "denovo"    = "De novo"
    )
  ) +
  theme_bw(base_size = 20) + 
  labs(
    x = "Genus", 
    y = "Normalized polyA Motif Occurrence" ###the titles also
  ) +
  theme(
    axis.title = element_text(size = 22, face = "bold"),            
    axis.text.x = element_text(angle = 45, hjust = 1, size = 18),  
    axis.text.y = element_text(size = 18),                         
    legend.title = element_text(size = 20, face = "bold"),         
    legend.text = element_text(size = 18),                         
    panel.grid.major.x = element_blank()
  )


###HGT 
###tree

hgt_data <- df[df$updated_final_type == 'HGT', ]
hgt_families <- data.frame(table(hgt_data$merged_families))

genes_per_species <- genes_per_species[match(tree$tip.label, genes_per_species$V2), ]
tree_data <- ggtree(tree)$data

tip_metadata_simple <- df %>%
  group_by(species, updated_final_type) %>%
  summarise(count = n(), .groups = 'drop') %>%
  pivot_wider(names_from = updated_final_type, values_from = count, values_fill = 0) %>%
  dplyr::rename(label = species)

inner_metadata <- df %>%
  group_by(recalculated_origin_forty, final_age_, updated_final_type) %>%
  summarise(count = n(), .groups = 'drop') %>%
  pivot_wider(names_from = updated_final_type, values_from = count, values_fill = 0) %>%
  dplyr::rename(label = recalculated_origin_forty)

tree_data <- tree_data %>%
  left_join(tip_metadata_simple, by = "label") %>%
  left_join(inner_metadata, by = "label", suffix = c("_tip", "_internal")) %>%
  mutate(hgt_value = coalesce(HGT_tip, HGT_internal))

tree_data$hgt_all_genes <- tree_data$hgt_value / genes_per_species$V1

clade_metadata <- read.csv('data/metadata/343taxa_speicies-name_clade-name_color-code.txt',sep ='\t')
clade_metadata <- clade_metadata[,c("old_speceis_names", "Major.clade")]
ordered_clades <- unique(clade_metadata$Major.clade)
num_clades <- length(ordered_clades)
clade_colors_vector <- colorRampPalette(brewer.pal(12, "Set3"))(num_clades)
clade_colors <- setNames(clade_colors_vector, ordered_clades)

tip_metadata_bar <- df %>%
  group_by(species, genus, updated_final_type) %>%
  summarise(count = n(), .groups = 'drop') %>%
  pivot_wider(names_from = updated_final_type, values_from = count, values_fill = 0)

tip_metadata_bar <- tip_metadata_bar %>%
  left_join(clade_metadata, by = c("species" = "old_speceis_names")) 

tip_metadata_bar <- tip_metadata_bar[match(tree$tip.label, tip_metadata_bar$species), ]
tip_metadata_bar$hgt_all_genes <- tip_metadata_bar$HGT / genes_per_species$V1
tip_metadata_bar <- tip_metadata_bar %>%
  mutate(Major.clade = factor(Major.clade, levels = ordered_clades))

p <- ggtree(tree, layout = "circular")

p <- p +
  geom_fruit(
    data = tip_metadata_bar,
    geom = geom_bar,
    mapping = aes(y = species, x = hgt_all_genes, fill = Major.clade), 
    stat = "identity",
    orientation = "y",
    width = 0.8, 
    offset = 0.05, 
    size = 0.4,
    grid.params = list(
      vline = TRUE, 
      color = "grey85",  
      size = 0.15,       
      linetype = "dashed"
    ),
    axis.params = list(
      axis = "x", 
      text.size = 2,     
      text.angle = 0, 
      hjust = 0.5,
      color = "grey50"   
    )
  ) +
  scale_fill_manual(
    values = clade_colors, 
    name = "Saccharomycotina Order" 
  )

p <- p %<+% tree_data +
  new_scale_color() + 
  geom_tree(aes(color = HGT_internal), size = 1) + 
  scale_color_gradient(
    low = "grey", 
    high = "black", 
    name = "HGT origination events",
    guide = guide_colorbar(position = "left") 
  ) +
  theme(
    legend.position = "right", 
    legend.justification = "top" 
  )


ggsave(
  "figures/figure_9a.svg", 
  plot = p,
  width = 8, 
  height = 6,
  device = svglite)


long_motifs_hgt <- long_motifs %>%
  left_join(dplyr::select(df, code_name, cleaned_species), by = c("sequence" = "code_name")) %>%
  filter(type %in% my_order) %>%
  mutate(raw_motif_count = normalized_occurrence * length_x) %>%
  group_by(type, motif) %>%
  summarise(
    total_motifs = sum(raw_motif_count, na.rm = TRUE),
    total_length = sum(length_x, na.rm = TRUE),
    p = total_motifs / total_length,
    se = sqrt((p * (1 - p)) / total_length),
    .groups = 'drop'
  ) %>%
  mutate(type = factor(type, levels = my_order))

my_colors_hgt <- c(
  "denovo" = "#FAD510FF",
  "HGT" = "sienna4"
)

ggplot(long_motifs_hgt, aes(x = type, y = p, fill = type)) +
  geom_col(color = "black", alpha = 0.8) +
  geom_errorbar(
    aes(ymin = p - se,
        ymax = p + se),
    width = 0.2,
    linewidth = 0.8
  ) +
  scale_fill_manual(
    values = my_colors_hgt,
    name = "Gene Type",
    labels = c("denovo" = "denovo", "HGT" = "HGT")
  ) +
  guides(fill = guide_legend(reverse = TRUE)) +
  facet_wrap(~ motif, scales = "free_y") +
  theme_bw(base_size = 16) +
  labs(
    title = NULL,
    x = NULL,
    y = NULL
  ) +
  theme(
    legend.position = c(0.85, 0.75),
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 11),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(face = "bold", size = 14)
  )

###supplemenary comparison of properties between hgt and de novo

properties_to_plot <- c("length_x", "gc", "isoelectric", "cost", "propensity", "diso_pct" ,"hydropathicity", "nuc_coverage")
length_col <- rlang::sym("length_x")
cost_col   <- rlang::sym("cost")
tm_col     <- rlang::sym("propensity")
properties_to_pivot <- c("length_x", "gc", "isoelectric", "propensity", "diso_pct", "hydropathicity", "nuc_coverage", "Normalized Cost")
my_colors <- c("denovo" = "#FAD510FF", "HGT" = "sienna4")

property_labels <- c(
  "Length (aa) (log scale)", "GC Content (%)", "Isoelectric Point (pH)", 
  "TM propensity", "Disorder Pct (%)", 
  "Hydropathicity (GRAVY)", "Nucleosomal coverage", "Normalized Cost"
)
names(property_labels) <- properties_to_pivot

df_long <- df %>%
  filter(updated_final_type %in% c("denovo", "HGT")) %>%
  mutate(across(all_of(properties_to_plot), as.numeric)) %>%
  mutate(
    `Normalized Cost` = !!cost_col / !!length_col
  ) %>%
  pivot_longer(
    cols = all_of(properties_to_pivot), 
    names_to = "Raw_Property",
    values_to = "Value"
  ) %>%
  mutate(Property = property_labels[Raw_Property]) %>%
  mutate(Property = factor(Property, levels = property_labels)) %>%
  filter(!is.na(Value)) %>%
  filter(is.finite(Value)) 

ggplot(df_long, aes(x = final_type, y = Value, fill = final_type)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  stat_compare_means(
    method = "wilcox.test",          
    label = "p.format",              
    label.x.npc = "center",          
    vjust = -0.5                     
  ) +
  scale_fill_manual(values = my_colors) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  facet_wrap(
    ~ Property, 
    scales = "free_y", 
    strip.position = "left"
  ) +
  labs(x = NULL, y = NULL, fill = "Gene Type") +
  theme_minimal() +
  theme(
    strip.placement = "outside", 
    strip.text.y.left = element_text(size = 12, color = "black", face = "bold", angle = 90),
    strip.background = element_blank(),
    axis.title = element_text(size = 18, color = "gray30"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.text.x = element_blank(),  
    axis.ticks.x = element_blank(),
    legend.position = "bottom",  
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 18)
  )

###just  for formatting
df_length <- df %>%
  filter(updated_final_type %in% c("denovo", "HGT")) %>%
  mutate(length_x = as.numeric(length_x)) %>%
  filter(!is.na(length_x) & is.finite(length_x))

ggplot(df_length, aes(x = updated_final_type, y = length_x, fill = updated_final_type)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  stat_compare_means(
    method = "wilcox.test",          
    label = "p.format",              
    label.x.npc = "center",          
    vjust = -0.5                     
  ) +
  scale_fill_manual(values = c("denovo" = "#FAD510FF", "HGT" = "sienna4")) +
  scale_y_log10(expand = expansion(mult = c(0.05, 0.15))) + 
  labs(
    x = NULL, 
    y = "Length (aa) (log scale)", 
    fill = "Gene Type"
  ) +
  theme_minimal() +
  theme(
    axis.title.y = element_text(size = 18, color = "black", face = 'bold'),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.text.x = element_blank(),  
    axis.ticks.x = element_blank(),
    legend.position = "bottom",  
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 18)
  )