# **Phylogeography branch**

## **This branch is for the Phylogeography figure of the discussion. The input data was sourced from the three Australian published studies of *Uncinaria stenocephala* benzimidazole amplicon metabarcoding. The first two studies represented Chapters 3 and 5 of the thesis. The third study is based on a DESS experimental control study of *Ancylostoma caninum* and *U. stenocephala.***

### **The coding (part 1):**

The input files are generated from the three studies' LabArchives repositories. The files were refined to eliminate all excess knowledge apart from the ASVs (Sequences), samples as columns and the status of resistance (Status). Only *U. stenocephala* ASVs were retained. The three files are merged using the following libraries and codes:

```         
library(dplyr)
library(tidyverse)
library(phytools)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(tidyr)
library(ape)
library(RColorBrewer)
library(ggplot2)
library(scatterpie)
library(ggtree
library(patchwork)
```

```         
DESS<- read.csv("C:/Users/tsto3616/thesis-drafts/DESS_hooks_BZ167.csv")

head(DESS)

Unamb<- read.csv("C:/Users/tsto3616/thesis-drafts/unambiguous_hooks_BZ167.csv")

head(Unamb)

widespread<- read.csv("C:/Users/tsto3616/thesis-drafts/widespread_hooks_BZ167.csv")

widespread

merged1<- DESS %>% full_join(Unamb, by= c("Sequences", "Status"))

merged1

merged_Uste <- merged1 %>% full_join(widespread, by=c("Sequences", "Status"))

merged_Uste

# Replace all NA values in a dataframe with 0
merged_Uste[is.na(merged_Uste)] <- 0

merged_Uste
```

The fasta file of all ASVs and a metadata csv was generated through the coding:

```         
# assign ASV labels
merged_Uste$ASV_ID <- paste0("ASV", seq_len(nrow(merged_Uste)))

# Convert sequences to DNAStringSet
asv_set <- DNAStringSet(merged_Uste$Sequences)
names(asv_set) <- merged_Uste$ASV_ID   # use ASV IDs as headers

# Write to FASTA file
writeXStringSet(asv_set, filepath = "C:/Users/tsto3616/thesis/Uste_ASVs.fasta")

write.csv(merged_Uste, "C:/Users/tsto3616/thesis/Uste_ASVs_metadata.csv")
```

### **Generation of the phylogenetic tree:**

The phylogenetic tree was generated in MEGA12 using default ClustrW and model selection for the ideal maximum likelihood model, which was reported to be K2 with invariable sites was reported to be best. The phylogenetic tree was cleaned in MEGA12 and imported into R as a nwk file to be analysed for the phylogeographical model.

### **The coding (part2):**

The metadata was manually manipulated in excel and had the location (with rough GPS coordinates) added. First data manipulation and loading the files:

```         
tree <- read.tree("Uste_ASVs_phylogeny.nwk")
loc <- read.csv("Uste_ASVs_metadata.csv")

tree$tip.label
loc

keep<- c("Samples", "Long", "Lat", "Status")

location<- loc %>% pivot_longer(cols = -all_of(keep), names_to = "ASVs", values_to = "Presence")

coords <- as.matrix(location[, c("Lat","Long")])
rownames(coords) <- location$ASVs
rownames(coords) <- location$Presence
rownames(coords) <- location$Samples
```

Then go and reformat the tree so that the tree can be geographically and phylogenetically clustered based on the longitude, latitude and the genetic distance using R. Adjusted posthoc for the number of clusters to make the most biological sense:

```         
# Compute pairwise patristic distances between tips
D <- cophenetic(tree)

# Cluster tips using hclust on the distance matrix
hc <- hclust(as.dist(D))

# Cut into N groups
groups <- cutree(hc, k = 7)

# Suppose groups is your clade assignment vector: names = tip labels, values = clade IDs
N <- length(unique(groups))
palette <- brewer.pal(N, "Set3")

# Map each tip to its clade colour
tip_cols <- setNames(palette[groups], names(groups))

# Plot tree with coloured tips
plot(tree, tip.color = tip_cols, cex = 0.7)

location <- location %>%
  mutate(clade = groups[ASVs])

# Aggregate coordinates by clade (centroid example)
clade_coords <- location %>%
  filter(Presence > 0) %>%
  group_by(clade) %>%
  summarise(Long = mean(Long), Lat = mean(Lat))
```

Now to create the pie-chart graph of East Coast of Australia and most of New Zealand, this code enables the demonstration of pie charts for each representative cluster of the phylogeny:

```         
anz <- ne_countries(scale = "medium",
                    country = c("Australia","New Zealand"),
                    returnclass = "sf")
                    
pie_data <- location %>%
  filter(Presence > 0) %>%
  group_by(Samples, Long, Lat, clade) %>%
  summarise(count = sum(Presence), .groups = "drop") %>%
  pivot_wider(names_from = clade,
              values_from = count,
              values_fill = 0,
              names_prefix = "Clade_")

pie_data <- pie_data %>%
  mutate(Long = jitter(Long, amount = 1),
         Lat  = jitter(Lat, amount = 1))

head(pie_data)

p_map <- ggplot(data = anz) +
  geom_sf(fill = "white", color = "black") +
  geom_scatterpie(
    data = pie_data,
    aes(x = Long, y = Lat, colour = Status),   # map outline color to Status
    cols = names(pie_data)[grepl("Clade_", names(pie_data))],
    pie_scale = 0.9, lwd = 0.2
  ) +
  coord_sf(xlim = c(140, 180), ylim = c(-45, -30)) +
  scale_fill_brewer(palette = "Set3") +
  scale_color_manual(values = c("Resistant" = "red", "Sensitive" = "black")) +
  theme_minimal() +
  theme(legend.position = "none")
```

Now we add our complete tree with tip labels as shapes to denote the F167Y ASVs (triangles):

```         
highlight_asvs <- c("ASV2", "ASV13", "ASV19", "ASV25", "ASV44")

tip_meta <- data.frame(
  ASVs = tree$tip.label,
  clade = groups[tree$tip.label],
  highlight_flag = tree$tip.label %in% highlight_asvs
)


# tip_meta must include: ASVs (matching tree$tip.label), clade, resistant_flag
p_tree <- ggtree(tree, layout = "rectangular") %<+% tip_meta +
  geom_tiplab(
    aes(color = factor(clade)),
    angle = 90,      # horizontal labels
    hjust = 0,      # left-justify relative to tip
    vjust = 1,      # anchor for downward orientation
    size = 2
  ) + coord_cartesian(clip = "off") +
  geom_tippoint(aes(shape = highlight_flag, color = factor(clade)), size = 2) +
  scale_color_manual(values = palette) +
  theme_tree() +
  theme(legend.position = "none") +
  coord_flip() +     # make the tree vertical
  scale_y_reverse() +   # flip vertically so tips face down
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.15)))

p_tree
```

And now we create a stacked plot with the phylogeny below the pie charts:

```         
combined <- p_map / p_tree

combined

ggsave("C:/Users/tsto3616/thesis-drafts/phylogeography.svg", plot = combined, width = 8, height = 6, dpi = 300)
```
