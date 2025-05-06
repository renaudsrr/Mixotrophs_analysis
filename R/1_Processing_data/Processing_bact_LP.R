# Script by SERRE Renaud
# Internship UQAM-GRIL under supervision of BEISNER Beatrix

library(dplyr)
library(openxlsx)
library(tidyverse)
library(gridExtra)
library(corrplot)
library(glue)
library(vegan)
library(grid) 
library(patchwork)
library(forcats)
library(cowplot)
library(viridis)
library(ggsci)

rm(list=ls()) 
graphics.off() 
cat("\014")

# Import data -------------------------------------------------------------

bact = read.csv("/Users/renaudsrr/Desktop/STAGE_MTL/MODELISATION/DATA/2_LP_NLA/Main_asvtable_taxonomy_rmcyano_rarefied.csv", check.names = FALSE) # avec methode de rarefaction relou
data_phyto = read.csv("/Users/renaudsrr/Desktop/STAGE_MTL/MODELISATION/DATA/2_LP_NLA/HarmonizedPhyto_LakePulseNLA2017_19092022.csv")
info_genus = read.csv("/Users/renaudsrr/Desktop/STAGE_MTL/MODELISATION/DATA/Mixotroph_strat/info_genus.csv")
data_python = read.csv("/Users/renaudsrr/Desktop/STAGE_MTL/MODELISATION/DATA/0_new_csv/data_python/df_LP_NLA.csv")

lakes_LP = data_phyto$lake_id[data_phyto$Survey=="Lake Pulse"]
data_python = data_python[data_python$lake_id %in% lakes_LP,] 
data_python = data_python %>%
  rename(rich_genus_phyto = rich_genus,
         H_phyto = shannon,
         one_minus_D_phyto = op_simpson,
         J_phyto = eveness_piel)

bact = bact[,-1]
bact = bact %>%
  # select(-c(asv_code,sequence,kingdom,clade,lineage,tribe)) %>%
  select(-c(asv_code,sequence,kingdom)) %>%
  rename(Phylum  = phylum,
         Class  = class, 
         Order  = order) # %>%
  # mutate(
  #   Phylum = ifelse(is.na(Phylum), "Unknown", Phylum),# pour eviter biais annotation
  #   Class = ifelse(is.na(Class), "Unknown", Class),
  #   Order = ifelse(is.na(Order), "Unknown", Order))

bact_long = bact %>%
  pivot_longer(
    cols = -c(Phylum, Class, Order,clade,lineage,tribe),
    names_to = "lake_id",
    values_to = "abundance"
  ) %>% select(lake_id,Phylum,Class,Order,clade,lineage,tribe,abundance)
dim(bact_long)[1]
bact_long = na.omit(bact_long)
dim(bact_long)[1]


# Phylum ------------------------------------------------------------------
summary_data_bact_phylum = bact_long %>%
  filter(abundance > 0) %>%
  group_by(lake_id, Phylum) %>%
  summarize(
    ab_phylum = sum(abundance),
    .groups = "drop"
  ) %>%
  group_by(lake_id) %>%
  mutate(ab_tot_phylum_bact = sum(ab_phylum)) %>%
  summarize(
    rich_phylum_bact = n_distinct(Phylum),
    ab_tot_phylum_bact = unique(ab_tot_phylum_bact),
    H_phylum_bact = -sum((ab_phylum / ab_tot_phylum_bact) * log(ab_phylum / ab_tot_phylum_bact)),
    one_minus_D_phylum_bact = 1 - sum((ab_phylum / ab_tot_phylum_bact)^2),
    .groups = "drop"
  ) %>%
  mutate(J_phylum_bact = H_phylum_bact / log(rich_phylum_bact))

data_python_bact_phylum = data_python %>%
  left_join(summary_data_bact_phylum,by="lake_id")


hist(data_python_bact_phylum$rich_phylum_bact)
hist(data_python_bact_phylum$H_phylum_bact)
hist(data_python_bact_phylum$one_minus_D_phylum_bact)
hist(data_python_bact_phylum$J_phylum_bact)

plot(x=data_python_bact_phylum$prev_Mixo,y=data_python_bact_phylum$rich_phylum_bact)
plot(x=data_python_bact_phylum$prev_Mixo,y=data_python_bact_phylum$H_phylum_bact)
plot(x=data_python_bact_phylum$prev_Mixo,y=data_python_bact_phylum$one_minus_D_phylum_bact)
plot(x=data_python_bact_phylum$prev_Mixo,y=data_python_bact_phylum$J_phylum_bact)

# Class -------------------------------------------------------------------

summary_data_bact_class = bact_long %>%
  filter(abundance > 0) %>%
  group_by(lake_id, Class) %>%
  summarize(
    ab_class = sum(abundance),
    .groups = "drop"
  ) %>%
  group_by(lake_id) %>%
  mutate(ab_tot_class_bact = sum(ab_class)) %>%
  summarize(
    rich_class_bact = n_distinct(Class),
    ab_tot_class_bact = unique(ab_tot_class_bact),
    H_class_bact = -sum((ab_class / ab_tot_class_bact) * log(ab_class / ab_tot_class_bact)),
    one_minus_D_class_bact = 1 - sum((ab_class / ab_tot_class_bact)^2),
    .groups = "drop"
  ) %>%
  mutate(J_class_bact = H_class_bact / log(rich_class_bact))

data_python_bact_class = data_python %>%
  left_join(summary_data_bact_class,by="lake_id")


hist(data_python_bact_class$rich_class_bact)
hist(data_python_bact_class$H_class_bact)
hist(data_python_bact_class$one_minus_D_class_bact)
hist(data_python_bact_class$J_class_bact)

plot(x=data_python_bact_class$prev_Mixo,y=data_python_bact_class$rich_class_bact)
plot(x=data_python_bact_class$prev_Mixo,y=data_python_bact_class$H_class_bact)
plot(x=data_python_bact_class$prev_Mixo,y=data_python_bact_class$one_minus_D_class_bact)
plot(x=data_python_bact_class$prev_Mixo,y=data_python_bact_class$J_class_bact)

# test -------------------------------------------------------------------

numeric_vars = data_python_bact_class %>%
  select(where(is.numeric)) %>%
  select(-prev_Mixo) %>%
  colnames()

for (var in numeric_vars) {
  gg = ggplot(data_python_bact_class, aes_string(x = "prev_Mixo", y = var)) +
    geom_point(alpha = 0.5) +
    geom_smooth(method = "gam", se = T, color = "red") + # ajuste un gam sur chaque var
    labs(x = "prev_Mixo",
         y = var) +
    theme_minimal()

  print(gg)
}


# ANALYSE TAXONOMIQUE -----------------------------------------------------

data_taxon_bact = bact_long %>%
  left_join(summary_data_bact_class,by="lake_id") %>%
  left_join(data_python[,c("lake_id","prev_Mixo")],by="lake_id") %>%
  select(lake_id,abundance,ab_tot_class_bact,Phylum,Class,Order,everything())

gg_H = ggplot(data_python_bact_class, aes(x = prev_Mixo, y = H_class_bact)) +
  geom_point(alpha = 0.3, color = "black") +
  theme_minimal()
gg_1_D = ggplot(data_python_bact_class, aes(x = prev_Mixo, y = one_minus_D_class_bact)) +
  geom_point(alpha = 0.3, color = "black") +
  theme_minimal()
gg_J = ggplot(data_python_bact_class, aes(x = prev_Mixo, y = J_class_bact)) +
  geom_point(alpha = 0.3, color = "black") +
  theme_minimal()

all_type_taxon = c("Phylum","Class","Order","clade","lineage","tribe")
palette_taxon = c(
  "#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e",
  "#e6ab02", "#a6761d", "#666666", "#8dd3c7", "#fb8072",  
  "#80b1d3", "#fdb462", "#b3de69", "#fccde5", "#bc80bd",  
  "#cab2d6", "#ffffb3", "#b15928",                        
  "black", "white" # Other et NA                                
)

bins = 10 # nb separation % prevalence mixo
intervals = 100 %/% bins  
data_taxon_bact = data_taxon_bact %>%
  mutate(bins = prev_Mixo %/% intervals + 1)

top_n_taxa = 15  # nombre de modalités à garder
i = 1


for (taxon in all_type_taxon) {
  
  data_taxon_grouped = data_taxon_bact %>%
    mutate(taxon_grouped = fct_lump_n(.data[[taxon]], n = top_n_taxa, other_level = "Other"))
  
  table_taxon = data_taxon_grouped %>%
    group_by(bin = bins, taxon_grouped) %>%
    summarise(
      freq = n(),
      ab_tot = sum(abundance, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      percent_freq = 100 * freq / ave(freq, bin, FUN = sum),
      percent_ab = 100 * ab_tot / ave(ab_tot, bin, FUN = sum)
    )
  
  # write.csv(table_taxon,paste0("/Users/renaudsrr/Desktop/STAGE_MTL/Scripts/Data/new_csv/table_taxons/",i,"_Table_taxons_",taxon,".csv"),row.names = F)
  
  richesse_taxon = data_taxon_grouped %>%
    group_by(bin = bins) %>%
    summarise(richesse = n_distinct(.data[[taxon]])) 
  
  gg_rich_taxon = ggplot(richesse_taxon, aes(x = bin, y = richesse)) +
    geom_line(color = "black") +
    theme_minimal() +
    labs(y = paste0("Richesse en ", taxon), x = "Bin de % Mixo")
  
  # recup legend
  gg_for_legend = ggplot(table_taxon, aes(x = bin, y = percent_ab, fill = taxon_grouped)) +
    geom_bar(stat = "identity") +
    labs(fill = taxon) +
    scale_fill_manual(values = palette_taxon) +
    theme_minimal(base_size = 9) +  
    theme(legend.position = "right",
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 12))
  legend_plot = cowplot::get_legend(gg_for_legend)
  
  # Freq apparition taxons
  gg_freq1 = ggplot(table_taxon, aes(x = bin, y = freq, fill = taxon_grouped)) +
    geom_bar(stat = "identity") +
    labs(y = "Fréquence apparition taxon") +
    scale_fill_manual(values = palette_taxon) +
    theme_minimal() +
    theme(legend.position = "none")
  
  # % apparition taxons
  gg_percent1 = ggplot(table_taxon, aes(x = bin, y = percent_freq, fill = taxon_grouped)) +
    geom_bar(stat = "identity") +
    labs(y = "% apparition taxon") +
    scale_fill_manual(values = palette_taxon) +
    theme_minimal() +
    theme(legend.position = "none")
  
  # Biomass cumulée taxons
  gg_freq2 = ggplot(table_taxon, aes(x = bin, y = ab_tot, fill = taxon_grouped)) +
    geom_bar(stat = "identity") +
    labs(y = "abundance cumulé") +
    scale_fill_manual(values = palette_taxon) +
    theme_minimal() +
    theme(legend.position = "none")
  
  # % Biomass taxon
  gg_percent2 = ggplot(table_taxon, aes(x = bin, y = percent_ab, fill = taxon_grouped)) +
    geom_bar(stat = "identity") +
    labs(y = "% abundance par taxon") +
    scale_fill_manual(values = palette_taxon) +
    theme_minimal() +
    theme(legend.position = "none")
  
  # 8 plots sans légende
  plots_no_legend = plot_grid(
    gg_freq2, gg_freq1, gg_H, gg_rich_taxon,
    gg_percent2, gg_percent1, gg_1_D, gg_J,
    ncol = 4
  )
  
  # Ajout legend
  final_plot = plot_grid(legend_plot, plots_no_legend, ncol = 2, rel_widths = c(0.1, 1.1))  
  filename = paste0("/Users/renaudsrr/Desktop/STAGE_MTL/Scripts/figures/panels/panels_bact/",i,"_Panel_", taxon, ".png")
  ggsave(filename, plot = final_plot, width = 20, height = 12, dpi = 300, bg = "white")  
  i = i+1
}
