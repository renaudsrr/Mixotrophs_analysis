# Script by SERRE Renaud
# Internship UQAM-GRIL under supervision of BEISNER Beatrix

rm(list = ls())
graphics.off()
cat("\014") 

#####################################################################################################
# Taxons majoritaire pour relation div~prevMixo ---------------------------------------------------------------------------------------
#####################################################################################################

# --- CHOIX DATASET A ANALYSER ---------------------------------------------------------------------------
choice = 3
if (choice ==1){source("/Users/renaudsrr/Desktop/STAGE_MTL/MODELISATION/CODE/R/1_Processing_data/1_Processing_CARBBAS.R")
} else if (choice ==2){source("/Users/renaudsrr/Desktop/STAGE_MTL/MODELISATION/CODE/R/1_Processing_data/2_Processing_LP_NLA.R")
} else if (choice ==3){source("/Users/renaudsrr/Desktop/STAGE_MTL/MODELISATION/CODE/R/1_Processing_data/3_Processing_IISD_ELA.R")}
choice = 3
# --------------------------------------------------------------------------------------------------------

library(dplyr)
library(tidyverse)
library(grid) 
library(forcats)
library(cowplot)

# ------------------------------------------------------------------------------

if (choice ==1){
  
}

if (choice ==2){
  data_taxon = data_phyto %>%
    rename("Genus" = "Target_taxon",
           "biovol_tot" = "total.biov") %>%
    left_join(taxo_info,by=c("Genus")) %>%
    rename(biovol_um3.mL = biovol)
    left_join(summary_data_phyto,by="lake_id") %>%
    select(-c(classic.group,Survey.x,Survey.y,biovol_tot)) %>%
    select(lake_id,biovol_um3.mL,biovol_mixo_um3.mL,biovol_tot_um3.mL,Empire,Kingdom,Phylum,Class,Order,Family,Genus,everything())
  
  data_taxon = data_taxon[-which(data_taxon$prev_Mixo==100),]
}

if (choice ==3){
  data_taxon = data_phyto %>%
    left_join(summary_data_phyto,by=c("lake_id","date_sampling")) 
}
# Prep plots --------------------------------------------------------------

# div_metrics
gg_H = ggplot(summary_data_phyto, aes(x = prev_Mixo, y = shannon)) +
  geom_point(alpha = 0.3, color = "black") +
  theme_minimal()
gg_1_D = ggplot(summary_data_phyto, aes(x = prev_Mixo, y = op_simpson)) +
  geom_point(alpha = 0.3, color = "black") +
  theme_minimal()
gg_J = ggplot(summary_data_phyto, aes(x = prev_Mixo, y = eveness_piel)) +
  geom_point(alpha = 0.3, color = "black") +
  theme_minimal()

all_type_taxon = c("Empire","Kingdom","Phylum",
                   "Class","Order","Family","Genus","strat")

palette_taxon = c(
  "#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e",
  "#e6ab02", "#a6761d", "#666666", "#8dd3c7", "#fb8072",  
  "#80b1d3", "#fdb462", "#b3de69", "#fccde5", "#bc80bd",  
  "#cab2d6", "#ffffb3", "#b15928",                        
  "black", "white" # Other et NA                                
)

bins = 10 # nb separation % prevalence mixo
intervals = 100 %/% bins  
data_taxon = data_taxon %>%
  mutate(bins = prev_Mixo %/% intervals + 1)

top_n_taxa = 15  # nombre de modalités à garder
i = 1

for (taxon in all_type_taxon) {
  
  data_taxon_grouped = data_taxon %>%
    mutate(taxon_grouped = fct_lump_n(.data[[taxon]], n = top_n_taxa, other_level = "Other"))
  
  table_taxon = data_taxon_grouped %>%
    group_by(bin = bins, taxon_grouped) %>%
    summarise(
      freq = n(),
      biovol_tot = sum(biovol_um3.mL, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      percent_freq = 100 * freq / ave(freq, bin, FUN = sum),
      percent_biovol = 100 * biovol_tot / ave(biovol_tot, bin, FUN = sum)
    )
  
  write.csv(table_taxon,paste0("/Users/renaudsrr/Desktop/STAGE_MTL/MODELISATION/DATA/0_new_csv/table_taxons/",i,"_Table_taxons_",taxon,".csv"),row.names = F)
  
  richesse_taxon = data_taxon_grouped %>%
    group_by(bin = bins) %>%
    summarise(richesse = n_distinct(.data[[taxon]])) 
  
  gg_rich_taxon = ggplot(richesse_taxon, aes(x = bin, y = richesse)) +
    geom_line(color = "black") +
    theme_minimal() +
    labs(y = paste0("Richesse en ", taxon), x = "Bin de % Mixo")
  
  # recup legend
  gg_for_legend = ggplot(table_taxon, aes(x = bin, y = percent_biovol, fill = taxon_grouped)) +
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
  gg_freq2 = ggplot(table_taxon, aes(x = bin, y = biovol_tot, fill = taxon_grouped)) +
    geom_bar(stat = "identity") +
    labs(y = "Biovol cumulé") +
    scale_fill_manual(values = palette_taxon) +
    theme_minimal() +
    theme(legend.position = "none")
  
  # % Biomass taxon
  gg_percent2 = ggplot(table_taxon, aes(x = bin, y = percent_biovol, fill = taxon_grouped)) +
    geom_bar(stat = "identity") +
    labs(y = "% biovol par taxon") +
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
  final_plot = plot_grid(legend_plot, plots_no_legend, ncol = 2, rel_widths = c(0.15, 1))
  
  if (choice==1){filename = paste0("/Users/renaudsrr/Desktop/STAGE_MTL/MODELISATION/figures/panels_CARBBAS/",i,"_Panel_", taxon, ".png")}
  if (choice==2){filename = paste0("/Users/renaudsrr/Desktop/STAGE_MTL/MODELISATION/figures/panels_LP_NLA/",i,"_Panel_", taxon, ".png")}
  if (choice==3){filename = paste0("/Users/renaudsrr/Desktop/STAGE_MTL/MODELISATION/figures/panels_IISD_ELA/",i,"_Panel_", taxon, ".png")}
  ggsave(filename, plot = final_plot, width = 18, height = 9, dpi = 600, bg = "white")
  
  i = i+1
}
