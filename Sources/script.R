################################################################################
#####################        Script R 10/12/2024         #######################
#####################     Developed by Maxime Naour      #######################
################################################################################

# Loading libraries
.libPaths(c("C:/Users/Maxime/AppData/Local/Programs/R/R-4.4.1/library", .libPaths()))
library(devtools)
library(BiocManager)
library(igraph); packageVersion("igraph")
library(ggrepel); packageVersion("ggrepel")
library(ggprism); packageVersion("ggprism")
library(readxl); packageVersion("readxl")
library(openxlsx); packageVersion("openxlsx")
library(writexl); packageVersion("writexl")
library(jsonlite); packageVersion("jsonlite")
library(microbiome); packageVersion("microbiome")
library(microbiomeutilities); packageVersion("microbiomeutilities")
library(dplyr); packageVersion("dplyr")
library(vegan); packageVersion("vegan")
library(ggplot2); packageVersion("ggplot2")
library(ggh4x); packageVersion("ggh4x")
library(ggforce); packageVersion("ggforce")
library(ggnewscale); packageVersion("ggnewscale")
library(phyloseq); packageVersion("phyloseq")
library(phyloseq.extended); packageVersion("phyloseq.extended")
library(biomformat); packageVersion("biomformat")
library(DESeq2); packageVersion("DESeq2")
library(tibble); packageVersion("tibble")  
library(plyr); packageVersion("plyr")
library(gcookbook); packageVersion("gcookbook")
library(tidyverse); packageVersion("tidyverse")
library(viridis); packageVersion("viridis")
library(gdata); packageVersion("gdata")
library(data.table); packageVersion("data.table")
library(ade4); packageVersion("ade4")
library(effsize); packageVersion("effsize")
library(ComplexHeatmap); packageVersion("ComplexHeatmap")
library(circlize); packageVersion("circlize")
library(RColorBrewer); packageVersion("RColorBrewer")
library(cluster); packageVersion("cluster")
library(factoextra); packageVersion("factoextra")
library(gridExtra); packageVersion("gridExtra")
library(magick); packageVersion("magick")
library(fpc); packageVersion("fpc")
library(mixOmics); packageVersion("mixOmics") # devtools::install_github("mixOmicsTeam/mixOmics")
library(ropls); packageVersion("ropls")
library(ggtree); packageVersion("ggtree")
library(ALDEx2); packageVersion("ALDEx2") #BiocManager::install("ALDEx2") & install.packages("MASS", lib = "C:/Users/Maxime/AppData/Local/R/win-library/4.2")

# Color palettes
ggsci.pals <- c("npg", "aaas", "nejm", "lancet", "jama", "jco", "ucscgb", "d3", 
                "locuszoom", "igv", "cosmic", "uchicago", "startrek", "tron", 
                "futurama", "rickandmorty", "simpsons", "flatui", "frontiers", "gsea", "material")

# Loop to count the number of colours in each palette
df <- data.frame()
for (pal in ggsci.pals) {
  n <- 1
  max_colors_found <- FALSE
  last_valid_colors <- NULL
  while(!max_colors_found) {
    pal_fun_name <- paste0("pal_", pal)
    if (!exists(pal_fun_name, where = asNamespace("ggsci"))) {
      message(pal, " does not exist in ggsci")
      break
    }
    pal_fun <- get(pal_fun_name, envir = asNamespace("ggsci"))
    color_generator <- pal_fun()
    colors <- color_generator(n)
    if (any(is.na(colors))) {
      num_colors <- sum(!is.na(colors))
      message(pal, " a un maximum de ", num_colors, " colours. The colours are :")
      print(last_valid_colors)
      max_colors_found <- TRUE
      df <- rbind(df, data.frame(palette = rep(pal, num_colors), color = last_valid_colors, order = 1:num_colors))
    } else {
      last_valid_colors <- colors
      n <- n + 1
    }
  }
}

# Display the colours in each palette
ggplot(df, aes(x = palette, y = order, fill = color)) +
  geom_tile() +
  scale_fill_identity() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "Palette", y = "Colour order", fill = "Colour") +
  guides(fill = "none")


# Importing data

# Move to the working directory
setwd("/Users/Maxime/Desktop/Projects/HB/")

# Import modules tables
func.df <- read.csv(file = "Data/atg16_tuft/atg16_meteor2_output_modules.tsv", 
                    header = TRUE, sep = "\t", stringsAsFactors = FALSE)
head(func.df)

# Import functional tables
kegg.df <- read.csv(file = "Data/atg16_tuft/atg16_meteor2_output_kegg_as_genes_sum.tsv", 
                header = TRUE, sep = "\t", stringsAsFactors = FALSE)
head(kegg.df)

mustard.df <- read.csv(file = "Data/atg16_tuft/atg16_meteor2_output_mustard_as_genes_sum.tsv", 
                    header = TRUE, sep = "\t", stringsAsFactors = FALSE)
head(mustard.df)

dbcan.df <- read.csv(file = "Data/atg16_tuft/atg16_meteor2_output_dbcan_as_genes_sum.tsv", 
                    header = TRUE, sep = "\t", stringsAsFactors = FALSE)
head(dbcan.df)


# TAXONOMY PART

# Move to the working directory
setwd("/Users/Maxime/Desktop/Projects/HB/")

# Import metadata
metadata <- read.csv("Data/atg16_tuft/metadata_atg16.csv", 
                     header = TRUE, sep = ",", stringsAsFactors = FALSE)
rownames(metadata) <- metadata$SampleID
head(metadata)

# Import report from METEOR2 analysis
report <- read.csv(file = "Data/atg16_tuft/atg16_meteor2_output_report.tsv", 
                   header = TRUE, sep = "\t", stringsAsFactors = FALSE)
head(report)

# Import tax table
tax.df <- read.csv(file = "Data/atg16_tuft/atg16_meteor2_output_msp_taxonomy.tsv", 
                header = TRUE, sep = "\t", stringsAsFactors = FALSE)
rownames(tax.df) <- tax.df$msp_name
tax.df$msp_name <- NULL
head(tax.df)

# Import MSP tax table
msp.df <- read.csv(file = "Data/atg16_tuft/atg16_meteor2_output_msp.tsv", 
                header = TRUE, sep = "\t", stringsAsFactors = FALSE)
rownames(msp.df) <- msp.df$msp_name
msp.df$msp_name <- NULL
head(msp.df)

# Import the Newick phylogenetic tree file
tree_file <- "Data/atg16_tuft/mm_5_0_gut_1252MSPs.nwk"
tree <- phyloseq::read_tree_greengenes(tree_file)

# Convert df to matrix
otu_mat <- as.matrix(msp.df)
tax_mat <- as.matrix(tax.df)

# Prepare matrices for generating phyloseq object
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE)
TAX <- tax_table(tax_mat)

# Create phyloseq object
ps <- phyloseq(OTU, TAX)
print(ps)

# Importing the metadata and phylogenetic tree file
sample_data(ps) <- metadata
phyloseq::phy_tree(ps) <- tree

# Display the first part of each component of object phyloseq
head(sample_names(ps))
head(sample_data(ps))
head(otu_table(ps))
head(tax_table(ps))
phy_tree(ps)

# Save the phyloseq (ps) object
saveRDS(ps, "Data/atg16_tuft/ps_atg16_coverage.rds")

# Processing report

# Move to the working directory
setwd("/Users/Maxime/Desktop/Projects/HB/")

# Lire l'objet ps
ps <- readRDS("Data/atg16_tuft/ps_atg16_coverage.rds")

# Extraction des noms d'échantillons
samples <- rownames(sample_data(ps))

# Import metadata
metadata <- read.csv("Data/atg16_tuft/metadata_atg16.csv", 
                     header = TRUE, sep = ",", stringsAsFactors = FALSE)
rownames(metadata) <- metadata$SampleID
head(metadata)

# Associer les données du tableau `report` aux échantillons dans `ps`
report_subset <- report %>%
  dplyr::filter(sample %in% samples) %>%
  mutate(raw = total_read_count, 
         mapped = mapped_read_count)

# Créer un dataframe combinant les informations de `sample_data` et `report`
df <- data.frame(
  samples = samples,
  raw = report_subset$raw,
  mapped = report_subset$mapped,
  ttt = sample_data(ps)$Condition
)

# Spécifier l'ordre des conditions
cond <- c("Atg16wt", "Atg16ko", "Pou2f3ko", "Atg16ko; Pou2f3ko")

# Modifier l'ordre des conditions dans le dataframe
df$ttt <- factor(df$ttt, levels = cond)
head(df$ttt)

# Préparer le dataframe
df_melt <- df %>%
  mutate(Sex = metadata[as.character(samples), "Sex"]) %>%
  pivot_longer(-c(samples, ttt, Sex), names_to = "step", values_to = "reads") %>%
  arrange(ttt, factor(Sex, levels = c("Female", "Male")), samples)

# Redéfinir les facteurs pour les échantillons et conditions
df_melt$samples <- factor(df_melt$samples, levels = unique(df_melt$samples))
df_melt$ttt <- factor(df_melt$ttt, levels = cond)
df_melt$step <- factor(df_melt$step, levels = c("raw", "mapped"))

# Associer des couleurs spécifiques aux conditions
cond.colors <- RColorBrewer::brewer.pal(length(cond), "Set1")
names(cond.colors) <- cond

# Identifier le nombre minimum et maximum de reads
min_reads <- min(df_melt$reads, na.rm = TRUE)
max_reads <- max(df_melt$reads, na.rm = TRUE)

# Ajuster la limite inférieure pour inclure la plus petite valeur
min_reads <- min_reads - (0.05 * (max_reads - min_reads))
if (min_reads < 0) min_reads <- 0 

# Associer des couleurs spécifiques pour le sexe
sex.colors <- c("Male" = "#7BCFE5", "Female" = "#FEB2B2")

# Matcher les conditions (ttt) du dataframe df_melt avec les échantillons du dataframe metadata 
# Cela ajoute une colonne ttt à metadata, permettant de lier chaque échantillon à sa condition expérimentale.
metadata$ttt <- df_melt$ttt[match(metadata$SampleID, df_melt$samples)]

# Générer le barplot
process.p <- ggplot(df_melt, aes(x = samples, y = reads, fill = step)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
  scale_fill_brewer(palette = "Set2", name = "Step") +
  ggnewscale::new_scale_fill() +
  geom_bar(aes(fill = ttt), stat = "identity", alpha = 0, show.legend = TRUE) +
  scale_fill_manual(name = "Condition", values = cond.colors, guide = guide_legend(override.aes = list(alpha = 1))) +
  ggnewscale::new_scale_fill() +
  geom_tile(
    data = metadata,
    aes(x = SampleID, y = -0.02 * max_reads, fill = Sex),
    height = 0.03 * max_reads,
    inherit.aes = FALSE
  ) +
  scale_fill_manual(name = "Sex", values = sex.colors, guide = guide_legend(title = "Sex")) +
  scale_y_continuous(
    limits = c(-0.03 * max(df_melt$reads), max(df_melt$reads)), 
    labels = function(y) ifelse(y > 0, y, ""),
    oob = scales::squish
  ) +
  facet_grid(.~ttt, scales = "free_x", space = "free_x") +
  scale_x_discrete(labels = function(x) {
    sapply(x, function(sample) {
      sprintf('<span style="color:%s;">%s</span>', cond.colors[as.character(df[df$samples == sample, "ttt"])], sample)
    })
  }) +
  coord_cartesian(ylim = c(-0.03 * max(df_melt$reads), max(df_melt$reads)), expand = FALSE) +
  theme_classic2() +
  theme(
    axis.text.x = ggtext::element_markdown(angle = 45, vjust = 1, hjust = 1, face = "bold", size = 12),
    axis.text.y = element_text(face = "bold", size = 12),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.title = element_text(face = "bold", size = 14),
    strip.text = element_blank(),
    plot.margin = margin(10, 10, 10, 10)
  ) +
  labs(
    x = "",
    y = "Reads count",
    title = "Read counts at each step"
  )

# Afficher le graphique
print(process.p)

# Enregistrer le graphique
ggsave("Figures/Taxonomic/processing/overview_processing_atg16.png",
       process.p, width = 14, height = 10, dpi = 300, bg = "white")

# Liste des configurations à générer
configs <- list(
  list(name = "by_cond", grouping_vars = c("ttt", "step"), facet = FALSE),
  list(name = "by_cond_and_sex", grouping_vars = c("ttt", "step", "Sex"), facet = TRUE)
)

# Boucle pour générer les graphiques selon les configurations
for (config in configs) {
  # Calcul des moyennes et écarts-types
  df_summary <- df_melt %>%
    group_by(across(all_of(config$grouping_vars))) %>%
    summarise(
      mean = mean(reads, na.rm = TRUE),
      sd = sd(reads, na.rm = TRUE),
      n = n(),
      .groups = 'drop'
    )
  
  # Identifier les limites de l'axe y
  min_reads <- min(df_summary$mean, na.rm = TRUE)
  max_reads <- max(df_summary$mean, na.rm = TRUE)
  
  # Ajuster la limite inférieure
  min_reads <- min_reads - (0.05 * (max_reads - min_reads))
  if (min_reads < 0) min_reads <- 0
  
  # Initialisation du barplot
  process.bp <- ggplot(df_summary, aes(x = ttt, y = mean, fill = step)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
    geom_errorbar(
      aes(ymin = mean - sd, ymax = mean + sd),
      width = 0.2,
      position = position_dodge(0.7)
    ) +
    labs(
      x = "Condition",
      y = "Reads count",
      fill = "Step",
      title = "Average read counts at different processing steps"
    ) +
    scale_fill_brewer(palette = "Set2", name = "Step") +
    scale_y_continuous(limits = c(min_reads, max_reads), oob = scales::squish) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, face = "bold", size = 12),
      axis.text.y = element_text(face = "bold", size = 12),
      legend.position = "right",
      legend.title = element_text(face = "bold", size = 12),
      legend.text = element_text(size = 12),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      axis.title = element_text(face = "bold", size = 14),
      strip.text = element_text(face = "bold", size = 14)
    )
  # Ajouter les facettes
  if (config$facet) {
    process.bp <- process.bp + facet_wrap(~Sex, nrow = 1)
  }
  
  # Afficher le graphique
  print(process.bp)
  
  # Enregistrer le graphique
  ggsave(
    filename = paste0("Figures/Taxonomic/processing/overview_processing_", config$name, ".png"),
    plot = process.bp, width = 10, height = 8, dpi = 300, bg = "white"
  )
}

# Rarefaction curves

# Move to the working directory
setwd("/Users/Maxime/Desktop/Projects/HB/")

# Lire l'objet ps
ps <- readRDS("Data/atg16_tuft/ps_atg16_coverage.rds")
head(ps@sam_data)

# Spécifier l'ordre des conditions
cond <- c("Atg16wt", "Atg16ko", "Pou2f3ko", "Atg16ko; Pou2f3ko")

# Sélectionner les individus selon leur Condition
ps.filt <- subset_samples(ps, Condition %in% cond)

# Modifier l'ordre des niveaux de facteur pour les Conditions
sample_data(ps.filt)$Condition <- factor(sample_data(ps.filt)$Condition, levels = cond)

# Get the total number of reads for each sample and sort the results
sort(sample_sums(ps.filt))

# Conversion en "pseudo-comptes" entiers
a <- 250
otu_table(ps.filt) <- round(otu_table(ps.filt) * a)
saveRDS(ps.filt, "Data/atg16_tuft/ps_atg16_coverage_alpha_corr.rds")

# Set the cutoff
threshold <- min(sample_sums(ps.filt))
print(threshold)

# Associer des couleurs spécifiques aux conditions
cond.colors <- RColorBrewer::brewer.pal(length(cond), "Set1")
names(cond.colors) <- cond

# Création du graphique de rarefaction
rare.plot <- ggrare(physeq = ps.filt, step = 100, color = "Condition", label = "SampleID", se = FALSE, parallel = TRUE) + 
  geom_line(linewidth = 0.7) + 
  scale_color_manual(name = "Condition", values = cond.colors) +
  geom_vline(xintercept = threshold, linetype = "dashed", color = "black") +
  labs(
    x = "Reads count", y = "MSP richness", fill = "Step", 
    title = "Rarefaction curves",
    caption = "Normalized data"
  ) +
  theme_light() +
  theme(
    axis.text.x = element_text(face = "bold", size = 12),
    axis.text.y = element_text(face = "bold", size = 12),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.title = element_text(face = "bold", size = 14), 
    plot.caption = element_text(face = "bold", size = 12)
  )
print(rare.plot)

# Enregistrer le plot
ggsave(rare.plot, file = paste0("Figures/Taxonomic/processing/rerafy_curves_atg16_cutoff_",threshold,".png"), 
       width = 12, height = 10, dpi = 300, bg = "white")


# Richesse et Alpha-diversity

# Move to the working directory
setwd("/Users/Maxime/Desktop/Projects/HB/")

# Lire l'objet ps
ps <- readRDS("Data/atg16_tuft/ps_atg16_coverage_alpha_corr.rds")
print(ps)
dim(otu_table(ps))
head(ps@otu_table)

# Spécifier l'ordre des conditions
cond <- c("Atg16wt", "Atg16ko", "Pou2f3ko", "Atg16ko; Pou2f3ko")

# Sélectionner les individus selon leur Condition
ps.filt <- subset_samples(ps, Condition %in% cond)

# Modifier l'ordre des niveaux de facteur pour les Conditions
sample_data(ps.filt)$Condition <- factor(sample_data(ps.filt)$Condition, levels = cond)

# Associer des couleurs spécifiques aux conditions
cond.colors <- RColorBrewer::brewer.pal(length(cond), "Set1")
names(cond.colors) <- cond

# Calcul des mesures de diversité
div.est <- estimate_richness(ps.filt, measures = c("Observed", "Shannon", "Simpson", "InvSimpson", "Fisher"))
print(div.est)

# Conversion en dataframe
div.df <- as.data.frame(div.est)
div.df$Sample <- rownames(div.df)
div.df <- merge(div.df, data.frame(sample_data(ps.filt)), by.x="Sample", by.y="row.names")

# Fonction pour générer les analyses + plots
rich.alpha <- function(df, measure, split_cond, title, colors) {
  # Calculer les tailles des échantillons
  sample_sizes <- df %>%
    dplyr::group_by(Condition) %>%
    dplyr::summarize(N = n(), .groups = "drop")
  
  # Ajouter les étiquettes avec les tailles des échantillons
  df <- df %>%
    dplyr::left_join(sample_sizes, by = "Condition") %>%
    dplyr::mutate(Condition_label = paste0(Condition, "\n(n=", N, ")"))
  
  # Créer l'ordre des niveaux pour Condition et Condition_label
  df$Condition <- factor(df$Condition, levels = cond)
  df$Condition_label <- factor(
    df$Condition_label,
    levels = paste0(cond, "\n(n=", sample_sizes$N[match(cond, sample_sizes$Condition)], ")")
  )
  
  # Comparaisons pour les tests statistiques
  comp <- combn(levels(df$Condition), 2, simplify = FALSE)
  
  # Comparaisons pour les tests statistiques
  if (!is.null(split_cond) && split_cond != "") {
    test_result <- df %>%
      group_by(!!sym(split_cond)) %>%
      summarise(
        p_value = kruskal.test(as.formula(paste(measure, "~ Condition")), data = cur_data())$p.value,
        .groups = "drop"
      )
    
    caption.text <- paste0(
      "Kruskal-Wallis p-values:\n",
      paste(test_result[[split_cond]], sprintf("%.3f", test_result$p_value), sep = ": ", collapse = "; ")
    )
  } else {
    test_result <- kruskal.test(as.formula(paste(measure, "~ Condition")), data = df)
    caption.text <- paste0(
      "Kruskal-Wallis p-value: ", sprintf("%.3f", test_result$p.value)
    )
  }
  
  # Créer le plot
  p <- ggplot(df, aes(x = Condition, y = !!rlang::sym(measure), fill = Condition)) +
    geom_violin(trim = TRUE, scale = "width", alpha = 0.7) +
    geom_boxplot(width = 0.2, outlier.shape = NA, color = "black") +
    geom_jitter(width = 0.2, size = 1.5, color = "black") +
    scale_x_discrete(labels = levels(df$Condition_label)) +
    labs(title = title, y = measure, x = "", caption = caption.text) +
    theme_bw() +
    theme(
      panel.grid.major.y = element_line(color = "grey80", size = 0.75),
      legend.position = "none",
      axis.title.y = element_text(size = 14, face = "bold"),
      axis.text.x = element_text(color = "black", size = 12, face = "bold", angle = 45, hjust = 1, vjust = 1), 
      axis.text.y = element_text(size = 12, face = "bold"), 
      plot.caption = element_text(size = 12, hjust = 0.5, vjust = 1)
    ) +
    scale_fill_manual(values = colors) +
    stat_compare_means(comparisons = comp, label = "p.format", method = "t.test") #p.format or p.signif
  
  # Appliquer le splitting si split_cond est spécifiés
  if (!is.null(split_cond) && split_cond != "") {
    p <- p + facet_grid(as.formula(paste(". ~", split_cond)))
  }
  
  return(p)
}

# Générer les plots
rich.plots <- lapply(c("Observed", "Shannon", "Simpson", "InvSimpson", "Fisher"), function(measure) {
  rich.alpha(div.df, measure, split_cond = "Sex", title = measure, colors = cond.colors)
})

# Combiner les plots
comb <- do.call(grid.arrange, c(rich.plots, ncol = 3))

# Sauvegarder les plots combinés
ggsave("Figures/Taxonomic/div/adiv_atg16_coverage.png", 
       comb, width = 16, height = 12, bg = "white", dpi = 300)

# Générer toutes les combinaisons possibles de comparaisons
all.comp <- combn(cond, 2, simplify = FALSE)

# Créer un dataframe pour les p-values
test_results <- data.frame(Measure = character(), Comparison = character(), P_value = numeric(), stringsAsFactors = FALSE)

# Appliquer la fonction mise à jour à tous les indices et enregistrer les p-values
for (measure in c("Observed", "Shannon", "Simpson", "InvSimpson", "Fisher")) {
  for (comp in all_comparisons) {
    test_res <- wilcox.test(
      div.df[[measure]][div.df$Condition == comp[1]],
      div.df[[measure]][div.df$Condition == comp[2]],
      exact = FALSE
    )
    test_results <- rbind(test_results, data.frame(
      Measure = measure,
      Comparison = paste(comp, collapse = " vs "),
      P_value = test_res$p.value
    ))
  }
}

# Sauvegarder les résultats des p-values dans un fichier CSV
write.csv(test_results, "Data/output/alpha_diversity_stats.csv", row.names = FALSE)


# LOG-TRANSFORMATION

# Move to the working directory
setwd("/Users/Maxime/Desktop/Projects/HB/")

# Lire l'objet ps
ps <- readRDS("Data/atg16_tuft/ps_atg16_coverage_alpha_corr.rds")
print(ps)
dim(otu_table(ps))
head(ps@otu_table)

# Spécifier l'ordre des conditions
cond <- c("Atg16wt", "Atg16ko", "Pou2f3ko", "Atg16ko; Pou2f3ko")

# Sélectionner les individus selon leur Condition
ps.filt <- subset_samples(ps, Condition %in% cond)

# Modifier l'ordre des niveaux de facteur pour les Conditions
sample_data(ps.filt)$Condition <- factor(sample_data(ps.filt)$Condition, levels = cond)

# Paramètres pour "denom": 
# all : Moyenne géométrique de l'abondance de toutes les caractéristiques (méthode par défaut de clr)
# iqlr : Caractéristiques entre les 1er et 3e quartiles de variance dans les valeurs clr.
# zero : Caractéristiques non nulles dans chaque groupe.
# lvha : Caractéristiques avec faible variance (quartile inférieur) et haute abondance relative (quartile supérieur).
# median : Médiane de l'abondance de toutes les caractéristiques.
# user : Vecteur personnalisé d'indices de lignes spécifiés par l'utilisateur.

# Paramètres pour "denom": 
features <- "all"

# Vecteur des tailles de mc.samples 
mc.samples <- c(64, 128, 256, seq(500, 5000, by = 500))

# Répertoire de sauvegarde des figures
output_dir <- "Figures/Taxonomic/processing/LOG_TRANS/"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Initialiser la liste des résultats
var_list <- list()

# Charger les données
dataset <- as.data.frame(otu_table(ps.filt))
rownames(dataset) <- rownames(otu_table(ps.filt))

# Préparation des métadonnées
metadata_filt <- data.frame(sample_data(ps.filt))

# Tester pour chaque référence et chaque valeur de mc.samples
for (ref in cond) {
  cat("Processing reference =", ref, "\n")
  
  # Mettre à jour le groupe de référence
  metadata_filt$Condition <- relevel(as.factor(metadata_filt$Condition), ref = ref)
  mm <- model.matrix(~ Condition, metadata_filt)
  
  # Initialiser les résultats pour cette combinaison
  var_df <- data.frame(MC_Samples = numeric(), Variance = numeric())
  
  # Tester chaque valeur de mc.samples
  for (mc_sample in mc.samples) {
    cat("Processing mc.samples =", mc_sample, "\n")
    
    # Appliquer ALDEx2
    aldex.clr <- aldex.clr(dataset, mm, mc.samples = mc_sample, denom = features, useMC = FALSE, verbose = TRUE)
    mc_instances <- getMonteCarloInstances(aldex.clr)
    
    # Calculer les moyennes CLR pour chaque instance Monte Carlo
    clr_results <- sapply(mc_instances, function(instance) {
      rowMeans(instance, na.rm = TRUE)
    })
    
    # Calculer la variance pour chaque fonctionnalité à travers les instances Monte Carlo
    feature_variances <- apply(clr_results, 1, var)
    var_total <- sum(feature_variances)
    
    # Stocker les résultats
    var_df <- rbind(var_df, data.frame(MC_Samples = mc_sample, Variance = var_total))
  }
  
  # Enregistrer les résultats dans la liste
  var_list[[paste0("ref_", ref)]] <- var_df
  
  # Générer et sauvegarder le plot
  plot <- ggplot(var_df, aes(x = MC_Samples, y = Variance)) +
    geom_point(color = "black") +
    geom_smooth(method = "loess", formula = y ~ x, span = 0.4, level = 0.95, color = "red3", fill = "lightpink") +
    theme_bw() +
    labs(
      x = "MC Samples",
      y = "Variance of CLR means across MC instances",
      title = paste("Evolution of variance for reference", ref)
    ) +
    theme(plot.title = element_text(hjust = 0.5))
  
  # Sauvegarder le graphique
  ggsave(
    filename = paste0(output_dir, "variance_evolution_ref_", ref, ".png"),
    plot = plot,
    width = 8,
    height = 6
  )
}

# Convertir les colonnes nécessaires en facteurs
metadata_filt$gp <- as.factor(metadata_filt$Condition)

# Créer une matrice de modèle pour inclure toutes les comparaisons
mm <- model.matrix(~ 0 + gp, metadata_filt)
colnames(mm) <- gsub("gp", "", colnames(mm))

# Affichage pour vérification
cat("Dimensions de la matrice de modèle :", dim(mm), "\n")
cat("Colonnes de la matrice de modèle :", colnames(mm), "\n")

# Exécution d'ALDEx2 pour log-transformer les données
clr.res <- aldex.clr(reads = dataset, 
                              conds = mm,
                              mc.samples = 1000, 
                              denom = features,
                              useMC = TRUE, 
                              verbose = TRUE)

# Obtenir les valeurs log-transformées par CLR 
mc.instances <- getMonteCarloInstances(clr.res)
df.clr <- sapply(names(mc.instances), function(sample) rowMeans(mc.instances[[sample]], na.rm = TRUE))

# Fonction pour exécuter ALDEx2 pour chaque Condition comme référence pour les effets
comp.gps <- function(ref, data, metadata, gp) {
  data <- as.data.frame(data)
  metadata$gp <- relevel(as.factor(metadata[[gp]]), ref = ref)
  mm <- model.matrix(~ gp, metadata)
  colnames(mm) <- gsub("^Condition(.*)$", "Condition_\1", colnames(mm))
  print(dim(mm))
  print(colnames(mm))
  aldex.clr <- aldex.clr(data, mm, mc.samples = 1000, denom = features, useMC = FALSE, verbose = TRUE)
  glm <- aldex.glm(aldex.clr, verbose = TRUE)
  glm.effect <- aldex.glm.effect(aldex.clr, useMC = FALSE, CI = 95, verbose = TRUE)
  df <- data.frame(glm, glm.effect, check.names = FALSE)
  return(df)
}

# Automation to change reference groups
all <- list()

# Itérer sur tous les groupes pour redéfinir le groupe de référence (sans le dernier)
for (ref in cond[-length(cond)]) {
 cat("Computing for Condition with reference", ref, "\n")
 all[[ref]] <- comp.gps(as.character(ref), dataset, metadata_filt, "Condition")
}

# Check the dataframes
for (name in names(all)) {
  cat("\nResult for", name, ":\n")
  print(head(all[[name]], n = 2L))
}

# Liste de dataframes filtrés
df.filt <- list()

# Supprimer les colonnes inutiles et renommer les colonnes en ajoutant le groupe de référence
for (df in names(all)) {
  subdf <- all[[df]][, !grepl("^Intercept::", colnames(all[[df]])), drop = FALSE]
  colnames(subdf) <- paste0("ref_", df, "_", colnames(subdf))
  df.filt[[df]] <- subdf
}

# Afficher les noms des colonnes pour vérifier leur format
for (n in names(df.filt)) {
  cat("\nNoms des colonnes dans le dataframe :", n, "\n")
  print(grep("gpAtg16wt", colnames(df.filt[[n]]), value = TRUE, invert = TRUE))
}


# Créer une liste pour stocker les comparaisons
comp <- list()

# Extraire les colonnes pour chaque comparaison basée sur cond
for (i in seq_along(cond)[-length(cond)]) {
  ref <- cond[i]
  other_groups <- cond[(i + 1):length(cond)]
  
  for (comp_group in other_groups) {
    # Construire le pattern pour sélectionner les colonnes correspondantes
    pattern <- paste0("^ref_", ref, "_gp", comp_group, "[:.]")
    
    # Extraire et combiner les colonnes correspondant au pattern
    extracted_cols <- lapply(df.filt, function(df) {
      df[, grepl(pattern, colnames(df)), drop = FALSE]
    })
    
    # Ajouter les colonnes extraites dans la liste
    combined_df <- do.call(cbind, extracted_cols)
    if (ncol(combined_df) > 0) {
      comp[[length(comp) + 1]] <- combined_df
    }
  }
}

# Consolider les comparaisons dans un seul tableau
df_diff <- do.call(cbind, comp)

# Nettoyer les noms de colonnes
colnames(df_diff) <- sub("^[^\\.]+\\.(ref_.*)", "\\1", colnames(df_diff))

# Afficher les noms des colonnes pour chaque groupe de référence
cat("\nNames of columns in the final dataframe for each reference group :\n")
for (ref in cond) {
  cat("\nReference :", ref, "\n")
  print(grep(paste0("^ref_", ref, "_"), colnames(df_diff), value = TRUE))
}

# Convertir df.clr en dataframe si nécessaire
df.clr <- as.data.frame(df.clr)

# Fusionner df.clr et df_diff
df_final <- cbind(df.clr, df_diff)

# Ajouter les row.names comme première colonne
df_final <- cbind(msp = rownames(df_final), df_final)

# Enregistrer le dataframe final
write_csv(df_final, "Data/output/df_final_aldex2_log_clr_glm.csv")
write_xlsx(df_final, "Data/output/df_final_aldex2_log_clr_glm.xlsx")

# Ecraser la table de comptage par celle log-transformée dans l'objet phyloseq
ps.clr <- ps.filt
otu_table(ps.clr) <- otu_table(as.matrix(df.clr), taxa_are_rows = TRUE)

# Sauvegarder le nouvel objet phyloseq log-transformé
saveRDS(ps.clr, file = "Data/atg16_tuft/ps_atg16_coverage_alpha_corr_clr.rds")


# Beta-diversity with CLR-transformed data

# Move to the working directory
setwd("/Users/Maxime/Desktop/Projects/HB/")

# Lire l'objet ps non transformée (pour philR)
ps <- readRDS("Data/atg16_tuft/ps_atg16_coverage_alpha_corr.rds")
print(ps)
dim(otu_table(ps))
head(ps@otu_table)

# Lire l'objet ps log-transformé
ps.clr <- readRDS("Data/atg16_tuft/ps_atg16_coverage_alpha_corr_clr.rds")
print(ps.clr)
dim(otu_table(ps.clr))
head(ps.clr@otu_table)

# Spécifier l'ordre des conditions
cond <- c("Atg16wt", "Atg16ko", "Pou2f3ko", "Atg16ko; Pou2f3ko")

# Sélectionner les individus selon leur Condition
ps.filt <- subset_samples(ps, Condition %in% cond)
ps.clr.filt <- subset_samples(ps.clr, Condition %in% cond)

# Modifier l'ordre des niveaux de facteur pour les Conditions
sample_data(ps.filt)$Condition <- factor(sample_data(ps.filt)$Condition, levels = cond)
sample_data(ps.clr.filt)$Condition <- factor(sample_data(ps.clr.filt)$Condition, levels = cond)

# Associer des couleurs spécifiques aux conditions
cond.colors <- RColorBrewer::brewer.pal(length(cond), "Set1")
names(cond.colors) <- cond

# Fonction pour générer les analyses de beta-diversité
beta.ord <- function(ps, method, feat = NULL, seed = NULL) {
  data <- as.data.frame(t(otu_table(ps)))
  metadata <- as(sample_data(ps), "data.frame")
  
  compatible_methods <- c("NMDS", "MDS", "RDA", "philR")
  if (!method %in% compatible_methods) {
    stop(paste("La méthode", method, "n'est pas compatible avec les données CLR transformées."))
  }
  
  scores_df <- NULL
  method_title <- ""
  axis_names <- NULL
  explained_var <- NULL
  
  if (method == "philR") {
    tree <- phy_tree(ps)
    otus <- otu_table(ps) + 1e-6
    philr_obj <- philr(t(otus), tree, part.weights = "enorm.x.gm.counts", ilr.weights = "blw.sqrt")
    dist_philr <- dist(philr_obj)
    ord <- cmdscale(dist_philr, eig = TRUE, k = 2)
    method_title <- "philR"
    scores_df <- as.data.frame(ord$points)
    colnames(scores_df) <- c("philR1", "philR2")
    axis_names <- c("philR1", "philR2")
    explained_var <- ord$eig[1:2] / sum(ord$eig) * 100
  } else if (method == "NMDS") {
    if (!is.null(seed)) set.seed(seed)
    dist_aitchison <- vegdist(data, method = "euclidean")
    ord <- metaMDS(dist_aitchison, trymax = 100)
    method_title <- "NMDS"
    scores_df <- as.data.frame(ord$points)
    colnames(scores_df) <- c("NMDS1", "NMDS2")
    axis_names <- c("NMDS1", "NMDS2")
  } else if (method == "MDS") {
    dist_aitchison <- vegdist(data, method = "euclidean")
    ord <- cmdscale(dist_aitchison, eig = TRUE, k = 2)
    method_title <- "MDS"
    scores_df <- as.data.frame(ord$points)
    colnames(scores_df) <- c("MDS1", "MDS2")
    axis_names <- c("MDS1", "MDS2")
    explained_var <- ord$eig[1:2] / sum(ord$eig) * 100
  } else if (method == "RDA") {
    if (is.null(feat)) {
      stop("For RDA, you must specify the metadata columns with the 'feat' argument.")
    }
    formula <- as.formula(paste("data ~", paste(feat, collapse = " + ")))
    ord <- rda(formula, data = metadata)
    method_title <- "RDA"
    scores_df <- as.data.frame(scores(ord, display = "sites"))
    colnames(scores_df) <- c("RDA1", "RDA2")
    axis_names <- c("RDA1", "RDA2")
    explained_var <- summary(ord)$cont$importance[2, 1:2] * 100
  }
  
  scores_df$Condition <- sample_data(ps)$Condition
  
  if (method != "philR") {
    dist_aitchison <- vegdist(data, method = "euclidean")
    perm_result <- vegan::adonis2(dist_aitchison ~ Condition, data = metadata, permutations = 999)
  } else {
    perm_result <- vegan::adonis2(dist_philr ~ Condition, data = metadata, permutations = 999)
  }
  p_value <- perm_result$`Pr(>F)`[1]
  stars <- ifelse(p_value < 0.001, "***", ifelse(p_value < 0.01, "**", ifelse(p_value < 0.05, "*", "ns")))
  
  if (!is.null(explained_var)) {
    labs_x <- paste(axis_names[1], "[", round(explained_var[1], 2), "%]")
    labs_y <- paste(axis_names[2], "[", round(explained_var[2], 2), "%]")
  } else {
    labs_x <- axis_names[1]
    labs_y <- axis_names[2]
  }
  
  p_beta <- ggplot(scores_df, aes_string(x = axis_names[1], y = axis_names[2], color = "Condition")) +
    geom_point(size = 3) +
    theme_classic2() +
    labs(title = paste("Beta-diversity -", method_title),
         subtitle = paste("Ordination using", method_title),
         x = labs_x, y = labs_y) +
    scale_color_manual(values = cond.colors) +
    scale_fill_manual(values = cond.colors) +
    stat_ellipse(aes(color = Condition, fill = Condition), geom = "polygon", 
                 type = "norm", alpha = 0.05, linetype = 1, linewidth = 0.3) +
    annotate("text", x = Inf, y = Inf, label = paste("p =", format(p_value, digits = 2), stars),
             hjust = 1.1, vjust = 1.5, size = 4, color = "red3")
  
  centers <- aggregate(cbind(scores_df[[axis_names[1]]], scores_df[[axis_names[2]]]) ~ Condition, data = scores_df, FUN = mean)
  colnames(centers) <- c("Condition", paste(axis_names[1], "center", sep = "."), paste(axis_names[2], "center", sep = "."))
  
  segment_data <- scores_df %>%
    dplyr::left_join(centers, by = "Condition") %>%
    mutate(
      xend = scores_df[[axis_names[1]]],
      yend = scores_df[[axis_names[2]]],
      x = .[[paste(axis_names[1], "center", sep = ".")]],
      y = .[[paste(axis_names[2], "center", sep = ".")]]
    )
  
  p_beta <- p_beta +
    geom_segment(data = segment_data, aes(x = x, y = y, xend = xend, yend = yend),
                 alpha = 0.5)
  
  return(p_beta)
}


# Exécuter la fonction pour réaliser les analyses de beta-diversité
methods <- c("philR", "NMDS", "MDS", "RDA")
features <- c("Sex", "Condition", "Caecum_mg")
plots <- list()

for (method in methods) {
  if (method == "philR") {
    plots[[method]] <- beta.ord(ps.filt, method)
  } else {
    plots[[method]] <- beta.ord(ps.clr.filt, method, feat = features, seed = if (method == "NMDS") 123 else NULL)
  }
}

# Afficher les plots
for (method in methods) {
  print(plots[[method]])
}

# Enregistrer les plots
png("Figures/Taxonomic/div/beta_diversity_atg16_clr.png", units = "cm", width = 30, height = 25, res = 300)

grid.arrange(
  grobs = plots,
  ncol = 2,
  top = "Beta-diversity Analysis"
)

dev.off()


# COMPOSITIONAL

# Move to the working directory
setwd("/Users/Maxime/Desktop/Projects/HB/")

# Lire l'objet ps log-transformé
ps.clr <- readRDS("Data/atg16_tuft/ps_atg16_coverage_alpha_corr_clr.rds")
print(ps.clr)
dim(otu_table(ps.clr))
head(ps.clr@otu_table)

# Spécifier l'ordre des conditions
cond <- c("Atg16wt", "Atg16ko", "Pou2f3ko", "Atg16ko; Pou2f3ko")

# Sélectionner les individus selon leur Condition
ps.clr.filt <- subset_samples(ps.clr, Condition %in% cond)

# Modifier l'ordre des niveaux de facteur pour les Conditions
sample_data(ps.clr.filt)$Condition <- factor(sample_data(ps.clr.filt)$Condition, levels = cond)

# Associer des couleurs spécifiques aux conditions
cond.colors <- RColorBrewer::brewer.pal(length(cond), "Set1")
names(cond.colors) <- cond

# Revenir aux "pseudo-comptages bruts" à partir des données CLR
exp_data <- round(exp(otu_table(ps.clr.filt)))

# Ecraser la table de comptage log-transformées par cette table de comptage brut
ps.clr.raw <- ps.clr.filt
ps.clr.raw@otu_table <- exp_data

# Sauvegarder l'objet phyloseq
saveRDS(ps.clr.raw, "Data/atg16_tuft/ps_atg16_coverage_alpha_corr_clr_raw.rds")


# Lire l'objet ps log-transformé avec pseudo-comptage
ps.clr.raw <- readRDS("Data/atg16_tuft/ps_atg16_coverage_alpha_corr_clr_raw.rds")
print(ps.clr.raw)
dim(otu_table(ps.clr.raw))
head(ps.clr.raw@otu_table)

# Spécifier l'ordre des conditions
cond <- c("Atg16wt", "Atg16ko", "Pou2f3ko", "Atg16ko; Pou2f3ko")

# Sélectionner les individus selon leur Condition
ps.clr.raw.filt <- subset_samples(ps.clr.raw, Condition %in% cond)

# Modifier l'ordre des niveaux de facteur pour les Conditions
sample_data(ps.clr.raw.filt)$Condition <- factor(sample_data(ps.clr.raw.filt)$Condition, levels = cond)

# Associer des couleurs spécifiques aux conditions
cond.colors <- RColorBrewer::brewer.pal(length(cond), "Set1")
names(cond.colors) <- cond

# Fonction pour générer les barplots de l'analyse compositionnelle
comp_taxa <- function(ps_data_filt, tax_level, gp.colors, output_dir) {
  # Calculer le nombre de phyla pour définir la valeur de top
  top <- length(unique(psmelt(tax_glom(ps_data_filt, "phylum"))$phylum))
  
  # Agréger les donnéess pour les niveaux riches
  glom_tax <- aggregate_top_taxa2(ps_data_filt, top = top-1, tax_level)
  df_tax <- psmelt(glom_tax)
  
  # Supprimer les lignes avec tax_level non identifié
  df_tax <- df_tax[!is.na(df_tax[[tax_level]]),]
  
  # Calculer les proportions et convertir en pourcentage
  df_tax <- df_tax %>%
    dplyr::group_by(Sample) %>%
    dplyr::mutate(Proportion = Abundance / sum(Abundance) * 100) %>%
    ungroup()
  
  # Trier les échantillons par Condition
  df_tax <- df_tax %>%
    dplyr::mutate(
      Condition = factor(Condition, levels = levels(sample_data(ps_data_filt)$Condition))
    ) %>%
    dplyr::arrange(Condition) %>%
    dplyr::mutate(Sample = factor(Sample, levels = unique(Sample)))
  
  # Calculer l'abondance moyenne par Condition
  df_tax_mean <- df_tax %>%
    dplyr::group_by(Condition, !!sym(tax_level)) %>%
    dplyr::summarize(Proportion = mean(Proportion, na.rm = TRUE), .groups = "drop")
  
  # Palette de couleurs
  nb.cols <- length(unique(df_tax[[tax_level]]))
  palette <- ggsci::pal_frontiers(alpha = 0.5)(nb.cols)
  palette[is.na(palette)] <- c("#f0a5b4", "#f4bde8")[1:sum(is.na(palette))]
  palette <- c(palette[palette != "#706F6F7F"], "#706F6F7F")
  
  # Créer le barplot empilé pour tous les échantillons
  tax_plot <- ggplot(df_tax, aes(x = Sample, y = Proportion, fill = !!sym(tax_level))) +
    geom_bar(stat = "identity", position = "stack") +
    theme_classic2() +
    scale_fill_manual(values = palette, name = tax_level) +
    ggnewscale::new_scale_fill() +
    geom_bar(aes(fill = Condition), stat = "identity", alpha = 0, show.legend = TRUE) +
    scale_fill_manual(
      name = "Condition",
      values = cond.colors,
      guide = guide_legend(override.aes = list(alpha = 1))
    ) +
    scale_x_discrete(labels = function(x) {
      sapply(x, function(sample) {
        condition <- as.character(sample_info[sample_info$Sample == sample, "Condition"])
        sprintf('<span style="color:%s;">%s</span>', cond.colors[condition], sample)
      })
    }) +
    labs(
      x = "",
      y = "Relative Abundance (%)",
      fill = tax_level,
      title = paste("Relative Abundance of", tax_level, "by Sample"),
      caption = "CLR dataset converted in pseudo-counts"
    ) +
    theme(
      axis.text.x = ggtext::element_markdown(angle = 90, vjust = 0.5, hjust = 1, size = 12),
      axis.text.y = element_text(size = 12, face = "bold"),
      panel.grid.major.y = element_line(color = "grey80", linewidth = 0.75),
      legend.position = "right",
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 12),
      plot.caption = element_text(size = 10, face = "bold"),
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold", color = "black")
    )
  
  # Enregistrer le barplot
  ggsave(
    filename = file.path(output_dir, paste0("relabundance_", tax_level, "_samples.png")),
    plot = tax_plot, width = 10, height = 8, bg = "white", dpi = 300
  )
  
  # Créer le barplot empilé pour l'abondance moyenne
  tax_plot_mean <- ggplot(df_tax_mean, aes(x = Condition, y = Proportion, fill = !!sym(tax_level))) +
    geom_bar(stat = "identity", position = "stack") +
    theme_classic2() +
    scale_fill_manual(values = palette, name = tax_level) +
    labs(
      x = "Condition",
      y = "Relative Abundance (%)",
      fill = tax_level,
      title = paste("Relative Abundance of", tax_level, "by Condition"),
      caption = "CLR dataset converted in pseudo-counts"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 12, face = "bold"),
      axis.text.y = element_text(size = 12, face = "bold"),
      panel.grid.major.y = element_line(color = "grey80", linewidth = 0.75),
      legend.position = "right",
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 12),
      plot.caption = element_text(size = 10, face = "bold"),
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold", color = "black")
    )
  
  # Enregistrer le barplot
  ggsave(
    filename = file.path(output_dir, paste0("relabundance_mean_", tax_level, "_conditions.png")),
    plot = tax_plot_mean, width = 10, height = 8, bg = "white", dpi = 300
  )
}


# Utiliser la fonction pour chaque niveau taxonomique
tax_levels <- colnames(ps.clr.raw.filt@tax_table)[2:length(colnames(ps.clr.raw.filt@tax_table))]
output <- "Figures/Taxonomic/composition"
dir.create(output, showWarnings = FALSE)

for (tax_level in tax_levels) {
  comp_taxa(ps_data = ps.clr.raw.filt, tax_level = tax_level, gp.colors = cond.colors, output_dir = output)
}


# Analyses différentielles

# Volcano & MA

# Move to the working directory
setwd("/Users/Maxime/Desktop/Projects/HB/")

# Import log-transformed table
df.clr <- read.csv("Data/output/df_final_aldex2_log_clr_glm.csv")
head(df.clr)

# Read the log-transformed ps object
ps.clr <- readRDS("Data/atg16_tuft/ps_atg16_coverage_alpha_corr_clr.rds")
print(ps.clr)
dim(otu_table(ps.clr))
head(ps.clr@otu_table)
head(ps.clr@tax_table)

# Spécifier l'ordre des conditions
cond <- c("Atg16wt", "Atg16ko", "Pou2f3ko", "Atg16ko; Pou2f3ko")

# Sélectionner les individus selon leur Condition
ps.filt <- subset_samples(ps.clr, Condition %in% cond)

# Modifier l'ordre des niveaux de facteur pour les Conditions
sample_data(ps.filt)$Condition <- factor(sample_data(ps.filt)$Condition, levels = cond)

# Conversion of tax_table into dataframe
tax_df <- as.data.frame(ps.filt@tax_table) %>%
  rownames_to_column(var = "msp")

# Merge taxonomies with clr-transformed abundance table
df_merged <- df.clr %>%
  dplyr::left_join(tax_df, by = "msp")

# Identification des colonnes ajoutées
add_cols <- setdiff(names(tax_df), "msp")

# Déplacement des colonnes ajoutées entre la 1ère et 2ème colonne de df.clr
df_merged <- df_merged %>%
  dplyr::relocate(all_of(add_cols), .after = 1)

# Fonction pour normaliser les noms de colonnes
norm_colnames <- function(colnames_vector) {
  colnames_vector %>%
    gsub("[[:space:]]+", "", .) %>%
    gsub("[.]+", "_", .) %>%
    gsub(";", "_", .)
}

# Remplacer ".." par "." dans les noms de colonnes de df_merged
colnames(df_merged) <- norm_colnames(colnames(df_merged))
head(df_merged)

# Définir les combinaisons
comb_define <- function(cond1, cond2, df) {
  # Construction des préfixes pour les deux directions
  prefix_1_vs_2 <- paste0("ref_", cond1, "_gp", cond2)
  prefix_2_vs_1 <- paste0("ref_", cond2, "_gp", cond1)
  
  # Identification des colonnes de significativité brute et corrigée pour chaque comparaison
  signif_raw_1_vs_2  <- paste0(prefix_1_vs_2, ".pval")
  signif_holm_1_vs_2 <- paste0(prefix_1_vs_2, ".pval.holm")
  
  signif_raw_2_vs_1  <- paste0(prefix_2_vs_1, ".pval")
  signif_holm_2_vs_1 <- paste0(prefix_2_vs_1, ".pval.holm")
  
  # Colonnes avec les mesures différentielles pour "cond1 vs cond2"
  diff_1_vs_2 <- c(
    diff_btw    = paste0(prefix_1_vs_2, ".diff.btw"),
    rab_all     = paste0(prefix_1_vs_2, ".rab.all"),
    effect      = paste0(prefix_1_vs_2, ".effect"),
    effect_low  = paste0(prefix_1_vs_2, ".effect.low"),
    effect_high = paste0(prefix_1_vs_2, ".effect.high"),
    dispersion  = paste0(prefix_1_vs_2, ".diff.win")
  )
  
  # Colonnes avec les mesures différentielles pour "cond2 vs cond1"
  diff_2_vs_1 <- c(
    diff_btw    = paste0(prefix_2_vs_1, ".diff.btw"),
    rab_all     = paste0(prefix_2_vs_1, ".rab.all"),
    effect      = paste0(prefix_2_vs_1, ".effect"),
    effect_low  = paste0(prefix_2_vs_1, ".effect.low"),
    effect_high = paste0(prefix_2_vs_1, ".effect.high"),
    dispersion  = paste0(prefix_2_vs_1, ".diff.win")
  )
  
  return(list(
    signif_raw_cond1_vs_cond2  = signif_raw_1_vs_2,
    signif_holm_cond1_vs_cond2 = signif_holm_1_vs_2,
    diff_cond1_vs_cond2        = diff_1_vs_2,
    
    signif_raw_cond2_vs_cond1  = signif_raw_2_vs_1,
    signif_holm_cond2_vs_cond1 = signif_holm_2_vs_1,
    diff_cond2_vs_cond1        = diff_2_vs_1
  ))
}

# Générer toutes les combinaisons à 2 conditions
comb <- combn(cond, 2) %>% t() %>% as.data.frame()
colnames(comb) <- c("Cond1", "Cond2")

# Pour chaque ligne, on applique la fonction comb_define
comb <- comb %>%
  mutate(Comparison = pmap(list(Cond1, Cond2), ~ comb_define(..1, ..2, df_merged)))

# Transformer en data frame long
comb <- comb %>%
  rowwise() %>%
  do({
    # Extraire les conditions et le résultat de comb_define
    cond1 <- .$Cond1
    cond2 <- .$Cond2
    comp  <- .$Comparison
    
    # Pour la direction Cond1 vs Cond2
    df1 <- tibble(
      Cond1 = cond1,
      Cond2 = cond2,
      direction = "Cond1_vs_Cond2",
      type = c("signif", "signif_holm", names(comp$diff_cond1_vs_cond2)),
      colname = c(comp$signif_raw_cond1_vs_cond2, comp$signif_holm_cond1_vs_cond2, comp$diff_cond1_vs_cond2)
    )
    
    # Pour la direction Cond2 vs Cond1
    df2 <- tibble(
      Cond1 = cond1,
      Cond2 = cond2,
      direction = "Cond2_vs_Cond1",
      type = c("signif", "signif_holm", names(comp$diff_cond2_vs_cond1)),
      colname = c(comp$signif_raw_cond2_vs_cond1, comp$signif_holm_cond2_vs_cond1, comp$diff_cond2_vs_cond1)
    )
    
    bind_rows(df1, df2)
  }) %>%
  ungroup()

# Normalisation des noms de colonnes
comb <- comb %>% mutate(colname = norm_colnames(colname))
print(comb)

# Vérifier si toutes les colonnes listées dans "comb$colname" existent dans df_merged
missing_cols <- setdiff(comb$colname, colnames(df_merged))

# Afficher les colonnes manquantes
if (length(missing_cols) == 0) {
  message("Toutes les colonnes existent dans df_merged.")
} else {
  message("Colonnes manquantes dans df_merged :")
  print(missing_cols)
}

# Supprimer les lignes du tableau long contenant des colonnes manquantes
df.diff <- comb %>% 
  filter(colname %in% colnames(df_merged))
print(df.diff)
print(df.diff$colname)

# FONCTION
perform_dfa <- function(mapping, data, output_prefix, plot_title,
                        signif_types = c("signif", "signif_holm"),
                        ntop_signif = 10,
                        ntop_rab = 5,
                        direction = "Cond1_vs_Cond2",
                        thresholds_neg = 1.5,
                        thresholds_pos = 1.5,
                        thresholds_pvalue = c(0.05, 0.01),
                        size_pt = c(2.5, 1),
                        size_lab = 3.5,
                        colors_pt = c("grey50", "grey70", "#6598b1", "#e2003f"),
                        colors_lab = c("#1D4995", "#a20101", "grey50"),
                        labels = NULL,
                        perform_volcano = TRUE,
                        perform_MA = TRUE,
                        perform_MW = TRUE,
                        dim_volcano = c(9, 7),
                        dim_MA = c(9, 7),
                        dim_MW = c(9, 7)) {
  
  # Filtrer le mapping selon la direction et fixer l'ordre d'apparition
  mapping_filt <- mapping %>%
    dplyr::filter(direction == !!direction) %>%
    dplyr::mutate(Cond1 = factor(Cond1, levels = unique(Cond1)),
                  Cond2 = factor(Cond2, levels = unique(Cond2))) %>%
    dplyr::arrange(Cond1, Cond2)
  
  comparisons <- mapping_filt %>%
    dplyr::group_by(Cond1, Cond2, direction) %>%
    dplyr::group_split()
  
  # Calcul des seuils
  log2fc_threshold_neg <- -log2(thresholds_neg)
  log2fc_threshold_pos <- log2(thresholds_pos)
  pvalue_threshold     <- -log10(thresholds_pvalue[1])
  pvalue_threshold2    <- -log10(thresholds_pvalue[2])
  
  # Extraire les métadonnées (colonnes dont le nom ne commence pas par "ref")
  metadata <- data %>%
    dplyr::select(msp, dplyr::everything()[!grepl("^ref", colnames(data))])
  
  # Boucle sur chaque comparaison
  for(comp in comparisons) {
    comp_id <- paste(unique(comp$Cond1), unique(comp$Cond2), sep = "_vs_")
    message("Comparaison: ", comp_id)
    
    # Groupes pour la légende
    ref_group <- as.character(unique(comp$Cond1))[1]
    comp_group <- as.character(unique(comp$Cond2))[1]
    
    # Extraire les colonnes numériques
    diff_col <- as.character((comp %>% dplyr::filter(type == "diff_btw") %>% dplyr::pull(colname))[1])
    rab_col  <- as.character((comp %>% dplyr::filter(type == "rab_all") %>% dplyr::pull(colname))[1])
    disp_col <- as.character((comp %>% dplyr::filter(type %in% c("dispersion", "diff_win")) %>% dplyr::pull(colname))[1])
    message("diff_col: ", diff_col)
    message("rab_col: ", rab_col)
    message("disp_col: ", disp_col)
    
    # Créer le data frame de base
    base_dp <- data %>%
      dplyr::select(msp,
                    diff_btw = !!sym(diff_col),
                    rab_all  = !!sym(rab_col),
                    diff_win = !!sym(disp_col))
    
    dp <- base_dp %>%
      dplyr::left_join(metadata, by = "msp")
    
    # Boucle sur chaque type de p-value
    for(pval_type in signif_types) {
      pval_col <- as.character((comp %>% dplyr::filter(type == pval_type) %>% dplyr::pull(colname))[1])
      if(is.na(pval_col) || !(pval_col %in% colnames(data))) {
        warning(paste("La colonne", pval_col, "n'existe pas dans data pour la comparaison", comp_id))
        next
      }
      
      dp <- dp %>%
        dplyr::left_join(data %>% dplyr::select(msp, pval_temp = !!sym(pval_col)), by = "msp") %>%
        dplyr::mutate(pval = as.numeric(pval_temp)) %>%
        dplyr::select(-pval_temp)
      
      # Agrégation par 'labels' si fourni
      dp_aggr <- dp
      if(!is.null(labels)) {
        if(!(labels %in% colnames(dp))) {
          warning("La métadonnée '", labels, "' n'existe pas dans data.")
        } else {
          dp_aggr <- dp %>%
            dplyr::group_by(!!sym(labels)) %>%
            dplyr::summarize(across(c(diff_btw, rab_all, diff_win, pval), ~ mean(.x, na.rm = TRUE)),
                             msp = dplyr::first(msp),
                             .groups = "drop")
        }
      }
      dp_use <- dp_aggr
      
      # Pour Volcano et MW : définir la couleur et la couleur du label
      dp_use <- dp_use %>% dplyr::mutate(
        color = factor(case_when(
          pval <= thresholds_pvalue[1] & diff_btw <= log2fc_threshold_neg ~ colors_pt[3],
          pval <= thresholds_pvalue[1] & diff_btw >= log2fc_threshold_pos ~ colors_pt[4],
          TRUE ~ colors_pt[2]
        ), levels = c(colors_pt[3], colors_pt[4], colors_pt[2])),
        label_color = case_when(
          pval <= thresholds_pvalue[1] & diff_btw <= log2fc_threshold_neg ~ colors_lab[1],
          pval <= thresholds_pvalue[1] & diff_btw >= log2fc_threshold_pos ~ colors_lab[2],
          TRUE ~ colors_lab[3]
        )
      )
      
      # Sélectionner les points significatifs (pval <= seuil et diff_btw extrêmes)
      sig_dp <- dp_use %>% 
        dplyr::filter(pval <= thresholds_pvalue[1] &
                        (diff_btw <= log2fc_threshold_neg | diff_btw >= log2fc_threshold_pos))
      if(nrow(sig_dp) == 0) {
        warning(paste("Aucune fonction significative trouvée pour", pval_type, "dans", comp_id, ". Aucun label affiché."))
        all_labels <- dp_use[0, ]
      } else {
        neg_top <- sig_dp %>% 
          dplyr::filter(diff_btw <= log2fc_threshold_neg) %>% 
          dplyr::arrange(diff_btw) %>% 
          dplyr::slice_head(n = ntop_signif)
        pos_top <- sig_dp %>% 
          dplyr::filter(diff_btw >= log2fc_threshold_pos) %>% 
          dplyr::arrange(desc(diff_btw)) %>% 
          dplyr::slice_head(n = ntop_signif)
        all_labels <- dplyr::distinct(dplyr::bind_rows(neg_top, pos_top), msp, .keep_all = TRUE)
      }
      
      # Pour le MA plot : sélectionner les ntop_rab points les plus abondants
      if(nrow(sig_dp) > 0) {
        ma_labels <- sig_dp %>% 
          dplyr::arrange(desc(abs(rab_all))) %>% 
          dplyr::slice_head(n = ntop_rab)
      } else {
        ma_labels <- dp_use[0, ]
      }
      
      # Créer l'union des points annotés 
      annotated_points <- dplyr::distinct(dplyr::bind_rows(all_labels, ma_labels), msp, .keep_all = TRUE)
      
      # Utiliser 'labels' si fourni, sinon 'msp'
      if(!is.null(labels)) {
        dp_use <- dp_use %>% dplyr::mutate(
          label = ifelse(msp %in% annotated_points$msp, as.character(!!sym(labels)), "")
        )
        annotated_points <- annotated_points %>% dplyr::mutate(
          label = as.character(!!sym(labels))
        )
      } else {
        dp_use <- dp_use %>% dplyr::mutate(
          label = ifelse(msp %in% annotated_points$msp, msp, "")
        )
      }
      
      # Définir la taille des points en fonction de l'annotation :
      # points annotés = size_pt[1] et autres = size_pt[2]
      dp_use <- dp_use %>% mutate(point_size = ifelse(msp %in% annotated_points$msp, size_pt[1], size_pt[2]))
      
      ### Volcano Plot
      if(perform_volcano) {
        volcano_plot <- ggplot(dp_use, aes(x = diff_btw, y = -log10(pval), color = color, size = point_size)) +
          geom_point(alpha = 0.8) +
          scale_color_manual(
            values = setNames(c(colors_pt[3], colors_pt[4], colors_pt[2]),
                              c(colors_pt[3], colors_pt[4], colors_pt[2])),
            breaks = c(colors_pt[3], colors_pt[4]),
            labels = c(ref_group, comp_group),
            name = "Higher abundant in:"
          ) +
          scale_size_identity() +
          guides(color = guide_legend(override.aes = list(size = size_pt[1]))) +
          new_scale_color() +
          geom_text_repel(data = annotated_points, 
                          aes(x = diff_btw, y = -log10(pval), label = label, color = label_color),
                          size = size_lab, max.overlaps = Inf, alpha = 1, show.legend = FALSE) +
          scale_color_manual(values = c(colors_lab[1], colors_lab[2], colors_lab[3])) +
          theme_prism() +
          labs(title = plot_title, 
               subtitle = paste("Volcano plot -", comp_id),
               x = expression("Median Log"[2]*" Difference"),
               y = expression("-"*log[10]*"(pvalue)"),
               caption = paste("Top", ntop_signif+ntop_rab, "CLR-transformed data")) +
          geom_hline(yintercept = pvalue_threshold, linetype = "dashed", color = "darkred") +
          geom_hline(yintercept = pvalue_threshold2, linetype = "dashed", color = "#1ecb88") +
          geom_vline(xintercept = c(log2fc_threshold_neg, log2fc_threshold_pos), linetype = "dashed", color = "darkred") +
          annotate("rect", xmin = -Inf, xmax = log2fc_threshold_neg, 
                   ymin = pvalue_threshold2, ymax = Inf, fill = "#efd2d2", alpha = 0.3) +
          annotate("rect", xmin = log2fc_threshold_pos, xmax = Inf, 
                   ymin = pvalue_threshold2, ymax = Inf, fill = "#efd2d2", alpha = 0.3) +
          annotate("rect", xmin = -Inf, xmax = log2fc_threshold_neg, 
                   ymin = pvalue_threshold, ymax = pvalue_threshold2, fill = "#b9ecd1", alpha = 0.3) +
          annotate("rect", xmin = log2fc_threshold_pos, xmax = Inf, 
                   ymin = pvalue_threshold, ymax = pvalue_threshold2, fill = "#b9ecd1", alpha = 0.3) +
          theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 14, face = "bold"),
                axis.text.y = element_text(size = 14, face = "bold"),
                axis.title.x = element_text(size = 14, face = "bold"),
                axis.title.y = element_text(size = 14, face = "bold"),
                legend.position = "right",
                legend.title = element_text(size = 12, face = "bold"),
                legend.text = element_text(size = 12),
                plot.title = element_text(hjust = 0.5, size = 16, face = "bold", color = "black"),
                plot.caption = element_text(size = 10, face = "bold"))
        print(volcano_plot)
        volcano_filename <- paste0(output_prefix, comp_id, "_", pval_type, "_", labels, "_volcano.png")
        ggsave(volcano_filename, plot = volcano_plot, units = "in", width = dim_volcano[1], height = dim_volcano[2], dpi = 300, bg = "white")
      }
      
      ### MA Plot
      if(perform_MA) {
        # Créer un jeu de données pour le MA plot en ajoutant 'ma_color'
        dp_ma <- dp_use %>% mutate(ma_color = case_when(
          diff_btw <= log2fc_threshold_neg ~ colors_pt[3],
          diff_btw >= log2fc_threshold_pos ~ colors_pt[4],
          TRUE ~ colors_pt[2]
        ))
        # Assurer que les points annotés possèdent également 'ma_color'
        annotated_points <- annotated_points %>% mutate(ma_color = case_when(
          diff_btw <= log2fc_threshold_neg ~ colors_pt[3],
          diff_btw >= log2fc_threshold_pos ~ colors_pt[4],
          TRUE ~ colors_pt[2]
        ))
        
        ma_plot <- ggplot(dp_ma, aes(x = abs(rab_all), y = diff_btw, color = ma_color, size = point_size)) +
          geom_point(alpha = 0.8) +
          scale_color_manual(
            values = setNames(c(colors_pt[3], colors_pt[4], colors_pt[2]),
                              c(colors_pt[3], colors_pt[4], colors_pt[2])),
            breaks = c(colors_pt[3], colors_pt[4]),
            labels = c(ref_group, comp_group),
            name = "Higher abundant in:"
          ) +
          scale_size_identity() +
          guides(color = guide_legend(override.aes = list(size = size_pt[1]))) +
          new_scale_color() +
          geom_text_repel(data = annotated_points, 
                          aes(x = abs(rab_all), y = diff_btw, label = label, color = label_color),
                          size = size_lab, max.overlaps = Inf, alpha = 1, show.legend = FALSE) +
          scale_color_identity() +
          theme_prism() +
          labs(title = plot_title,
               subtitle = paste("MA Plot -", comp_id),
               x = expression("Median Log"[2]*" relative abundance"),
               y = expression("Median Log"[2]*" Difference"),
               caption = paste("Top", ntop_signif+ntop_rab, "CLR-transformed data"),
               color = "Significance\n(p < 0.05)") +
          geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
          geom_hline(yintercept = c(log2fc_threshold_neg, log2fc_threshold_pos), linetype = "dashed", color = "darkred") +
          theme(axis.text.x = element_text(size = 14, face = "bold"),
                axis.text.y = element_text(size = 14, face = "bold"),
                axis.title.x = element_text(size = 14, face = "bold"),
                axis.title.y = element_text(size = 14, face = "bold"),
                legend.position = "right",
                legend.title = element_text(size = 12, face = "bold"),
                legend.text = element_text(size = 12),
                plot.title = element_text(hjust = 0.5, size = 16, face = "bold", color = "black"),
                plot.caption = element_text(size = 10, face = "bold"))
        print(ma_plot)
        ma_filename <- paste0(output_prefix, comp_id, "_", pval_type, "_", labels, "_MA.png")
        ggsave(ma_filename, plot = ma_plot, units = "in", width = dim_MA[1], height = dim_MA[2], dpi = 300, bg = "white")
      }
      
      ### MW Plot
      if(perform_MW) {
        if(!"diff_win" %in% colnames(dp_use)) {
          warning("La colonne 'diff_win' est introuvable pour la comparaison ", comp_id, "... MW plot ignoré!")
        } else {
          mw_plot <- ggplot(dp_use, aes(x = diff_win, y = diff_btw, color = color, size = point_size)) +
            geom_point(alpha = 0.8) +
            scale_color_manual(
              values = setNames(c(colors_pt[3], colors_pt[4], colors_pt[2]),
                                c(colors_pt[3], colors_pt[4], colors_pt[2])),
              breaks = c(colors_pt[3], colors_pt[4]),
              labels = c(ref_group, comp_group),
              name = "Higher abundant in:"
            ) +
            scale_size_identity() +
            guides(color = guide_legend(override.aes = list(size = size_pt[1]))) +
            new_scale_color() +
            geom_text_repel(data = annotated_points, 
                            aes(x = diff_win, y = diff_btw, label = label, color = label_color),
                            size = size_lab, max.overlaps = Inf, alpha = 1, show.legend = FALSE) +
            scale_color_identity() +
            theme_prism() +
            labs(title = plot_title, 
                 subtitle = paste("MW Plot -", comp_id),
                 x = expression("Median Log"[2]*" Dispersion"),
                 y = expression("Median Log"[2]*" Difference"),
                 caption = paste("Top", ntop_signif+ntop_rab, "CLR-transformed data")) +
            geom_hline(yintercept = c(log2fc_threshold_neg, log2fc_threshold_pos), linetype = "dashed", color = "darkred") +
            geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
            geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
            geom_abline(intercept = 0, slope = -1, linetype = "dashed", color = "grey50") +
            theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 14, face = "bold"),
                  axis.text.y = element_text(size = 14, face = "bold"),
                  axis.title.x = element_text(size = 14, face = "bold"),
                  axis.title.y = element_text(size = 14, face = "bold"),
                  legend.position = "right",
                  legend.title = element_text(size = 12, face = "bold"),
                  legend.text = element_text(size = 12),
                  plot.title = element_text(hjust = 0.5, size = 16, face = "bold", color = "black"),
                  plot.caption = element_text(size = 10, face = "bold"))
          print(mw_plot)
          mw_filename <- paste0(output_prefix, comp_id, "_", pval_type, "_", labels, "_MW.png")
          ggsave(mw_filename, plot = mw_plot, units = "in", width = dim_MW[1], height = dim_MW[2], dpi = 300, bg = "white")
        }
      }
    }
  }
}


# Utilisation de la fonction :
# - 'df.diff' est le tableau de correspondance (mapping)
# - 'df_merged' contient les valeurs (data)
signif.cols <- c("signif", "signif_holm")

perform_dfa(
  mapping = df.diff,
  data = df_merged,
  output_prefix = "Figures/Taxonomic/DAA/",
  plot_title = "Differential Abundance Analysis",
  signif_types = signif.cols,
  ntop_signif = 15,
  ntop_rab = 5,
  labels = "species",
  direction = "Cond1_vs_Cond2",
  thresholds_neg = 1.5,
  thresholds_pos = 1.5,
  thresholds_pvalue = c(0.05, 0.01),
  size_pt = c(2.5, 1),
  size_lab = 3.5,
  colors_pt = c("grey50", "grey70", "#6598b1", "#e2003f"),
  colors_lab = c("#1D4995", "#a20101", "grey50"),
  perform_volcano = TRUE,
  perform_MA = TRUE,
  perform_MW = TRUE,
  dim_volcano = c(9, 7),
  dim_MA = c(9, 7),
  dim_MW = c(9, 7)
)

