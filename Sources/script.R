################################################################################
#####################        Script R 10/12/2024         #######################
#####################     Developed by Maxime Naour      #######################
################################################################################

# Loading libraries
.libPaths(c("C:/Users/Maxime/AppData/Local/Programs/R/R-4.4.1/library", .libPaths()))
library(devtools)
library(BiocManager)
library(igraph); packageVersion("igraph")
library(readxl); packageVersion("readxl")
library(openxlsx); packageVersion("openxlsx")
library(writexl); packageVersion("writexl")
library(jsonlite); packageVersion("jsonlite")
library(microbiome); packageVersion("microbiome")
library(microbiomeutilities); packageVersion("microbiomeutilities")
library(dplyr); packageVersion("dplyr")
library(vegan); packageVersion("vegan")
library(ggplot2); packageVersion("ggplot2")
library(ggforce); packageVersion("ggforce")
library(ggnewscale); packageVersion("ggnewscale")
library(phyloseq); packageVersion("phyloseq")
library(phyloseq.extended); packageVersion("phyloseq.extended")
#library(biomformat); packageVersion("biomformat")
library(DESeq2); packageVersion("DESeq2")
#library(tibble); packageVersion("tibble")  
#library(plyr); packageVersion("plyr")
library(gcookbook); packageVersion("gcookbook")
library(tidyverse); packageVersion("tidyverse")
library(viridis); packageVersion("viridis")
library(gdata); packageVersion("gdata")
#library(data.table); packageVersion("data.table")
#library(ade4); packageVersion("ade4")
library(effsize); packageVersion("effsize")
library(ComplexHeatmap); packageVersion("ComplexHeatmap")
library(circlize); packageVersion("circlize")
#library(RColorBrewer); packageVersion("RColorBrewer")
#library(cluster); packageVersion("cluster")
library(factoextra); packageVersion("factoextra")
library(gridExtra); packageVersion("gridExtra")
library(magick); packageVersion("magick")
library(fpc); packageVersion("fpc")
library(mixOmics); packageVersion("mixOmics") # devtools::install_github("mixOmicsTeam/mixOmics")
library(ropls); packageVersion("ropls")
library(ggtree); packageVersion("ggtree")
library(ALDEx2, lib.loc = "C:/Users/Maxime/AppData/Local/R/win-library/4.2"); packageVersion("ALDEx2") #BiocManager::install("ALDEx2") & install.packages("MASS", lib = "C:/Users/Maxime/AppData/Local/R/win-library/4.2")

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

# Préparer le dataframe pour ggplot
df_melt <- df %>%
  pivot_longer(-c(samples, ttt), names_to = "step", values_to = "reads") %>%
  arrange(ttt, samples)

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

# Générer le barplot
process.p <- ggplot(df_melt, aes(x = samples, y = reads, fill = step)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
  scale_fill_brewer(palette = "Set2", name = "Step") +
  ggnewscale::new_scale_fill() +
  geom_bar(aes(fill = ttt), stat = "identity", alpha = 0, show.legend = TRUE) +
  labs(
    x = "", y = "Reads count",
    fill = "Step",
    title = "Read counts at different processing steps",
    color = "Condition"
  ) +
  scale_fill_manual(
    name = "Condition",
    values = cond.colors,
    guide = guide_legend(override.aes = list(alpha = 1))
  ) +
  scale_y_continuous(limits = c(min_reads, max_reads), oob = scales::squish) +
  scale_x_discrete(labels = function(x) {
    sapply(x, function(sample) {
      sprintf('<span style="color:%s;">%s</span>', cond.colors[as.character(df[df$samples == sample, "ttt"])], sample)
    })
  }) +
  theme_classic() +
  theme(
    axis.text.x = ggtext::element_markdown(angle = 45, vjust = 1, hjust = 1, face = "bold", size = 12),
    axis.text.y = element_text(face = "bold", size = 12),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.title = element_text(face = "bold", size = 14)
  )
print(process.p)

# Enregistrer le graphique
ggsave("Figures/Taxonomic/processing/overview_processing_atg16.png", 
       process.p, width = 14, height = 10, dpi = 300, bg = "white")

# Calcul des moyennes et écarts-types par Condition
df_summary <- df_melt %>%
  group_by(ttt, step) %>%
  summarise(mean = mean(reads, na.rm = TRUE),
            sd = sd(reads, na.rm = TRUE),
            n = n(), .groups = 'drop')

# Identifier le nombre minimum et maximum des moyennes pour ajuster l'axe y
min_reads <- min(df_melt$reads, na.rm = TRUE)
max_reads <- max(df_melt$reads, na.rm = TRUE)

# Ajuster la limite inférieure pour inclure la plus petite valeur
min_reads <- min_reads - (0.05 * (max_reads - min_reads))
if (min_reads < 0) min_reads <- 0

# Création du barplot par condition
process.bp <- ggplot(df_summary, aes(x = ttt, y = mean, fill = step)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), 
                width = 0.2, 
                position = position_dodge(0.7)) +
  labs(
    x = "", y = "Reads count", fill = "Step", 
    title = "Average read counts at different processing steps by condition"
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
    axis.title = element_text(face = "bold", size = 14)
  )
print(process.bp)

# Enregistrer le graphique
ggsave("Figures/Taxonomic/processing/overview_processing_atg16_summary.png", 
       process.bp, width = 10, height = 8, dpi = 300, bg = "white")


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

# Fonction pour les graphiques avec annotations
rich.alpha <- function(df, measure, title, colors) {
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
  comp <- list(
    c(levels(df$Condition_label)[1], levels(df$Condition_label)[2]),
    c(levels(df$Condition_label)[1], levels(df$Condition_label)[3]),
    c(levels(df$Condition_label)[2], levels(df$Condition_label)[4]),
    c(levels(df$Condition_label)[3], levels(df$Condition_label)[4])
  )

  test_result <- kruskal.test(as.formula(paste(measure, "~ Condition")), data = df)
  caption_text <- paste0(
    "\nKruskal-Wallis p-value: ", sprintf("%.3f", test_result$p.value)
  )
  p <- ggplot(df, aes(x = Condition_label, y = !!rlang::sym(measure), fill = Condition)) +
    geom_violin(trim = TRUE, scale = "width", alpha = 0.7) +
    geom_boxplot(width = 0.2, outlier.shape = NA, color = "black") +
    geom_jitter(width = 0.2, size = 1.5, color = "black") +
    labs(title = title, y = measure, x = "", caption = caption_text) +
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
  return(p)
}

# Générer les plots
rich.plots <- lapply(c("Observed", "Shannon", "Simpson", "InvSimpson", "Fisher"), function(measure) {
  rich.alpha(div.df, measure, measure, cond.colors)
})

# Combiner les plots
comb <- do.call(grid.arrange, c(rich.plots, ncol = 3))

# Sauvegarder les plots combinés
ggsave("Figures/Taxonomic/div/alpha_diversity_atg16_coverage.png", 
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
metadata$gp <- as.factor(metadata$Condition)

# Créer une matrice de modèle pour inclure toutes les comparaisons
mm <- model.matrix(~ 0 + gp, metadata)
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
