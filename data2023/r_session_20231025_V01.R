##### 2.1 Loading libraries and microbiome data #####

# load ggplot2 library (graphics)
library(ggplot2)

# loading phyloseq library (microbiome analysis)
library(phyloseq)

# OTU data
InputBiomFile <- "~/data2023/mbtmicrobiome2023.biom"

# Samples' metadata
InputMapFile <- "~/data2023/sample-metadata-2023.tsv"

# prepare phyloseq object by loading both files
BiomData <- import_biom(InputBiomFile, parseFunction = parse_taxonomy_greengenes)
SampleData <- import_qiime_sample_data(InputMapFile)

# create phyloseq object by merging OTU and sample data
ExperimentPhyloseqObject <- merge_phyloseq(BiomData, SampleData)

# checking the features of our original microbiome data
ExperimentPhyloseqObject

# create a temporary phyloseq object for working
psTemp <- ExperimentPhyloseqObject

# subset samples
# Prune OTUs with low abundances from all samples 
psTemp <- prune_taxa(taxa_sums(psTemp) > 100, psTemp)

# Prune samples with no metadata
psTemp <-  subset_samples(psTemp, Gender != "NA") 

# checking the features of filtered our microbiome data
psTemp

# Q: What are the differences between the "ExperimentPhyloseqObject" and "psTemp" objects?


##### 2.2 Sample ordination #####

# Calculate distance and ordination
iDist <- distance(psTemp, method="bray")
iMDS  <- ordinate(psTemp, distance=iDist)

# plot sample ordination, e.g. by Gender
plot_ordination(psTemp, iMDS, color="Gender")
# Q: Are there any clear separations between the gender groups? Why/Why not?

# repeat the ordination plot using other metadata information (e.g. cat yes/no)
plot_ordination(psTemp, iMDS, color="Cat") + 
  geom_text(aes(label=SampleID), vjust = -1)

# Q: Are there any clear separations between the plotted groups? Why/Why not?


##### 2.3 Microbial communities #####

### Abundances per sample
plot_bar(psTemp, "SampleID", fill="Phylum")

# Rarefaction to an even depth
ps.rarefied <- rarefy_even_depth(psTemp)

# Remove lines
ps.rarefied.glom <- tax_glom(ps.rarefied, "Phylum")

# Plot abundances
plot_bar(ps.rarefied.glom, "SampleID", fill="Phylum")
# Q: What are the differences between the abundance plots before and after rarefaction?

### Abundances per category, e.g. "Gender"
# Merge samples by the selected category
mergedGP <- merge_samples(psTemp, "Gender")

# Rarefaction to an even depth
ps.rarefied <- rarefy_even_depth(mergedGP)

# Remove lines
ps.rarefied.glom <- tax_glom(ps.rarefied, "Phylum")

# Plot abundances for the example category "Gender"
plot_bar(ps.rarefied.glom, fill="Phylum")


##### 2.4 Diversity #####

# Observe the Shannon diversity on the individual samples
plot_richness(psTemp, x = "SampleID", measures = c("Shannon")) 

# Repeat the same analyses at the Gender level
plot_richness(psTemp, x = "Gender", color = "Gender", measures = c("Shannon")) 

# Assign our plot to a variable
Our_Richness_plot <- plot_richness(psTemp, x = "Gender", color = "Gender", measures = c("Shannon")) 

# Improving our plot e.g. by adding proper labels
Our_Richness_plot + geom_boxplot(data = Our_Richness_plot$data, aes(x = Gender, y = value, color = Gender), alpha = 0.1) + # boxplot
  labs(title = "Richness (Shannon alpha diversity)", subtitle = "MBT Class 2023") + # title and subtitle
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) # x-axis labels: 0 degree rotation, 0.5 horizontal position
