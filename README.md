# UMR MBT Microbiome Praktikum

Scripts and tutorials for analyzing microbiome data. Lab practice for the lecture: <br>[Moderne molekulare und Hochdurchsatz-Technologien in der medizinischen Grundlagenforschung](https://lsf.uni-rostock.de/qisserver/rds?state=verpublish&status=init&vmfile=no&moduleCall=webInfo&publishConfFile=webInfo&publishSubDir=veranstaltung&navigationPosition=lectures%2Csearch&breadcrumb=searchLectures&topitem=lectures&subitem=search&veranstaltung.veranstid=132915)

---

## Description

These sessions will cover the use of a variety of software tools needed for the analysis of microbiome data, from the handling of the Illumina sequencing data, to the processing of 16S rRNA amplicon data. During these sessions, the students will be able to:

* Evaluate the quality of an Illumina sequencing run, including data filtering;
* Carry out assemblies of 16S rRNA amplicons;
* Assign Illumina runs to OTUs from 16S rRNA databases;
* Learn the basics on the R programming language and environment; and 
* State how to manipulate microbiome data including count tables and sample metadata. 

##### Venues and Dates

<!--

- Session 1, Setup: [Schillingallee 69a, 70 - SR1, Med. Th. Institut](https://lsf.uni-rostock.de/qisserver/rds?state=verpublish&status=init&vmfile=no&moduleCall=webInfo&publishConfFile=webInfoRaum&publishSubDir=raum&keep=y&raum.rgid=2114); 29.9.2021, 13 hrs. 	
- Sessions 2-3, Data analysis: [IBIMA](https://ibima.med.uni-rostock.de) Computerraum ([Ernst-Heydemann-Str. 8](https://goo.gl/maps/JGDWhPDLHxG2), 3. Etage, Nr. 3016); 13.10 and 20.10.2021, 10:00 - 11:30 hrs.

-->

- Sessions 1-2: [IBIMA](https://ibima.med.uni-rostock.de) Computerraum ([Ernst-Heydemann-Str. 8](https://goo.gl/maps/JGDWhPDLHxG2), 3. Etage, Nr. 3016); 12.10 and 19.10.2022, 10:00 - 11:30 hrs.

##### Presentations

- Slides: [`MBTPraktikum.2021.V02.pdf`](https://drive.google.com/file/d/1llSavTsSPvWGhzqsqVciQ8zzlLMfkPm_/view?usp=sharing)

##### Software

- All required software, packages and data are accessible through our the virtual [Binder](https://mybinder.org/v2/gh/barrantesisrael/mbt.microbiome.2021/main?urlpath=rstudio) environment



---

## Session 1

Launch our interactive course here: [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/barrantesisrael/mbt.microbiome.2021/main?urlpath=rstudio)

##### 1.1 Introduction to the command line for Bioinformatics and Quality control with the [FASTQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) tool


```bash
# try basic Linux/OSX commands, such as
$ ls
$ cd data2022

# observe the first ten lines -what do you notice?
$ head Platz1_R1.head.fastq

# count the total number of lines with the command below -how many READS are in this FASTQ?
# hint: each read in FASTQ format consists of four lines
$ wc -l Platz1_R1.head.fastq

# find specific words, e.g "AATATT"
$ grep "AATATT" Platz1_R1.head.fastq

# combine commands. e.g.
$ grep "AATATT" Platz1_R1.head.fastq | wc -l

# quality check with fastqc
$ fastqc --quiet Platz1_R1.head.fastq
# click to HTML output and select "View in Web Browser"
# a new tab opens with the fastqc results
# repeat the analyses with the second sequencing pair, e.g. Platz1_R2.head.fastq
```


##### 1.2 Amplicon assembly with [pandaseq](https://github.com/neufeld/pandaseq)

Continue on the `Terminal` from the Binder session

```bash
# Current path: ~/data2022
# pandaseq options
$ pandaseq -h

# amplicon assembly with pandaseq
$ pandaseq -f Platz1_R1.head.fastq -r Platz1_R2.head.fastq -w Platz1.fa -g log.txt

# observe the first ten lines of the FASTA output
$ head Platz1.fa

# how many sequences are in this FASTA?
# hint: count the total number of HEADER line symbols (">") with any of the commands below 
$ grep ">" Platz1.fa | wc -l
$ grep -c ">" Platz1.fa

# Q: what is the rate of FASTA sequences vs FASTQ reads?
# and what does this tell us about our sequencing and assembly quality and efficiency?
```

##### 1.3 OTU Assignment with [kraken2](https://ccb.jhu.edu/software/kraken2/) against the 16S [Greengenes](https://greengenes.secondgenome.com/) database


```bash
# Current path: ~/data2022
# Note: Run this step together with your neighbor/colleague to avoid saturation of the storage server

# Run your samples against Greengenes with your own data
# the example here is the with the Platz10 run
# replace the FASTQ filenames with your own FASTQ names
$ kraken2 --db 16S_Greengenes_k2db --use-names --output output.txt --report report.txt --paired Platz10_R1.head.fastq Platz10_R2.head.fastq

# inspect your individual results (file: report.txt) within the RStudio window
# Q: What are the most predominant genera in your personal Illumina runs?
```

##### For Session 2:

- Aside of the regular attendance list, prepare a table with the following data: Platz-number, Gender (F/M), smoking (yes/no), pets (yes/no). We will visualize if these common confounders have a relationship with the microbiome data.
- The whole OTU table can be visually inspected here: [`mbtmicrobiome2022.tsv`](https://github.com/barrantesisrael/mbt.microbiome.2021/blob/main/data2022/mbtmicrobiome2022.tsv)


---

## Session 2

##### 2.1 Loading libraries and microbiome data 

```r
# load ggplot2 library (graphics)
library(ggplot2, quietly = TRUE)

# loading phyloseq library (microbiome analysis)
library(phyloseq, quietly = TRUE)

# OTU data
InputBiomFile <- "~/data2022/mbtmicrobiome2022.biom"

# Samples' metadata
InputMapFile <- "~/data2022/sample-metadata-2022.tsv"

# prepare phyloseq object by loading both files
BiomData <- import_biom(InputBiomFile, parseFunction = parse_taxonomy_greengenes)
SampleData <- import_qiime_sample_data(InputMapFile)

# create phyloseq object by merging OTU and sample data
ExperimentPhyloseqObject <- merge_phyloseq(BiomData, SampleData)

# create a temporary phyloseq object for working
psTemp <- ExperimentPhyloseqObject

# checking the features of our microbiome data
psTemp

# subset samples
# Prune OTUs with low abundances from all samples 
psTemp <- prune_taxa(taxa_sums(psTemp) > 100, psTemp)

# Prune samples with no metadata
psTemp <-  subset_samples(psTemp, Geschlecht != "ND") 
```

##### 2.2 Sample ordination

```r
# Calculate distance and ordination
iDist <- distance(psTemp, method="bray")
iMDS  <- ordinate(psTemp, distance=iDist)

# plot sample ordination
plot_ordination(psTemp, iMDS, color="Geschlecht")

# plot sample ordination, including labels
plot_ordination(psTemp, iMDS, color="Geschlecht") + 
  geom_text(aes(label=X.SampleID), hjust=0, vjust=0)
  
# repeat the ordination plot, using diet and smoking habits information
plot_ordination(psTemp, iMDS, color="Raucher") + 
  geom_text(aes(label=X.SampleID), vjust = -1) +
  stat_ellipse() # using default ellipse
```

##### 2.3 Microbial communities

```r
# Plot abundances
plot_bar(psTemp, "X.SampleID", fill="Phylum")

# Rarefaction to an even depth
ps.rarefied <- rarefy_even_depth(psTemp)

# Remove lines
ps.rarefied.glom <- tax_glom(ps.rarefied, "Phylum")

# Plot abundances
plot_bar(ps.rarefied.glom, "X.SampleID", fill="Phylum")

# Separate according to metadata
plot_bar(ps.rarefied.glom, "X.SampleID", fill="Phylum", facet_grid="Geschlecht")
plot_bar(ps.rarefied.glom, "X.SampleID", fill="Phylum", facet_grid="Geschlecht~Raucher")

### Merge samples by a category, e.g. "Raucher"
mergedGP <- merge_samples(psTemp, "Raucher")

# Rarefaction to an even depth
ps.rarefied <- rarefy_even_depth(mergedGP)

# Remove lines
ps.rarefied.glom <- tax_glom(ps.rarefied, "Phylum")

# Plot abundances for the example category "Raucher"
plot_bar(ps.rarefied.glom, fill="Phylum")
```


---

## Bibliography

* Afgan E et al. (2018). The Galaxy platform for accessible, reproducible and collaborative biomedical analyses: 2018 update. _Nucleic Acids Research_, 46(W1), W537-W544. https://doi.org/10.1093/nar/gky379
* Callahan BJ et al. (2016). Bioconductor workflow for microbiome data analysis: from raw reads to community analyses. _F1000Research_, 5, 1492. https://doi.org/10.12688/f1000research.8986.1
* Cole JR et al. (2014). Ribosomal Database Project: data and tools for high throughput rRNA analysis. _Nucleic Acids Research_, 42(D1), D633–D642. https://doi.org/10.1093/nar/gkt1244
* Comeau AM et al. (2017). Microbiome Helper: a Custom and Streamlined Workflow for Microbiome Research. _mSystems_, 2(1): e00127-16; DOI:10.1128/mSystems.00127-16
* Gentleman RC et al. (2004). Bioconductor: open software development for computational biology and bioinformatics. _Genome Biology_, 5(10), R80. https://doi.org/10.1186/gb-2004-5-10-r80
* Goodrich JKK et al. (2014). Conducting a Microbiome Study. _Cell_, 158(2), 250–262. https://doi.org/10.1016/j.cell.2014.06.037
* Knight R et al. (2018). Best practices for analysing microbiomes. _Nature Reviews Microbiology_, 16, 410–422. https://doi.org/10.1038/s41579-018-0029-9
* Kozich JJ, Westcott SL, Baxter NT, Highlander SK, Schloss PD. (2013) Development of a dual-index sequencing strategy and curation pipeline for analyzing amplicon sequence data on the MiSeq Illumina sequencing platform. _Applied and Environmental Microbiology_, 79(17):5112-20. http://doi.org/10.1128/AEM.01043-13
* Lu, J., Rincon, N., Wood, D.E. et al. Metagenome analysis using the Kraken software suite. Nat Protoc (2022). https://doi.org/10.1038/s41596-022-00738-y
* Masella AP et al. (2012). PANDAseq: Paired-end assembler for illumina sequences. _BMC Bioinformatics_, 13(1), 31. https://doi.org/10.1186/1471-2105-13-31
* McDonald D et al. (2012). The Biological Observation Matrix (BIOM) format or: how I learned to stop worrying and love the ome-ome. _GigaScience_, 1(1), 7. https://doi.org/10.1186/2047-217X-1-7
* McMurdie PJ and Holmes S. (2013). Phyloseq: An R Package for Reproducible Interactive Analysis and Graphics of Microbiome Census Data. _PLoS ONE_, 8(4), e61217. https://doi.org/10.1371/journal.pone.0061217

---

## Contact

Dr. rer. nat. Israel Barrantes <br>
Junior Research Group Translational Bioinformatics (head)<br>
Institute for Biostatistics and Informatics in Medicine and Ageing Research, Office 3017<br>
Rostock University Medical Center<br>
Ernst-Heydemann-Str. 8<br>
18057 Rostock, Germany<br>

Email: israel.barrantes[bei]uni-rostock.de

---
Last update 2022/09/30


