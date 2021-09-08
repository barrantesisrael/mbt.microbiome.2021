# mbt.microbiome.2021
Launch the MBT Binder here: [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/barrantesisrael/mbt.microbiome.2021/main?urlpath=rstudio)


Alternatively: https://mybinder.org/v2/gh/barrantesisrael/mbt.microbiome.2021/main?urlpath=rstudio

---

### Session 1 

Sequencing files: Quality control, assembly and OTU assignment

```bash
# observe the first ten lines -what do you notice?
$ head paired1.fastq

# count the total number of lines with the command below -how many READS are in this FASTQ?
$ wc -l paired1.fastq

# quality check with fastqc
$ fastqc --quiet paired1.fastq
# click to HTML output and select "View in Web Browser"
# a new tab opens with the fastqc results

# amplicon assembly with pandaseq
$ pandaseq -f paired1.fastq -r paired2.fastq -w output.fa -g log.txt

# observe the first ten lines of the FASTQ output
$ head output.fa

# count the total number of HEADER lines with the command below -how many FRAGMENTS are in this FASTA?
$ grep -c ">" output.fa

# copy first 10 lines from the FASTA and paste it into the RDP server
# create a new script and copy and paste the code below
```



---
