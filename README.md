# Taxonomy Plot Generator for Metagenomic Samples
This GitHub project provides a tool to generate taxonomy abundance plots for metagenomic samples using the common output format (BIOM file) from popular taxonomic assigners such as Kraken, Bracken, and more. With this tool, you can visualize the composition of microbial communities at different taxonomic levels, including Phylum, Class, Order, Family, Genus, or Species.

## Requirements
- R v3.6

### Libraries
- BiocManager v1.30.22
- phyloseq v1.30.0
- ggplot2 v3.4.4
- RColorBrewer v1.1.3
- patchwork v1.1.3
- pals v1.8
- optparse v1.7.3
- dplyr v1.1.3

### Installation
```R
install.packages("BiocManager", version = "1.30.22")
BiocManager::install("phyloseq", version = "1.30.0")
install.packages("ggplot2", version = "3.4.4")
install.packages("RColorBrewer", version = "1.1.3")
install.packages("patchwork", version = "1.1.3")
install.packages("pals", version = "1.8")
install.packages("optparse", version = "1.7.3")
install.packages("dplyr", version = "1.1.3")
```


## Usage Options

The script accepts the following command-line options:

- `-i INPUT, --input`: Path to the BIOM file.
- `-o OUTPUT, --output`: Output directory where the plots will be saved.
- `-l LEVEL, --level`: Taxonomy level at which the plot will be generated (Phylum, Class, Order, Family, Genus, or Species).
- `-t TOP, --top`: Specify the number of top microbes to plot (according to the selected taxonomy level).
- `-c COLORS, --colors`: Specify the number of colors to use.
- `-h, --help`: Show this help message.

Make sure to specify the required options and provide values for each when running the script. For example:

```bash
Rscript taxonomy_abundance_plot.R -i PATH/TO/BIOM -o OUPUT_DIR/ --level Order --top 29 --colors 20
```