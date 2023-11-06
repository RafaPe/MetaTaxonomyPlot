#!/usr/bin/Rscript

library("phyloseq")
library("ggplot2")
library("RColorBrewer")
library("patchwork")
library("pals")
library("glue")
library("optparse")
# library("dplyr")
suppressMessages(library(dplyr))

get_most_abundant <- function(df, num, tax_rank)
{
    # An auxiliary df is created to perform grouping by taxonomic levels and calculate the total abundance of their taxa.
    # The df is then sorted based on the total abundance of the taxa.
    df_aux <- df[, c(tax_rank, "Abundance")]
    df_aux <- df_aux %>%
        group_by(!!sym(tax_rank)) %>%
            summarize(Total_Abundance = sum(Abundance)) %>% arrange(desc(Total_Abundance))

    # Select top most abundant
    most_abundant <- head(df_aux, n = num)

    # The orignal df is filtered by most abundant taxa selected before
    filtered_df <- df[df[[tax_rank]] %in% most_abundant[[tax_rank]], ]
    # The not most abundant taxa are stored in unfilterd_df
    unfiltered_df <- df[!(df[[tax_rank]] %in% most_abundant[[tax_rank]]), ]

    # The unfiltered df is groupped by sample name in the unfiltered_df and sum abundances
    df_aux3 <- unfiltered_df[, c("Sample", "Abundance")]
    df_aux4 <- df_aux3 %>%
        group_by(Sample) %>%
        summarize(Total_Abundance = sum(Abundance)) %>% arrange(desc(Total_Abundance))
      
    # We add a new entry to the filtered_df for each sample with its abundance sum
    for (sample in df_aux4$Sample)
    {
        new_row <- data.frame(OTU = 0, Sample = sample, Abundance = df_aux4$Total_Abundance[df_aux4$Sample == sample], Id = sample)
        num_cols <- ncol(filtered_df)

        if (num_cols > 4)
        {
            new_row <- cbind(new_row, matrix('Other', nrow = 1, ncol = num_cols - 4))
        }

        colnames(new_row) <- colnames(filtered_df)
        filtered_final <- rbind(filtered_df, new_row)
    }

    # The filtered df with the new entries is returned
    return (filtered_final)
}

plot_taxonomy <- function(merged_metagenomes, tax_level, n_top, n_colors, relative = FALSE)
{
    # The phyloseq object is trimmed to remove data prefixes
    merged_metagenomes@tax_table@.Data <- substring(merged_metagenomes@tax_table@.Data, 4)

    # The names of the columns are adjusted
    colnames(merged_metagenomes@tax_table@.Data)<- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

    # Only the bacterial taxons are kept
    merged_metagenomes <- subset_taxa(merged_metagenomes, Kingdom == "Bacteria")

    # The entries where the specified taxonomic level was not identified are dropped
    if (tax_level == "Phylum") 
    {
        merged_metagenomes <- subset_taxa(merged_metagenomes, Phylum != "")
    } 
    else if (tax_level == "Class") 
    {
        merged_metagenomes <- subset_taxa(merged_metagenomes, Class != "")
    } 
    else if (tax_level == "Order") 
    {
        merged_metagenomes <- subset_taxa(merged_metagenomes, Order != "")
    } 
    else if (tax_level == "Family") 
    {
        merged_metagenomes <- subset_taxa(merged_metagenomes, Family != "")
    } 
    else if (tax_level == "Genus") 
    {
        merged_metagenomes <- subset_taxa(merged_metagenomes, Genus != "")
    } 
    else if (tax_level == "Species") 
    {
        merged_metagenomes <- subset_taxa(merged_metagenomes, Species != "")
    }

    #If the relative abundance plot is required, we normalice the abundance.
    if (!relative)
    {
        glom <- tax_glom(physeq = merged_metagenomes, taxrank = tax_level)
        df <- psmelt(glom)
        if (tax_level == 'Species')
        {
            df$Species <- paste(df$Genus, df$Species, sep = " ") #The genus and species name is concatenated
        }
    }
    else
    {
        percentages <- transform_sample_counts(merged_metagenomes, function(x) x*100 / sum(x) )
        glom <- tax_glom(percentages, taxrank = tax_level)
        df <- psmelt(glom)
        if (tax_level == 'Species')
        {
            df$Species <- paste(df$Genus, df$Species, sep = " ") #The genus and species name is concatenated
        }
    }

    # The most abundant taxa is selected
    df <- get_most_abundant(df, n_top, tax_level)

    # The color palette is generated
    colors<- colorRampPalette(alphabet(n_colors)) (length(unique(df[[tax_level]])))

    # We generate the stacked bar plot
    plot <- ggplot(data= df, aes(x=Sample, y=Abundance, fill=!!sym(tax_level)))+geom_bar(aes(), stat="identity", position="stack")+scale_fill_manual(values = colors)

    return (plot)
}

# Definition of arguments to call from the command line
option_list <- list(
    make_option(
        c("--input", "-i"),
        type = "character",
        default=NULL,
        help = "Path to BIOM file"
    ),
    make_option(
        c("--output", "-o"),
        type = "character",
        default=NULL,
        help = "Output directory where the plots will be saved"
    ),
    make_option(
        c("--level", "-l"),
        type = "character",
        default="Phylum",
        help = "Taxonomy level at which the plot will be generated (Phylum, Class, Order, Family, Genus or Species)"
    ),
    make_option(
        c("--top", "-t"),
        type = "numeric",
        default=9,
        help = "Specify number of top species to plot (according to selected taxonomy level)"
    ),
    make_option(
        c("--colors", "-c"),
        type = "numeric",
        default=7,
        help = "Specify number of colors to use"
    )
)

#Argument parser
opt_parser <- OptionParser(usage = "Usage: script.R [options]", option_list = option_list)
args <- commandArgs(trailingOnly = TRUE)
opt <- parse_args(opt_parser, args = args)

#Throwing error in case no input was provided
if (is.null(opt$input))
{
    stop("Input not specified, use -i flag to define the path to your BIOM file")
}


# MAIN FUNCTION
merged_metagenomes <- import_biom(opt$input)

abs <- plot_taxonomy(merged_metagenomes, tax_level = opt$level, n_top =  opt$top, n_colors = opt$colors)
rel <- plot_taxonomy(merged_metagenomes, tax_level = opt$level, n_top =  opt$top, n_colors = opt$colors, relative = TRUE)

ggsave(paste(opt$output, "AbsoluteAbundance.png"), plot = abs, width = 4500, height = 2950, units = "px", dpi = 400)
ggsave(paste(opt$output, "RelativeAbundance.png"), plot = rel, width = 4500, height = 2950, units = "px", dpi = 400)