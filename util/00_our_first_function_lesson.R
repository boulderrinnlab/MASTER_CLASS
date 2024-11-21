#! Function Lesson !!

# This function has two parameters, and we have to let R know (x & y)
# function is inside { }, then we need a 'return' to get the answer
fun <- function(x, y) {
  ans <- x * y
  return(ans)
}
#let's try it 
fun(2,-4)
# Note that the object ans doesn't exist in our environment. 
# ans is a local variable that is erased after used. global variables are in our environment.
# It's good to remember the "scope" of vars in a function don't really exist outside the function.
# Finally note a function can only return one thing!

# Now it's important to document a function !
# Below is the common description style note " #' " annotation for function


#' A function to multiply two numbers
#'
#' @description 
#' This function will multiply the input values of X and Y
#' 
#' @param x one number you'd like to multiply
#' @param  y the other number you'd like to multiply
fun <- function(x, y) {
  ans <- x * y
  return(ans)
}

# Awesome we can now source this function and it's well documented !
# Let's do the same for our function to plot a avg and sd TPM values:

#' FUNCTION to plot TPM values for genes of interest
#'
#' @description 
#' This function will take a list of genes and plot the avg and sd
#' TPM values as a line plot.
#' 
#' @param x list of genes 
#' @param  y TPM values dataframe - requires WT_0_1 format


plot_gene_expression <- function(genes_of_interest, df) {
  # Now we are going to do all the steps inside the function {}
  df_plot_genes <- df[rownames(df) %in% genes_of_interest, ]
  # Move rownames to a column
  df_plot_genes <- rownames_to_column(df_plot_genes, var = "gene_id")
  # Pivot longer
  df_plot_genes <- df_plot_genes %>%
    pivot_longer(cols = -gene_id,
                 names_to = 'Time_Replicate',
                 values_to = 'TPM')
  # Separate sample name into time and replicate
  df_plot_genes <- df_plot_genes %>%
    separate(
      Time_Replicate,
      into = c('WT', 'Time', 'Replicate'),
      sep = '_',
      extra = "merge"
    ) %>%
    select(-WT)
  # Convert Time to numeric values
  df_plot_genes$Time <- as.numeric(df_plot_genes$Time)
  # Group by Time and calculate mean and std of TPM
  df_plot_genes <- df_plot_genes %>%
    group_by(gene_id, Time) %>%
    summarise(mean_TPM = mean(TPM), sd_TPM = sd(TPM))
  
  # Plot the line plot with error bars
  ggplot(df_plot_genes, aes(x = Time, y = mean_TPM)) +
    geom_line() +
    geom_errorbar(aes(ymin = mean_TPM - sd_TPM, ymax = mean_TPM + sd_TPM), width = 0.2) +
    labs(x = 'Time (hours)', y = 'Average TPM') +
    ggtitle('Average TPM across Time Points') +
    facet_wrap( ~ gene_id)
}

# Example usage:
# genes_of_interest <- c(
  #"ENSMUSG00000094125.1",
  #"ENSMUSG00000079247.2",
  #"ENSMUSG00000112027.1",
  #"ENSMUSG00000066632.3",
  #"ENSMUSG00000000031.16"
#)

# plot_gene_expression(genes_of_interest, TPM_filtered)

# LFCplot.R

#' LFCplot: A Function to Plot Log2 Fold Changes Over Time
#'
#' This function takes a list of genes of interest and DESeq2 results to produce a
#' faceted line plot showing log2 fold changes over a time course. The time points
#' are extracted from the DESeq2 result names, and each gene is plotted separately
#' with its own y-axis scale.
#'
#' @param genes_of_interest A data frame with a column named `gene` containing the gene IDs of interest.
#' @param res_df A data frame containing DESeq2 results with the following columns:
#'   \describe{
#'     \item{gene_id}{Gene identifiers matching those in `genes_of_interest`.}
#'     \item{baseMean}{The mean of normalized counts for all samples.}
#'     \item{log2FoldChange}{The log2 fold change between the conditions.}
#'     \item{lfcSE}{The standard error of the log2 fold change.}
#'     \item{stat}{The test statistic.}
#'     \item{pvalue}{The p-value for the statistical test.}
#'     \item{padj}{The adjusted p-value (Benjamini-Hochberg correction).}
#'     \item{gene_name}{The name of the gene (optional).}
#'     \item{result_name}{The name of the comparison, which should include the time point information.}
#'   }
#'
#' @return A ggplot object with a faceted line plot of log2 fold changes over time for each gene.
#'
#' @examples
#' # Define a set of genes of interest
#' genes_of_interest <- data.frame(gene = c("ENSMUSG00000094125.1", 
#'                                          "ENSMUSG00000079247.2", 
#'                                          "ENSMUSG00000112027.1", 
#'                                          "ENSMUSG00000066632.3", 
#'                                          "ENSMUSG00000000031.16"))
#'
#' # Use the function to plot log2 fold changes
#' LFCplot(genes_of_interest, res_df)
#'
#' @export
LFCplot <- function(genes_of_interest, res_df) {
  library(dplyr)
  library(ggplot2)
  
  # Filter DESeq2 results to include only the genes of interest
  df_lfc_plot <- res_df[res_df$gene_id %in% genes_of_interest$gene, ]
  
  # Extract the time point from the 'result_name' and convert to numeric
  df_lfc_plot <- df_lfc_plot %>%
    mutate(time_point = gsub("time_point_|_vs_0", "", result_name)) %>%
    mutate(time_point = as.numeric(time_point))
  
  # Sort the data by gene and time point
  df_lfc_plot <- df_lfc_plot %>%
    arrange(gene_id, time_point)
  
  # Generate the line plot with facet_wrap for each gene
  plot <- ggplot(df_lfc_plot, aes(x = time_point, y = log2FoldChange, group = gene_id)) +
    geom_line() +
    geom_point() +
    facet_wrap(~ gene_id, scales = "free_y") +  # Facet by gene
    labs(title = "Log2 Fold Change Over Time",
         x = "Time Point (hours)",
         y = "log2 Fold Change") +
    theme_minimal()
  
  return(plot)
}


# ATACSEQ IMPORT PEAKS

#' FUNCTION to import peak files from ATACseq into List of GRanges objects
#'
#' @description 
#' Any peak file can be loaded in and made into a GRanges from one file 
#' to thouusands of peak files. This function will load and get into GRange format.
#' 
#' @param peak_file, file path to the peak files - generic name "consensus_file_path".  
#' @param  regex, # Warning custom to application # there is a regex to get "sample name" for input into function 
#' 
#' 
import_peaks <- function(consensus_file_path) {
  
  # List all MACS2 broadPeak files in the directory
  peak_files <- list.files(consensus_file_path, full.names = TRUE, pattern = ".broadPeak")
  
  # Extract sample names directly from file names
  sample_names <- sapply(peak_files, function(file) {
    str_extract(file, "(KO|WT)_control_\\d+")
  })
  
  # Initialize an empty list to store filtered GRanges objects
  peak_list <- list()
  
  # Import, filter, and store GRanges objects
  for (i in seq_along(peak_files)) {
    # Import peaks as GRanges
    peaks <- rtracklayer::import(peak_files[i])
    
    # Filter to keep only seqnames starting with "chr"
    filtered_peaks <- peaks[grepl("^chr", as.character(seqnames(peaks)))]
    
    # Add filtered peaks to the list with a sample-specific name
    peak_list[[sample_names[i]]] <- filtered_peaks
  }
  
  return(peak_list)
}





#' CREATE CONSENSUS PEAKS
#' this function will take multiple replicate .broadPeak files (also narrow)
#' find peaks that overlap in all the replicates. 
#' @description 
#' input set of chipseq replicate peak files
#' this function then creates one merged file peaks in all samples
#' @param sample_name
#' This will be extracted with names(GR_list) in the lapply at end of fun
#' You will need a "dbps" or some object for the lapply that has the 
#' name of each dbp in the named GRanges list
#' 
#' @param peak_list
#' Named list of GRanges for each chipseq replicate
#' peak_list can be generated using import_peaks function above

consensus_from_reduced <- function(dbp, peak_list) {
  dbp_peaks <- peak_list[grepl(as.character(dbp), names(peak_list))]
  suppressWarnings(all_peaks <- GenomicRanges::reduce(unlist(as(dbp_peaks, "GRangesList"))))
  all_peaks <- all_peaks[grepl("chr", seqnames(all_peaks))]
  
  # peak_exists <- lapply(dbp_peaks, function(x) {
  #   as.numeric(countOverlaps(all_peaks, x) > 0))
  # }) %>%
  # bind_rows() OR bind_cols()
  peak_exists <- matrix(NA, nrow = length(all_peaks), ncol = length(dbp_peaks))
  for(i in 1:length(dbp_peaks)) {
    suppressWarnings(peak_exists[,i] <- as.numeric(countOverlaps(all_peaks, dbp_peaks[[i]]) > 0))
  }
  
  # filter to consensus requiring peaks to be in all replicates
  dbp_consensus <- all_peaks[rowSums(peak_exists) == ncol(peak_exists)]
  # Required only two replicates == dbp_consensus <- all_peaks[rowSums(peak_exists) > 1]
  return(dbp_consensus)
}


