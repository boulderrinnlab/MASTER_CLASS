#! Function Lesson !!

###################################
# Our first function - multiplication
###################################

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


###################################
# PLOT TPM values 
###################################

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



###################################
# Log fold change plot from Deseq2 values
###################################

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



###################################
# ATACSEQ IMPORT PEAKS
###################################

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


###################################
# FIND COMMON PEAKS
###################################

#' FUNCTION to identify common peaks across multiple ATAC-seq samples
#'
#' @description 
#' This function processes a list of GRanges objects, where each GRanges represents 
#' ATAC-seq peaks from a single sample. It identifies genomic regions (peaks) 
#' that are shared across all samples by iteratively calculating overlaps. 
#' The resulting common peaks are labeled with unique identifiers for downstream 
#' analysis or visualization in tools like IGV.
#'
#' @param gr_list A list of GRanges objects. Each GRanges represents the peaks from 
#'                a single ATAC-seq sample. The list should have names corresponding 
#'                to the sample identifiers (e.g., "sample1", "sample2", etc.).
#' 
#' @return A GRanges object containing the genomic intervals (peaks) that are common 
#'         across all input GRanges objects. Each interval is labeled with a unique 
#'         identifier in the metadata column `name`, using the format "common_peak_<number>".
#'
#' @details 
#' 1. The function validates that the input is a list of GRanges objects.
#' 2. It iteratively identifies overlaps across all GRanges objects in the list.
#' 3. Only intervals shared across all samples are retained at each step.
#' 4. The resulting peaks are assigned unique names for easy tracking.
#'
#' @examples
#' # Example list of GRanges objects
#' gr_list <- list(
#'   sample1 = GRanges(seqnames = "chr1", ranges = IRanges(start = c(1, 50), end = c(20, 70))),
#'   sample2 = GRanges(seqnames = "chr1", ranges = IRanges(start = c(10, 60), end = c(30, 80))),
#'   sample3 = GRanges(seqnames = "chr1", ranges = IRanges(start = c(15, 55), end = c(25, 75)))
#' )
#'
#' # Find common peaks
#' common_peaks <- find_common_peaks(gr_list)
#'
#' # View results
#' print(common_peaks)
#'
#' @export
find_common_peaks <- function(gr_list) {
  # Validate input
  if (!is.list(gr_list) || !all(sapply(gr_list, inherits, "GRanges"))) {
    stop("Input must be a list of GRanges objects.")
  }
  
  # Start with the first GRanges object
  common_peaks <- gr_list[[1]]
  
  # Iteratively find overlaps across all GRanges objects
  for (i in 2:length(gr_list)) {
    current_gr <- gr_list[[i]]
    
    # Find overlaps
    overlaps <- findOverlaps(common_peaks, current_gr)
    
    # Subset to overlapping regions
    common_peaks <- subsetByOverlaps(common_peaks, current_gr)
  }
  
  # Assign custom names to the common peaks
  mcols(common_peaks)$name <- paste0("common_peak_", seq_along(common_peaks))
  
  return(common_peaks)
}


###################################
# FIND My PEAKS
###################################


#' Find Unique Peaks in a GRanges Object
#'
#' This function identifies peaks that are unique to one GRanges object 
#' (`original_peaks`) by excluding peaks that overlap with another 
#' GRanges object (`common_peaks`). It is designed for genomic peak 
#' analyses, such as identifying condition-specific peaks in comparative 
#' datasets.
#'
#' @param common_peaks A GRanges object representing the common peaks 
#' across multiple conditions or datasets.
#' @param original_peaks A GRanges object representing the peaks for a 
#' specific condition or dataset.
#' 
#' @return A GRanges object containing peaks that are unique to 
#' `original_peaks`, i.e., those that do not overlap with any peaks 
#' in `common_peaks`.
#'
#' @examples
#' # Load required package
#' library(GenomicRanges)
#'
#' # Example GRanges objects
#' common_peaks <- GRanges(seqnames = "chr1", ranges = IRanges(start = c(1, 50, 100), end = c(10, 60, 110)))
#' original_peaks <- GRanges(seqnames = "chr1", ranges = IRanges(start = c(5, 70, 120), end = c(15, 80, 130)))
#'
#' # Identify unique peaks
#' unique_peaks <- find_my_peaks(common_peaks, original_peaks)
#'
#' # Inspect results
#' unique_peaks
#'
#' @export
find_my_peaks <- function(common_peaks, original_peaks) {
  # Find overlaps
  overlaps <- findOverlaps(original_peaks, common_peaks)
  
  # Identify peaks in original_peaks that are not in common_peaks
  unique_peaks <- original_peaks[-queryHits(overlaps)]
  
  # Return unique peaks
  return(unique_peaks)
}
# create find_my_peaks function
find_my_peaks <- function(common_peaks, original_peaks) {
  # Find overlaps
  overlaps <- findOverlaps(original_peaks, common_peaks)
  
  # Identify peaks in original_peaks that are not in common_peaks
  unique_peaks <- original_peaks[-queryHits(overlaps)]
  
  # return peaks in condition selected
  return(unique_peaks)
}




