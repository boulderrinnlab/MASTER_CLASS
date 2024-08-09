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
