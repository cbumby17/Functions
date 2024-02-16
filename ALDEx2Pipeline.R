#### ALDEx2 Pipeline (This needs to be wrapped in a single function) ####
aldex2pipeline <- function(reads, conds){
  reads <- as.data.frame(reads)
  conds <- as.character(conds)
  
  x.aldex <- ALDEx2::aldex(reads, conds, mc.samples = 128, test = "t", effect = TRUE, 
                   include.sample.summary = FALSE, denom = "all", verbose = FALSE, 
                   paired.test = FALSE, gamma = NULL)
  
  #Check to see whether the geometric mean is appropriate for our situation. 
  #These should be centered around zero. If they aren't, revisit ALDEx2 tutorial.
  hist(x.aldex$diff.btw)
  hist(x.aldex$effect)
  
  #The left panel is an Bland-Altman or MA plot that shows the relationship between (relative) Abundance and Difference.
  MAplot <- ALDEx2::aldex.plot(x.aldex, type="MA", test="welch", xlab="Log-ratio abundance",
                       ylab="Difference", main='Bland-Altman plot')
  #The middle panel is an effect plot that shows the relationship between Difference and Dispersion; the lines are equal difference and dispersion.
  MWplot <- ALDEx2::aldex.plot(x.aldex, type="MW", test="welch", xlab="Dispersion",
                       ylab="Difference", main='Effect plot')
  #The right hand plot is a volcano plot; the lines represent a posterior predictive p-value of 0.001 and 1.5 fold difference. 
  volcano <- ALDEx2::aldex.plot(x.aldex, type="volcano", test="welch", xlab="Difference",
                        ylab="-1(log10(q))", main='Volcano plot') 
  #In all plots features that are not significant are in grey or black. Features that are statistically significant are in red. 
  #The log-ratio abundance axis is the clr value for the feature.
  aldexplots <- plot_grid(MAplot, MWplot, volcano)
  
  #make table of significant taxa
  sig_aldex <- x.aldex %>%
    rownames_to_column(var = "OTU") %>%
    filter(we.ep < 0.05) %>%
    arrange(effect, we.ep) %>%
    dplyr::select(OTU, diff.btw, diff.win, effect, we.ep, we.eBH)
  sig_aldex 
  
  #Let's make a fancier volcano plot
  aldexres <- x.aldex %>%
    mutate(Significant = if_else(we.ep < 0.05,TRUE, FALSE)) %>%
    mutate(Taxon = as.character(rownames(x.aldex))) %>%
    mutate(TaxonToPrint = if_else(we.ep < 0.05, Taxon, "")) #only provide a label to signifcant results
  
  #Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2fc respectively positive or negative)
  ##The down and up-regulated taxa and to be reversed so the comparison was KO compared to WT.
  aldexres$diffexpressed <- "NO"
  #if log2Foldchange > 0.6 and pvalue < 0.05, set as "DOWN"
  aldexres$diffexpressed[aldexres$diff.btw > 0.6 & aldexres$we.ep < 0.05] <- "DOWN"
  #if log2Foldchange < -0.6 and pvalue < 0.05, set as "UP"
  aldexres$diffexpressed[aldexres$diff.btw < -0.6 & aldexres$we.ep < 0.05] <- "UP"
  
  #below, -diff.btw was used becase we needed to reverse the default comparison
  vplot <- ggplot2::ggplot(data = aldexres, aes(x = -diff.btw, y = -log10(we.ep), col = diffexpressed, label = TaxonToPrint)) +
    geom_point(size = 2) +
    scale_color_manual(values = c("#00AFBB", "grey", "#FF0000"), # to set the colours of our variable
                       labels = c("Decreased", "Not significant", "Increased")) + # to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)
    labs(color = "Differential Abundance",
         x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value")) +
    geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
    geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
    ggtitle('Differentially abundant taxa in Low vs High Inflammatory Response') + # Plot title 
    geom_text_repel(max.overlaps = Inf, show.legend = FALSE) + # To show all labels 
    theme_bw()
  
  return(list(aldexplots, vplot))
  
}
