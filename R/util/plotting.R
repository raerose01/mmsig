## Plotting functions for mmsig

rotatedAxisElementText = function(angle,position='x'){
    angle     = angle[1]; 
    position  = position[1]
    positions = list(x=0,y=90,top=180,right=270)
    if(!position %in% names(positions))
        stop(sprintf("'position' must be one of [%s]",paste(names(positions),collapse=", ")),call.=FALSE)
    if(!is.numeric(angle))
        stop("'angle' must be numeric",call.=FALSE)
    rads  = (angle - positions[[ position ]])*pi/180
    hjust = 0.5*(1 - sin(rads))
    vjust = 0.5*(1 + cos(rads))
    element_text(angle=angle,vjust=vjust,hjust=hjust)
}

scale_fill_sigs <- function(...){
  #ggplot2:::manual_scale('fill', 
  #                       values = setNames(c(RColorBrewer::brewer.pal(8, "Dark2"), "lightblue", "#FB8072"),
  #                                         c("SBS1", "SBS2", "SBS5", "SBS8", "SBS9", "SBS13", "SBS18", "SBS-MM1", "SBS35", "SBS84")), 
  #                       ...)
  
  ggplot2:::manual_scale('fill', 
                         values = setNames(c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928'),
                                           c("Signature.1", "Signature.2", "Signature.4", "Signature.5", "Signature.6", "Signature.7", "Signature.13", "Signature.8", "Signature.10", "Signature.30", "Signature.11", "Signature.9")), 
                         ...)
}

plot_signatures = function(sig,
                           sig_order = c("Signature.1","Signature.2","Signature.4","Signature.5","Signature.13"),
                           samples = FALSE){
  "
  Plotting function for mmSig package
  sig_order: order of signatures to plot in stacked bar chart, from top to bottom
  samples: boolean whether or not to plot sample names on the x axis
  "
  sigPlot <- sig %>%
      rownames_to_column(var = "sample") %>%
      dplyr::select(-mutations) %>%
      melt(id.var = "sample", variable.name = "Signature", value.name = "prop") %>%
      mutate(Signature = factor(Signature, levels = sig_order)) %>%
      ggplot(aes(sample, prop, fill = Signature)) +
      geom_col(width=1)+
      scale_fill_sigs()+
      scale_y_continuous(expand = c(0,0))+
      theme_bw()+
      labs(x = "Sample",
           y = "Relative contribution",
           fill = "Signature")+
      theme(text = element_text(size = 10, color = "black"),
            axis.text = element_text(color = "black"),
            axis.text.y = element_text(size = 12),
            axis.title = element_text(size = 14),
            panel.grid = element_blank(),
            panel.border = element_blank(),
            legend.position = "right",
            legend.text = element_text(size = 12),
            legend.title = element_text(size = 14),
            axis.title.x = element_blank())
  
  if(samples){
    sigPlot +
          theme(axis.text.x = rotatedAxisElementText(90, "top"))
  } else{
    sigPlot +
          theme(axis.text.x = element_blank(),
                axis.ticks.x = element_blank())
  }
}

bootSigsPlot <- function(mutSigsSummary){
    ggplot(mutSigsSummary, aes(signature, estimate, fill = signature))+
        geom_bar(position="dodge", stat="identity")+
        geom_errorbar(data = mutSigsSummary, mapping = aes(x = signature, ymin = CI025, ymax = CI975))+
        facet_grid(~sample)+
        theme(text = element_text(size = 12),
              axis.text.x = rotatedAxisElementText(90, 'top'),
              axis.title.x = element_blank(),
              strip.background = element_blank(),
              legend.position = 'none')+
        scale_fill_sigs()+
        labs(y = 'Relative contribution')
}
