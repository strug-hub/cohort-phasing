library(ggplot2)
library(ggpubr)
library(svglite)

prefix <- "./data/"

#plot corrdinate parameters
shift <-  1000
x_values <- c(0, 10000, 20000, 30000, Inf)

lr <- read.csv(paste(prefix, "longrangerv2.2.2_union_with_platinumgenome_16x2000nt_bin_size.tsv", sep=""), header = TRUE, sep= '\t')
ont <- read.csv(paste(prefix, "nanopre_16x2000nt_bin_size.tsv", sep=""), header = TRUE, sep= '\t')
ccs <- read.csv(paste(prefix, "pacbio_ccs_16x2000nt_bin_size.tsv", sep=""), header = TRUE, sep= '\t')
clr <- read.csv(paste(prefix, "pacbio_clr_16x2000nt_bin_size.tsv", sep=""), header = TRUE, sep= '\t')

lr_error <- read.csv(paste(prefix, "longrangerv2.2.2_union_with_platinumgenome_16x2000nt_error_rate.tsv", sep=""), header = TRUE, sep= '\t')
ont_error <- read.csv(paste(prefix, "nanopre_16x2000nt_error_rate.tsv", sep=""), header = TRUE, sep= '\t')
ccs_error <- read.csv(paste(prefix, "pacbio_ccs_16x2000nt_error_rate.tsv", sep=""), header = TRUE, sep= '\t')
clr_error <- read.csv(paste(prefix, "pacbio_clr_16x2000nt_error_rate.tsv", sep=""), header = TRUE, sep= '\t')

f_small <- function(bins, error, title, error_line, variant_count) {
  bins$phase_status <- factor(bins$phase_status, levels=c("unphased", "phased_incorrect", "phased_correct"))
  p <- ggplot() +
    geom_bar(data = bins, aes(x=(bin_start + shift)/1000, y=percent, fill=phase_status), stat="identity" ) + 
    scale_fill_manual(values=c("grey69","tomato", "turquoise3"), labels = c("Unphased", "Incorrectly Phased", "Correctly Phased")) +
    geom_segment(aes(x = 0, xend = Inf, y = 0, yend = 0),
                 arrow = arrow(length = unit(.2, 'cm'))) +
    labs(x="Distance between heterozygous variants (KB)\n", y = "Phase Percentage\n", title = title, fill = "Phase Status", color = "") + theme_bw() + 
    theme(
      axis.title.x = element_text(size = 18),
      axis.text.x = element_text(size = 16),
      axis.title.y = element_text(size = 18),
      axis.text.y = element_text(size = 16),
      legend.title = element_text(size=18),
      legend.text=element_text(size=18),
      plot.title = element_text(size=20))
  if (error_line) {
    p <- p + 
      geom_line(data = error, aes(x = (bin_start + shift)/1000, y = error_rate * (10/8), colour = "Percent of Incorrectly Phased Variant"), size = 2, alpha = 0.75) + 
      scale_color_manual(values="purple", labels = "Percent of Incorrectly Phased Variant") +
      scale_y_continuous(sec.axis = sec_axis(~.* (8/10), name = "Phasing Error Rate\n")) 
  }
  if (variant_count) {
    p <- p +
      geom_text(data=error,aes(x=(bin_start + shift)/1000,y=0.01,label=  sprintf("%s / %s", total_phased, total_vars)),hjust=0, size = 5, angle = 90, alpha = 0.6)
  }
  return(p)
}

lr_plot <- f_small(lr, lr_error, "LongRanger 10X Genomics", FALSE, FALSE)
ont_plot <- f_small(ont, ont_error, "Nanopore", FALSE, FALSE)
ccs_plot <- f_small(ccs, ccs_error, "PacBio Consensus Circular Sequence (CCS)", FALSE, FALSE)
clr_plot <- f_small(clr, clr_error, "PacBio Continuous Long Reads (CLR)", FALSE, FALSE)

p <- ggarrange(lr_plot, ont_plot, ccs_plot, clr_plot, common.legend = TRUE, legend="right")
q <- annotate_figure(p, top = text_grob("Probability of Phase Status vs Distance Between Heterozygous Variants\n", face = "bold", size = 25))

ggsave(file="30kb_bin_errors.svg", plot=q, width=20, height=15)
