library(ggplot2)
library(gridExtra)
library(ggpubr)
library(svglite)

dir <- "data/"

#plots the phasing stats of different combinations of technology 
compare <- read.csv(paste(dir, "avg_switchflip_allchr_allcombos.tsv_parsed.txt", sep=""), header = TRUE, sep= '\t')
compare_pos <- read.csv(paste(dir, "avg_switchflip_allchr_allcombos.tsv_pos.txt", sep=""), header = TRUE, sep= '\t')

stats <- read.csv(paste(dir, "stats_summary_allchr_allcombos.tsv_numdataset.tsv", sep=""), header = TRUE, sep= '\t')
phase_blocks_pos <- read.csv(paste(dir, "stats_summary_allchr_allcombos.tsv_phase_blocks.txt", sep=""), header = TRUE, sep= '\t')
perc_phased_var_pos <- read.csv(paste(dir, "stats_summary_allchr_allcombos.tsv_phased.txt", sep=""), header = TRUE, sep= '\t')
n50_pos <- read.csv(paste(dir, "stats_summary_allchr_allcombos.tsv_n50.txt", sep=""), header = TRUE, sep= '\t')

(phase_errors <- ggplot(compare, aes(x=factor(num_datasets), y=switch_flip_rate)) + 
    #geom_boxplot(aes(color = factor(num_datasets)), color = c("#F8766D","#B79F00", "#00BA38", "#00BFC4")) +
    stat_summary(fun = mean, geom = "errorbar", 
                 aes(ymax = ..y.., ymin = ..y.., group = factor(num_datasets)),
                 width = 0.75, linetype = "solid", colour = "red", size = 2) +
    geom_point(data = compare_pos, aes(x = x, y=y, colour=factor(dataset), shape = factor(dataset)), 
               show.legend = TRUE, size = 4) +
    scale_shape_manual(values=c(15, 16, 17, 18, 8)) +
    labs(x="Number of Technologies in Combination\n", y="Average Switch and Flip Error Rates", title="Phasing Errors", 
         color ="Technologies", shape = "Technologies" ) +
    theme(
      axis.title.x = element_text(size = 18),
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 14),
      axis.title.y = element_text(size = 18),
      plot.title = element_text(size=20),
      legend.text=element_text(size=18),
      legend.title = element_text(size=20)))


(num_block <- ggplot(data=stats,aes(x=as.factor(num_datasets), y=blocks)) +
    #geom_boxplot(data=stats,aes(x=as.factor(num_datasets), y=blocks), outlier.shape = NA, show.legend = FALSE) +
    stat_summary(fun = mean, geom = "errorbar", 
                 aes(ymax = ..y.., ymin = ..y.., group = factor(num_datasets)),
                 width = 0.75, linetype = "solid", colour = "red", size = 2) +
    geom_point(data=phase_blocks_pos, aes(x = x, y=y, colour=factor(dataset), shape = factor(dataset)), show.legend = TRUE, size = 4) + 
    scale_shape_manual(values=c(15, 16, 17, 18, 8)) +
    labs(x="Number of Technologies in Combination\n", y="Number of Phase Blocks", title = "Number of Phase Blocks",
         color ="Technologies", shape = "Technologies") +
    theme(
      axis.title.x = element_text(size = 18),
      axis.text.x = element_text(size = 16),
      axis.text.y = element_text(size = 16),
      axis.title.y = element_text(size = 18),
      plot.title = element_text(size=20),
      legend.text=element_text(size=18),
      legend.title = element_text(size=20)))

stats[22] <- stats[22]/1000000
n50_pos$y <- n50_pos$y/1000000
(n50 <- ggplot(data=stats,aes(x=as.factor(num_datasets), y=block_n50)) +
    #geom_boxplot(data=stats,aes(x=as.factor(num_datasets), y=blocks), outlier.shape = NA, show.legend = FALSE) +
    stat_summary(fun = mean, geom = "errorbar", 
                 aes(ymax = ..y.., ymin = ..y.., group = factor(num_datasets)),
                 width = 0.75, linetype = "solid", colour = "red", size = 2) +
    geom_point(data=n50_pos, aes(x = x, y=y, colour=factor(dataset), shape = factor(dataset)), show.legend = TRUE, size = 4) + 
    scale_shape_manual(values=c(15, 16, 17, 18, 8)) +
    labs(x="Number of Technologies in Combination\n", y="Phase Blocks N50 (MB)", title = "Phase Blocks N50",
         color ="Technologies", shape = "Technologies") +
    theme(
      axis.title.x = element_text(size = 18),
      axis.text.x = element_text(size = 16),
      axis.text.y = element_text(size = 16),
      axis.title.y = element_text(size = 18),
      plot.title = element_text(size=20),
      legend.text=element_text(size=18),
      legend.title = element_text(size=20)))

(percent_phased <- ggplot(data=stats,aes(x=as.factor(num_datasets), y=(phased/heterozygous_variants))) +
    #geom_boxplot(data=stats,aes(x=as.factor(num_datasets), y=unphased), outlier.shape = NA, show.legend = FALSE) +
    stat_summary(fun.y = mean, geom = "errorbar", 
                 aes(ymax = ..y.., ymin = ..y.., group = factor(num_datasets)),
                 width = 0.75, linetype = "solid", colour = "red", size = 2) +
    geom_point(data=perc_phased_var_pos, aes(x = x, y=y, colour=factor(dataset), shape = factor(dataset)), 
               show.legend = TRUE, size = 4) + 
    scale_shape_manual(values=c(15, 16, 17, 18, 8)) +
    labs(x="Number of Technologies in Combination\n", y="Percentage of Unphased Variants", title = "Percentage of Phased Variants",
         color ="Technologies", shape = "Technologies") +
    theme(
      axis.title.x = element_text(size =18),
      axis.text.x = element_text(size = 16),
      axis.text.y = element_text(size = 16),
      axis.title.y = element_text(size = 18),
      plot.title = element_text(size=20),
      legend.text=element_text(size=18),
      legend.title = element_text(size=20)))

p <- ggarrange(phase_errors, percent_phased, num_block, n50, common.legend = TRUE, legend="bottom", ncol = 4)
(q <- annotate_figure(p, top = text_grob("Combination of Phasing Technologies Statistics\n", face = "bold", size = 25)))

ggsave(file="combination_stats_compare.svg", plot=q, width=20, height=10)
