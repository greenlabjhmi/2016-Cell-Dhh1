library(ggplot2)

Dat <- read.table("Data/RPKM-All-Datasets.txt", header = TRUE)

pdf("Plots/TAIvGC.pdf")
ggplot(data = Dat) + stat_binhex(aes(x = TAI, y = GC), bins = 50) + scale_x_continuous(limits = c(0,1)) + 
                     scale_y_continuous(limits = c(0,1)) + theme_bw() + theme(panel.grid = element_blank())
dev.off()
