
library(scProcess)
library(ggpubr)
library('lib.r')

raw = readRDS('raw.rds')



library(conos)
library(cacoa)
library(ggplot2)
library(magrittr)
library(dropestr)
library(pbapply)

# define threshold

threshold.totalUMI <- 700
threshold.MT.ratio <- 0.2
threshold.doubletScores <- 0.4

cms  <- raw
# total UMI
lapply(names(cms[inter]), function(n)
  dropestr::PlotCellsNumberHist(colSums(cms[[n]]), estimate.cells.number=T, show.legend=F) +
    geom_vline(aes(xintercept=log10(threshold.totalUMI)))) %>%
  cowplot::plot_grid(plotlist=., ncol=3, labels=names(cms))



cms <- names(cms[inter]) %>% setNames(., .) %>%
  pblapply(function(n) cms[[n]] %>% .[, colSums(.) >= threshold.totalUMI])

# MT gene ratio

mit_frac_per_dataset <- pblapply(cms, function(cm) 
  colSums(cm[grep("MT-", rownames(cm)), ]) / colSums(cm))

lapply(mit_frac_per_dataset,summary)

lapply(mit_frac_per_dataset, function(fr)
  qplot(fr[fr < 0.5], xlab="Mit. fraction", ylab="#Cells", xlim=c(-0.01, 0.5), bins=30) +
  geom_vline(xintercept = threshold.MT.ratio, col ='red')
  ) %>%
  cowplot::plot_grid(plotlist=., ncol=4, labels=names(cms))


cms <- pbmapply(function(cm, frac) cm[, frac < threshold.MT.ratio], cms, mit_frac_per_dataset) %>%
  setNames(names(cms))


# doublest score

doublet_info <- sccore:::plapply(cms, function(x) { get.scrublet.scores(x) }, n.cores = 10, mc.preschedule = T, progress = TRUE)



doublet_info.all <- Reduce(c, doublet_info)

index <- doublet_info.all > threshold.doubletScores
doublet_info.all[index]='doublets'
doublet_info.all[!index]='singlets'

table(doublet_info.all)

singlets.cell <- names(doublet_info.all[doublet_info.all=='singlets'])

cms_filt <- lapply(cms, function(x) {
  x[!grepl("^MT-", rownames(x)),intersect(colnames(x),singlets.cell)]
})



saveRDS(cms_filt,'raw.filter.rds',sep=''))








