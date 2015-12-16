library(ggplot2)

res  <- read.csv('d_tissue.csv')
p1 <- ggplot(data = res, aes( x = nb_probes)) + geom_histogram()+ scale_x_discrete(breaks=c(1:12))+ facet_wrap(~ID.name, scale = "free_x") + theme_bw() 

# res$ID1  <- as.character(res$ID)
# ind  <- which(duplicated(res$pair_probes))
# res$ID1[ind]  <- as.character(res$pair_probes[ind])
ind  <- !duplicated(res$pair_probes)
res1  <- res[ind, ]



m <- ggplot(res1, aes(x=phase_abs_diff)) + geom_histogram( fill = "blue")+ labs(y = "# of probe pairs", x = "Phase difference") + theme_bw()
ggsave(paste0('d_phase_diff_arabidopsis.pdf'), m, width = 9, height = 7)
