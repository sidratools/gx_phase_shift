rm(list = ls())
library(plyr)
# bat  <- read.csv('bat_cor_quantile_mat_Output.csv')
# bat  <- read.csv('liver_cor_quantile_mat_Output.csv')
# prefix  <- 'bat'
# prefix  <- 'iwat'
prefix  <- 'iwat'
bat  <- read.csv(paste0(prefix, '_cor_quantile_mat_Output.csv', collapse = ''))
# bat  <- read.csv('')

bat$Phase_class[bat$Phase_class==20]  <- 6

# ann  <- read.table('Mouse430A_2_Master_Annotation_File.csv', sep = '\t', header = T)
load('Mouse430A_2_Master_Annotation_File.RData')

ii  <- bat$P_value <= 0.1
bat  <- bat[ii, ]
bat  <- droplevels(bat)

imat  <- match(bat$ProbeID,ann$Probe.Set.ID)

ann1  <-  subset(ann, select = c(Probe.Set.ID, Gene.Symbol)) 
# bat$Probe1 <- ann$Probe.Set.ID[imat]
bat$Gene.Symbol  <- ann$Gene.Symbol[imat]
# bat  <- head(bat, 1000)
bat1  <-  subset(bat, select = c(ProbeID, Phase_class, P_value, Gene.Symbol))



pro  <- function(df){
	x  <- df$Phase_class
	probes  <- df$ProbeID
	pvalues  <- df$P_value

	if (length(x) > 1) {
		mat_prob  <-  combn(probes, 2)
		mat_prob  <- matrix(as.character(mat_prob),  nc = dim(mat_prob)[2])
		y <-  apply(t(combn(x, 2)), 1, function(u){abs(diff(u))})
	 	prob_diff <-  apply(mat_prob, 2, paste, collapse = '-') 
# 		prob_diff  <-  1
# 		browser()

		} else
		{
			y  <- NA
			prob_diff  <- NA
		}
	data.frame(phase_abs_diff = y, pair_probes = as.character(prob_diff), nb_probes = length(x))
# 	data.frame(Y = y)
}

# 	ifelse(length(x) > 1, y <-  apply(t(combn(x, 2)), 1, function(u){abs(diff(u))}), y <- NA)
# 	data.frame(Y = y)
# }


idx.rm  <- grep('---', bat$Gene.Symbol)

bat2  <- bat1[-idx.rm, ]
# bat2  <- bat1
bat2$ProbeID  <- as.character(bat2$ProbeID)
bat3  <-  ddply(bat2, .(Gene.Symbol), .fun = pro )
# browser()

bat4  <- na.omit(bat3)


v1  <- unlist(lapply(as.list(2:11), function(x) {dim(combn(x, 2))[2]}))

dplt  <- data.frame(x = 2:11, y = v1)
library(ggplot2)
p_combn  <- ggplot(dplt) + geom_line(aes(x = x, y = y)) +  geom_point(aes(x = x, y = y)) + labs( x = "# of probes in Gene", y = "# of possible combinations") + theme_bw()

m <- ggplot(bat4, aes(x=phase_abs_diff)) + geom_histogram( fill = "blue")+ labs(y = "# of probe pairs", x = "Phase difference") + theme_bw()
ggsave(paste0('phase_diff_',  prefix,'.pdf'), m, width = 9, height = 7)

m1 <- ggplot(bat4, aes(x=phase_abs_diff)) + geom_histogram( fill = "blue") + facet_wrap(~nb_probes) + labs(y = "# of probe pairs", x = "Phase difference") +theme_bw()
ggsave(paste0('phase_diff_per_nb_probes_',  prefix,'.pdf'), m1, width = 9, height = 7)

write.csv(bat4, paste0(prefix, '_analysis.csv'))

ind0 <- with(bat4, phase_abs_diff != 0)

abs_dif  <- bat4$phase_abs_diff[ind0]

length(abs_dif)

ks  <- ks.test(abs_dif, "punif", 1, 6)



