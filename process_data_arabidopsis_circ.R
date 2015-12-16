rm(list = ls())
library(plyr)
library(ggplot2)
# bat  <- read.csv('bat_cor_quantile_mat_Output.csv')
# bat  <- read.csv('liver_cor_quantile_mat_Output.csv')
# prefix  <- 'bat'
# prefix  <- 'iwat'

v.prefix1  <- c("RefSeq.Transcript.ID", "AGI", "Gene.Symbol", "Entrez.Gene")
v.prefix  <- c('d', 'w')

vv  <- expand.grid(v.prefix, v.prefix1)

for (i in c(1:dim(vv)[1])) {
	print(i)
	     prefix  <- vv[i, 1]
	     prefix1  <- as.character(vv[i, 2])

bat  <- read.csv(paste0( 'Arabidopsis_', prefix, '_Results.csv', collapse = ''))
# bat  <- read.csv('')

bat$Phase[bat$Phase == 20]  <- 6

# ann  <- read.csv('ATH1-121501.na24.annot.csv')
load('ATH1-121501.na24.annot.RData')

# load('Mouse430A_2_Master_Annotation_File.RData')

ii  <- bat$P_value <= 1
bat  <- bat[ii, ]
bat  <- droplevels(bat)

imat  <- match(bat$ProbeID,ann$Probe.Set.ID)

# ann1  <-  subset(ann, select = c(Probe.Set.ID, AGI)) 
ann1  <-  subset(ann, select = c("Probe.Set.ID", prefix1)) 
# bat$Probe1 <- ann$Probe.Set.ID[imat]
# bat$Entrez.Gene  <- ann$Entrez.Gene[imat]
bat[, prefix1]  <- ann[, prefix1][imat]
# bat  <- head(bat, 1000)
bat1  <-  subset(bat, select = c("ProbeID", "Phase", "P_value", prefix1))



pro  <- function(df){
	x  <- df$Phase
	probes  <- df$ProbeID
	pvalues  <- df$P_value

	if (length(x) > 1) {
		mat_prob  <-  combn(probes, 2)
		mat_prob  <- matrix(as.character(mat_prob),  nc = dim(mat_prob)[2])

		confid <-  combn(as.character(pvalues), 2)
		confid <- matrix(as.character(confid),  nc = dim(confid)[2])
		y <-  apply(t(combn(x, 2)), 1, function(u){abs(diff(u))})
	 	prob_diff <-  apply(mat_prob, 2, paste, collapse = '-') 
	 	conf2 <-  apply(confid, 2, paste, collapse = '-') 
# 		prob_diff  <-  1
# 		browser()

		} else
		{
			y  <- NA
			prob_diff  <- NA
			conf2 <- NA
		}
	data.frame(phase_abs_diff = y, pair_probes = as.character(prob_diff), confidence = as.character(conf2), nb_probes = length(x))
# 	data.frame(Y = y)
}

# 	ifelse(length(x) > 1, y <-  apply(t(combn(x, 2)), 1, function(u){abs(diff(u))}), y <- NA)
# 	data.frame(Y = y)
# }


idx.rm  <- grep('---', bat[, prefix1])

 bat2  <- bat1[-idx.rm, ]
# bat2  <- bat1
bat2$ProbeID  <- as.character(bat2$ProbeID)
bat3  <-  ddply(bat2, as.name(prefix1), .fun = pro)
# browser()

bat4  <- na.omit(bat3)


v1  <- unlist(lapply(as.list(2:11), function(x) {dim(combn(x, 2))[2]}))

dplt  <- data.frame(x = 2:11, y = v1)
p_combn  <- ggplot(dplt) + geom_line(aes(x = x, y = y)) +  geom_point(aes(x = x, y = y)) + labs( x = "# of probes in Gene", y = "# of possible combinations") + theme_bw()

bat
# m <- ggplot(bat4, aes(x=phase_abs_diff)) + geom_histogram( fill = "blue")+ labs(y = "# of probe pairs", x = "Phase difference") + theme_bw()
# ggsave(paste0('phase_diff_',  prefix,'_', prefix1,'_arabidopsis.pdf'), m, width = 9, height = 7)
# 
# m1 <- ggplot(bat4, aes(x=phase_abs_diff)) + geom_histogram( fill = "blue") + facet_wrap(~nb_probes) + labs(y = "# of probe pairs", x = "Phase difference") +theme_bw()
# ggsave(paste0('phase_diff_per_nb_probes_',  prefix,'_', prefix1,'_arabidopsis.pdf'), m1, width = 9, height = 7)

ind.dup  <- !duplicated(bat4$pair_probes)
bat5  <- bat4[ind.dup, ]
bat5$ID.name  <- prefix1
names(bat5)[1]  <- "ID"

write.csv(bat5, paste0(prefix,'_', prefix1, '_analysis_arabidopsis.csv'))

ind0 <- with(bat5, phase_abs_diff != 0)

abs_dif  <- bat5$phase_abs_diff[ind0]

length(abs_dif)

ks.test(abs_dif, "punif")
}



# d1  <-  read.csv('d_arabidopsis_v1.csv')
# ind.dup  <- !duplicated(d1$pair_probes)
# d1  <- d1[ind.dup, ]
# w1  <-  read.csv('w_arabidopsis_v1.csv')
# ind.dup  <- !duplicated(w1$pair_probes)
# w1  <- w1[ind.dup, ]

# pd1 <- ggplot(d1, aes(x=phase_abs_diff)) + geom_histogram( fill = "blue")+ labs(y = "# of probe pairs", x = "Phase difference") + theme_bw()
# ggsave(paste0('phase_diff_per_nb_probes_d_arabidopsis.pdf'), pd1, width = 9, height = 7)
# pw1 <- ggplot(w1, aes(x=phase_abs_diff)) + geom_histogram( fill = "blue")+ labs(y = "# of probe pairs", x = "Phase difference") + theme_bw()
# ggsave(paste0('phase_diff_per_nb_probes_w_arabidopsis.pdf'), pw1, width = 9, height = 7)

