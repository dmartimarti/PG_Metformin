## test script to plot growth curves from the pangenomic library with different 
## metformin concentrations (0, 50 , 100, 200)

# libraries
library(tidyverse)
library(readxl)
library(ComplexHeatmap)
library(circlize)
library(ggrepel)
library(forcats)
library(FactoMineR) 
library(factoextra)
library(here)
library(viridis)


# session options
options(width = 220)



### TIME SERIES

# Get timeseries data
time.data = read_csv('Output/Timeseries.csv', quote = "\"") %>%
	filter(Data == '595nm_f') %>%
	select(-Strain, -Pattern, -Reader) %>%
	rename(Strain = Media)  %>%
	gather(Time_s, OD, matches('\\d')) %>%
	filter(!is.na(OD)) %>% # Remove empty values if there are missmatches
	mutate(Time_s = as.numeric(Time_s),
		   Time_h = Time_s/3600,
		   Row = str_match_all(Well,'[:digit:]{1,}'), #Get plate row names from well name. 
		   Col = str_match_all(Well,'[:alpha:]{1,}'), #Get plate column names from well name
		   Row = factor(Row, levels = 1:12), #Make them categorical variable with set order them
		   Col = factor(Col, levels = LETTERS[1:8]),
		   Strain = as.factor(Strain),
		   Metformin_mM = as.factor(Metformin_mM),
		   Plate = as.factor(Plate)) %>%
	select(-File, -Data)


# get factor sort order depending on the wells
lvls = naturalsort::naturalsort(unique(time.data$Well))


# test for plot
cosa = time.data %>% filter(Strain == 'NT12098')

cosa %>% 
	ggplot(aes(x = Time_h, y = OD, colour = Metformin_mM)) +
		geom_line()


tsum = time.data %>%
	group_by(Strain, Plate, Metformin_mM, Well, Time_h) %>%
	summarise(Mean = mean(OD), SD = sd(OD), SE = SD/sqrt(length(OD))) %>%
	ungroup


plate = 1
tsum %>%
	filter(Plate == plate) %>%
	mutate(Well = factor(Well, levels = lvls)) %>%
	ggplot( aes(x = Time_h, y = Mean, fill = Metformin_mM, color = Metformin_mM)) +
	geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD), color = NA, alpha = 0.2) +
	geom_line(size = 1) +
	scale_x_continuous(breaks = seq(0, 24, by = 6)) +
	labs(
		 title = paste('Plate', plate, sep = ' '),
		 x = 'Time, h',
		 y = 'O.D.') +
	labs(fill = "Metformin_mM") +
	facet_wrap(vars(Well), ncol = 12) +
	scale_colour_viridis(discrete = TRUE) +
	scale_fill_viridis(discrete = TRUE) +
	theme(panel.grid.major = element_blank(),
	      panel.grid.minor = element_blank(),
	      panel.background = element_rect(fill = "white", colour = "grey50")) 


quartz.save(file = here('exploration', 'growth_curves_plate1.pdf'),
    type = 'pdf', dpi = 300, height = 8, width = 12)









#################
### load AUCs ###
#################



data = read_csv('Output/Summary.csv', quote = "\"") %>%
	rename(AUC_raw = `595nm_f_AUC`,
		   Sample = Strain,
		   Strain = Media) %>%
	mutate(Row = str_match_all(Well,'[:digit:]{1,}'), #Get plate row names from well name. 
		   Col = str_match_all(Well,'[:alpha:]{1,}'), #Get plate column names from well name
		   Row = factor(Row, levels = 1:12), #Make them categorical variable with set order them
		   Col = factor(Col, levels = LETTERS[1:8]),
		   Strain = as.factor(Strain),
		   Metformin_mM = as.factor(Metformin_mM),
		   Plate = as.factor(Plate)) %>%
	select(Strain, Well, Plate, Metformin_mM, Replicate, AUC_raw) %>%
	group_by(Metformin_mM, Plate, Well) %>%
	arrange(Metformin_mM, Strain, Well) %>%
	ungroup


data.wide = data %>% 
	pivot_wider(names_from = Metformin_mM, values_from = AUC_raw, names_prefix = 'met_') 


# v1
v1 = data.wide %>%
	filter(met_0 > 0.5) %>%
	mutate(x = met_50/met_0,
		   y = met_100/met_0,
		   z = met_200/met_0) 

# v1$x = normalize(v1$x)
# v1$y = normalize(v1$y)
# v1$z = normalize(v1$z)


v1 %>%
	ggtern(aes(x,y,z, colour = met_0, size = met_0)) + geom_point(alpha = 0.6)



#v2

v2 = data.wide %>%
	mutate( x = met_50/met_0,
		 	y = (met_100 - (met_50 - met_100))/met_0,
		 	z = (met_200 - (met_50 - met_200))/met_0
		 	)

v2 %>%
	filter(met_0 > 0.5) %>%
	ggtern(aes(x,y,z)) + geom_point()


normalize <- function(x) {
return ((x - min(x)) / (max(x) - min(x)))
}







v2 = data.wide %>%
	filter(met_0 > 0.5) %>%
	mutate( x = met_50/met_0,
		 	y = (met_100 - (met_50 - met_100))/met_0,
		 	z = (met_200 - (met_50 - met_200))/met_0
		 	)


v2$x = normalize(v2$x)
v2$y = normalize(v2$y)
v2$z = normalize(v2$z)

v2 %>%
	ggtern(aes(x,y,z)) + geom_point() 


v2 %>%
	filter(Strain == 'NT12129') %>%
	ggtern(aes(x,y,z)) + geom_point()








####################################
### analysis of the ASKA library ###
####################################

getwd = c("/Users/dmarti14/Documents/MRC_Postdoc/Projects/Simren/Metformin_resistance/ASKA/biolog/Output")




### TIME SERIES

# Get timeseries data
time.data = read_csv('Timeseries.csv', quote = "\"") %>%
	filter(Data == '595nm_f') %>%
	select(-Strain, -Pattern, -Reader, -Replicate) %>%
	rename(Strain = Sheet1, 
		   Metformin_mM = Sheet2,
		   Replicate = Sheet3)  %>%
	gather(Time_s, OD, matches('\\d')) %>%
	filter(!is.na(Strain)) %>% # Remove empty values if there are missmatches
	mutate(Time_s = as.numeric(Time_s),
		   Time_h = Time_s/3600,
		   Row = str_match_all(Well,'[:digit:]{1,}'), #Get plate row names from well name. 
		   Col = str_match_all(Well,'[:alpha:]{1,}'), #Get plate column names from well name
		   Row = factor(Row, levels = 1:12), #Make them categorical variable with set order them
		   Col = factor(Col, levels = LETTERS[1:8]),
		   Strain = as.factor(Strain),
		   Metformin_mM = as.factor(Metformin_mM)) %>%
	select(-File, -Data)


# get factor sort order depending on the wells
lvls = naturalsort::naturalsort(unique(time.data$Well))


# test for plot
cosa = time.data %>% filter(Strain == 'artM', Replicate == 1)

cosa %>% 
	ggplot(aes(x = Time_h, y = OD, colour = Metformin_mM)) +
		geom_line()


tsum = time.data %>%
	group_by(Strain, Metformin_mM, Well, Time_h) %>%
	summarise(Mean = mean(OD), SD = sd(OD), SE = SD/sqrt(length(OD))) %>%
	ungroup



tsum %>%
	mutate(Well = factor(Well, levels = lvls)) %>%
	ggplot( aes(x = Time_h, y = Mean, fill = Metformin_mM, color = Metformin_mM)) +
	geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD), color = NA, alpha = 0.2) +
	geom_line(size = 1) +
	scale_x_continuous(breaks = seq(0, 24, by = 6)) +
	labs(x = 'Time, h',
		 y = 'O.D.') +
	labs(fill = "Metformin_mM") +
	facet_wrap(vars(Well), ncol = 12) +
	scale_colour_viridis(discrete = TRUE) +
	scale_fill_viridis(discrete = TRUE) +
	theme(panel.grid.major = element_blank(),
	      panel.grid.minor = element_blank(),
	      panel.background = element_rect(fill = "white", colour = "grey50")) 


quartz.save(file = here('exploration', 'growth_curves_plate1.pdf'),
    type = 'pdf', dpi = 300, height = 8, width = 12)






### load AUCs


data = read_csv('Output/Summary.csv', quote = "\"") %>%
	select(-Strain, -Pattern, -Reader, -Replicate) %>%
	rename(AUC_raw = `595nm_f_AUC`,
		   Strain = Sheet1, 
		   Metformin_mM = Sheet2,
		   Replicate = Sheet3) %>%
	filter(!is.na(Strain)) %>%
	mutate(Row = str_match_all(Well,'[:digit:]{1,}'), #Get plate row names from well name. 
		   Col = str_match_all(Well,'[:alpha:]{1,}'), #Get plate column names from well name
		   Row = factor(Row, levels = 1:12), #Make them categorical variable with set order them
		   Col = factor(Col, levels = LETTERS[1:8]),
		   Strain = as.factor(Strain),
		   Metformin_mM = as.factor(Metformin_mM)) %>%
	select(Strain, Well, Metformin_mM, Replicate, AUC_raw) 

data.sum = data %>%
	group_by(Strain, Metformin_mM, Well) %>%
	summarise(Mean = mean(AUC_raw, na.rm = TRUE), 
			  SD = sd(AUC_raw, na.rm = TRUE), 
			  SE = SD/sqrt(length(AUC_raw))) %>%
	ungroup %>%
	arrange(desc(Mean))


# plot joint results
data.sum %>% 
	# filter(Metformin_mM != 200) %>%
	ggplot(aes(x = reorder(Strain, -Mean), y = Mean, colour = Metformin_mM)) +
	geom_point(size = 3) +
	geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD)) +
	labs(x = 'Strain',
		 y = 'AUC') + 
	# facet_wrap(~Metformin_mM, ncol = 1, scales = 'free_y') +
	theme(panel.grid.major = element_line(colour = "grey90"),
	      panel.background = element_rect(fill = "white", colour = "black"),
	      axis.text.x = element_text(angle = 45, hjust = 1))

quartz.save(file = here('exploration', 'AUC_by_gene_join_ASKA.pdf'),
    type = 'pdf', dpi = 300, height = 8, width = 12)




# plot joint results
met = 200
data.sum %>% 
	filter(Metformin_mM == met) %>%
	ggplot(aes(x = reorder(Strain, -Mean), y = Mean, colour = Metformin_mM)) +
	geom_point(size = 3) +
	labs(x = 'Strain',
		 y = 'AUC') +
	geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD)) +
	# facet_wrap(~Metformin_mM, ncol = 1, scales = 'free_y') +
	theme(panel.grid.major = element_line(colour = "grey90"),
	      panel.background = element_rect(fill = "white", colour = "black"),
	      axis.text.x = element_text(angle = 45, hjust = 1))

quartz.save(file = here('exploration', paste0('AUC_by_gene_',met,'_ASKA.pdf')),
    type = 'pdf', dpi = 300, height = 5, width = 11)



data.sum %>%
	select(Strain, Metformin_mM, Mean) %>%
	pivot_wider(names_from = Metformin_mM, names_prefix = "M", values_from = Mean) %>%
	mutate(Norm50  = (M50/M0),
		   Norm100 = (M100/M0),
		   Norm200 = (M200/M0)) %>%
	ggplot(aes(x = Norm50, y = Norm100)) +
	geom_point(aes(size = Norm200 * 2, fill = Norm200), shape = 21) +
	scale_size(range = c(3, 8)) +
	xlim(0.5, 0.85) +
	ylim(0.3, 0.85) +
	labs(x = 'Norm AUC at 50 mM',
		 y = 'Norm AUC at 100 mM') +
	scale_fill_viridis(name = 'Norm AUC at\n200 mM', option = "magma") +
	geom_text_repel(aes(label = Strain), size = 3.5) +
	theme_classic() +
	guides(size = FALSE)
 
quartz.save(file = here('exploration', 'Scatter_norm_AUC_ASKA.pdf'),
    type = 'pdf', dpi = 300, height = 8, width = 10)

# write gene list table for enrichment
write.table(unique(data.sum[,1]), file = 'gene_names_ASKA.txt', quote = F, row.names = F, col.names = F)










####################################
### analysis of the KEIO library ###
####################################

getwd = c("/Users/dmarti14/Documents/MRC_Postdoc/Projects/Simren/Metformin_resistance/KEIO/biolog/Output")




### TIME SERIES

# Get timeseries data
time.data = read_csv('Output/Timeseries.csv', quote = "\"") %>%
	filter(Data == '595nm_f') %>%
	select(-Strain, -Pattern, -Reader, -Replicate) %>%
	rename(Strain = pattern, 
		   Metformin_mM = conc,
		   Replicate = rep)  %>%
	gather(Time_s, OD, matches('\\d')) %>%
	filter(!is.na(Strain)) %>% # Remove empty values if there are missmatches
	mutate(Time_s = as.numeric(Time_s),
		   Time_h = Time_s/3600,
		   Row = str_match_all(Well,'[:digit:]{1,}'), #Get plate row names from well name. 
		   Col = str_match_all(Well,'[:alpha:]{1,}'), #Get plate column names from well name
		   Row = factor(Row, levels = 1:12), #Make them categorical variable with set order them
		   Col = factor(Col, levels = LETTERS[1:8]),
		   Strain = as.factor(Strain),
		   Metformin_mM = as.factor(Metformin_mM)) %>%
	select(-File, -Data)


# get factor sort order depending on the wells
lvls = naturalsort::naturalsort(unique(time.data$Well))


# test for plot
cosa = time.data %>% filter(Strain == 'envZ', Replicate == 1)
cosa %>% 
	ggplot(aes(x = Time_h, y = OD, colour = Metformin_mM)) +
		geom_line()


tsum = time.data %>%
	group_by(Strain, Metformin_mM, Well, Time_h) %>%
	summarise(Mean = mean(OD), SD = sd(OD), SE = SD/sqrt(length(OD))) %>%
	ungroup



tsum %>%
	mutate(Well = factor(Well, levels = lvls)) %>%
	ggplot( aes(x = Time_h, y = Mean, fill = Metformin_mM, color = Metformin_mM)) +
	geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD), color = NA, alpha = 0.2) +
	geom_line(size = 1) +
	scale_x_continuous(breaks = seq(0, 24, by = 6)) +
	labs(x = 'Time, h',
		 y = 'O.D.') +
	labs(fill = "Metformin_mM") +
	facet_wrap(vars(Strain), ncol = 12) +
	scale_colour_viridis(discrete = TRUE) +
	scale_fill_viridis(discrete = TRUE) +
	theme(panel.grid.major = element_blank(),
	      panel.grid.minor = element_blank(),
	      panel.background = element_rect(fill = "white", colour = "grey50")) 

quartz.save(file = here('exploration', 'growth_curves_KEIO.pdf'),
    type = 'pdf', dpi = 300, height = 8, width = 12)




### load AUCs


data = read_csv('Output/Summary.csv', quote = "\"") %>%
	select(-Strain, -Pattern, -Reader, -Replicate) %>%
	rename(AUC_raw = `595nm_f_AUC`,
		   Strain = pattern, 
		   Metformin_mM = conc,
		   Replicate = rep) %>%
	filter(!is.na(Strain)) %>%
	mutate(Row = str_match_all(Well,'[:digit:]{1,}'), #Get plate row names from well name. 
		   Col = str_match_all(Well,'[:alpha:]{1,}'), #Get plate column names from well name
		   Row = factor(Row, levels = 1:12), #Make them categorical variable with set order them
		   Col = factor(Col, levels = LETTERS[1:8]),
		   Strain = as.factor(Strain),
		   Metformin_mM = as.factor(Metformin_mM)) %>%
	select(Strain, Well, Metformin_mM, Replicate, AUC_raw) 

data.sum = data %>%
	group_by(Strain, Metformin_mM, Well) %>%
	summarise(Mean = mean(AUC_raw, na.rm = TRUE), 
			  SD = sd(AUC_raw, na.rm = TRUE), 
			  SE = SD/sqrt(length(AUC_raw))) %>%
	ungroup %>%
	arrange(desc(Mean))


# plot joint results
data.sum %>% 
	# filter(Metformin_mM != 200) %>%
	ggplot(aes(x = reorder(Strain, -Mean), y = Mean, colour = Metformin_mM)) +
	geom_point(size = 3) +
	geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD)) +
	labs(x = 'Strain',
		 y = 'AUC') + 
	# facet_wrap(~Metformin_mM, ncol = 1, scales = 'free_y') +
	theme(panel.grid.major = element_line(colour = "grey90"),
	      panel.background = element_rect(fill = "white", colour = "black"),
	      axis.text.x = element_text(angle = 45, hjust = 1))

quartz.save(file = here('exploration', 'AUC_by_gene_join_KEIO.pdf'),
    type = 'pdf', dpi = 300, height = 8, width = 12)




# plot joint results
met = 200
data.sum %>% 
	filter(Metformin_mM == met) %>%
	ggplot(aes(x = reorder(Strain, -Mean), y = Mean, colour = Metformin_mM)) +
	geom_point(size = 3) +
	labs(x = 'Strain',
		 y = 'AUC') +
	geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD)) +
	# facet_wrap(~Metformin_mM, ncol = 1, scales = 'free_y') +
	theme(panel.grid.major = element_line(colour = "grey90"),
	      panel.background = element_rect(fill = "white", colour = "black"),
	      axis.text.x = element_text(angle = 45, hjust = 1))

quartz.save(file = here('exploration', paste0('AUC_by_gene_',met,'_KEIO.pdf')),
    type = 'pdf', dpi = 300, height = 5, width = 11)



data.sum %>%
	select(Strain, Metformin_mM, Mean) %>%
	pivot_wider(names_from = Metformin_mM, names_prefix = "M", values_from = Mean) %>%
	mutate(Norm50  = (M50/M0),
		   Norm100 = (M100/M0),
		   Norm200 = (M200/M0)) %>%
	ggplot(aes(x = Norm50, y = Norm100)) +
	geom_point(aes(size = Norm200, fill = Norm200), shape = 21) +
	xlim(0.4, 0.85) +
	ylim(0.4, 0.85) +
	labs(x = 'Norm AUC at 50 mM',
		 y = 'Norm AUC at 100 mM') +
	scale_fill_viridis(name = 'Norm AUC at\n200 mM', option = "magma") +
	geom_text_repel(aes(label = Strain), size = 3.5) +
	theme_classic() +
	scale_size(range = c(3, 8)) +
	guides(size = FALSE)
 
quartz.save(file = here('exploration', 'Scatter_norm_AUC_KEIO.pdf'),
    type = 'pdf', dpi = 300, height = 8, width = 10)


# write table with gene names
write.table(unique(data.sum[,1]), file = 'gene_names_KEIO.txt', quote = F, row.names = F, col.names = F)










