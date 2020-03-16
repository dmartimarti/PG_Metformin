#!/usr/bin/env Rscript

# libraries
suppressPackageStartupMessages(require(optparse));
suppressMessages(suppressWarnings(library(tidyverse)));
suppressMessages(library(readxl));
suppressMessages(suppressWarnings(library(here)));
suppressMessages(suppressWarnings(library(viridis)));

# option list
option_list = list(
	make_option(c("-i", "--input"), default = NULL, 
		help = "input file, timeseries data from Biospa"),
	make_option(c("-o", "--output"), default = NULL, 
		help = "output folder to save plots, it will be created if it doesn't exist")
); 


opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

# creates the directory
dir.create(opt$output, showWarnings = FALSE);



# Get timeseries data
time.data = read_csv(opt$input, quote = "\"") %>%
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
	select(-File, -Data);


# get factor sort order depending on the wells
lvls = naturalsort::naturalsort(unique(time.data$Well));




tsum = time.data %>%
	group_by(Strain, Plate, Metformin_mM, Well, Time_h) %>%
	summarise(Mean = mean(OD), SD = sd(OD), SE = SD/sqrt(length(OD))) %>%
	ungroup;


plate = as.integer(unique(tsum$Plate));

for (i in plate) {

	tsum %>%
		filter(Plate == i) %>%
		mutate(Well = factor(Well, levels = lvls)) %>%
		ggplot( aes(x = Time_h, y = Mean, fill = Metformin_mM, color = Metformin_mM)) +
		geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD), color = NA, alpha = 0.2) +
		geom_line(size = 1) +
		scale_x_continuous(breaks = seq(0, 24, by = 6)) +
		labs(
			 title = paste('Plate', i, sep = ' '),
			 x = 'Time, h',
			 y = 'O.D.') +
		ylim(c(0,1.4)) +
		labs(fill = "Metformin_mM") +
		facet_wrap(vars(Well), ncol = 12) +
		scale_colour_viridis(discrete = TRUE) +
		scale_fill_viridis(discrete = TRUE) +
		theme(panel.grid.major = element_blank(),
		      panel.grid.minor = element_blank(),
		      panel.background = element_rect(fill = "white", colour = "grey50")) ;
	
	ggsave(here(opt$output, paste0('growth_curves_plate',i,'.pdf')), width = 120, height = 80, units = 'mm', scale = 2, device = 'pdf')
}


