# libraries
library(tidyverse)
library(readxl)
# library(ComplexHeatmap)
# library(circlize)
# library(ggrepel)
# library(PFun)
# library(forcats)
# library(FactoMineR) # for PCA
# library(factoextra) # for PCA


#  "/Users/dmarti14/Documents/MRC_Postdoc/Projects/Simren/biolog/test_1"


# Get timeseries data
# It will take a while
data.b.ts = read_csv('Output/Timeseries.csv',quote = "\"") %>%
  filter(Data == '595nm_f') %>% #Select only the relevant data to speed it up
  gather(Time_s, OD, matches('\\d')) %>%
  filter(!is.na(OD)) %>%
  filter(!is.na(Media)) %>%
  mutate(Strain = as.character(Strain),
         Time_s = as.numeric(Time_s),
         Time_h = Time_s/3600,
         Row = str_match_all(Well,'[:digit:]{1,}'), #Get plate row names from well name. 
         Col = str_match_all(Well,'[:alpha:]{1,}'), #Get plate column names from well name
         Row = factor(Row, levels = 1:12), #Make them categorical variable with set order them
         Col = factor(Col, levels = LETTERS[1:8]),
         Strain = as.factor(Strain),
         Met = as.factor(Met),
         Well = as.factor(Well)) %>%
  select(Data:Media, Replicate = Replicate_y, Time_s:Col, -Replicate_x)




tsum = data.b.ts %>%
  group_by(Strain, Met, Time_h) %>%
  summarise(Mean = mean(OD), SD = sd(OD), SE = SD/sqrt(length(OD))) %>%
  ungroup 



tsum %>%
  ggplot(aes(x = Time_h, y = Mean, fill = Strain, color = Strain)) +
  geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD), color = NA, alpha = 0.2) +
  geom_line() +
  scale_x_continuous(breaks = seq(0, 24, by = 12)) +
  ylab("OD") +
  xlab("Time, h") +
  #facet_wrap(~MetaboliteU,ncol = 2)+
  facet_grid((~Met)) +
  theme_light()

quartz.save(file = here('summary', 'Growth_curves.pdf'),
	type = 'pdf', dpi = 300, height = 4, width = 10)










