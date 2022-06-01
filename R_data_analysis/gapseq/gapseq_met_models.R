
# libraries ---------------------------------------------------------------

library(tidyverse)
library(here)
library(cowplot)
library(sybil)
library(BacArena)

# sybil::SYBIL_SETTINGS("SOLVER","cplexAPI") # (optional)


# functions -------------------------------------------------------------


getMetaboliteProduction <- function(mod) {
  require(sybil)
  require(data.table)
  
  # MTF
  sol.mtf <- optimizeProb(mod, algorithm = "mtf")
  dt.mtf  <- data.table(ex = mod@react_id,
                        mtf.flux = sol.mtf@fluxdist@fluxes[1:mod@react_num])
  dt.mtf.tmp <- copy(dt.mtf[grepl("^EX_cpd[0-9]+_e0", ex)])
  
  # FVA
  sol.fv <- fluxVar(mod, react = mod@react_id[grep("^EX_cpd[0-9]+_e0", mod@react_id)])
  
  dt <- data.table(ex       = rep(mod@react_id[grep("^EX_cpd[0-9]+_e0", mod@react_id)],2),
                   rxn.name = rep(mod@react_name[grep("^EX_cpd[0-9]+_e0", mod@react_id)],2),
                   dir      = c(rep("l",length(grep("^EX_cpd[0-9]+_e0", mod@react_id))),
                                rep("u",length(grep("^EX_cpd[0-9]+_e0", mod@react_id)))),
                   fv       = sol.fv@lp_obj)
  dt <- dcast(dt, ex + rxn.name ~ dir, value.var = "fv")[(u>1e-6 & l >= 0)]
  
  dt <- merge(dt, dt.mtf, by = "ex")
  
  return(dt[order(-mtf.flux)])
}


# read metadata -------------------------


library(readxl)
metadata = read_excel("~/Documents/MRC_postdoc/Pangenomic/metadata/MAIN_metadata.xlsx")




# get the files names
files_list = list.files(path = ".", pattern = ".RDS")

gnm_names = list.files(path = ".", pattern = ".RDS") %>% 
  str_sub(start = 1, end = -5)

# filter metadata for the bugs we want to analyse: AUS and ECOREF
meta_filt = metadata %>% 
  mutate(Genome = str_sub(fasta, start = 1, end = -7), .before = ID) %>% 
  filter(Discard == 'No') %>% 
  filter(Origin %in% c('AUS', 'ECOREF')) %>% 
  filter(Genome %in% gnm_names) %>% 
  distinct(Genome, .keep_all = T) 



# read models -------------------------------------------------------------


b1 <- readRDS("NT12060_237.RDS") # 
mg <- readRDS("NT12001_189.RDS") # for MG1655

getMetaboliteProduction(b1)[1:10]
getMetaboliteProduction(mg)[1:10]


b1_met = getMetaboliteProduction(b1)



# sybil sims --------------------------------------------------------------

### 1. read model ####

b1 <- readRDS("NT12060_237.RDS") 

### 2. run the FBA ####

optL = optimizeProb(b1, 
                    algorithm = 'fba', 
                    retOptSol= F)

### 3. get the obj function ####

optL$obj








# BacArena ----------------------------------------------------------------



# Small fix to D/L-Lactate secretion (*) and model names
# bl <- rmReact(bl, react = "EX_cpd00221_e0")
b1@mod_desc <- "E. coli 12060"
mg@mod_desc <- "E. coli MG1655"

# Construct the organism objects for BacArena simulations
eure <- Bac(b1)
bilo <- Bac(mg)

# Construct the arena size 10x10 grid cells
arena <- Arena(n = 10, m = 10)

# For each organism, populate randomly 2 grid cells in the Arena as 
# 'starter culture'
arena <- addOrg(arena, eure, amount = 2)
arena <- addOrg(arena, bilo, amount = 2)

# add substrates to arena
arena_subs <- fread("../media/LBmed.csv") # same as gapfill medium
arena_subs[, ex.rxn := paste0("EX_", compounds, "_e0")]

arena <- addSubs(arena, smax = arena_subs$maxFlux, 
                 mediac = arena_subs$ex.rxn, unit = "mM", addAnyway = T)
# Remove acetate from initial substrate list to see effect of Cross-Feeding
# arena <- rmSubs(arena, mediac = "EX_cpd00029_e0") 


# Simulation for 13 time steps
CF_sim <- simEnv(arena, time=10, sec_obj = "mtf", 
                 diff_par = T, cl_size = 6)

# Plot levels of Acetate, Buyrate, and Lactate as well as growth
par(mfrow=c(1,2))
plotCurves2(CF_sim,legendpos = "topleft",
            subs = c("cpd00211_e0","cpd00029_e0","cpd00159_e0", "cpd02799"),
            dict = list(cpd00011_e0 = "CO2", 
                        cpd00029_e0 = "Acetate", 
                        cpd00013_e0 = "Lactate",
                        cpd02799_e0 = 'Cosa'))


par(mfrow=c(1,2))
plotCurves2(CF_sim,legendpos = "topleft",
            subs = c("cpd02799_c0"),
            dict = list(cpd02799_c0 = 'Cosa'))

