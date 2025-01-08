library(tidyverse)
library(openxlsx)
#library(Taxonstand)
library(zoo)

# 

## Import, format and correct data
flora0 <- readr::read_delim("../rawdata/inserimento_flora_2024.csv", delim=";", ) %>% 
  select(-Plot) %>% 
  #mutate(Date=lubridate::as_date("Year", format=c("%d/%m/%y"))) %>% 
  mutate(Year=strptime(Year, format("%d/%m/%Y"))) %>% 
  mutate(Year=lubridate::year(Year)) %>% # PROBLEMA ----------------------------------------
  relocate(Quadrat, .before=Year) %>% 
  mutate(Type=factor(Type)) %>% 
  rename(Layer=Type, Species_original = Species) %>% 
  mutate(Treatment=str_extract(Quadrat, pattern="G|20|40")) %>% 
  mutate(Treatment=fct_recode(factor(Treatment, 
                                     levels=c("G", "20", "40")), 
                              "Gap"="G", 
                              "Margin"="20", 
                              "Interior"="40")) %>% 
  mutate(Series=str_extract(Quadrat, "W|E|N|S")) %>% 
  mutate(Series=fct_collapse(factor(Series), 
                             "E"=c("W", "E"))) %>% 
  
  ## Assign most likely species name when unsure, based on species observed in that plot over the years
  ## NOT TO CHANGE ON ACCESS DB
  mutate(Species_original=str_replace(Species_original, 
                                      pattern="Arum italicum italicum", 
                                      replacement="Arum maculatum")) %>% 
  mutate(Species_original=replace(Species_original, 
                                  list=Quadrat=="GW" & Year==2022 & 
                                    Species_original=="Polystichum lonchitis", 
                                  values="Polystichum aculeatum")) %>% 
  mutate(Species_original=replace(Species_original, 
                                  list=Quadrat=="GW" & Year==2022 & 
                                    Species_original=="Stachys sylvatica", 
                                  values="Galeopsis tetrahit")) %>%
  mutate(Species_original=str_replace(Species_original, 
                                      pattern=" gr\\. murorum", 
                                      replacement=" murorum")) %>% 
  mutate(Species_original=str_replace(Species_original, 
                                      pattern="Small fern", 
                                      replacement="Cystopteris fragilis")) %>%
  mutate(Species_original=str_replace(Species_original, 
                                      pattern="Dryopteris", 
                                      replacement="Dryopteris dilatata")) %>%
  mutate(Species_original=str_replace(Species_original, 
                                      pattern="Vicia sepium|Lathyrus$", 
                                      replacement="Lathyrus pratensis")) %>% 
  mutate(Species_original=str_replace(Species_original, 
                                      pattern="Silene italica italica|Silene$", 
                                      replacement="Silene italica")) %>%
  mutate(Species_original=str_replace(Species_original, 
                                      pattern="Plantula foglie opposte dentate", 
                                      replacement="Stachys sylvatica")) %>% 
  mutate(Species_original=str_replace(Species_original, 
                                      pattern="Labiata$", 
                                      replacement="Stachys sylvatica")) %>% 
  mutate(Species_original=str_replace(Species_original, 
                                      pattern="Acer$", 
                                      replacement="Acer pseudoplatanus")) %>% 
  mutate(Species_original=str_replace(Species_original, 
                                      pattern="Senecio rupestris|Senecio vulgaris", 
                                      replacement="Senecio vulgaris")) %>% 
  mutate(Species_original=str_replace(Species_original, 
                                      pattern="Potentilla|Potentilla micrantha", 
                                      replacement="Fragaria vesca")) %>% 
  mutate(Species_original=str_replace(Species_original, 
                                      pattern="Asarum europaeum", 
                                      replacement="Arum maculatum")) %>% 
  mutate(Species_original=str_replace(Species_original, 
                                      pattern="Arabis", 
                                      replacement="Arabis alpina")) %>% 
  mutate(Species_original=str_replace(Species_original, 
                                      pattern="Carex$", 
                                      replacement="Carex divulsa")) %>% 
  mutate(Species_original=str_replace(Species_original, 
                                      pattern="Epipactis atrorubens", 
                                      replacement="Epipactis microphylla")) %>% 
  mutate(Species_original=replace(Species_original, 
                                  list=Quadrat=="E40" & Species_original=="Cephalanthera rubra", 
                                  values="Epipactis microphylla")) %>% 
  mutate(Species_original=str_replace(Species_original, 
                                      pattern="Senecio gr. Murorum", 
                                      replacement="Hieracium murorum")) %>% 
  mutate(Species_original=str_replace(Species_original, 
                                      pattern="Corydalis cava cava|Anemone apennina apennina", 
                                      replacement="Adoxa moschatellina")) 
  

### Standardizzazione tassonomica

checklist0 <- flora0 %>% 
  distinct(Species_original) %>% 
  arrange(Species_original) %>% 
  pull(Species_original)


#tpl_output <- TPL(checklist0)
#save(tpl_output, file="../intermediate_steps/tpl_output.RData")

load(file="../intermediate_steps/tpl_output.RData")

checklist <- data.frame(
    species_original=tpl_output$Taxon,
    new_genus=tpl_output$New.Genus, 
    new_species=tpl_output$New.Species, 
    taxonomic_status_2023=tpl_output$New.Taxonomic.status, 
    family=tpl_output$Family) %>% 
  mutate(species_resolved=species_original) %>% 
  mutate(new_species=replace_na(new_species, "sp.")) %>% 
  mutate(species_tmp=paste(new_genus, new_species)) %>% 
  mutate(species_resolved=ifelse(taxonomic_status_2023!="", species_tmp, species_resolved)) %>% 
  mutate(species_resolved=ifelse(new_genus==species_resolved, paste(species_resolved, "sp."), species_resolved)) %>% 
  select(-species_tmp)

flora <- flora0 %>% 
  left_join(checklist %>% 
              select(species_original, species_resolved), 
            by=c("Species_original"="species_original")) %>% 
  relocate(species_resolved, .before=Cover) %>% 
  rename(Species_resolved=species_resolved) %>% 
  ## manually resolve last four taxa
  mutate(Species_resolved=replace(Species_resolved, 
                                  list = Species_original=="Asplenium adiantum-nigrum adiantum-nigrum", 
                                  values="Asplenium adiantum-nigrum")) %>% 
  mutate(Species_resolved=replace(Species_resolved, 
                                  list = Species_original=="Festuca altissima", 
                                  values="Festuca drymeia")) %>% 
  mutate(Species_resolved=replace(Species_resolved, 
                                  list = Species_original=="Lamium garganicum longiflorum", 
                                  values="Lamium garganicum")) %>% 
  mutate(Species_resolved=replace(Species_resolved, 
                                  list = Species_original=="Hordelymus europaeus", 
                                  values="Hordelymus europaeus"))





### Import plot-level data
header0 <- read_delim("../rawdata/inserimento_stazionali_2024.csv", delim=";") %>% 
  select(-Plot, Date=Day) %>% 
  mutate(Date=lubridate::dmy(Date) ) %>% 
  mutate(Year=lubridate::year(Date)) %>% 
  relocate(Quadrat, .before=1)  %>% 
  rowwise() %>% 
  mutate(Openness = mean(c_across(N_Emis:W_Emis), na.rm=TRUE)) %>% 
  ungroup() %>% 
  select(Quadrat, Year, Openness) %>% 
  arrange(Quadrat, Year) %>% 
  mutate(Treatment=str_extract(Quadrat, pattern="G|20|40")) %>% 
  mutate(Treatment=fct_recode(factor(Treatment, 
                                     levels=c("G", "20", "40")), 
                              "Gap"="G", 
                              "Margin"="20", 
                              "Interior"="40")) %>% 
  mutate(Series=str_extract(Quadrat, "W|E|N|S")) %>% 
  mutate(Series=fct_collapse(factor(Series), 
                             "E"=c("W", "E")))



