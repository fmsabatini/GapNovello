library(tidyverse)
library(openxlsx)
library(Taxonstand)

## Import, format and correct data
flora0 <- readr::read_delim("../rawdata/inserimento_flora.csv", delim=";", ) %>% 
  select(-Plot) %>% 
  mutate(Date=lubridate::as_date(Year, format=c("%d/%m/%y") )) %>% 
  mutate(Year=lubridate::year(Date)) %>% 
  relocate(Quadrat, .before=Year) %>% 
  mutate(Type=factor(Type)) %>% 
  mutate(Species=replace(Species, 
                         list=Species=="Vincetoxium", 
                         values="Vincetoxicum hirundinaria")) %>% 
  mutate(Species=str_replace(Species, pattern=" gr. ", replacement=" ")) %>% 
  rename(Layer=Type, Species_original = Species) %>% 
  ## Correct some typos in raw data
  mutate(Quadrat=replace(Quadrat, list=Quadrat=="S40\r\nS40", values="S40")) %>% 
  mutate(Layer=replace(Layer, 
                       list= (Quadrat=="N40") & 
                             (Species_original=="Fagus sylvatica sylvatica") & 
                             (Year==2016), 
                       values= "O")) %>% 
  mutate(Layer=replace(Layer, 
                       list= (Quadrat=="S20") & 
                         (Species_original=="Fagus sylvatica sylvatica") & 
                         (Year==2012), 
                       values= "O")) %>% 
  mutate(Quadrat=factor(Quadrat)) %>% 
  mutate(Species_original=str_replace(Species_original, 
                                      pattern="Arum italicum italicum", 
                                      replacement="Arum maculatum"))

### Standardizzazione tassonomica

checklist0 <- flora0 %>% 
  distinct(Species_original) %>% 
  arrange(Species_original) %>% 
  pull(Species_original)


# tpl_output <- TPL(checklist0)
# save(tpl_output, file="../intermediate_steps/tpl_output.RData")

load(file="../intermediate_steps/tpl_output.RData")

checklist <- data.frame(
    species_original=tpl_output$Taxon,
    new_genus=tpl_output$New.Genus, 
    new_species=tpl_output$New.Species, 
    taxonomic_status_2023=tpl_output$New.Taxonomic.status) %>% 
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
  rename(Species_resolved=species_resolved)


### Graph of species abundances with time (by plot)
ggplot(data=flora %>% 
         filter(Layer=="U") %>% 
         filter(Quadrat=="N20") %>% 
         mutate(Species_resolved=factor(Species_resolved))) + 
  geom_line(aes(x=Year, y=Cover, group=Species_resolved, col=Species_resolved)) #+ 
  #scale_color_brewer(type="div")
  #ggplot2::facet_wrap(~Quadrat, nrow = 3, ncol=3)
  
  





