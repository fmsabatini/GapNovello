library(tidyverse)
library(patchwork)
library(zoo)
library(nlme)
#library(lme4)


load("../intermediate_steps/GapNovello_Data.RData")


### Ancillary functions ####
## Calculate rolling mean of change
ch <- function(x){x[length(x)]-x[1]}

## Pool species cover across consecutive years
pool <- function(x){ifelse(sum(!is.na(x))>=2, max(x, na.rm=T), NA)} # to fix omissions
commission <- function(x){ifelse(sum(!is.na(x))==1, 0, NA)} # to fix commissions

##  # Combine cover accounting for layers
combine.cover <- function(x){
  while (length(x)>1){
    x[2] <- x[1]+(100-x[1])*x[2]/100
    x <- x[-1]
  }
  return(x)
}


### Figure 2 - Study area ####
library(sf)
library(rnaturalearth)
library(ggspatial)
library(basemaps)


europe <- ne_countries(scale = "large", returnclass = "sf", continent="europe")   

gap_coordinates <- c(13.511308, 42.500128)
mybbox <- sf::st_bbox(europe)
mybbox[[1]] <- gap_coordinates[1]-.1
mybbox[[3]] <- gap_coordinates[1]+.1
mybbox[[2]] <- gap_coordinates[2]-.1
mybbox[[4]] <- gap_coordinates[2]+.1
mybbox <- mybbox %>% 
  sf::st_as_sfc()


## Plot
ylabs <- seq(35, 50, by=2)
xlabs <- seq(0, 30, by=2)

left <- ggplot() +
  theme_minimal() +
  geom_sf(data=europe, fill=gray(0.8), col=gray(0.6)) + 
  geom_sf(data=mybbox, aes(label="a"), 
          col="black", 
          fill=NA, 
          linewidth=.8) + 
  coord_sf(xlim=c(7, 19), ylim=c(36, 47)) + 
  scale_y_continuous(breaks=ylabs) + 
  scale_x_continuous(breaks=xlabs) + 
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1)) + 
  ggspatial::annotation_scale()

## zoom in 
gap_coordinates_sf <- data.frame(x=gap_coordinates[1],
                              y=gap_coordinates[2]) %>% 
  st_as_sf(coords=c(1:2), crs=4326) %>% 
  st_transform(crs = 32633) 

plot_coordinates <- data.frame(
  Quadrat <-  c("S40","N40","E40","GW","GN","GS","N20","E20","S20"),
  x <- st_coordinates(gap_coordinates_sf)[1] + c(0, 0, +40, -5, 0, 0, 0, 20, 0),
  y <- st_coordinates(gap_coordinates_sf)[2] + c(-40, 40, 0, 0, 5, -5, +20, 0, -20))   %>% 
  st_as_sf(coords=2:3, crs=32633) %>% 
  st_buffer(dist = 2.4, endCapStyle = "SQUARE") 


mybbox_small <- st_bbox(gap_coordinates_sf)
mybbox_small[1] <- mybbox_small[1]-50
mybbox_small[2] <- mybbox_small[2]-50
mybbox_small[3] <- mybbox_small[3]+50
mybbox_small[4] <- mybbox_small[4]+50
mybbox_small <- mybbox_small %>% 
  st_as_sfc(crs=32633)



## import georeferences screenshot from ORTOFOTO abruzzo 2018
myrast <- terra::rast("../figure_tables/ORTOFOTO_ABRUZZO_gapNovello_georef.tif") %>% 
  terra::project("EPSG:32633")

right <- ggplot() +
  geom_sf(data=mybbox_small, fill=NA) + 
  terrainr::geom_spatial_rgb(data=myrast, 
                   mapping = aes(x = x,
                                 y = y,
                                 r = red,
                                 g = green,
                                 b = blue)
  ) + 
  geom_sf(data=plot_coordinates, fill="deepskyblue", alpha=0.5, col="deepskyblue", linewidth=1.5) + 
    theme_minimal() +
  coord_sf(xlim=st_bbox(mybbox_small)[c(1,3)], 
           ylim=st_bbox(mybbox_small)[c(2,4)]) + 
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1), 
        axis.text = element_blank(), 
        axis.title = element_blank()) + 
  ggspatial::annotation_scale(aes(text_col="white", text_cex=2), location="bl") + 
  ggspatial::annotation_north_arrow(height=unit(1, "cm"), 
                                    width=unit(0.5, "cm"), 
                                    location="tl")

(both <- left + right)

ggsave(plot = both, 
       file.path("../figure_tables/", "Figure2_Studyarea.png"), 
       dpi=300, bg="white", 
       type="cairo", width = 8, height=4, units = "in")



### RQ1 - Account for pseudoturnover  ####
## Pool species entries using a 3-year moving window
flora_pooled <- flora %>% 
  filter(Layer=="U") %>% 
  dplyr::select(Year, Quadrat, Species_resolved, Cover) %>% 
  complete(Year, nesting(Quadrat, Species_resolved)) %>% 
  arrange(Quadrat, Species_resolved, Year) %>% 
  group_by(Quadrat, Species_resolved) %>%
  #round one control for omission error
  mutate(Cover_roll=rollapply(Cover, 3, pool, align="center", fill=NA)) %>% 
  # in the first and last year of the time series I keep the observed values
  mutate(Cover_roll2=coalesce(Cover, Cover_roll)) %>% 
  # round two - control for commission error
  mutate(Cover_roll3=rollapply(Cover_roll2, 3, commission, align="center", fill=NA)) %>% 
  mutate(Cover_roll4=coalesce(Cover_roll3, Cover_roll2)) %>% 
  mutate(Cover_roll5=replace(Cover_roll4, list=is.na(Cover_roll4), values=0)) %>% 
  dplyr::select(Year:Cover, Cover_roll=Cover_roll5) %>% 
  ungroup()




## pool species across consecutive years to reduce omission error
########## pool species cover across n consecutive years
flora_wide_pooled <- flora_pooled %>% 
  dplyr::select(-Cover) %>% 
  pivot_wider(names_from=Species_resolved, values_from = "Cover_roll", values_fill = 0) %>% 
  unite("Quadrat_year", Quadrat, Year) %>% 
  arrange(Quadrat_year) %>% 
  #View()
  column_to_rownames("Quadrat_year")

### Same for matrix T
# Convert species x plot level to wide format
flora_wide <- flora %>% 
  filter(Layer=="U") %>% 
  ## flatten vegetation layers
  #group_by(Quadrat, Species_resolved, Year) %>% 
  #summarize(Cover=combine.cover(Cover)) %>% 
  dplyr::select(Quadrat, Year, Species_resolved, Cover) %>% 
  pivot_wider(names_from=Species_resolved, values_from = "Cover", values_fill = 0) %>% 
  unite("Quadrat_year", Quadrat, Year) %>% 
  arrange(Quadrat_year) %>% 
  column_to_rownames("Quadrat_year")


### Approach to decompose the species x site matrix into
### TOT variation = Real turnover + Pseudoturnover matrix
### TOT variation == Flora_wide
### Real turnover == flora_wide_pooled
### Pseudoturnover == Flora_wide - flora_wide_pooled

tot <- vegan::decostand(flora_wide, "hellinger")
real <- vegan::decostand(flora_wide_pooled, "hellinger")
pseudo <- tot-real

#ss.tot <- adespatial::beta.div(Y=tot, method = "euclidean")$beta
#ss.real <- adespatial::beta.div(Y=real, method = "euclidean")$beta
#ss.pseudo <- adespatial::beta.div(Y=pseudo, method = "euclidean")$beta

(ss.tot <- sum((t(tot)-rowMeans(t(tot)))^2))
(ss.real <- sum((t(real)-rowMeans(t(real)))^2))
(ss.pseudo <- sum((t(pseudo)-rowMeans(t(pseudo)))^2))



## double check values
# Calculate mean abundance for each site for the total matrix
site_means_total <- colMeans(tot)
# Calculate deviations from site means for the total matrix
deviations_total <- tot - outer(rep(1, nrow(tot)), site_means_total, FUN = "*")
# Calculate variance for the total matrix
variance_total <- sum(deviations_total^2)
# Display result
cat("Variance of Total Matrix (A + B):", variance_total, "\n")


# Assuming your initial matrix is the sum of matrices real + pseudo
# Calculate mean abundance for each site for A and B
site_means_real <- colMeans(real)
site_means_pseudo <- colMeans(pseudo)

# Calculate deviations from site means for A and B
deviations_real <- real - outer(rep(1, nrow(real)), site_means_real, FUN = "*")
deviations_pseudo <- pseudo - outer(rep(1, nrow(pseudo)), 
                                    site_means_pseudo, FUN = "*")

# Calculate variances for A and B
variance_real <- sum(deviations_real^2)
variance_pseudo <- sum(deviations_pseudo^2)

# Calculate covariance between A and B
covariance_real_pseudo <- sum(deviations_real * deviations_pseudo)

# Calculate total variance for the combined matrix
total_variance <- variance_real + variance_pseudo + 2 * covariance_real_pseudo

# Display results
cat("Variance of A:", variance_real, "\n")
cat("Variance of B:", variance_pseudo, "\n")
cat("Covariance between A and B:", covariance_real_pseudo, "\n")
cat("Total Variance (A + B + 2AB):", total_variance, "\n")

variance_real + variance_pseudo + 2*covariance_real_pseudo



real_dist <- vegan::vegdist(real, method="euclidean") ##already transformed into hellinger
mydata <- data.frame(qyear=names(real_dist)) %>%
  separate(qyear, into=c("Quadrat", "Year"), sep="_") %>% 
  left_join(flora %>% 
              distinct(Quadrat, Treatment), 
            by="Quadrat")
SS.deco.real <- vegan::adonis2(formula= real_dist ~ Treatment*Year + Quadrat, data=mydata, strata=mydata$Quadrat, by="term")



tot_dist <- vegan::vegdist(tot, method="euclidean")
mydata <- data.frame(qyear=names(tot_dist)) %>% 
  separate(qyear, into=c("Quadrat", "Year"), sep="_") %>% 
  left_join(flora %>% 
              distinct(Quadrat, Treatment), 
            by="Quadrat")
SS.deco.tot <- vegan::adonis2(formula= tot_dist ~ Treatment*Year + Quadrat, data=mydata, strata=mydata$Quadrat, by="term")


#### Account for therophytes and spring ephemerals ####
thero_list <- checklist_elle %>% 
  filter(GF=="T") %>% 
  pull(Species_resolved)

spring_list <- checklist_elle %>% 
  filter(Spring_ephemeral) %>% 
  pull(Species_resolved)

### Exclude therophytes from moving window
flora_wide_nothero <- flora_pooled %>% 
  rowwise() %>% 
  mutate(isthero=Species_resolved %in% thero_list) %>% 
  mutate(Cover_roll=replace(Cover_roll, 
                            list=isthero,
                            values=Cover)) %>% 
  ungroup() %>% 
  dplyr::select(-Cover, -isthero) %>% 
  mutate(Cover_roll=replace(Cover_roll, list=is.na(Cover_roll), values=0)) %>%
  pivot_wider(names_from=Species_resolved, values_from = "Cover_roll", values_fill = 0) %>% 
  unite("Quadrat_year", Quadrat, Year) %>% 
  arrange(Quadrat_year) %>% 
  #View()
  column_to_rownames("Quadrat_year")

nothero <- vegan::decostand(flora_wide_nothero, "hellinger")

(ss.nothero <- adespatial::beta.div(Y=nothero, method = "euclidean")$beta)


nothero_dist <- vegan::vegdist(nothero, method="euclidean") ##already transformed into hellinger
mydata <- data.frame(qyear=names(nothero_dist)) %>% 
  separate(qyear, into=c("Quadrat", "Year"), sep="_") %>% 
  left_join(flora %>% 
              distinct(Quadrat, Treatment), 
            by="Quadrat")
SS.deco.nothero <- vegan::adonis2(formula= real_dist ~ Treatment*Year + Quadrat, data=mydata, strata=mydata$Quadrat, by="term")

variance_nothero <- ss.nothero[1]


#### Apply moving window exclusively to spring ephemerals ####
flora_wide_nospring <- flora_pooled %>% 
  rowwise() %>% 
  mutate(isspring=Species_resolved %in% spring_list) %>% 
  mutate(Cover_roll=replace(Cover_roll, 
                            list=isspring,
                            values=Cover)) %>% 
  ungroup() %>% 
  dplyr::select(-Cover, -isspring) %>% 
  mutate(Cover_roll=replace(Cover_roll, list=is.na(Cover_roll), values=0)) %>%
  pivot_wider(names_from=Species_resolved, values_from = "Cover_roll", values_fill = 0) %>% 
  unite("Quadrat_year", Quadrat, Year) %>% 
  arrange(Quadrat_year) %>% 
  #View()
  column_to_rownames("Quadrat_year")

nospring <- vegan::decostand(flora_wide_nospring, "hellinger")

(ss.nospring <- adespatial::beta.div(Y=nospring, method = "euclidean")$beta)


nospring_dist <- vegan::vegdist(nospring, method="euclidean") ##already transformed into hellinger
mydata <- data.frame(qyear=names(nospring_dist)) %>% 
  separate(qyear, into=c("Quadrat", "Year"), sep="_") %>% 
  left_join(flora %>% 
              distinct(Quadrat, Treatment), 
            by="Quadrat")
SS.deco.nospring <- vegan::adonis2(formula= real_dist ~ Treatment*Year + Quadrat, data=mydata, strata=mydata$Quadrat, by="term")

variance_nospring <- ss.nospring[1]





parted <- data.frame(
  Fraction=c(rownames(SS.deco.tot), "Pseudoturnover"),
  Share_tot=c(SS.deco.tot$SumOfSqs, 0),
  Share_real=c(SS.deco.real$SumOfSqs, total_variance-variance_real), 
  Share_nothero=c(SS.deco.nothero$SumOfSqs, total_variance-variance_nothero),
  Share_nospring=c(SS.deco.nospring$SumOfSqs, total_variance-variance_nospring))  %>% 
  filter(Fraction != "Total") %>% 
  bind_cols({.} %>% 
              summarize_at(.vars = vars(Share_tot:Share_nospring),
                           .funs = list("sum"=~sum(.))) ) %>% 
  mutate(Share_tot_perc=round(Share_tot/Share_tot_sum*100,2),
         Share_real_perc=round(Share_real/Share_real_sum*100,2), 
         Share_nothero_perc=round(Share_nothero/Share_nothero_sum*100,2), 
         Share_nospring_perc=round(Share_nospring/Share_nospring_sum*100,2), 
  ) 


parted_long <- parted %>% 
  dplyr::select(-ends_with("_sum")) %>% 
  pivot_longer(-Fraction, names_to="Variable", values_to="Explained Variation (%)") %>% 
  mutate(Variable=str_remove(Variable, pattern="^Share_")) %>% 
  separate(Variable, into=c("model", "Abs_Perc")) %>% 
  mutate(Abs_Perc=replace_na(Abs_Perc, "abs")) %>% 
  mutate(Fraction=factor(Fraction, 
                         levels=c("Treatment", 
                                  "Treatment:Year", 
                                  "Year",
                                  "Residual",
                                  "Quadrat", 
                                  "Pseudoturnover",
                                  "Total"))) %>% 
  mutate(Fraction=fct_recode(Fraction,
                              "Residual (Quadrat x Year)" = "Residual")) %>% 
  mutate(model=factor(model, levels=c( "tot",  "nothero", "nospring", "real"), labels=c("Total","No T", "SE only", "Real"))) %>% 
  arrange(Fraction)

## format table
library(flextable)

flextable(parted_long %>% 
               filter(Abs_Perc=="perc") %>% 
               select(-Abs_Perc) %>% 
               pivot_wider(names_from = model, values_from=`Explained Variation (%)`)) %>% 
  theme_zebra() %>% 
  save_as_docx(path = "../figure_tables/tableS1.docx")


#### Figure 3 ####

ggplot(data=parted_long %>% 
         filter(Abs_Perc=="perc") %>% 
         filter(Fraction!="Total") %>% 
         filter(model %in% c("Total", "Real")) %>% 
         group_by(model) %>% 
         mutate(label_y=cumsum(`Explained Variation (%)`))) +
  geom_col(aes(x=model, y=`Explained Variation (%)`, fill=Fraction, group=model), width=0.5) +
  geom_text(aes(x=model, y=label_y, label = ifelse(round(`Explained Variation (%)`,1) != 0, 
                                                   round(`Explained Variation (%)`,1),
                                                   NA), group=model), 
            hjust = 1.5, 
            colour = "white") +
  theme_minimal() + 
  scale_fill_brewer(palette = "Paired") + 
  #guides(fill = guide_legend(reverse=TRUE)) + 
  scale_x_discrete(name=NULL) + 
  coord_flip() + 
  theme(legend.position = "bottom", 
        legend.title = element_blank()) 
#axis.text.y = element_blank())

ggsave(last_plot(), file="../figure_tables/Figure2_ExplainedVariation_pseudoturnover.png", 
       width = 7, height=3, units="in",dpi = 300, bg = "white")


#### Figure S2 ####

ggplot(data=parted_long %>% 
         filter(Abs_Perc=="perc") %>% 
         filter(Fraction!="Total") %>% 
         #filter(model %in% c("Total", "Real"))
         group_by(model) %>% 
         mutate(label_y=cumsum(`Explained Variation (%)`))
) +
  geom_col(aes(x=model, y=`Explained Variation (%)`, fill=Fraction, group=model), width=0.5) +
  geom_text(aes(x=model, y=label_y, label = round(`Explained Variation (%)`,1), group=model), 
            hjust = 1.5, 
            colour = "white") +
  theme_minimal() + 
  scale_fill_brewer(palette = "Paired") + 
  #guides(fill = guide_legend(reverse=TRUE)) + 
  scale_x_discrete(name=NULL) + 
  coord_flip() + 
  theme(legend.position = "bottom", 
        legend.title = element_blank()) 
#axis.text.y = element_blank())

ggsave(last_plot(), file="../figure_tables/AppendixS2_ExplainedVariation_alternatives.png", 
       width = 7, height=3, units="in",dpi = 300, bg = "white")



### RQ2: Rate of Change ####
#### Graphs with raw data ####
### Graph of species richness with time by plot
flora_pooled2 <- flora_pooled %>% 
  filter(Cover_roll!=0) %>% 
  left_join(flora %>% 
              distinct(Quadrat, Treatment, Series), 
            by="Quadrat") %>% 
  # change below here for graphs with T data
  dplyr::select(-Cover) #%>% 
  #filter(!is.na(Cover)) %>% 
  #mutate(Cover_roll = Cover)
  


#### Figure 4 ####
(topleft <- ggplot(data=flora_pooled2 %>% 
                     group_by(Quadrat, Treatment, Series, Year) %>% 
                     #distinct(Species_resolved, .keep_all=T) %>% 
                     summarize(SR=n()) ) + 
   geom_line(aes(x=Year, y=SR, group=Quadrat, col=Treatment, linetype=Series)) + 
   scale_x_continuous(breaks = seq(2012, 2024, by=2)) + 
   scale_y_continuous(name="Species Richness") + 
   scale_color_brewer(palette = "Dark2", direction = -1) + 
   theme_bw()
)

### Graph of herb-layer total cover with time by plot
(topright <- ggplot(data=flora_pooled2 %>% 
                      #filter(Layer!="O") %>% 
                      group_by(Quadrat, Treatment, Series, Year) %>% 
                      summarize(Tot_cover=sum(Cover_roll)) %>%
                      arrange(Quadrat, Year)) + 
    geom_line(aes(x=Year, y=Tot_cover, group=Quadrat, col=Treatment, linetype=Series)) + 
    scale_x_continuous(breaks = seq(2012, 2024, by=2)) + 
    scale_y_continuous(name="Total Herb Layer Cover (%)") +
    scale_color_brewer(palette = "Dark2", direction = -1) + 
    theme_bw()
)

mylegend <- ggpubr::get_legend(topright + 
                                 theme(legend.position = "bottom"))

#### Helical Graphs ####
# See: https://www.sciencedirect.com/science/article/pii/S157495412200406X?casa_token=_7VEA2W0D6kAAAAA:SEYaqqXXmvrEoWgMfLFP18S4vSB_iaDXCnxPbz6Gp1yvKnnXbTdUMO_D_o2_S_UeSsGi0S2G
flora_sr <- flora_pooled2 %>% 
  #filter(Layer!="O") %>% 
  group_by(Quadrat, Treatment, Series, Year) %>% 
  summarize(SR=n(), Tot_cover=sum(Cover_roll)) %>%
  arrange(Quadrat, Year) %>% 
  mutate(SR_roll=rollapply(SR, 3, mean, align="right", fill=NA)) %>% 
  mutate(SRch_roll=rollapply(SR, 3, ch, align="right", fill=NA)) %>% 
  mutate(Cov_roll=rollapply(Tot_cover, 3, mean, align="right", fill=NA)) %>% 
  mutate(Covch_roll=rollapply(Tot_cover, 3, ch, align="right", fill=NA)) %>% 
  mutate(row_id=row_number()) %>% 
  # Join dissimilarity values between year and year +1
  #  left_join(flora_beta2 %>% 
  #               filter(Year1== Year2-1),
  #            by=c("Quadrat", "Year"="Year2", "Treatment")) %>% 
  filter(!is.na(SR_roll)) %>% 
  # smoothed rolling mean of SR, change in SR, and sorensen dissimilarity
  mutate(SR_roll_lo=predict(loess(SR_roll~row_id, span=0.75))) %>%
  mutate(SRch_roll_lo=predict(loess(SRch_roll~row_id, span=0.75))) %>% 
  mutate(Cov_roll_lo=predict(loess(Cov_roll~row_id, span=0.75))) %>% 
  mutate(Covch_roll_lo=predict(loess(Covch_roll~row_id, span=0.75))) #%>% 
#  mutate(sor_roll_lo=predict(loess(beta_sor~row_id)))


## Graph of smoothed rolling mean of change in SR vs SR
(bottomleft <- ggplot(data=flora_sr, 
                      aes(x=SRch_roll_lo, y=SR_roll_lo, 
                          group=Quadrat, col=Treatment, linetype=Series)) + 
    geom_path(arrow = arrow(type="closed", length=unit(0.1, "inches"))) + 
    geom_text(label=flora_sr$Year, size=3, hjust=-0.2, show.legend=F,
              col="black", check_overlap = T, alpha=1, fontface="bold") +
    scale_color_brewer(palette = "Dark2", direction = -1) + 
    scale_y_continuous(name="Species Richness") +
    scale_x_continuous(name="Yearly Change in Richness", expand = expansion(add = 1)) +
    theme_bw()
)

## Graph of smoothed rolling mean of change in Cover vs Cover
(bottomright <- ggplot(data=flora_sr, 
                       aes(x=Covch_roll_lo, y=Cov_roll_lo, 
                           group=Quadrat, col=Treatment, linetype=Series)) + 
    geom_path(arrow = arrow(type="closed", length=unit(0.1, "inches"))) + 
    geom_text(label=flora_sr$Year, size=3, hjust=-0.2, show.legend=F,
              col="black", check_overlap = T, alpha=1, fontface="bold") + 
    scale_color_brewer(palette = "Dark2", direction = -1) + 
    scale_y_continuous(name="Total Cover (%)") +
    scale_x_continuous(name="Yearly Change Cover (%)", expand = expansion(add = 10)) +
    theme_bw()
)


Figure4 <- ((topleft + theme(legend.position="none")  |
               topright+ theme(legend.position="none")) / 
              (bottomleft + theme(legend.position="none") |
                 bottomright + theme(legend.position="none"))) / 
  mylegend + plot_layout(height=c(6,6,1)) + plot_annotation(tag_levels = list(c("a", "b", 'c', "d", "")))

ggsave(Figure4, file="../figure_tables/Figure4_SR_Cover_Change.png", 
       width = 9, height=8, units="in",dpi = 300)
#### end of Figure 3 

# calculate summaries for text
flora_pooled2 %>% 
  count(Quadrat, Year, Treatment) %>% 
  group_by(Treatment, Year) %>% 
  summarize(mean(n)) %>% 
  View()

flora_pooled2 %>% 
  count(Quadrat, Year, Treatment) %>% 
  rename(SR=n) %>% 
  #group_by(Treatment, Year) %>% 
  #summarize(SR=mean(SR)) %>% 
  group_by(Quadrat, Treatment) %>% 
  mutate(DeltaSR=SR-lag(SR)) %>%
  group_by(Treatment) %>% 
  arrange(desc(abs(DeltaSR))) %>% 
  slice_head(n=5) %>% 
  View()

#same for cover
flora_pooled2 %>% 
  group_by(Quadrat, Treatment, Year) %>% 
  summarize(Cover=sum(Cover_roll)) %>% 
  group_by(Treatment, Year) %>% 
  summarize(mean(Cover)) %>% 
  View()


flora_pooled2 %>% 
  group_by(Quadrat, Treatment, Year) %>% 
  summarize(Cover=sum(Cover_roll)) %>% 
  group_by(Treatment, Year) %>% 
  summarize(Cover=mean(Cover)) %>% 
  mutate(DeltaCov=Cover-lag(Cover)) %>%
  #group_by(Treatment) %>% 
  #arrange(desc(abs(DeltaSR))) %>% 
  #slice_head(n=5) %>% 
  View()


### NMDS ####
# Convert species x plot level to wide format
flora_wide <- flora_pooled2 %>% 
  ## filter(Layer=="U") %>% 
  ## flatten vegetation layers
  #group_by(Quadrat, Species_resolved, Year) %>% 
  #summarize(Cover=combine.cover(Cover)) %>% 
  dplyr::select(Quadrat, Year, Species_resolved, Cover_roll) %>% 
  pivot_wider(names_from=Species_resolved, values_from = "Cover_roll", values_fill = 0) %>% 
  unite("Quadrat_year", Quadrat, Year) %>% 
  arrange(Quadrat_year) %>% 
  column_to_rownames("Quadrat_year")
# Compute dissimilarity matrix
flora_dist <- vegan::vegdist((flora_wide), method="bray")

# Compute NMDS
set.seed(1000)
flora_mds <- vegan::metaMDS(flora_dist, k=2, try = 100, trymax=100)

# Passively project species  

flora_envfit0 <- vegan::envfit(flora_mds, flora_wide)
flora_envfit <- data.frame(flora_envfit0$vectors$arrows) %>%  
   bind_cols(data.frame(r2=flora_envfit0$vectors$r, p=flora_envfit0$vectors$pvals)) %>% 
   rownames_to_column("Species") %>% 
   filter(p<0.05) %>% 
   filter(r2>0.33) %>% 
   arrange(desc(r2))


signspecies <- flora_envfit$Species


##### Figure 5 ####
# Convert data to data.frame and add ancillary values
flora_mds_scores <- flora_mds$points %>% 
  as.data.frame() %>% 
  rownames_to_column("Quadrat_year") %>% 
  as_tibble() %>% 
  separate(Quadrat_year, into=c("Quadrat", "Year"), sep="_") %>% 
  left_join(flora %>% 
              distinct(Quadrat, Treatment, Series), 
            by="Quadrat") %>% 
  arrange(Quadrat, Year) %>% 
  rename(NMS1=MDS1, NMS2=MDS2)



Figure5 <- ggplot(data=flora_mds_scores, aes(x=NMS1, y=NMS2)) +
  #  geom_point(aes(x=MDS1, y=MDS2, col=Quadrat)) + 
  geom_path(aes(group=Quadrat, col=Treatment, linetype=Series), arrow=arrow(angle = 15, length = unit(0.07, "inches"))) + 
  geom_point(data=flora_mds_scores %>% 
               filter(Year==2012), 
                      aes(group=Quadrat, col=Treatment), show_guide = FALSE) +
  scale_color_brewer(palette = "Dark2", direction = -1) + 
  geom_text(label=flora_mds_scores$Year, size=3, hjust=-0.2, show.legend =F,
            col= "black",
            check_overlap = F, 
            alpha=rep(c(1,0,0,1,0,0,1,0,0,0,1)/3, length.out=nrow(flora_mds_scores)), 
            fontface="bold") + 
#  geom_point(data=flora_envfit, 
#             aes(x=NMDS1, y=NMDS2), pch=17) + 
#  ggrepel::geom_text_repel(data=flora_envfit, 
#               aes(x=NMDS1, y=NMDS2, label=Species), max.overlaps = 12) + 
  theme_bw() + 
  theme(panel.grid = element_blank())


## Alternative to Figure 5 to show loads of envfit species - not shown in papaer
ggplot(data=flora_mds_scores, aes(x=NMS1, y=NMS2)) +
  #  geom_point(aes(x=MDS1, y=MDS2, col=Quadrat)) + 
  geom_path(aes(group=Quadrat, col=Treatment, linetype=Series), arrow=arrow(angle = 15, length = unit(0.07, "inches")), alpha=0.2) + 
  geom_point(data=flora_mds_scores %>% 
               filter(Year==2012), 
             aes(group=Quadrat, col=Treatment), show_guide = FALSE, alpha=0.2) +
  scale_color_brewer(palette = "Dark2", direction = -1) + 
  geom_text(label=flora_mds_scores$Year, size=3, hjust=-0.2, show.legend =F,
            col= "black",
            check_overlap = F, 
            alpha=rep(c(1,0,0,1,0,0,1,0,0,0,1)/3, length.out=nrow(flora_mds_scores)), 
            fontface="bold") + 
  geom_point(data=flora_envfit, 
               aes(x=NMDS1, y=NMDS2), pch=17) + 
    ggrepel::geom_text_repel(data=flora_envfit, 
                 aes(x=NMDS1, y=NMDS2, label=Species), max.overlaps = 12) + 
  theme_bw() + 
  theme(panel.grid = element_blank())



ggsave(Figure5, file="../figure_tables/Figure5_NMDS.png", 
       width = 5, height=4, units="in",dpi = 300)





### GLMMs SR, Cover####

### Ancillary function to make graphs of predictions with C.I.
myplot_lme <- function(mydata, response_var, pred_data, mymod){
  response_var <- sym(response_var)
  mydata <- mydata %>% 
    mutate(response_var=!!response_var)
  
  ggplot(data=pred_data) +
    geom_point(data=mydata, 
               aes(x=jitter(LogDeltaYear), y=response_var, 
                   colour=Treatment), alpha=0.7) +
    geom_line(data=mydata, 
              aes(x=LogDeltaYear, y=predict(mymod), 
                  group=Quadrat, color=Treatment), alpha=0.3) +
    geom_line(aes(y=pred, x=LogDeltaYear, color=Treatment)) +
    #geom_ribbon(data=pred_data, aes(x=DeltaYear,ymin=pred-2*SE2,ymax=pred+2*SE2, fill=Treatment),alpha=0.2) +
    geom_ribbon(aes(x=LogDeltaYear,ymin=pred-2*SE,ymax=pred+2*SE, fill=Treatment),alpha=0.2) +
    scale_x_continuous(name="Year", 
                       trans="exp", 
                       breaks = log(seq(0,12, by=2)+1), 
                       labels = seq(2012, 2024, by=2 )) + 
    scale_y_continuous(name=response_var) + 
    scale_alpha_continuous(name=NULL) +
    theme_bw() + 
    scale_color_brewer(palette = "Dark2", direction = -1) + 
    scale_fill_brewer(palette = "Dark2", direction = -1)
}


#### Data Preparation ####
floraSR <- flora_pooled2 %>% 
  #filter(Layer=="U") %>% 
  group_by(Quadrat, Treatment, Series, Year) %>% 
  summarize(SR=n(), Tot_cover=sum(Cover_roll)) %>% 
  left_join({.} %>% 
              ungroup() %>% 
              filter(Year==2012) %>% 
              select(Quadrat, SR0=SR, Tot_cover0=Tot_cover), 
            by="Quadrat") %>% 
  mutate(DeltaSR=SR-SR0, 
         DeltaCov=Tot_cover-Tot_cover0) %>% 
  left_join(header0 %>% 
              select(Quadrat, Year, Openness), 
            by=c("Quadrat", "Year")) %>% 
  group_by(Quadrat, Treatment, Series) %>% 
  #mutate(OpennessT_1=lag(Openness)) %>% 
  mutate(DeltaYear=Year-2012, LogDeltaYear=log(DeltaYear+1)) %>% 
  ungroup() %>% 
  filter(complete.cases(.)) # %>% 
  ## Deprecated code to calculate species gains and losses
#  left_join(flora %>% 
#              filter(Layer=="U") %>% 
#              group_by(Quadrat, Treatment, Series, Year) %>% 
#              arrange(Species_resolved) %>% 
#              summarise(species_list = list(Species_resolved)) %>% 
#              mutate(species_list_lag = lag(species_list), 
#                     Gains = map2_int(species_list, 
#                                      species_list_lag, 
#                                      ~length(setdiff(.x, .y))), 
#                     Losses = map2_int(species_list, 
#                                       species_list_lag, 
#                                       ~length(setdiff(.y, .x)))) %>% 
#              filter(Year != 2012) %>% 
#              mutate(DeltaYear=Year-2012, LogDeltaYear=log(DeltaYear+1)) %>% 
#              ungroup() %>% 
#              dplyr::select(Quadrat, Year, Gains, Losses), 
#            by=c("Quadrat", "Year"))#


#### GLMMs change of SR with time - Model selection ####
mod1 <- lme(DeltaSR ~ Treatment * LogDeltaYear, data=floraSR, random=~1|Quadrat)
mod2 <- lme(DeltaSR ~ Treatment * LogDeltaYear +0, 
            data=floraSR, 
            random=(~1|Quadrat), 
            correlation = corAR1(form = ~1 | Quadrat))
mod3 <- lme(DeltaSR ~ Treatment * DeltaYear, data=floraSR, 
            random=(~1+Year|Quadrat), 
            correlation = corAR1(form = ~1 | Quadrat))
mod4 <- lme(DeltaSR ~ Treatment * I(DeltaYear)^2, data=floraSR, 
            random=(~1|Quadrat), 
            correlation = corAR1(form = ~1 | Quadrat))
mod5 <- lme(DeltaSR ~ Treatment * DeltaYear + Openness, data=floraSR, 
            random=(~1 |Quadrat), 
            correlation = corAR1(form = ~1 | Quadrat))

MuMIn::AICc(mod1, mod2, mod3, mod4, mod5)
## Best model is mod2: Treatment*DeltaYear with corAR
## Quadratic term not needed
## Openness not needed

summary(mod1)
summary(mod2)

### Check model Residuals 
# standardized residuals versus fitted values by gender
plot(mod2, resid(., type = "p") ~ fitted(.) | Treatment, abline = 0)
plot(mod2, resid(., type = "p") ~ DeltaYear | Treatment, abline = 0)
# box-plots of residuals by Subject
plot(mod2, Quadrat ~ resid(.))
# observed versus fitted values by Subject
plot(mod2, DeltaSR ~ fitted(.) | Quadrat, abline = c(0,1))


### predict response
newdat <- expand.grid(Treatment=unique(floraSR$Treatment),
                      LogDeltaYear=seq(from=min(floraSR$LogDeltaYear),
                                       to=max(floraSR$LogDeltaYear), 
                                       length.out=101))
newdat$pred <- predict(mod2, newdata=newdat,level=0)

## code to predict response at round intervals
newdat <- expand.grid(Treatment=unique(floraSR$Treatment),
                     LogDeltaYear=log(0:13+1))
newdat$pred <- predict(mod2, newdata=newdat,level=0)
newdat %>% 
 mutate(DeltaYear=exp(LogDeltaYear)-1) %>% 
 arrange(Treatment, LogDeltaYear) %>% 
 group_by(Treatment) %>% 
 mutate(diff = pred - lag(pred, n=1L)) %>% 
 View()



#create design matrix
Designmat <- model.matrix(eval(eval(mod2$call$fixed)[-2]), newdat[-ncol(newdat)])

#compute standard error for predictions
predvar <- diag(Designmat %*% mod2$varFix %*% t(Designmat))
newdat$SE <- sqrt(predvar) 
newdat$SE2 <- sqrt(predvar+mod2$sigma^2)

#### Plot predicted responses of DeltaSR 
left_SR <- myplot_lme(mydata = floraSR, "DeltaSR", newdat, mod2)



#### GLMMs change of cover with time - Model selection ####
mod_cov1 <- lme(DeltaCov ~ Treatment * LogDeltaYear, data=floraSR, random=~1|Quadrat)
mod_cov2 <- lme(DeltaCov ~ Treatment * LogDeltaYear +0, data=floraSR, 
                random=(~1|Quadrat), 
                correlation = corAR1(form = ~1 | Quadrat))
mod_cov3 <- lme(DeltaCov ~ Treatment * LogDeltaYear + Openness, data=floraSR, 
                random=(~1|Quadrat), 
                correlation = corAR1(form = ~1 | Quadrat))
MuMIn::AICc(mod_cov1, 
            mod_cov2, 
            mod_cov3)

summary(mod_cov2)

## Check model residuals
# standardized residuals versus fitted values by gender
plot(mod_cov2, resid(., type = "p") ~ fitted(.) | Treatment, abline = 0)
plot(mod_cov2, resid(., type = "p") ~ DeltaYear | Treatment, abline = 0)
# box-plots of residuals by Subject
plot(mod_cov2, Quadrat ~ resid(.))
# observed versus fitted values by Subject
plot(mod_cov2, DeltaCov ~ fitted(.) | Quadrat, abline = c(0,1))


## Predict new data to plot response
newdat2 <- newdat
newdat2$pred <- predict(mod_cov2, newdata=newdat2,level=0)

#create design matrix
Designmat <- model.matrix(eval(eval(mod_cov2$call$fixed)[-2]), newdat[-ncol(newdat2)])

#compute standard error for predictions
predvar <- diag(Designmat %*% mod_cov2$varFix %*% t(Designmat))
newdat2$SE <- sqrt(predvar) 
newdat2$SE2 <- sqrt(predvar+mod_cov2$sigma^2)


## code to predict response at round intervals

newdat <- expand.grid(Treatment=unique(floraSR$Treatment),
                      LogDeltaYear=log(0:13+1))
newdat$pred <- predict(mod_cov2, newdata=newdat,level=0)
newdat %>% 
  mutate(DeltaYear=exp(LogDeltaYear)-1) %>% 
  arrange(Treatment, LogDeltaYear) %>% 
  group_by(Treatment) %>% 
  mutate(diff = pred - lag(pred, n=1L)) %>% 
  View()


#### Plot predicted responses of DeltaCov 
right_Cov <- myplot_lme(mydata = floraSR, "DeltaCov", newdat2, mod_cov2)



##### Figure 6 ####
(left_SR + 
    theme(legend.position = "none") + 
    scale_y_continuous(name="Delta Species Richness")) | 
  (right_Cov + scale_y_continuous(name="Delta Cover")) + 
  plot_annotation(tag_levels = "a")

ggsave(filename="../figure_tables/Figure6_SR_lme.png", width=8, height=4, units = "in", dpi=300)




### same but based on matrix T
floraSR_T <- flora %>% 
  filter(Layer=="U") %>% 
  group_by(Quadrat, Treatment, Series, Year) %>% 
  summarize(SR=n(), Tot_cover=sum(Cover)) %>% 
  left_join({.} %>% 
              ungroup() %>% 
              filter(Year==2012) %>% 
              select(Quadrat, SR0=SR, Tot_cover0=Tot_cover), 
            by="Quadrat") %>% 
  mutate(DeltaSR=SR-SR0, 
         DeltaCov=Tot_cover-Tot_cover0) %>% 
  left_join(header0 %>% 
              select(Quadrat, Year, Openness), 
            by=c("Quadrat", "Year")) %>% 
  group_by(Quadrat, Treatment, Series) %>% 
  #mutate(OpennessT_1=lag(Openness)) %>% 
  mutate(DeltaYear=Year-2012, LogDeltaYear=log(DeltaYear+1)) %>% 
  ungroup() %>% 
  filter(complete.cases(.))


modT <- nlme::lme(DeltaSR ~ Treatment * LogDeltaYear + 0, 
            data=floraSR_T, 
            random=(~1|Quadrat), 
            correlation = corAR1(form = ~1 | Quadrat))

mod_covT <- lme(DeltaCov ~ Treatment * LogDeltaYear +0, 
                data=floraSR_T, 
                random=(~1|Quadrat), 
                correlation = corAR1(form = ~1 | Quadrat))


# library(flextable)
# library(modelsummary)
# modelsummary(list("SR (R Matrix)" = mod2, 
#                   "SR (T Matrix)" = modT, 
#                   "Cover (R Matrix)" = mod_cov2, 
#                   "Cover (T Matrix)" = mod_covT), 
#              fmt = 2,  
#              estimate  = c("{estimate}{stars}"))
# 
# 




## #### GLMMS change of Gains
## floraSRc <- floraSR %>% 
##   filter(Year != 2012)
## 
## mod1_gain <- lme(Gains ~ Treatment * LogDeltaYear, data=floraSRc, random=~1|Quadrat)
## plot(mod1_gain, resid(., type = "p") ~ fitted(.) | Treatment, abline = 0)
## plot(mod1_gain, resid(., type = "p") ~ DeltaYear | Treatment, abline = 0)
## # box-plots of residuals by Subject
## plot(mod1_gain, Quadrat ~ resid(.))
## # observed versus fitted values by Subject
## plot(mod1_gain, Gains ~ fitted(.) | Quadrat, abline = c(0,1))
## summary(mod1_gain)
## 
## 
## 
## ### predict response
## newdat <- expand.grid(Treatment=unique(floraSRc$Treatment),
##                       LogDeltaYear=seq(from=min(floraSRc$LogDeltaYear),
##                                        to=max(floraSRc$LogDeltaYear), 
##                                        length.out=101))
## newdat$pred <- predict(mod1_gain, newdata=newdat,level=0)
## 
## 
## #create design matrix
## Designmat <- model.matrix(eval(eval(mod1_gain$call$fixed)[-2]), newdat[-ncol(newdat)])
## 
## #compute standard error for predictions
## predvar <- diag(Designmat %*% mod1_gain$varFix %*% t(Designmat))
## newdat$SE <- sqrt(predvar) 
## newdat$SE2 <- sqrt(predvar+mod1_gain$sigma^2)
## 
## #### Plot predicted responses of Gains
## (ggmodgains <- myplot_lme(floraSRc, "Gains", newdat, mod1_gain))
## 
## 
## 
## 
## #### GLMMS change of Losses
## mod1_loss <- lme(Losses ~ Treatment * LogDeltaYear, data=floraSRc, random=~1|Quadrat)
## plot(mod1_loss, resid(., type = "p") ~ fitted(.) | Treatment, abline = 0)
## plot(mod1_loss, resid(., type = "p") ~ DeltaYear | Treatment, abline = 0)
## # box-plots of residuals by Subject
## plot(mod1_loss, Quadrat ~ resid(.))
## # observed versus fitted values by Subject
## plot(mod1_loss, Losses ~ fitted(.) | Quadrat, abline = c(0,1))
## summary(mod1_loss)
## 
## 
## 
## ### predict response
## newdat <- expand.grid(Treatment=unique(floraSRc$Treatment),
##                       LogDeltaYear=seq(from=min(floraSRc$LogDeltaYear),
##                                        to=max(floraSRc$LogDeltaYear), 
##                                        length.out=101))
## newdat$pred <- predict(mod1_loss, newdata=newdat,level=0)
## 
## 
## #create design matrix
## Designmat <- model.matrix(eval(eval(mod1_loss$call$fixed)[-2]), newdat[-ncol(newdat)])
## 
## #compute standard error for predictions
## predvar <- diag(Designmat %*% mod1_loss$varFix %*% t(Designmat))
## newdat$SE <- sqrt(predvar) 
## newdat$SE2 <- sqrt(predvar+mod1_loss$sigma^2)
## 
## #### Plot predicted responses of gains
## (ggmodlosses <- myplot_lme(floraSRc, "Losses", newdat, mod1_loss))
## 
## 
## 
## 
## ##### Figure xxx gain/losses
## gggains <- ggplot(data=floraSRc) + 
##   geom_boxplot(aes(y=Gains, x=Treatment, col=Treatment, group=Treatment), alpha=0.7) + 
##   scale_color_brewer(palette = "Dark2", direction = -1) + 
##   scale_fill_brewer(palette = "Dark2", direction = -1) +
##   theme_bw() + 
##   theme(legend.position = "none")
## 
## gglosses <- ggplot(data=floraSRc) + 
##   geom_boxplot(aes(y=Losses, x=Treatment, col=Treatment, group=Treatment), alpha=0.7) +    
##   theme_bw() + 
##   scale_color_brewer(palette = "Dark2", direction = -1) + 
##   scale_fill_brewer(palette = "Dark2", direction = -1) + 
##   theme(legend.position = "none")
## 
## (gggains + ggmodgains) / (gglosses + (ggmodlosses  + theme(legend.position = "none")))
## 
## ggsave(filename="../figure_tables/FigureXXX_GainLosses_lme.png", width=8, height=6, units = "in", dpi=300)
## 
## 







