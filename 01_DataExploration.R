library(tidyverse)
library(betapart)
library(patchwork)
library(wesanderson)

### Ancillary functions ####

## Calculate rolling mean of change
ch <- function(x){x[3]-x[1]}

# Combine cover accounting for layers
combine.cover <- function(x){
  while (length(x)>1){
    x[2] <- x[1]+(100-x[1])*x[2]/100
    x <- x[-1]
  }
  return(x)
}



### General Description ####

### Graph of species abundances with time (by plot)
ggplot(data=flora %>% 
         filter(Layer=="U") %>% 
         filter(Quadrat=="GN") %>% 
         mutate(Species_resolved=factor(Species_resolved)) %>% 
         left_join({.} %>% 
                     count(Species_resolved), 
                   by=c("Species_resolved")) %>% 
         group_by(Species_resolved, Year) %>% 
         filter(n>2) %>% 
         summarize(Summed_cover=sum(Cover))) + 
  geom_line(aes(x=Year, y=Summed_cover)) + #, group=Species_resolved, col=Species_resolved)) +
  #scale_color_brewer(type="div")
  ggplot2::facet_wrap(~Species_resolved, ncol=5) + 
  scale_x_continuous(breaks = seq(2012, 2022, by=2)) + 
  theme_bw()
  

### Graph of species richness with time by plot
(topleft <- ggplot(data=flora %>% 
         group_by(Quadrat, Treatment, Series, Year) %>% 
           #distinct(Species_resolved, .keep_all=T) %>% 
         summarize(SR=n()) ) + 
  geom_line(aes(x=Year, y=SR, group=Quadrat, col=Treatment, linetype=Series)) + 
  scale_x_continuous(breaks = seq(2012, 2022, by=2)) + 
  scale_y_continuous(name="Species Richness") + 
  scale_color_brewer(palette = "Dark2", direction = -1) + 
  theme_bw()
)

### Graph of herb-layer total cover with time by plot
(topright <- ggplot(data=flora %>% 
                    filter(Layer!="O") %>% 
                    group_by(Quadrat, Treatment, Series, Year) %>% 
                    summarize(Tot_cover=sum(Cover)) %>%
                    arrange(Quadrat, Year)) + 
  geom_line(aes(x=Year, y=Tot_cover, group=Quadrat, col=Treatment, linetype=Series)) + 
  scale_x_continuous(breaks = seq(2012, 2022, by=2)) + 
  scale_y_continuous(name="Total Herb Layer Cover (%)") +
  scale_color_brewer(palette = "Dark2", direction = -1) + 
  theme_bw()
)

mylegend <- ggpubr::get_legend(topright + 
                     theme(legend.position = "bottom"))

# ## Graph of canopy openness per plot
# ggplot(data=header0) +
#   #geom_line(aes(x=Year, y=Openness,group=Quadrat, col=Treatment, linetype=Series )) + 
#   geom_line(aes(x=Year, y=Openness, group=Quadrat, col=Treatment, linetype=Series )) + 
#   scale_x_continuous(breaks = seq(2012, 2022, by=2)) + 
#   theme_bw()








### RQ1: Rate of Change ####

### Helical Graphs
# See: https://www.sciencedirect.com/science/article/pii/S157495412200406X?casa_token=_7VEA2W0D6kAAAAA:SEYaqqXXmvrEoWgMfLFP18S4vSB_iaDXCnxPbz6Gp1yvKnnXbTdUMO_D_o2_S_UeSsGi0S2G

flora_sr <- flora %>% 
  filter(Layer!="O") %>% 
  group_by(Quadrat, Treatment, Series, Year) %>% 
  summarize(SR=n(), Tot_cover=sum(Cover)) %>%
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
    scale_x_continuous(name="Yearly Change in Richness") +
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
    scale_x_continuous(name="Yearly Change Cover (%)") +
  theme_bw()
)

Figure2 <- ((topleft + theme(legend.position="none")  |
               topright+ theme(legend.position="none")) / 
              (bottomleft + theme(legend.position="none") |
                 bottomright + theme(legend.position="none"))) / 
  mylegend + plot_layout(height=c(6,6,1)) + plot_annotation(tag_levels = list(c("a", "b", 'c', "d", "")))

ggsave(Figure2, file="../figure_tables/Figure2_SR_Cover_Change.png", 
       width = 9, height=8, units="in",dpi = 300)


## Graph of smoothed rolling mean of change in openness vs openness

header_roll <- header0 %>% 
  group_by(Quadrat, Treatment, Series) %>% 
  mutate(Openness_roll=rollapply(Openness, 3, mean, align="right", fill=NA)) %>% 
  mutate(Opennessch_roll=rollapply(Openness, 3, ch, align="right", fill=NA)) %>% 
  filter(!is.na(Openness_roll)) %>% 
  mutate(row_id=row_number()) %>% 
  mutate(Openness_roll_lo=predict(loess(Openness_roll~row_id))) %>%
  mutate(Opennessch_roll_lo=predict(loess(Opennessch_roll~row_id)))


ggplot(data=header_roll, 
       aes(x=Opennessch_roll_lo, y=Openness_roll_lo, 
           group=Quadrat, col=Treatment, linetype=Series)) + 
  geom_path(arrow = arrow(type="closed", length=unit(0.1, "inches"))) + 
  geom_text(label=flora_sr$Year, size=4, hjust=-0.2, show.legend=F,
            col="black", check_overlap = T, alpha=1, fontface="bold") + 
  theme_bw()






# ## Graph of smoothed rolling mean of beta_sorensen SR vs SR
# ggplot(data=flora_sr, 
#       aes(x=sor_roll_lo, y=SR_roll_lo, 
#           group=Quadrat, col=Treatment, linetype=Series)) + 
#  geom_path(arrow = arrow(type="closed", length=unit(0.1, "inches"))) + 
#  geom_text(label=flora_sr$Year, size=4, hjust=-0.2, show.legend=F,
#            col="black", check_overlap = T, alpha=1, fontface="bold") + 
#  theme_bw()



### Proportional  effective species turnover index ####
# Beta_year = 1 - S_year/S_accumulated

flora_eff <- flora %>% 
  filter(Layer!="O") %>% 
  group_by(Quadrat, Treatment, Series, Year) %>% 
  summarize(SR=n()) %>%
  arrange(Quadrat, Year) %>% 
  left_join(flora %>% 
              filter(Layer!="O") %>% 
              group_by(Quadrat, Treatment, Series) %>% 
              distinct(Species_resolved, .keep_all=T) %>% 
              summarize(SR_acc=n()), 
            by=c("Quadrat", "Treatment", "Series")) %>% 
  mutate(Beta_year=(1-SR/SR_acc))

ggplot(data=flora_eff, 
       aes(x=Year, y=Beta_year, col=Treatment)) + 
  geom_point()  + 
  #geom_smooth(method=lm, formula = y ~ splines::bs(x, 3), se = FALSE) +
  geom_smooth(aes(fill=Treatment), alpha=1/7, method=lm, formula = y ~ x, se = T) +
  scale_color_brewer(palette = "Dark2", direction = -1) + 
  scale_fill_brewer(palette = "Dark2", direction = -1) + 
  scale_x_continuous(breaks = seq(2012, 2022, by=2)) + 
  theme_bw()


### Graph of beta eff with time by plot
ggplot(data=flora_eff) + 
  geom_line(aes(x=Year, y=Beta_year, group=Quadrat, col=Treatment, linetype=Series)) + 
  scale_x_continuous(breaks = seq(2012, 2022, by=2)) + 
  scale_color_brewer(palette = "Dark2", direction = -1) + 
  theme_bw()



### Specaccum by site ####
## After how many years am I confident I sampled most of the species?

plot_index <- sapply(rownames(flora_wide), function(x){strsplit(x, split = "_")[[1]][[1]]})

specaccum_df0 <- NULL
for(i in unique(plot_index)){
  tmp_flora <- flora_wide[which(plot_index==i),]
  tmp_flora <- tmp_flora[, which(colSums(tmp_flora)>0)] #delete empty species
  rownames(tmp_flora) <- data.frame(QuadratYear=rownames(tmp_flora)) %>% 
    separate(QuadratYear, into=c("Quadrat", "Year"), sep="_") %>% 
    pull(Year)
  specaccum_df0 <- bind_rows(specaccum_df0, 
                              specaccum(tmp_flora, method="collector")[[4]])
}
specaccum_df <- specaccum_df0 %>% 
  mutate(Quadrat=unique(plot_index)) %>% 
  left_join(flora %>% 
              distinct(Quadrat, Treatment, Series), 
            by="Quadrat") %>% 
  pivot_longer(cols=-c(Quadrat, Treatment, Series), 
               values_to="Cumulative_richness", 
               names_to="Year") %>% 
  mutate(Year=as.numeric(Year)) %>% 
  left_join({.} %>% 
              dplyr::select(Quadrat, Cumulative_richness) %>% 
              group_by(Quadrat) %>% 
              summarize(Total_SR=max(Cumulative_richness)), 
            by="Quadrat") %>% 
  mutate(Sampled_richness=Cumulative_richness/Total_SR*100)

ggplot(data=specaccum_df) + 
  geom_line(aes(x=Year, y=Sampled_richness, 
                group=Quadrat, col=Treatment, linetype=Series)) + 
  scale_x_continuous(breaks = seq(2012, 2022, by=2)) + 
  scale_color_brewer(palette = "Dark2", direction = -1) + 
  theme_bw() + 
  facet_wrap(.~Treatment, scales = "free_y")







### Modelling SR change from t0 ####
floraSR <- flora %>% 
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




ggplot(data=floraSR) +
  geom_line(aes(x=DeltaYear, y=DeltaSR, group=Quadrat, col=Treatment, linetype=Series)) + 
  scale_x_continuous(breaks = seq(0,12, by=2)) + 
  theme_bw()

ggplot(data=floraSR) +
  geom_point(aes(x=Openness, y=SR, group=Quadrat, col=Treatment, linetype=Series)) + 
  theme_bw()

write_csv(floraSR, file = "../intermediate_steps/floraSR.csv")



library(nlme)
#library(effects)

mod1 <- lme(DeltaSR ~ Treatment * LogDeltaYear, data=floraSR, random=~1|Quadrat)
summary(mod1)
mod2 <- lme(DeltaSR ~ Treatment * LogDeltaYear +0, data=floraSR, 
            random=(~1|Quadrat), 
            correlation = corAR1(form = ~1 | Quadrat))
summary(mod2)
MuMIn::AICc(mod1, mod2)

# set.seed(1999)
# floraSRj <- floraSR %>% 
#   mutate(DeltaSR=jitter(DeltaSR), DeltaYear=jitter(DeltaYear))
# mod2a <- lme(DeltaSR ~ Treatment : DeltaYear + Year +0, data=floraSRj, 
#             random=(~1|Quadrat), 
#             correlation = corAR1(form = ~1 | Quadrat))
# mod2b <- lme(DeltaSR ~ Treatment : DeltaYear +0, data=floraSRj, 
#              random=(~1|Quadrat), 
#              correlation = corAR1(form = ~1 | Quadrat))
# mod2c <- lme(DeltaSR ~ Treatment : DeltaYear +0, data=floraSRj, 
#              random=(~1|Quadrat))
# 
# MuMIn::r.squaredGLMM(mod2)
# MuMIn::r.squaredGLMM(mod2a)
# MuMIn::r.squaredGLMM(mod2b)
# MuMIn::r.squaredGLMM(mod2c)


summary(mod2)
mod3 <- lme(DeltaSR ~ Treatment * DeltaYear, data=floraSR, 
            random=(~1+Year|Quadrat), 
            correlation = corAR1(form = ~1 | Quadrat))


MuMIn::AICc(mod1, mod2, mod3)
## Phi is extremely low!

mod4 <- lme(DeltaSR ~ Treatment * I(DeltaYear)^2, data=floraSR, 
            random=(~1|Quadrat), 
            correlation = corAR1(form = ~1 | Quadrat))
summary(mod4)
## Quadratic term not needed

mod5 <- lme(DeltaSR ~ Treatment * DeltaYear + Openness, data=floraSR, 
            random=(~1 |Quadrat), 
            correlation = corAR1(form = ~1 | Quadrat))
summary(mod5)
## Openness not needed
MuMIn::AICc(mod1, mod2, mod3, mod4, mod5)
## Best model is mod1: Treatment*DeltaYear -- But I feel corAR is needed being longitudinal data


# standardized residuals versus fitted values by gender
plot(mod2, resid(., type = "p") ~ fitted(.) | Treatment, abline = 0)
plot(mod2, resid(., type = "p") ~ DeltaYear | Treatment, abline = 0)
# box-plots of residuals by Subject
plot(mod2, Quadrat ~ resid(.))
# observed versus fitted values by Subject
plot(mod2, DeltaSR ~ fitted(.) | Quadrat, abline = c(0,1))


newdat <- expand.grid(Treatment=unique(floraSR$Treatment),
                      LogDeltaYear=seq(from=min(floraSR$LogDeltaYear),
                                    to=max(floraSR$LogDeltaYear), 
                                    length.out=101))

#predict response

newdat$pred <- predict(mod2, newdata=newdat,level=0)
#create design matrix
Designmat <- model.matrix(eval(eval(mod2$call$fixed)[-2]), newdat[-ncol(newdat)])

#compute standard error for predictions
predvar <- diag(Designmat %*% mod2$varFix %*% t(Designmat))
newdat$SE <- sqrt(predvar) 
newdat$SE2 <- sqrt(predvar+mod2$sigma^2)



#left_SR <- 
myplot_lme <- function(response_var, pred_data, mymod){
  response_var <- sym(response_var)
  floraSR <- floraSR %>% 
    mutate(response_var=!!response_var)
  print(floraSR)
  
  ggplot(data=pred_data) +
    geom_point(data=floraSR, aes(x=LogDeltaYear, y=response_var, colour=Treatment)) +
    geom_line(data=floraSR, aes(x=LogDeltaYear, y=predict(mymod), group=Quadrat, color=Treatment), alpha=0.3) +
    geom_line(aes(y=pred, x=LogDeltaYear, color=Treatment)) +
    #geom_ribbon(data=pred_data, aes(x=DeltaYear,ymin=pred-2*SE2,ymax=pred+2*SE2, fill=Treatment),alpha=0.2) +
    geom_ribbon(aes(x=LogDeltaYear,ymin=pred-2*SE,ymax=pred+2*SE, fill=Treatment),alpha=0.2) +
    scale_x_continuous(name="Year", 
                       trans="exp", 
                       breaks = log(seq(0,12, by=2)+1), 
                       labels = seq(2012, 2024, by=2 )) + 
    scale_y_continuous(name=response_var) + 
    theme_bw() + 
    scale_color_brewer(palette = "Dark2", direction = -1) + 
    scale_fill_brewer(palette = "Dark2", direction = -1)
}

left_SR <- myplot_lme("DeltaSR", newdat, mod2)


## Model total cover

mod_cov0 <- lme(DeltaCov ~ Treatment * DeltaYear, data=floraSR, random=~1|Quadrat)
mod_cov1 <- lme(DeltaCov ~ Treatment * LogDeltaYear, data=floraSR, random=~1|Quadrat)
mod_cov2 <- lme(DeltaCov ~ Treatment * LogDeltaYear +0, data=floraSR, 
            random=(~1|Quadrat), 
            correlation = corAR1(form = ~1 | Quadrat))
summary(mod_cov2)
MuMIn::AICc(mod_cov0, mod_cov1, mod_cov2)

# standardized residuals versus fitted values by gender
plot(mod_cov2, resid(., type = "p") ~ fitted(.) | Treatment, abline = 0)
plot(mod_cov2, resid(., type = "p") ~ DeltaYear | Treatment, abline = 0)
# box-plots of residuals by Subject
plot(mod_cov2, Quadrat ~ resid(.))
# observed versus fitted values by Subject
plot(mod_cov2, DeltaCov ~ fitted(.) | Quadrat, abline = c(0,1))

newdat2 <- newdat
newdat2$pred <- predict(mod_cov2, newdata=newdat2,level=0)
#create design matrix
Designmat <- model.matrix(eval(eval(mod_cov2$call$fixed)[-2]), newdat[-ncol(newdat2)])

#compute standard error for predictions
predvar <- diag(Designmat %*% mod_cov2$varFix %*% t(Designmat))
newdat2$SE <- sqrt(predvar) 
newdat2$SE2 <- sqrt(predvar+mod_cov2$sigma^2)

right_Cov <- myplot_lme("DeltaCov", newdat2, mod_cov2)


(left_SR + 
  theme(legend.position = "none") + 
  scale_y_continuous(name="Delta Species Richness")) | 
(right_Cov + scale_y_continuous(name="Delta Cover")) + 
plot_annotation(tag_levels = "a")

ggsave(filename="../figure_tables/Figure3_SR_lme.png", width=8, height=4, units = "in", dpi=300)









### NMDS   ####

# Convert species x plot level to wide format
flora_wide <- flora %>% 
  filter(Layer=="U") %>% 
  ## flatten vegetation layers
  group_by(Quadrat, Species_resolved, Year) %>% 
  summarize(Cover=combine.cover(Cover)) %>% 
  pivot_wider(names_from=Species_resolved, values_from = "Cover", values_fill = 0) %>% 
  unite("Quadrat_year", Quadrat, Year) %>% 
  arrange(Quadrat_year) %>% 
  column_to_rownames("Quadrat_year")

# Compute dissimilarity matrix
flora_dist <- vegan::vegdist((flora_wide), method="bray")

# Compute NMDS
flora_mds <- vegan::metaMDS(flora_dist, k=2, )

# Graph of NMDS
# Convert data to data.frame and add ancillary values
flora_mds_scores <- flora_mds$points %>% 
  as.data.frame() %>% 
  rownames_to_column("Quadrat_year") %>% 
  as_tibble() %>% 
  separate(Quadrat_year, into=c("Quadrat", "Year"), sep="_") %>% 
  left_join(flora %>% 
              distinct(Quadrat, Treatment, Series), 
            by="Quadrat") %>% 
  arrange(Quadrat, Year)

ggplot(data=flora_mds_scores) +
#  geom_point(aes(x=MDS1, y=MDS2, col=Quadrat)) + 
  geom_path(aes(x=MDS1, y=MDS2, group=Quadrat, col=Treatment, linetype=Series), arrow=arrow(angle = 15, length = unit(0.15, "inches"))) + 
  theme_bw()



#### Test of Shimadzu et al. 2015 framework and functions

### Look at temporal turnover between t and t0, as decomposed 
### in its composition (D1) and overall abundance (D2) components

DD <- function(x, ref.t=1, zero.rm=FALSE){
  lmb <- apply(x, 1, sum)
  D2 <- log(lmb/lmb[ref.t])
  x.p <- t(apply(x, 1, function(z)z/sum(z)))
  Pt <- x.p[ref.t,]
  if(zero.rm==FALSE){
    D1 <- -t(apply(x.p, 1, function(z)ifelse(Pt==0, 0, log(Pt/z)))) %*% Pt
  }else{
    D1 <- -t(apply(x.p, 1, function(z)ifelse(Pt==0|z==0, 0, log(Pt/z)))) %*% Pt
  }
  D <- D1 + D2
  data.frame(D, D1, D2)
}

quadrat_list <- rownames(flora_wide) %>% 
  as_tibble() %>% 
  separate(value, into=c("Quadrat", "Year"), sep = "_")
dd <- NULL

for(i in unique(quadrat_list$Quadrat)){
  aa <- flora_wide[which(quadrat_list$Quadrat==i),]
  aa <- aa[,which(colSums(aa)>0)]
  
  bb <- DD(aa, zero.rm=T) %>% 
    rownames_to_column("Quadrat_year") %>% 
    separate(Quadrat_year, into=c("Quadrat", "Year"), sep = "_") %>% 
    mutate(Year=as.numeric(Year)) %>% 
    pivot_longer(D:D2, values_to = "D", names_to="Component")
  dd <- bind_rows(dd, bb)
}
dd <- dd %>% 
  left_join(flora %>% 
              distinct(Quadrat, Treatment, Series), 
            by=c("Quadrat"))

ggplot(data=dd, 
       aes(x=Year, y=D, col=Treatment, fill=Treatment, #linetype=Series, 
           )) +
  geom_line(aes(group=Quadrat)) + 
  #geom_smooth(alpha=0.2, se = T) + 
  theme_bw() + 
  facet_grid(Component~.) + 
  scale_x_continuous(breaks = seq(2012, 2022, by=2)) +
  scale_color_brewer(palette = "Dark2", direction = -1) + 
  scale_fill_brewer(palette = "Dark2", direction = -1)
#  scale_color_manual(values = wes_palette("Moonrise3")) + 
#  scale_fill_manual(values = wes_palette("Moonrise3")) 








# ## Dispersion from stand centroid
# library(vegan)
# plot_index <- sapply(rownames(flora_wide), function(x){strsplit(x, split = "_")[[1]][[1]]})
# 
# 
# flora_dist <- vegan::vegdist(flora_wide, method="bray")
# 
# flora_betadisper <- vegan::betadisper(flora_dist, group=plot_index, type="centroid")
# anova(flora_betadisper) # no overall difference among plots
# plot(flora_betadisper)
# betadist <- flora_betadisper$distances
# 
# 
# ## Variation partitioning of forest plots
# 
# flora_dist <- vegan::vegdist(flora_wide[str_detect(rownames(flora_wide), "20_|40_|^G"),], 
#                              method="jaccard")
# plot_index <- sapply(rownames(as.matrix(flora_dist)), 
#                      function(x){strsplit(x, split = "_")[[1]][[1]]})
# plot_index <- sapply(rownames(as.matrix(flora_dist)), 
#                      function(x){str_extract(x, "20_|40_|^G")})
# 
# 
# flora_betadisper <- vegan::betadisper(flora_dist, group=plot_index, type="centroid")
# anova(flora_betadisper)
# plot(flora_betadisper)
# betadist <- flora_betadisper$distances
# 
# k <- 3 #number of groups
# n <- 33 #observations per group
# N <- k*n
# GrandMean <- mean(betadist)
# GroupMeans <- tapply(betadist, plot_index, "mean")
# (SSB <- sum (n*(GroupMeans-GrandMean)^2))
# (SSE <- sum ((betadist - rep(GroupMeans, each=n))^2))
# #flora_ss <- (flora_wide - colMeans(flora_wide))^2
# (SST <- sum( (betadist-GrandMean)^2))
# 
# (MSB <- SSB/(k-1))
# (MSE <- SSE/(N-k))
# 
# (MSadd <- (MSB-MSE)/n)
# (Var_add <- MSadd/(MSE+MSadd) * 100)
# 
# #Var_flora <- SST / (ncol(flora_wide)-1)
# 
# betadist_df <- betadist %>% 
#   as.data.frame() %>% 
#   rename("betadist"=1) %>% 
#   rownames_to_column("Quadrat_year") %>% 
#   as_tibble() %>% 
#   separate(Quadrat_year, into=c("Quadrat", "Year"), sep="_") %>% 
#   mutate(Treatment=str_extract(Quadrat, "20|40|^G")) %>% 
#   as.data.frame()
# 
# #By Treatment
# mod1 <- lm(betadist ~ Treatment*Year, data=betadist_df)
# aa <- (aov(betadist ~ Treatment*Year, data=betadist_df))
# summary(aa)
# bb <- summary(aa)[[1]][[3]]
# MSadd <- (bb[1:3]-bb[4])/n
# (Var_add <- MSadd/(MSE+MSadd) * 100)
# 
# 
# #By year
# aa <- (aov(betadist ~ Year, data=betadist_df))
# summary(aa)
# bb <- summary(aa)[[1]][[3]]
# MSadd <- (bb[1]-bb[2])/n
# (Var_add <- MSadd/(MSE+MSadd) * 100)


### Decompose the SS of the species x plot matrix as in Legendre & De Caceres
### Try to decompose total variation as in a Random effect Anova

plot_index <- sapply(rownames(flora_wide), function(x){strsplit(x, split = "_")[[1]][[1]]})
flora_dist <- vegan::vegdist(flora_wide, method="hellinger")
mydata <- data.frame(qyear=names(flora_dist)) %>% 
  separate(qyear, into=c("Quadrat", "Year"), sep="_") %>% 
  left_join(flora %>% 
              distinct(Quadrat, Treatment), 
            by="Quadrat")
adonis2(formula= flora_dist ~ Treatment*Year + Quadrat, data=mydata, strata=mydata$Quadrat)







### Partitioning beta diversity into turnover and nestedness
# Obj: test to what extent temporal change between t and t+1 
# is due to nestedness and turnover, and assess whether this 
# proportion 1) varies b/w gap and forest interior plots and 
# 2) varies with time

# I expect variation is mostly due to nestedness in forest interior
# and tunrover in gaps. I also expect that the first years in the 
# forest gap are dominated by nestedness, and turnover increases
# after colonization (which takes a few years)




flora_betapart <- betapart::beta.pair((flora_wide>0)*1)

flora_beta2 <- flora_betapart$beta.sne %>% 
  as.matrix() %>% 
  as.data.frame() %>% 
  rownames_to_column("Quadrat_year") %>% 
  pivot_longer(-Quadrat_year, names_to = "Quadrat_year2", values_to="beta_sne") %>% 
  left_join(flora_betapart$beta.sim %>% 
              as.matrix() %>% 
              as.data.frame() %>% 
              rownames_to_column("Quadrat_year") %>% 
              pivot_longer(-Quadrat_year, names_to = "Quadrat_year2", values_to="beta_sim"), 
            by=c("Quadrat_year", "Quadrat_year2")) %>% 
  left_join(flora_betapart$beta.sor %>% 
              as.matrix() %>% 
              as.data.frame() %>% 
              rownames_to_column("Quadrat_year") %>% 
              pivot_longer(-Quadrat_year, names_to = "Quadrat_year2", values_to="beta_sor"), 
            by=c("Quadrat_year", "Quadrat_year2")) %>% 
  separate(Quadrat_year, into=c("Quadrat1", "Year1")) %>% 
  separate(Quadrat_year2, into=c("Quadrat2", "Year2")) %>% 
  filter(Quadrat1==Quadrat2) %>% 
  select(-Quadrat2, Quadrat=Quadrat1) %>% 
  mutate_at(.vars=vars(starts_with("Year")), 
            .funs = list(~as.numeric(.))) %>% 
  #filter(Year1== Year2-1) %>% 
  mutate(Treatment=str_extract(Quadrat, pattern="G|20|40")) %>% 
  mutate(Treatment=fct_recode(factor(Treatment, 
                                     levels=c("G", "20", "40")), 
                              "Gap"="G", 
                              "Margin"="20", 
                              "Interior"="40"))


#ggplot(data=flora_beta2 %>% 
#         pivot_longer(beta_sne:beta_sor, names_to="beta_component", values_to#="beta")) + 
#  geom_density(aes(x=beta, group=Treatment, col=Treatment)) + 
#  theme_bw() + 
#  facet_grid(beta_component~.)

ggplot(data=flora_beta2 %>% 
         filter(Year1== Year2-1) %>% 
         pivot_longer(beta_sne:beta_sor, names_to="beta_component", values_to="beta")) + 
  geom_boxplot(aes(y=beta, group=Treatment, col=Treatment)) + 
  theme_bw() + 
  facet_grid(beta_component~.)

ggplot(data=flora_beta2 %>% 
         filter(Year1== Year2-1),# %>% 
         #mutate(Prop_turnover=beta_sim/beta_sor*100) %>% 
         #pivot_longer(Prop_turnover, names_to="beta_component", values_to="beta"
  aes(x=Year2, y=beta_sor, col=Treatment, fill=Treatment)) + 
  geom_point() + 
  geom_smooth(alpha=0.2, method="lm") + 
  theme_bw() + 
  #facet_grid(beta_component~.) + 
  scale_x_continuous(breaks = seq(2012, 2022, by=2)) + 
  scale_color_brewer(palette = "Dark2", direction = -1) + 
  scale_fill_brewer(palette = "Dark2", direction = -1)


## Partitioning change between nestedness and turnover betwene t and t0
ggplot(data=flora_beta2 %>% 
         filter(Year1==2012) %>%
         pivot_longer(beta_sne:beta_sor, 
                      names_to="beta_component", values_to="Beta Diversity") %>% 
         mutate(beta_component=factor(beta_component, 
                                      levels=c("beta_sim", "beta_sne", "beta_sor"), 
                                      labels=c("Turnover", "Nestedness", "Total"))) %>% 
         mutate(beta_component=fct_rev(beta_component)) %>% 
         dplyr::rename(Year=Year2) %>% 
         left_join(flora %>% 
                     distinct(Quadrat, Treatment, Series), 
                   by=c("Quadrat", "Treatment")), 
       aes(y=`Beta Diversity`, x=Year, col=Treatment, fill=Treatment)) + 
  geom_point() + 
  geom_smooth(method="lm", se=T, alpha=0.2) +
  theme_bw() + 
  facet_grid(beta_component~.) + 
  #scale_color_discrete(name="Beta Diversity\nComponent") + 
  scale_x_continuous(breaks = seq(2012, 2022, by=2)) + 
  scale_color_brewer(palette = "Dark2", direction = -1) + 
  scale_fill_brewer(palette = "Dark2", direction = -1)


##Plot temporal decay

ggplot(data=flora_beta2 %>% 
         filter(Year1<Year2) %>% 
         mutate(Delta_Year=Year2-Year1) %>% 
         pivot_longer(beta_sne:beta_sor, names_to="beta_component", values_to="beta"), 
       aes(x=Delta_Year, y=beta, alpha=0.3)) + 
  geom_point() + 
  geom_smooth() + 
  theme_bw() + 
  facet_grid(Treatment~beta_component) + 
  scale_x_continuous(breaks = seq(1,10, by=2))

