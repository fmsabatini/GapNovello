library(tidyverse)
library(betapart)
library(patchwork)
library(wesanderson)
library(vegan)
library(ggrepel)

### Ancillary functions ####
## Calculate rolling mean of change
ch <- function(x){x[length(x)]-x[1]}

## Pool species cover across consecutive years
pool <- function(x){ifelse(sum(!is.na(x))>=2, max(x, na.rm=T), NA)}


# Combine cover accounting for layers
combine.cover <- function(x){
  while (length(x)>1){
    x[2] <- x[1]+(100-x[1])*x[2]/100
    x <- x[-1]
  }
  return(x)
}

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

plot_index <- sapply(rownames(flora_wide), function(x){strsplit(x, split = "_")[[1]][[1]]})


## Pool species entries using a 3-year moving window
flora_pooled <- flora %>% 
  filter(Layer=="U") %>% 
  dplyr::select(Year, Quadrat, Species_resolved, Cover) %>% 
  complete(Year, nesting(Quadrat, Species_resolved)) %>% 
  arrange(Quadrat, Species_resolved, Year) %>% 
  group_by(Quadrat, Species_resolved) %>%
  mutate(Cover_roll=rollapply(Cover, 3, pool, align="center", fill=NA)) %>% 
  mutate(Cover_roll=coalesce(Cover, Cover_roll)) %>% 
  mutate(Cover_roll=replace(Cover_roll, list=is.na(Cover_roll), values=0)) %>% 
  dplyr::select(-Cover) %>% 
  ungroup()



### General Description ####
### Graph of species abundances with time (by plot)
ggplot(data=flora_pooled %>% 
         #filter(Layer=="U") %>% 
         filter(Quadrat=="GN") %>% 
         mutate(Species_resolved=factor(Species_resolved)) %>% 
         left_join({.} %>% 
                     count(Species_resolved), 
                   by=c("Species_resolved")) %>% 
         group_by(Species_resolved, Year) %>% 
         filter(n>2) %>% 
         summarize(Summed_cover=sum(Cover_roll))) + 
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

flora_eff <- flora_pooled %>% 
  group_by(Quadrat, Year) %>% 
  filter(Cover_roll>0) %>% 
  summarize(SR=n()) %>%
  arrange(Quadrat, Year) %>% 
  left_join(flora_pooled %>% 
              filter(Cover_roll>0) %>% 
              group_by(Quadrat) %>% 
              distinct(Species_resolved, .keep_all=T) %>% 
              summarize(SR_acc=n()), 
            by="Quadrat") %>% 
  left_join(flora %>% 
              distinct(Quadrat, Series, Treatment), 
            by="Quadrat") %>% 
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
# Set contrasts for year variable (adjust according to your data)

mod1 <- lme(DeltaSR ~ Treatment * LogDeltaYear, data=floraSR, random=~1|Quadrat)
summary(mod1)
mod2 <- lme(DeltaSR ~ Treatment * LogDeltaYear + 0, data=floraSR, 
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

## code to predict response at round intervals

#newdat <- expand.grid(Treatment=unique(floraSR$Treatment),
#                      LogDeltaYear=log(0:11+1))
#newdat$pred <- predict(mod2, newdata=newdat,level=0)
#newdat %>% 
#  mutate(DeltaYear=exp(LogDeltaYear)-1) %>% 
#  arrange(Treatment, LogDeltaYear) %>% 
#  group_by(Treatment) %>% 
#  mutate(diff = pred - lag(pred, n=1L))



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

## code to predict response at round intervals

newdat <- expand.grid(Treatment=unique(floraSR$Treatment),
                      LogDeltaYear=log(0:12+1))
newdat$pred <- predict(mod_cov2, newdata=newdat,level=0)
newdat %>% 
  mutate(DeltaYear=exp(LogDeltaYear)-1) %>% 
  arrange(Treatment, LogDeltaYear) %>% 
  group_by(Treatment) %>% 
  mutate(diff = pred - lag(pred, n=1L))



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





###
library(lme4)
mod1b <- lmer(DeltaSR ~ Treatment * LogDeltaYear + (1|Quadrat), data=floraSR)
partR2(mod1b, 
       partvars = c("Treatment", 
                    "LogDeltaYear", 
                    "Treatment:LogDeltaYear"), 
       R2_type="marginal")


##variation decomposition to account for interactions as suggested in Stoffel, Fig 2 C
mod1c <- lmer(DeltaSR ~Treatment * LogDeltaYear + (1|Quadrat), data=floraSR)
part1 <- partR2(mod1c, partvars = c("Treatment:LogDeltaYear"), nboot=100)
mod2c <- lmer(DeltaSR ~Treatment + LogDeltaYear + (1|Quadrat), data=floraSR)
part2 <- partR2(mod2c, partvars = c("Treatment", "LogDeltaYear"), nboot=100)
bb <- mergeR2(part1, part2)
forestplot(bb)
forestplot(bb, "IR2")



### NMDS   ####
# Compute dissimilarity matrix
flora_dist <- vegan::vegdist((flora_wide), method="hellinger")

# Compute NMDS
set.seed(1000)
flora_mds <- vegan::metaMDS(flora_dist, k=2, try = 100, trymax=100)

# Passively project species
flora_envfit0 <- envfit(flora_mds, flora_wide)
flora_envfit <- data.frame(flora_envfit0$vectors$arrows) %>%  
  bind_cols(data.frame(r2=flora_envfit0$vectors$r, p=flora_envfit0$vectors$pvals)) %>% 
  rownames_to_column("Species") %>% 
  filter(p<0.05) %>% 
  filter(r2>0.25) %>% 
  arrange(desc(r2))


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

ggplot(data=flora_mds_scores, aes(x=MDS1, y=MDS2)) +
#  geom_point(aes(x=MDS1, y=MDS2, col=Quadrat)) + 
  geom_path(aes(group=Quadrat, col=Treatment, linetype=Series), arrow=arrow(angle = 15, length = unit(0.15, "inches"))) + 
  scale_color_brewer(palette = "Dark2", direction = -1) + 
  geom_text(label=flora_mds_scores$Year, size=3, hjust=-0.2, show.legend=F,
            col= "black",
            check_overlap = F, 
            alpha=rep(c(1,0,0,1,0,0,1,0,0,0,1)/3, length.out=nrow(flora_mds_scores)), 
            fontface="bold") + 
  geom_point(data=flora_envfit, 
                  aes(x=NMDS1, y=NMDS2), pch=17) + 
#  geom_text_repel(data=flora_envfit, 
#             aes(x=NMDS1, y=NMDS2, label=Species), max.overlaps = 12) + 
  theme_bw() + 
  theme(panel.grid = element_blank())





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
SS.deco <- adonis2(formula= flora_dist ~ Treatment*Year + Quadrat, data=mydata, strata=mydata$Quadrat)

## Calculate LCBD
aa <- adespatial::beta.div(Y=flora_wide, method = "hellinger")

bb <- data.frame(
  Quadrat_year=names(aa$LCBD), 
  LCBD=aa$LCBD,
  p.LCBD=aa$p.LCBD) %>% 
  as_tibble() %>% 
  separate("Quadrat_year", into=c("Quadrat", "Year"), sep="_") %>% 
  mutate(Treatment=str_extract(Quadrat, pattern="G|20|40")) %>% 
  mutate(Treatment=fct_recode(factor(Treatment, 
                                     levels=c("G", "20", "40")), 
                              "Gap"="G", 
                              "Margin"="20", 
                              "Interior"="40")) %>% 
  arrange(Treatment)

ggplot(data=bb, 
       aes(x=Year, y=Quadrat, size=LCBD, col=p.LCBD>0.05)) + 
  geom_point()


ggplot(data=bb, 
       aes(x=Quadrat, y=LCBD, fill=Year, )) + 
  geom_bar(stat="identity")






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

# ggplot(data=flora_beta2 %>% 
#          filter(Year1== Year2-1) %>% 
#          pivot_longer(beta_sne:beta_sor, names_to="beta_component", values_to="beta")) + 
#   geom_boxplot(aes(y=beta, group=Treatment, col=Treatment)) + 
#   theme_bw() + 
#   facet_grid(beta_component~.)


## Partitioning change between nestedness and turnover betwene t and t0
ll <- ggplot(data=flora_beta2 %>% 
         filter(Year1==2012) %>%
         pivot_longer(beta_sne:beta_sor, 
                      names_to="beta_component", values_to="Beta Diversity") %>% 
         mutate(beta_component=factor(beta_component, 
                                      levels=c("beta_sim", "beta_sne", "beta_sor"), 
                                      labels=c("Turnover", "Nestedness", "Total"))) %>% 
         mutate(beta_component=fct_rev(beta_component)) %>% 
        filter(beta_component=="Total") %>% 
         dplyr::rename(Year=Year2) %>% 
         left_join(flora %>% 
                     distinct(Quadrat, Treatment, Series), 
                   by=c("Quadrat", "Treatment")), 
       aes(y=`Beta Diversity`, x=Year, col=Treatment, fill=Treatment)) + 
  geom_point() + 
  #geom_smooth(method="lm", se=T, alpha=0.2) +
  geom_smooth() + 
  theme_bw() + 
  facet_grid(Treatment ~.) + 
  #scale_color_discrete(name="Beta Diversity\nComponent") + 
  scale_x_continuous(breaks = seq(2012, 2022, by=2)) + 
  scale_color_brewer(palette = "Dark2", direction = -1) + 
  scale_fill_brewer(palette = "Dark2", direction = -1)

rr <- ggplot(data=flora_beta2 %>% 
               filter(Year1== Year2-1),# %>% 
             #mutate(Prop_turnover=beta_sim/beta_sor*100) %>% 
             #pivot_longer(Prop_turnover, names_to="beta_component", values_to="beta"
             aes(x=Year2, y=beta_sor, col=Treatment, fill=Treatment)) + 
  geom_point() + 
  #geom_smooth(alpha=0.2, method="lm") + 
  geom_smooth(alpha=0.2) +
  theme_bw() + 
  facet_grid(Treatment~.) + 
  scale_x_continuous(breaks = seq(2012, 2022, by=2)) + 
  scale_color_brewer(palette = "Dark2", direction = -1) + 
  scale_fill_brewer(palette = "Dark2", direction = -1)



ll +theme(legend.position="null") | rr



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



#### Life-history traits ####
elle <- read_csv("../rawdata/Ellenberg_Pignatti_out.csv")
checklist_elle <- checklist %>% 
  mutate(pignatti_match=species_resolved) %>% 
  mutate(pignatti_match=factor(pignatti_match)) %>% 
  mutate(pignatti_match=fct_recode(pignatti_match,
                                   `Festuca drymeia`="Festuca drymeja",
                                     `Hieracium sylvaticum`="Hieracium murorum", 
                                     `Mycelis muralis`="Lactuca muralis", 
                                     `Lamiastrum galeobdolon`="Lamium galeobdolon", 
                                     `Pulmonaria officinalis`="Pulmonaria apennina",
                                     `Rubus hirtus`="Rubus proiectus")) %>% 
  left_join(elle, by=c("pignatti_match"="Species")) %>% 
  mutate(GF=replace(GF, 
                    list=species_resolved %in% c("Hieracium murorum", "Trifolium pratense", "Fragaria vesca", "Veronica montana", "Veronica officinalis"), 
                    values="H")) %>% 
  mutate(LF=replace(LF, 
                    list=species_resolved %in% c("Hieracium murorum", "Trifolium pratense"), 
                    values="Scap"))  %>% 
  mutate(LF=replace(LF, 
                  list=species_resolved %in% c("Veronica montana", "Veronica officinalis"),
                  values="Rept")) %>% 
  dplyr::select(-c(1:5, 7)) %>% 
  rename(Species_resolved=species_resolved) %>% 
  distinct() %>% 
  mutate(Spring_ephemeral=F) %>% 
  mutate(Spring_ephemeral=replace(Spring_ephemeral, 
                                  list=Species_resolved %in% c("Adoxa moschatellina", "Corydalis cava", 
                                                               "Anemone apennina", "Cardamine bulbifera", 
                                                               "Anemome ranunculoides", "Cardamine enneaphyllos"), 
                                  values=TRUE))

### Decompose the SS of the species x plot matrix after excluding therophytes and spring ephemerals

## SS of the total matrix
plot_index <- sapply(rownames(flora_wide), function(x){strsplit(x, split = "_")[[1]][[1]]})
flora_dist <- vegan::vegdist(flora_wide, method="hellinger")
mydata <- data.frame(qyear=names(flora_dist)) %>% 
  separate(qyear, into=c("Quadrat", "Year"), sep="_") %>% 
  left_join(flora %>% 
              distinct(Quadrat, Treatment), 
            by="Quadrat")
SS.deco0 <- adonis2(formula= flora_dist ~ Treatment*Year + Quadrat, data=mydata, strata=mydata$Quadrat)


# exclude terophytes and spring ephemerals
flora_dist <- vegan::vegdist(flora_wide[,which(colnames(flora_wide) %in% 
                                                 (checklist_elle %>% 
                                                    filter(Spring_ephemeral==F & GF!="T") %>% 
                                                    pull(Species_resolved)))], 
                             method="hellinger")
mydata <- data.frame(qyear=names(flora_dist)) %>% 
  separate(qyear, into=c("Quadrat", "Year"), sep="_") %>% 
  left_join(flora %>% 
              distinct(Quadrat, Treatment), 
            by="Quadrat")
SS.deco2 <- adonis2(formula= flora_dist ~ Treatment*Year + Quadrat, data=mydata, strata=mydata$Quadrat)

#exclude all species never having more than 0.5 cover
flora_dist <- vegan::vegdist(flora_wide[,apply(flora_wide, "max", MARGIN=2)>0.5], 
                             method="hellinger")
mydata <- data.frame(qyear=names(flora_dist)) %>% 
  separate(qyear, into=c("Quadrat", "Year"), sep="_") %>% 
  left_join(flora %>% 
              distinct(Quadrat, Treatment), 
            by="Quadrat")
SS.deco3 <- adonis2(formula= flora_dist ~ Treatment*Year + Quadrat, data=mydata, strata=mydata$Quadrat)

## pool species across consecutive years to reduce omission error
##########à pool species cover across n consecutive years
flora_wide_pooled <- flora_pooled %>% 
  pivot_wider(names_from=Species_resolved, values_from = "Cover_roll", values_fill = 0) %>% 
  unite("Quadrat_year", Quadrat, Year) %>% 
  arrange(Quadrat_year) %>% 
  #View()
  column_to_rownames("Quadrat_year")

flora_dist <- vegan::vegdist(flora_wide_pooled, method="hellinger")

mydata <- data.frame(qyear=names(flora_dist)) %>% 
  separate(qyear, into=c("Quadrat", "Year"), sep="_") %>% 
  left_join(flora %>% 
              distinct(Quadrat, Treatment), 
            by="Quadrat")
SS.deco4 <- adonis2(formula= flora_dist ~ Treatment*Year + Quadrat, data=mydata, strata=mydata$Quadrat)

### decompose variation only for NON gap plots
toexclude <- which(mydata$Treatment == "Gap")
mydata_nongap <- mydata %>% 
  filter(Treatment != "Gap")

SS.deco_nongap <- adonis2(formula= as.dist(as.matrix(flora_dist)[-toexclude, -toexclude]) ~
                            Treatment*Year + Quadrat, data=mydata_nongap, 
                          strata=mydata_nongap$Quadrat)


parted <- data.frame(
  Fraction=rownames(SS.deco0),
  Share_raw=SS.deco$SumOfSqs,
  Share_pooled=SS.deco4$SumOfSqs, 
  Share_nongap=SS.deco_nongap$SumOfSqs)  %>% 
  bind_cols({.} %>% 
              summarize_at(.vars = vars(Share_raw:Share_nongap),
                           .funs = list("sum"=~sum(.)/2))
  ) %>% 
  mutate(Share_raw_perc=round(Share_raw/Share_raw_sum*100,2),
         Share_pooled_perc=round(Share_pooled/Share_pooled_sum*100,2), 
         Share_nongap_perc=round(Share_nongap/Share_nongap_sum*100,2)) 


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
                                  "Quadrat", 
                                  "Residual", 
                                  "Total"))) %>% 
  mutate(model=factor(model, levels=c("pooled", "nongap", "raw"))) %>% 
  arrange(Fraction)

ggplot(data=parted_long %>% 
         filter(Abs_Perc=="perc") %>% 
         filter(Fraction!="Total") %>% 
         filter(model=="raw")) +
  geom_col(aes(x=model, y=`Explained Variation (%)`, fill=Fraction, group=model), width=0.5) +
  theme_minimal() + 
  scale_fill_brewer(palette = "Dark2") + 
  #guides(fill = guide_legend(reverse=TRUE)) + 
  scale_x_discrete(name=NULL) + 
  coord_flip() + 
  theme(legend.position = "bottom", 
        legend.title = element_blank(), 
        axis.text.y = element_blank())

ggsave(last_plot(), file="../figure_tables/ExplainedVariation_rawonly.png", 
       width = 8, height=2, units="in",dpi = 300, bg="white")



ggplot(data=parted_long %>% 
         filter(Abs_Perc=="perc") %>% 
         filter(Fraction!="Total") %>% 
         filter(model!="nongap")) +
  geom_col(aes(x=model, y=`Explained Variation (%)`, fill=Fraction, group=model), width=0.5) +
  theme_minimal() + 
  scale_fill_brewer(palette = "Dark2") + 
  #guides(fill = guide_legend(reverse=TRUE)) + 
  scale_x_discrete(name=NULL) + 
  coord_flip() + 
  theme(legend.position = "bottom", 
        legend.title = element_blank())

ggsave(last_plot(), file="../figure_tables/ExplainedVariation.png", 
       width = 8, height=2, units="in",dpi = 300)


### True Diversity approach - to account for undersampling
Qjac <- function(x, q){
  source("D:/Nextcloud/MyFunctions/div.dec.R")
    x <- flora_wide
  nr <- nrow(x) #plots
  nc <- ncol(x) #species
  distij <- matrix(NA, nrow=nr, ncol=nr, dimnames = list(rownames(flora_wide), rownames(flora_wide)))
  diag(distij) <- 0
  
  for(i in 1:(nr-1)){
    for(j in (i+1):nr){
      xij <- flora_wide[c(i,j),]
      distij[i,j] <- distij[j,i] <- 1-((2/div.dec(xij, q=q, w="even")[[2]])-1)
    }
  }
  return(as.dist(distij))
}

dist0 <- Qjac(flora_wide, q=0)
dist2 <- Qjac(flora_wide, q=2)



flora_betaq <- dist0 %>% 
  as.matrix() %>% 
  as.data.frame() %>% 
  rownames_to_column("Quadrat_year") %>% 
  pivot_longer(-Quadrat_year, names_to = "Quadrat_year2", values_to="beta0") %>% 
  left_join(dist2 %>% 
              as.matrix() %>% 
              as.data.frame() %>% 
              rownames_to_column("Quadrat_year") %>% 
              pivot_longer(-Quadrat_year, names_to = "Quadrat_year2", values_to="beta2"), 
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
                              "Interior"="40")) %>% 
  pivot_longer(cols=beta0:beta2, names_to="q", values_to="beta")


### beta diversity between t and t0
(ll <- ggplot(data=flora_betaq %>% 
                filter(Year1==2012) %>%
                dplyr::rename(Year=Year2),  
              aes(y=beta, x=Year, col=Treatment, fill=Treatment)) + 
    geom_point() + 
    #geom_smooth(method="lm", se=T, alpha=0.2) +
    geom_smooth() + 
    theme_bw() + 
    facet_grid(q~.) + 
    #scale_color_discrete(name="Beta Diversity\nComponent") + 
    scale_x_continuous(breaks = seq(2012, 2023, by=2)) + 
    scale_color_brewer(palette = "Dark2", direction = -1) + 
    scale_fill_brewer(palette = "Dark2", direction = -1)
)


### Showing beta diversity between t and t-1
(rr <- ggplot(data=flora_betaq %>% 
               filter(Year1== Year2-1),# %>% 
             aes(x=Year2, y=beta, col=Treatment, fill=Treatment)) + 
  geom_point() + 
  #geom_smooth(alpha=0.2, method="lm") + 
  geom_smooth(alpha=0.2, method="lm") +
  theme_bw() + 
  facet_grid(q~.) + 
  scale_x_continuous(name="Year", breaks = seq(2012, 2023, by=2)) + 
  scale_color_brewer(palette = "Dark2", direction = -1) + 
  scale_fill_brewer(palette = "Dark2", direction = -1)
)


ll + theme(legend.position="none") | rr

mydata <- data.frame(qyear=names(dist0)) %>% 
  separate(qyear, into=c("Quadrat", "Year"), sep="_") %>% 
  left_join(flora %>% 
              distinct(Quadrat, Treatment), 
            by="Quadrat")
SS.deco0 <- adonis2(formula= dist0 ~ Treatment*Year + Quadrat, data=mydata, strata=mydata$Quadrat)
SS.deco2 <- adonis2(formula= dist2 ~ Treatment*Year + Quadrat, data=mydata, strata=mydata$Quadrat)
adonis2(formula= vegdist(flora_wide, method="jaccard") ~ Treatment*Year + Quadrat, data=mydata, strata=mydata$Quadrat)


parted_q <- data.frame(
  Fraction=rownames(SS.deco0),
  Share0=SS.deco0$SumOfSqs,
  Share2=SS.deco2$SumOfSqs)  %>% 
  bind_cols({.} %>% 
              summarize_at(.vars = vars(Share0:Share2),
                           .funs = list("sum"=~sum(.)/2))
  ) %>% 
  mutate(Share0_perc=round(Share0/Share0_sum*100,2),
         Share2_perc=round(Share2/Share2_sum*100,2))


  
parted_long_q <- parted_q %>% 
  dplyr::select(-ends_with("_sum")) %>% 
  pivot_longer(-Fraction, names_to="Variable", values_to="Explained Variation (%)") %>% 
  mutate(Variable=str_remove(Variable, pattern="^Share_")) %>% 
  separate(Variable, into=c("q", "Abs_Perc"), sep = "_") %>% 
  mutate(Abs_Perc=replace_na(Abs_Perc, "abs")) %>% 
  mutate(q=str_extract(q, pattern="[0-9]")) %>% 
  mutate(Fraction=factor(Fraction, 
                         levels=c("Treatment", 
                                  "Treatment:Year", 
                                  "Year", 
                                  "Quadrat", 
                                  "Residual", 
                                  "Total"))) %>% 
  mutate(q=factor(q, levels=c("2", "0"), labels=c("q = 2", "q = 0"))) %>% 
  arrange(Fraction)

ggplot(data=parted_long_q %>% 
         filter(Abs_Perc=="perc") %>% 
         filter(Fraction!="Total") ) +
  geom_col(aes(x=q, y=`Explained Variation (%)`, fill=Fraction, group=q), width=0.5) +
  theme_minimal() + 
  scale_fill_brewer(palette = "Dark2") + 
  #guides(fill = guide_legend(reverse=TRUE)) + 
  scale_x_discrete(name=NULL) + 
  coord_flip() + 
  theme(legend.position = "bottom", 
        legend.title = element_blank())

ggsave(last_plot(), file="../figure_tables/ExplainedVariation_Q.png", 
       width = 8, height=2, units="in",dpi = 300)
