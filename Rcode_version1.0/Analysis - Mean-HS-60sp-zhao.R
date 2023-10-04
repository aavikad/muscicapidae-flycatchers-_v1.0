## PGLS analysis for 60 Muscicapidae Flycatchers species ZHAO tree ##
### Ecology and not phylogeny influences sensitivity to climate change in Muscicapidae Flycatchers in Eastern Himalayan and Indo-Burman hotspots
## Aavika Dhanda, Michał T. Jezierski1, Tim Coulson, Sonya Clegg ##

library(ggpubr)
library(dplyr)
library(cli)
library(nlme)
library(ape)
library(MuMIn)
library(rlang)
library(tidyverse)
library(scales)
library(vctrs)
library(phylosignal)
library(phylobase)
library(ggtree)
library(BiocManager)

## Step 1: Load in the tree and data
zhao_tree <- read.tree('zhao_60consensus.nex')  ## Zhao: https://www.sciencedirect.com/science/article/pii/S1055790322002597
ssi<- read.csv("ssi.csv")

## PRESENT CLIMATIC CONDITIONS ### ====
## Step2a:run  phylogenetic model with bio8 as a response variable =====
presentvars<- read.csv('bio_mean_all_species60.csv') ## only considering variables with most contribution to present maxent models-bio8 ##
presentvars<- cbind(presentvars, ssi$s.)
colnames(presentvars)[colnames(presentvars) == "ssi$s."] <- "S."
View(presentvars)

### checking model residuals ##
mod1<- lm(bio8_contribution ~ S., data=presentvars)
res1<- resid(mod1)
plot(fitted(mod1), res1)
abline(0,0)
qqnorm(res1)
qqline(res1) #near normal

rownames(presentvars) <- presentvars$Species
zhao_8 <- subset(presentvars, select = bio8_contribution)
zhao_tree_match <- keep.tip(zhao_tree, tip = rownames(zhao_8))
mat_zhao <- vcv(zhao_tree_match, corr = T)
tip.heights <- diag(mat_zhao)
fit_zhao_1 <- gls(bio8_contribution ~ S.,
                    correlation=corBrownian(phy=zhao_tree_match, form = ~Species), weights = varFixed(~tip.heights), data=presentvars)
summary(fit_zhao_1) #p-value< 0.05 significant

## Step2b:run  phylogenetic model with elevation as a response variable =====
mod2<- lm(elev_contribution ~ S., data=presentvars)
res2<- resid(mod2)
plot(fitted(mod2), res2)
abline(0,0)
qqnorm(res2)
qqline(res2) #near normal

zhao_elev <- subset(presentvars, select = elev_contribution)
zhao_tree_match2 <- keep.tip(zhao_tree, tip = rownames(zhao_elev))
mat_zhao <- vcv(zhao_tree_match2, corr = T)
tip.heights <- diag(mat_zhao)
fit_zhao_2 <- gls(elev_contribution ~ S.,
                  correlation=corBrownian(phy=zhao_tree_match2, form = ~Species), weights = varFixed(~tip.heights), data= presentvars)
summary(fit_zhao_2) #p-value> 0.05 not significant

### Explore phylogenetic signal in bio8 (since it resulted in significant model)
zhao_tree_match <- keep.tip(zhao_tree, tip = rownames(zhao_8))
zhao_tree_ph4 <- as(zhao_tree_match, 'phylo4')
zhao_ph4d_hs <- phylobase::phylo4d(zhao_tree_ph4, zhao_8)
### Global phylogenetic signal
phyloSignal(zhao_ph4d_hs)
### Global phylogenetic correlogram
hscor_zhao_hs <- phyloCorrelogram(zhao_ph4d_hs)
plot(hscor_zhao_hs)

local.i_zhao_ph4d_hs <- lipaMoran(zhao_ph4d_hs,prox.phylo = 'nNodes', as.p4d = T)
points.col_zhao_hs <- lipaMoran(zhao_ph4d_hs, prox.phylo = 'nNodes')$p.value
points.col_zhao_hs <- ifelse(points.col_zhao_hs < 0.05, 'red', 'black')
dotplot_zhao_ph4d_hs<- dotplot.phylo4d(local.i_zhao_ph4d_hs, dot.col = points.col_zhao_hs,
                                       trait.bg.col = "white",
                                       tip.labels=sp_names$scientific_name,
                                       data.xlim=c(-6,6))


### FUTURE CLIMATIC CONDITIONS ### =====
### checking model residuals ##
futurevars<- read.csv("future60.csv")
futurevars<- cbind(futurevars, ssi$s.)
colnames(futurevars)[colnames(futurevars) == "ssi$s."] <- "S."
View(futurevars)

mod3<- lm(tmax_contribution126~ S., data=futurevars)
res3<- resid(mod3)
plot(fitted(mod3), res3)
abline(0,0)
qqnorm(res3)
qqline(res3)# near normal

mod4<- lm(tmin_contribution126~ S., data=futurevars)
res4<- resid(mod4)
plot(fitted(mod4), res4)
abline(0,0)
qqnorm(res4)
qqline(res4) # not normal

mod5<- lm(tmax_contribution245~ S., data=futurevars)
res5<- resid(mod5)
plot(fitted(mod5), res5)
abline(0,0)
qqnorm(res5)
qqline(res5) #near normal

mod6<- lm(tmin_contribution245~ S., data=futurevars)
res6<- resid(mod6)
plot(fitted(mod6), res6)
abline(0,0)
qqnorm(res6)
qqline(res6) #not normal

mod7<- lm(tmax_contribution370~ S., data=futurevars)
res7<- resid(mod7)
plot(fitted(mod7), res7)
abline(0,0)
qqnorm(res7)
qqline(res7) #near normal

mod8<- lm(tmin_contribution370~ S., data=futurevars)
res8<- resid(mod8)
plot(fitted(mod8), res8)
abline(0,0)
qqnorm(res8)
qqline(res8) # not normal

mod9<- lm(tmax_contribution585~ S., data=futurevars)
res9<- resid(mod9)
plot(fitted(mod9), res9)
abline(0,0)
qqnorm(res9)
qqline(res9) #near normal

mod10<- lm(tmin_contribution585~ S., data=futurevars)
res10<- resid(mod10)
plot(fitted(mod10), res10)
abline(0,0)
qqnorm(res10)
qqline(res10) # not normal


### taking logs of variables which were not normal ##
future_log<- data.frame(log(futurevars$tmin_contribution126+1),
                        log(futurevars$tmin_contribution245+1),
                        log(futurevars$tmin_contribution370+1),
                        log(futurevars$tmin_contribution585+1))

future_log<- cbind(futurevars$species, future_log, futurevars$tmax_contribution126,
                   futurevars$tmax_contribution245, futurevars$tmax_contribution370,
                   futurevars$tmax_contribution585, futurevars$S.)

### data frame with tmax_contribution and tmin_contribution under different scenarios ###

names(future_log)[1]<- "species"
names(future_log)[2]<- "tmin_contribution126"
names(future_log)[3]<- "tmin_contribution245"
names(future_log)[4]<- "tmin_contribution370"
names(future_log)[5]<- "tmin_contribution585"
names(future_log)[6]<- "tmax_contribution126"
names(future_log)[7]<- "tmax_contribution245"
names(future_log)[8]<- "tmax_contribution370"
names(future_log)[9]<- "tmax_contribution585"
names(future_log)[10]<- "S."
View(future_log)

mod4a<- lm(tmin_contribution126~ S., data=future_log)
res4a<- resid(mod4a)
plot(fitted(mod4a), res4a)
abline(0,0)
qqnorm(res4a)
qqline(res4a) # not normal

mod6a<- lm(tmin_contribution245~ S., data=future_log)
res6a<- resid(mod6a)
plot(fitted(mod6a), res6a)
abline(0,0)
qqnorm(res6a)
qqline(res6a) # not normal

mod7a<- lm(tmin_contribution370~ S., data=future_log)
res7a<- resid(mod7a)
plot(fitted(mod7a), res7a)
abline(0,0)
qqnorm(res7a)
qqline(res7a) #not normal

mod8b<- lm(tmin_contribution585~ S., data=future_log)
res8b<- resid(mod8b)
plot(fitted(mod8b), res8b)
abline(0,0)
qqnorm(res8b)
qqline(res8b) # not normal

######  Step 3: run  phylogenetic model with tmax as a response variable =====
####### ssp 126 ###### 

rownames(future_log) <- future_log$species
zhao_tmax126 <- subset(future_log, select = tmax_contribution126)
zhao_tree_match3 <- keep.tip(zhao_tree, tip = rownames(zhao_tmax126))
mat_zhao <- vcv(zhao_tree_match3, corr = T)
tip.heights <- diag(mat_zhao)
fit_zhao_3 <- gls(tmax_contribution126 ~ S.,
                  correlation=corBrownian(phy=zhao_tree_match3, form = ~species), weights = varFixed(~tip.heights), data= future_log)
summary(fit_zhao_3) #p-value< 0.05 significant

### Explore phylogenetic signal in max_temp126 (since it resulted in a significant model) 
zhao_tree_ph3 <- as(zhao_tree_match3, 'phylo4')
zhao_ph4d_hs3 <- phylobase::phylo4d(zhao_tree_ph3, zhao_tmax126) 
### Global phylogenetic signal
phyloSignal(zhao_ph4d_hs3)
### Global phylogenetic correlogram
hscor_zhao_hs3 <- phyloCorrelogram(zhao_ph4d_hs3)
plot(hscor_zhao_hs3)

local.i_zhao_ph4d_hs3 <- lipaMoran(zhao_ph4d_hs3,prox.phylo = 'nNodes', as.p4d = T)
points.col_zhao_hs3 <- lipaMoran(zhao_ph4d_hs3, prox.phylo = 'nNodes')$p.value
points.col_zhao_hs3 <- ifelse(points.col_zhao_hs3 < 0.05, 'red', 'black')
dotplot_zhao_ph4d_hs3<- dotplot.phylo4d(local.i_zhao_ph4d_hs3, dot.col = points.col_zhao_hs3,
                                        trait.bg.col = "white",
                                        tip.labels=sp_names$scientific_name,
                                        data.xlim=c(-6,6))

####### ssp 245 ######
zhao_tmax245 <- subset(future_log, select = tmax_contribution245)
zhao_tree_match4 <- keep.tip(zhao_tree, tip = rownames(zhao_tmax245))
mat_zhao <- vcv(zhao_tree_match4, corr = T)
tip.heights <- diag(mat_zhao)

fit_zhao_4 <- gls(tmax_contribution245 ~ S.,
                  correlation=corBrownian(phy=zhao_tree_match4, form = ~species), weights = varFixed(~tip.heights), data= future_log)
summary(fit_zhao_4) #p-value <0.05 significant

## Explore phylogenetic signal in max_temp245 (since it resulted in significant model)
zhao_tree_ph4 <- as(zhao_tree_match4, 'phylo4')
zhao_ph4d_hs4 <- phylobase::phylo4d(zhao_tree_ph4, zhao_tmax245)
### Global phylogenetic signal
phyloSignal(zhao_ph4d_hs4)
### Global phylogenetic correlogram
hscor_zhao_hs4 <- phyloCorrelogram(zhao_ph4d_hs4)
plot(hscor_zhao_hs4)

local.i_zhao_ph4d_hs4 <- lipaMoran(zhao_ph4d_hs4,prox.phylo = 'nNodes', as.p4d = T)
points.col_zhao_hs4 <- lipaMoran(zhao_ph4d_hs4, prox.phylo = 'nNodes')$p.value
points.col_zhao_hs4 <- ifelse(points.col_zhao_hs4 < 0.05, 'red', 'black')
dotplot_zhao_ph4d_hs4<- dotplot.phylo4d(local.i_zhao_ph4d_hs4, dot.col = points.col_zhao_hs4,
                                        trait.bg.col = "white",
                                        tip.labels=sp_names$scientific_name,
                                        data.xlim=c(-6,6))

####### ssp 370 ######
zhao_tmax370 <- subset(future_log, select = tmax_contribution370)
zhao_tree_match5 <- keep.tip(zhao_tree, tip = rownames(zhao_tmax370))
mat_zhao <- vcv(zhao_tree_match5, corr = T)
tip.heights <- diag(mat_zhao)

fit_zhao_5 <- gls(tmax_contribution370 ~ S.,
                  correlation=corBrownian(phy=zhao_tree_match5, form = ~species), weights = varFixed(~tip.heights), data= future_log)
summary(fit_zhao_5) #p-value <0.05 significant

## Explore phylogenetic signal in max_temp370 (since it resulted in a significant model)
zhao_tree_ph5 <- as(zhao_tree_match5, 'phylo4')
zhao_ph4d_hs5 <- phylobase::phylo4d(zhao_tree_ph5, zhao_tmax370)
### Global phylogenetic signal
phyloSignal(zhao_ph4d_hs5)
### Global phylogenetic correlogram
hscor_zhao_hs5 <- phyloCorrelogram(zhao_ph4d_hs5)
plot(hscor_zhao_hs5)

local.i_zhao_ph4d_hs5 <- lipaMoran(zhao_ph4d_hs5,prox.phylo = 'nNodes', as.p4d = T)
points.col_zhao_hs5 <- lipaMoran(zhao_ph4d_hs5, prox.phylo = 'nNodes')$p.value
points.col_zhao_hs5 <- ifelse(points.col_zhao_hs5 < 0.05, 'red', 'black')
dotplot_zhao_ph4d_hs5<- dotplot.phylo4d(local.i_zhao_ph4d_hs5, dot.col = points.col_zhao_hs5,
                                        trait.bg.col = "white",
                                        tip.labels=sp_names$scientific_name,
                                        data.xlim=c(-6,6))

####### ssp 585 ######
zhao_tmax585 <- subset(future_log, select = tmax_contribution585)
zhao_tree_match6 <- keep.tip(zhao_tree, tip = rownames(zhao_tmax585))
mat_zhao <- vcv(zhao_tree_match6, corr = T)
tip.heights <- diag(mat_zhao)

fit_zhao_6 <- gls(tmax_contribution585 ~ S.,
                   correlation=corBrownian(phy=zhao_tree_match6, form = ~species), weights = varFixed(~tip.heights), data= future_log)
summary(fit_zhao_6) #p-value <0.05 significant

## Explore phylogenetic signal in tmax_contribution585 (since it resulted in a significant model)
zhao_tree_ph6 <- as(zhao_tree_match6, 'phylo4')
zhao_ph4d_hs6 <- phylobase::phylo4d(zhao_tree_ph6, zhao_tmax585)
### Global phylogenetic signal
phyloSignal(zhao_ph4d_hs6)
### Global phylogenetic correlogram
hscor_zhao_hs6 <- phyloCorrelogram(zhao_ph4d_hs6)
plot(hscor_zhao_hs6)

local.i_zhao_ph4d_hs6 <- lipaMoran(zhao_ph4d_hs6,prox.phylo = 'nNodes', as.p4d = T)
points.col_zhao_hs6 <- lipaMoran(zhao_ph4d_hs6, prox.phylo = 'nNodes')$p.value
points.col_zhao_hs6 <- ifelse(points.col_zhao_hs6 < 0.05, 'red', 'black')
dotplot_zhao_ph4d_hs6<- dotplot.phylo4d(local.i_zhao_ph4d_hs6, dot.col = points.col_zhao_hs6,
                                        trait.bg.col = "white",
                                        tip.labels=sp_names$scientific_name,
                                        data.xlim=c(-6,6))

###### Step 4 -- Make phylogenetic tree with the distribution of bio8 contribution ========

sci_names_zhao <-gsub('_', ' ', zhao_tree_match$tip.label) # new scientific names without '_' 
zhao_tree_match$tip.label <- sci_names_zhao # replace scientific names with _ to without
presentvars$Species<- gsub('_', ' ', presentvars$Species)
rownames(presentvars)<- presentvars$Species
View(presentvars)

zhao_tree_match$tip.label[match(presentvars$Species,zhao_tree_match$tip.label)]

f<- ggtree(zhao_tree_match, layout= "circular", size= 1.0, branch.length = 'none', aes(color = bio8_contribution) , 
       continuous= 'colour') %<+% presentvars +scale_color_viridis_c() +
  geom_tiplab(size=4, hjust = -0.2)+
  theme(legend.position="left") +
  xlim(-15, 30)
f


########### Step 5-- Make phylogenetic trees with distribution of tmax contribution ===

sci_names_zhao <-gsub('_', ' ', zhao_tree_match$tip.label) # new scientific names without '_' 
zhao_tree_match$tip.label <- sci_names_zhao # replace scientific names with _ to without
futurevars<- read.csv("future60.csv")
futurevars$species<- gsub('_', ' ', futurevars$species)
rownames(futurevars)<- futurevars$species
View(futurevars)

zhao_tree_match$tip.label[match(futurevars$species,zhao_tree_match$tip.label)]


g <- ggtree(zhao_tree_match, layout= "circular", size= 1.0, branch.length = 'none', aes(color = tmax_contribution126) , 
            continuous= 'colour') %<+% futurevars +scale_color_viridis_c() +
  geom_tiplab(size=4, hjust = -0.2)+
  theme(legend.position="left") +
  xlim(-15, 30)
g

h <- ggtree(zhao_tree_match, layout= "circular", size= 1.0, branch.length = 'none', aes(color = tmax_contribution245) , 
            continuous= 'colour') %<+% futurevars +scale_color_viridis_c() +
  geom_tiplab(size=4, hjust = -0.2)+
  theme(legend.position="left") +
  xlim(-15, 30)
h

i <- ggtree(zhao_tree_match, layout= "circular", size= 1.0, branch.length = 'none', aes(color = tmax_contribution370) , 
            continuous= 'colour') %<+% futurevars +scale_color_viridis_c() +
  geom_tiplab(size=4, hjust = -0.2)+
  theme(legend.spacing= unit(3,"cm"),
        legend.position="left") +
  xlim(-15, 30)
i

j <- ggtree(zhao_tree_match, layout= "circular", size= 1.0, branch.length = 'none', aes(color = tmax_contribution585) , 
            continuous= 'colour') %<+% futurevars +scale_color_viridis_c() +
  geom_tiplab(size=4, hjust = -0.2)+
  theme(legend.spacing= unit(3,"cm"),
        legend.position="left") +
  xlim(-15, 30)
j