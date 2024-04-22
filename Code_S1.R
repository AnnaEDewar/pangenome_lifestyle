## ----global-options, include=FALSE-------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo=FALSE, warning=FALSE, message=FALSE)
options(tinytex.compile.min_times = 3) # Allows correct page numbers in contents page


## ----packages, include=FALSE-------------------------------------------------------------------------------------------
library(tidyverse)
library(MCMCglmm)
library(ape)
library(phylotools)
library(phytools)
library(phylopath)
library(ggsignif)
library(patchwork)
library(treedataverse)
library(phylolm)
library(kableExtra)
library(FactoMineR)
library(factoextra)


## ----data & tree set-up, include=FALSE---------------------------------------------------------------------------------
#### Data 
pangenome_lifestyles <- read.csv("pangenome_lifestyles_all.csv", header=T)

#### Tree 
dataTree<-read.nexus("PanX_tree.nex")
is.ultrametric(dataTree) #SHould say 'TRUE'

## Format tree for bayestraits
is.rooted(dataTree)
dataTree


## check ALL species in data set 'pangenome_lifestyles' are in tree (should be 126)
pangenome_lifestyles$Species[which((pangenome_lifestyles$Species %in% dataTree$tip.label)==TRUE)]


#### MCMCglmm set-up

## Set priors
# Uninformative prior for when 1 random effect
prior <- list(R=list(V = 1, nu = 0.002), G=list(G1=list(V=1,nu=0.002))) 

# Uninformative prior for when 2 random effects
prior2 <- list(R=list(V = 1, nu = 0.002), G=list(G1=list(V=1,nu=0.002), G2=list(V=1,nu=0.002))) 

## Order lifestyle traits in order of 'least' to 'most' variable 
pangenome_lifestyles$Host_or_free <- factor(pangenome_lifestyles$Host_or_free, levels=c("Host","Both","Free","Unknown"))
pangenome_lifestyles$Obligate_facultative <- factor(pangenome_lifestyles$Obligate_facultative, levels=c("Obligate","Facultative","Unknown",""))
pangenome_lifestyles$Intra_or_extracellular <- factor(pangenome_lifestyles$Intra_or_extracellular, levels=c("Intracellular","Both","Extracellular", "Unknown", ""))
pangenome_lifestyles$Effect_on_host <- factor(pangenome_lifestyles$Effect_on_host, levels=c("Pathogen","Both","Mutualist","Unknown", ""))
pangenome_lifestyles$Motility <- factor(pangenome_lifestyles$Motility, levels=c("Non-motile","Both","Motile"))
pangenome_lifestyles$Category_primary_env <- factor(pangenome_lifestyles$Category_primary_env, levels=c("Host","Free","Unknown"))

## Make data set a 'data.frame' (to avoid any errors from data format)
pangenome_lifestyles$pangenome_fluidity <- pangenome_lifestyles$genome_fluidity

pangenome_lifestyles <- as.data.frame(pangenome_lifestyles)



## ----Plotting & printing functions, include=FALSE----------------------------------------------------------------------
# Function to plot 0s in axes nicely
prettyZero <- function(l){
  max.decimals = max(nchar(str_extract(l, "\\.[0-9]+")), na.rm = T)-1
  lnew = formatC(l, replace.zero = T, zero.print = "0",
                 digits = max.decimals, format = "f", preserve.width=T)
  return(lnew)
}

lm_summary_table <- function(model) {
  coeffs <- coef(summary(model))
  df <- c(summary(model)$df[1], summary(model)$df[2])
  
  # Add stars for significant p-values
  p_values <- coeffs[, 4]
  stars <- ifelse(p_values < 0.001, "***",
                  ifelse(p_values < 0.01, "**",
                         ifelse(p_values < 0.05, "*", "")))
  
  table_data <- data.frame(
    "Estimate" = format(coeffs[, 1], digits = 4, scientific = FALSE),
    "Std.Error" = format(coeffs[, 2], digits = 4, scientific = FALSE),
    "t value" = format(coeffs[, 3], digits = 4, scientific = FALSE),
    "p-value" = signif(coeffs[, 4], digits = 4),
    "signif." = stars
    
  )
  
  return(table_data)
  
}

summary_mcmc_glmm <- function(mcmc_model) {
  summary <- summary(mcmc_model)
  summary_subset <- as.data.frame(summary$solutions)
  
  p_values <- summary_subset$pMCMC
  stars <- ifelse(p_values < 0.001, "***",
                  ifelse(p_values < 0.01, "**",
                         ifelse(p_values < 0.05, "*",
                                ifelse(p_values < 0.1, ".",""))))
  summary_subset$signif. <- stars
  
  # Extract numeric columns
  num_cols <- sapply(summary_subset, is.numeric)
  numeric_table <- summary_subset[, num_cols]
  
  # Format numeric columns to 4 significant figures
  formatted_table <- format(numeric_table, scientific = FALSE, digits = 4)
  
  # Replace original numeric columns with formatted columns
  summary_subset[, num_cols] <- formatted_table
  
  return(summary_subset)
}



## ----Host/free vs pangenome fluidity lm, echo=FALSE--------------------------------------------------------------------
pangenome_lifestyles_host_free_no_unknown <- pangenome_lifestyles %>%
  filter(Host_or_free != "Unknown")

pangenome_lifestyles_host_free_no_unknown <- as.data.frame(pangenome_lifestyles_host_free_no_unknown)

summary_lm_1 <- lm_summary_table(lm(pangenome_fluidity~Host_or_free,data=pangenome_lifestyles_host_free_no_unknown))


knitr::kable(summary_lm_1, 
             caption= "Results from a linear model with pangenome fluidity as the response variable and host-associated, free-living or both as the explanatory variable.", digits=4) %>% kable_styling(latex_options = "HOLD_position")


## ----Host/free vs genome fluidity MCMCglmm, echo=TRUE, cache=TRUE------------------------------------------------------
## Label tree nodes with numbers
# 'dataTree' is our phylogeny as an ultrametric tree in nexus format
dataTreeNode <- makeNodeLabel(dataTree, method = "number")

## Make matrix of tree called 'Ainv'
INtree <- inverseA(dataTreeNode, nodes="TIPS") # Converts phylogeny into a covariance matrix
Ainv <- INtree$Ainv #Extracts covariance values

# Weakly informative prior 
prior <- list(R=list(V = 1, nu = 0.002), G=list(G1=list(V=1,nu=0.002))) 

mcmc_model_1 <- MCMCglmm(pangenome_fluidity ~ Host_or_free, random=~Species, 
                         data=pangenome_lifestyles_host_free_no_unknown,
                         nitt=50000, 
                         prior=prior, ginverse = list(Species = Ainv),verbose = FALSE)


## ----Host/free vs pangenome fluidity MCMCglmm output, echo=FALSE-------------------------------------------------------
summary_mcmc_model_1 <- summary_mcmc_glmm(mcmc_model_1)

mFixed <- mean(mcmc_model_1$Sol[,2]) * mcmc_model_1$X[, 2] + 
  mean(mcmc_model_1$Sol[,3]) * mcmc_model_1$X[, 3] 

mVarF<- var(mFixed)

# Fixed effect variance / (fixed effect variance + random effect variance)

fixed_r2_model_1 <- mVarF/(mVarF+sum(apply(mcmc_model_1$VCV,2,mean)))

#random effect variance/(random effect variance + residual variance)

random_r2_model_1 <- (sum(apply(mcmc_model_1$VCV,2,mean)[-2]))/((sum(apply(mcmc_model_1$VCV,2,mean)))+mVarF)

#Conditionall R2 (total variance explained by the model)

conditional_r2_model_1 <- ((sum(apply(mcmc_model_1$VCV,2,mean)[-2]))+mVarF)/((sum(apply(mcmc_model_1$VCV,2,mean)))+mVarF)

R2_table_1 <- matrix(c(fixed_r2_model_1,random_r2_model_1,conditional_r2_model_1), ncol=1)
R2_table_1 <- format(R2_table_1, digits=4)
rownames(R2_table_1) <- c("Fixed effect", "Random effect", "Total model")
colnames(R2_table_1) <- c("R-squared value")

knitr::kable(list(summary_mcmc_model_1,R2_table_1),
             caption="Results from the above MCMCglmm with pangenome fluidity as the response variable, whether a species is host-associated, free-living or both as the fixed effect, and phylogeny as a random effect.", digits=4) %>% kable_styling(latex_options = "HOLD_position")



## ----Host/Free/Both graph, echo=FALSE, fig.cap="Host-reliance correlates with pangenome fluidity.", fig.width=4, fig.height=3----
ggplot(pangenome_lifestyles_host_free_no_unknown, aes(x=Host_or_free,y=pangenome_fluidity)) +
  geom_jitter(width=0.25, aes(colour=Host_or_free)) +
  labs(y="Pangenome\nfluidity", x="Host- or free-living") +
  stat_summary(aes(y = pangenome_fluidity,group=1), fun=mean, geom="crossbar",width=0.4,colour="black") +  
  scale_y_continuous(limits=c(0,0.48),breaks=c(0,0.1,0.2,0.3,0.4), labels=prettyZero, expand=c(0,0)) +
  theme_classic() +
  scale_colour_manual(values=c("#CC79A7", "#56B4E9", "#009E73")) +
  geom_signif(comparisons=list(c("Both", "Host")), annotation="***",
              y_position = 0.42, tip_length = 0, vjust=0.4) +
  geom_signif(comparisons=list(c("Free", "Host")), annotation="**",
              y_position = 0.44, tip_length = 0, vjust=0.4) +
  theme(plot.title = element_text(size=15, hjust=0.5), axis.title = element_text(size=13),
        axis.text = element_text(size=12), axis.title.y=element_text(angle=360,vjust=0.5),
        legend.position="none")



## ----Primary env data--------------------------------------------------------------------------------------------------
pangenome_lifestyles_primary <- pangenome_lifestyles_host_free_no_unknown %>%
  select(Species, Primary_environment, Category_primary_env)

knitr::kable(head(pangenome_lifestyles_primary), 
             caption= "Example of species' primary environments, and their 'Host' or 'Free' catgeories.", digits=4)%>% kable_styling(latex_options = "HOLD_position")


## ----Primary environment lm, echo=FALSE--------------------------------------------------------------------------------
pangenome_lifestyles_primary_env_no_unknown <- pangenome_lifestyles %>%
  filter(Category_primary_env != "Unknown")

pangenome_lifestyles_primary_env_no_unknown <- as.data.frame(pangenome_lifestyles_primary_env_no_unknown)

summary_lm_2 <- lm_summary_table(lm(pangenome_fluidity~Category_primary_env,data=pangenome_lifestyles_primary_env_no_unknown))

knitr::kable(summary_lm_2, 
             caption= "Results from a linear model with pangenome fluidity as the response variable and whether a species' primary environment is host-associated or free-living as the explanatory variable.", digits=4)%>% kable_styling(latex_options = "HOLD_position")


## ----Primary env MCMCglmm, echo=FALSE, cache=TRUE----------------------------------------------------------------------
mcmc_model_2 <- MCMCglmm(pangenome_fluidity ~ Category_primary_env, random=~Species, data=pangenome_lifestyles_primary_env_no_unknown,
                         nitt=50000,
                         prior=prior, ginverse = list(Species = Ainv),verbose = FALSE)

# R squared
mFixed <- mean(mcmc_model_2$Sol[,2]) * mcmc_model_2$X[, 2]

mVarF<- var(mFixed)

# Fixed effect variance / (fixed effect variance + random effect variance)

fixed_r2_model_2 <- mVarF/(mVarF+sum(apply(mcmc_model_2$VCV,2,mean)))

#random effect variance/(random effect variance + residual variance)

random_r2_model_2 <- (sum(apply(mcmc_model_2$VCV,2,mean)[-2]))/((sum(apply(mcmc_model_2$VCV,2,mean)))+mVarF)

#Conditionall R2 (total variance explained by the model)

conditional_r2_model_2 <- ((sum(apply(mcmc_model_2$VCV,2,mean)[-2]))+mVarF)/((sum(apply(mcmc_model_2$VCV,2,mean)))+mVarF)


R2_table_2 <- matrix(c(fixed_r2_model_2,random_r2_model_2,conditional_r2_model_2), ncol=1)
R2_table_2 <- format(R2_table_2, digits=4)
rownames(R2_table_2) <- c("Fixed effect", "Random effect", "Total model")
colnames(R2_table_2) <- c("R-squared value")

knitr::kable(list(summary_mcmc_glmm(mcmc_model_2),R2_table_2), caption="Results from a MCMCglmm with pangenome fluidity as the response variable, whether a species' primary environment is host-associated or free-living as the fixed effect, and phylogeny as a random effect.", digits=4)%>% kable_styling(latex_options = "HOLD_position")



## ----primary env graph, fig.cap="Whether a species' primary environment is host-associated or free-living correlates with pangenome fluidity.", fig.width=4.5, fig.height=3----
ggplot(pangenome_lifestyles_primary_env_no_unknown, aes(x=Category_primary_env,y=pangenome_fluidity)) +
  geom_jitter(width=0.25, aes(colour=Category_primary_env)) +
  labs(y="Pangenome\nfluidity", x="Primary environment:\nHost- or Free-living") +
  stat_summary(aes(y = pangenome_fluidity,group=1), fun=mean, geom="crossbar",width=0.4,colour="black") +  
  scale_y_continuous(limits=c(0,0.48),breaks=c(0,0.1,0.2,0.3,0.4), labels=prettyZero, expand=c(0,0)) +
  theme_classic() +
  scale_colour_manual(values=c("#CC79A7", "#009E73")) +
  geom_signif(comparisons=list(c("Free", "Host")), annotation="*",
              y_position = 0.44, tip_length = 0, vjust=0.4) +
  theme(plot.title = element_text(size=15, hjust=0.5), axis.title = element_text(size=13),
        axis.text = element_text(size=12), axis.title.y=element_text(angle=360,vjust=0.5),
        legend.position="none")




## ----Host or free 4 lm, echo=FALSE-------------------------------------------------------------------------------------
pangenome_lifestyles_host_free_4 <- pangenome_lifestyles %>%
  mutate(Host_free_4 = case_when(Host_or_free=="Host" ~ "Host",
                                 Host_or_free=="Free" ~ "Free",
                                 Host_or_free=="Both" & Category_primary_env == "Free" ~ "Mostly_free",
                                 Host_or_free=="Both" & Category_primary_env == "Host" ~ "Mostly_host",
                                 Host_or_free=="Unknown"~"Unknown"))

# Reorder factor
pangenome_lifestyles_host_free_4$Host_free_4 <- factor(pangenome_lifestyles_host_free_4$Host_free_4, 
                                                       levels=c("Host","Mostly_host","Mostly_free", "Free","Unknown"))

pangenome_lifestyles_host_free_4_no_unknown <- pangenome_lifestyles_host_free_4 %>%
  filter(Host_or_free != "Unknown")

summary_lm_3 <- lm_summary_table(lm(pangenome_fluidity~Host_free_4, data=pangenome_lifestyles_host_free_4_no_unknown))

knitr::kable(summary_lm_3, 
             caption= "Results from a linear model with pangenome fluidity as the response variable and whether a species' is always host-associated, Mostly host-associated, mostly free-living, or always free-living as the explanatory variable.", digits=4)%>% kable_styling(latex_options = "HOLD_position")


## ----Host or free 4 MCMCglmm, echo=FALSE, cache=TRUE-------------------------------------------------------------------
mcmc_model_3 <- MCMCglmm(pangenome_fluidity ~ Host_free_4, random=~Species, data=pangenome_lifestyles_host_free_4_no_unknown,
                         nitt=50000,
                         prior=prior, ginverse = list(Species = Ainv),verbose = FALSE)

# R squared
mFixed <- mean(mcmc_model_3$Sol[,2]) * mcmc_model_3$X[, 2] +
  mean(mcmc_model_3$Sol[,3]) * mcmc_model_3$X[, 3] +
  mean(mcmc_model_3$Sol[,4]) * mcmc_model_3$X[, 4]

mVarF<- var(mFixed)

# Fixed effect variance / (fixed effect variance + random effect variance)
fixed_r2_model_3 <- mVarF/(mVarF+sum(apply(mcmc_model_3$VCV,2,mean)))

#random effect variance/(random effect variance + residual variance)
random_r2_model_3 <- (sum(apply(mcmc_model_3$VCV,2,mean)[-2]))/((sum(apply(mcmc_model_3$VCV,2,mean)))+mVarF)

#Conditionall R2 (total variance explained by the model)
conditional_r2_model_3 <- ((sum(apply(mcmc_model_3$VCV,2,mean)[-2]))+mVarF)/((sum(apply(mcmc_model_3$VCV,2,mean)))+mVarF)

R2_table_3 <- matrix(c(fixed_r2_model_3,random_r2_model_3,conditional_r2_model_3), ncol=1)
R2_table_3 <- format(R2_table_3, digits=4)
rownames(R2_table_3) <- c("Fixed effect", "Random effect", "Total model")
colnames(R2_table_3) <- c("R-squared value")

knitr::kable(list(summary_mcmc_glmm(mcmc_model_3),R2_table_3), caption="Results from a MCMCglmm with pangenome fluidity as the response variable, whether a species' is always host-associated, mostly host-associated, mostly free-living or always free-living as the fixed effect, and phylogeny as a random effect.", digits=4)%>% kable_styling(latex_options = "HOLD_position")


## ----Host free 4 graph, fig.cap="Whether a species' is host-associated, free-living, or a mixture correlates with pangenome fluidity.", fig.width=5.5, fig.height=3----
ggplot(pangenome_lifestyles_host_free_4_no_unknown, aes(x=Host_free_4, y=pangenome_fluidity)) +
  geom_jitter(width=0.25, aes(colour=Host_free_4)) +
  labs(y="Pangenome \nFluidity", x="Host- or free-living") +
  stat_summary(aes(y = pangenome_fluidity,group=1), fun=mean, geom="crossbar",width=0.4,colour="black") +
  scale_y_continuous(limits=c(0,0.50),breaks=c(0,0.1,0.2,0.3,0.4), labels=prettyZero, expand=c(0,0)) +
  theme_classic() +
  scale_x_discrete(labels=c("Host", "Mostly host", "Mostly free", "Free")) +
  scale_colour_manual(values=c("#CC79A7", "#569AE9", "#56C8E9", "#009E73")) +
  geom_signif(comparisons=list(c("Mostly_host", "Host")), annotation="***",y_position = 0.405, tip_length = 0, vjust=0.4) +
  geom_signif(comparisons=list(c("Mostly_free", "Host")), annotation="***",y_position = 0.43, tip_length = 0, vjust=0.4) +
  geom_signif(comparisons=list(c("Free", "Host")), annotation="**",y_position = 0.455, tip_length = 0, vjust=0.4) +
  theme(plot.title = element_text(size=15, hjust=0.5), axis.title = element_text(size=13),
        axis.text = element_text(size=12), axis.title.y=element_text(angle=360,vjust=0.5),
        legend.position = "none")



## ----Ob/fac vs pangenome fluidity lm, echo=FALSE-----------------------------------------------------------------------
### Obligate/Faculative
# (all host/both where Ob/Fac != Unknown)
# N= 115
pangenome_lifestyles_ob_fac <- pangenome_lifestyles %>%
  filter(Host_or_free != "Unknown") %>% filter(Host_or_free != "Free") %>%
  filter(Obligate_facultative != "Unknown")

pangenome_lifestyles_ob_fac <- as.data.frame(pangenome_lifestyles_ob_fac)

summary_lm_4 <- lm_summary_table(lm(pangenome_fluidity~Obligate_facultative, data=pangenome_lifestyles_ob_fac))

knitr::kable(summary_lm_4, 
             caption= "Results from a linear model with pangenome fluidity as the response variable and whether a species' is obligately or facultatively host-reliant as the explanatory variable.", digits=4)%>% kable_styling(latex_options = "HOLD_position")


## ----Ob/fac vs pangenome fluidity MCMCglmm, echo=FALSE, cache=TRUE-----------------------------------------------------
## MCMCglmm model
mcmc_model_4 <- MCMCglmm(pangenome_fluidity ~ Obligate_facultative, random=~Species, data=pangenome_lifestyles_ob_fac,
                         nitt=50000,
                         prior=prior, ginverse = list(Species = Ainv),verbose = FALSE)

# R squared
mFixed <- mean(mcmc_model_4$Sol[,2]) * mcmc_model_4$X[, 2]
mVarF<- var(mFixed)

# Fixed effect variance / (fixed effect variance + random effect variance)

fixed_r2_model_4 <- mVarF/(mVarF+sum(apply(mcmc_model_4$VCV,2,mean)))
#random effect variance/(random effect variance + residual variance)

random_r2_model_4 <- (sum(apply(mcmc_model_4$VCV,2,mean)[-2]))/((sum(apply(mcmc_model_4$VCV,2,mean)))+mVarF)
#Conditional R2 (total variance explained by the model)

conditional_r2_model_4 <- ((sum(apply(mcmc_model_4$VCV,2,mean)[-2]))+mVarF)/((sum(apply(mcmc_model_4$VCV,2,mean)))+mVarF)

R2_table_4 <- matrix(c(fixed_r2_model_4,random_r2_model_4,conditional_r2_model_4), ncol=1)
R2_table_4 <- format(R2_table_4, digits=4)
rownames(R2_table_4) <- c("Fixed effect", "Random effect", "Total model")
colnames(R2_table_4) <- c("R-squared value")

knitr::kable(list(summary_mcmc_glmm(mcmc_model_4),R2_table_4), caption="Results from a MCMCglmm with pangenome fluidity as the response variable, whether a species' is obligately or facultatively host-reliant as the fixed effect, and phylogeny as a random effect.", digits=4)%>% kable_styling(latex_options = "HOLD_position")



## ----Ob/fac vs pangenome fluidity figure, fig.cap="Obligately host-reliant species have a lower pangenome fluidity than facultatively host-reliant species", fig.width=4.5, fig.height=3----

ggplot(pangenome_lifestyles_ob_fac, aes(x=Obligate_facultative,y=pangenome_fluidity)) +
  geom_jitter(width=0.25, aes(colour=Obligate_facultative)) +
  labs(y="Pangenome\nfluidity", x="Host reliance") +
  stat_summary(aes(y = pangenome_fluidity,group=1), fun=mean, geom="crossbar",width=0.4,colour="black") + 
  geom_signif(comparisons=list(c("Obligate", "Facultative")), annotation="***",y_position = 0.42, tip_length = 0, vjust=0.4) +
  scale_y_continuous(limits=c(0,0.50),breaks=c(0,0.1,0.2,0.3,0.4), labels=prettyZero, expand=c(0,0)) +
  theme_classic() +
  scale_colour_manual(values=c("#CC79A7", "#009E73")) +
  theme(plot.title = element_text(size=15, hjust=0.5), axis.title = element_text(size=13),
        axis.text = element_text(size=12), axis.title.y=element_text(angle=360,vjust=0.5),
        legend.position = "none")



## ----extra/intra vs pangenome fluidity lm, echo=FALSE------------------------------------------------------------------
pangenome_lifestyles_intra_extra <- pangenome_lifestyles %>%
  filter(Host_or_free != "Unknown") %>% filter(Host_or_free != "Free") %>%
  filter(Intra_or_extracellular != "Unknown")

summary_lm_5 <- lm_summary_table(lm(pangenome_fluidity~Intra_or_extracellular, data=pangenome_lifestyles_intra_extra))

knitr::kable(summary_lm_5, 
             caption= "Results from a linear model with pangenome fluidity as the response variable and whether a species lives inside or outside host cells as the explanatory variable.", digits=4)%>% kable_styling(latex_options = "HOLD_position")


## ----extra/intra vs pangenome fluidity MCMCglmm, echo=FALSE, cache=TRUE------------------------------------------------
## MCMCglmm model
mcmc_model_5 <- MCMCglmm(pangenome_fluidity ~ Intra_or_extracellular, random=~Species, data=pangenome_lifestyles_intra_extra,
                         nitt=50000,
                         prior=prior, ginverse = list(Species = Ainv),verbose = FALSE)

# R squared
mFixed <- mean(mcmc_model_5$Sol[,2]) * mcmc_model_5$X[, 2] +
  mean(mcmc_model_5$Sol[,3]) * mcmc_model_5$X[, 3]

mVarF<- var(mFixed)

# Fixed effect variance / (fixed effect variance + random effect variance)
fixed_r2_model_5 <- mVarF/(mVarF+sum(apply(mcmc_model_5$VCV,2,mean)))

#random effect variance/(random effect variance + residual variance)
random_r2_model_5 <- (sum(apply(mcmc_model_5$VCV,2,mean)[-2]))/((sum(apply(mcmc_model_5$VCV,2,mean)))+mVarF)

#Conditional R2 (total variance explained by the model)
conditional_r2_model_5 <- ((sum(apply(mcmc_model_5$VCV,2,mean)[-2]))+mVarF)/((sum(apply(mcmc_model_5$VCV,2,mean)))+mVarF)

R2_table_5 <- matrix(c(fixed_r2_model_5,random_r2_model_5,conditional_r2_model_5), ncol=1)
R2_table_5 <- format(R2_table_5, digits=4)
rownames(R2_table_5) <- c("Fixed effect", "Random effect", "Total model")
colnames(R2_table_5) <- c("R-squared value")

knitr::kable(list(summary_mcmc_glmm(mcmc_model_5),R2_table_5), caption="Results from a MCMCglmm with pangenome fluidity as the response variable, whether a species lives inside or outside host cells as the fixed effect, and phylogeny as a random effect.", digits=4)%>% kable_styling(latex_options = "HOLD_position")


## ----extra/intra vs pangenome fluidity MCMCglmm - extracellular as intercept, echo=FALSE, cache=TRUE-------------------
# Reorder factor 
pangenome_lifestyles_intra_extra_reordered <- pangenome_lifestyles_intra_extra

pangenome_lifestyles_intra_extra_reordered$Intra_or_extracellular <- factor(pangenome_lifestyles_intra_extra_reordered$Intra_or_extracellular, levels=c("Intracellular","Both","Extracellular", "Unknown"))

## MCMCglmm model
mcmc_model_5b <- MCMCglmm(pangenome_fluidity ~ Intra_or_extracellular, random=~Species, data=pangenome_lifestyles_intra_extra_reordered,
                         nitt=50000,
                         prior=prior, ginverse = list(Species = Ainv),verbose = FALSE)

# R squared
mFixed <- mean(mcmc_model_5b$Sol[,2]) * mcmc_model_5b$X[, 2] +
  mean(mcmc_model_5b$Sol[,3]) * mcmc_model_5b$X[, 3]

mVarF<- var(mFixed)

# Fixed effect variance / (fixed effect variance + random effect variance)
fixed_r2_model_5b <- mVarF/(mVarF+sum(apply(mcmc_model_5b$VCV,2,mean)))

#random effect variance/(random effect variance + residual variance)
random_r2_model_5b <- (sum(apply(mcmc_model_5b$VCV,2,mean)[-2]))/((sum(apply(mcmc_model_5b$VCV,2,mean)))+mVarF)

#Conditional R2 (total variance explained by the model)
conditional_r2_model_5b <- ((sum(apply(mcmc_model_5b$VCV,2,mean)[-2]))+mVarF)/((sum(apply(mcmc_model_5b$VCV,2,mean)))+mVarF)

R2_table_5b <- matrix(c(fixed_r2_model_5b,random_r2_model_5b,conditional_r2_model_5b), ncol=1)
R2_table_5b <- format(R2_table_5, digits=4)
rownames(R2_table_5b) <- c("Fixed effect", "Random effect", "Total model")
colnames(R2_table_5b) <- c("R-squared value")

knitr::kable(list(summary_mcmc_glmm(mcmc_model_5b),R2_table_5b), caption="Results from a MCMCglmm with pangenome fluidity as the response variable, whether a species lives inside or outside host cells as the fixed effect, and phylogeny as a random effect; instead with extracellular as the intercept.", digits=4)%>% kable_styling(latex_options = "HOLD_position")


## ----extra/intra vs pangenome fluidity graph, fig.cap="Intracellular species have a lower pangenome fluidity than extracellular species", fig.width=4.5, fig.height=3----
ggplot(pangenome_lifestyles_intra_extra, aes(x=Intra_or_extracellular,y=pangenome_fluidity)) +
  geom_jitter(width=0.25, aes(colour=Intra_or_extracellular)) +
  labs(y="Pangenome \nfluidity", x="Host location") +
  stat_summary(aes(y = pangenome_fluidity,group=1), fun=mean, geom="crossbar",width=0.4,colour="black") +  
  geom_signif(comparisons=list(c("Extracellular", "Intracellular")), annotation="*",y_position = 0.43, tip_length = 0, vjust=0.4) +
  geom_signif(comparisons=list(c("Extracellular", "Both")), annotation="***",y_position = 0.41, tip_length = 0, vjust=0.4) +
  scale_y_continuous(limits=c(0,0.48),breaks=c(0,0.1,0.2,0.3,0.4), labels=prettyZero, expand=c(0,0)) +
  theme_classic() +
  scale_colour_manual(values=c("#CC79A7", "#56B4E9", "#009E73")) +
  theme(plot.title = element_text(size=15, hjust=0.5), axis.title = element_text(size=13),
        axis.text = element_text(size=12), axis.title.y=element_text(angle=360,vjust=0.5),
        legend.position = "none")



## ----eff_host vs pangenome fluidity lm---------------------------------------------------------------------------------
pangenome_lifestyles_effect_host <- pangenome_lifestyles %>%
  filter(Host_or_free != "Unknown") %>% filter(Host_or_free != "Free") %>%
  filter(Effect_on_host != "Unknown")

summary_lm_6 <- lm_summary_table(lm(genome_fluidity~Effect_on_host, data=pangenome_lifestyles_effect_host))

knitr::kable(summary_lm_6, 
             caption= "Results from a linear model with pangenome fluidity as the response variable and whether a species' effect on its host as the explanatory variable.")%>% kable_styling(latex_options = "HOLD_position")



## ----eff_host vs pangenome fluidity mcmcglmm, cache=TRUE---------------------------------------------------------------
mcmc_model_6 <- MCMCglmm(pangenome_fluidity ~ Effect_on_host, random=~Species, data=pangenome_lifestyles_effect_host,
                         nitt=50000,
                         prior=prior, ginverse = list(Species = Ainv),verbose = FALSE)

# R squared
mFixed <- mean(mcmc_model_6$Sol[,2]) * mcmc_model_6$X[, 2] +
  mean(mcmc_model_6$Sol[,3]) * mcmc_model_6$X[, 3]

mVarF<- var(mFixed)

# Fixed effect variance / (fixed effect variance + random effect variance)
fixed_r2_model_6 <- mVarF/(mVarF+sum(apply(mcmc_model_6$VCV,2,mean)))

#random effect variance/(random effect variance + residual variance)
random_r2_model_6 <- (sum(apply(mcmc_model_6$VCV,2,mean)[-2]))/((sum(apply(mcmc_model_6$VCV,2,mean)))+mVarF)

#Conditional R2 (total variance explained by the model)
conditional_r2_model_6 <- ((sum(apply(mcmc_model_6$VCV,2,mean)[-2]))+mVarF)/((sum(apply(mcmc_model_6$VCV,2,mean)))+mVarF)

R2_table_6 <- matrix(c(fixed_r2_model_6,random_r2_model_6,conditional_r2_model_6), ncol=1)
R2_table_6 <- format(R2_table_6, digits=4)
rownames(R2_table_6) <- c("Fixed effect", "Random effect", "Total model")
colnames(R2_table_6) <- c("R-squared value")

knitr::kable(list(summary_mcmc_glmm(mcmc_model_6),R2_table_6), caption="Results from a MCMCglmm with pangenome fluidity as the response variable, whether a species has a pathogenic or mutualistic effect on its host(s) as the fixed effect, and phylogeny as a random effect.", digits=4)%>% kable_styling(latex_options = "HOLD_position")



## ----effect on host vs pangenome fluidity graph, fig.cap="Non-motile species have a lower pangenome fluidity than motile species", fig.width=4.5, fig.height=3----

ggplot(pangenome_lifestyles_effect_host, aes(x=Effect_on_host,y=pangenome_fluidity)) +
  geom_jitter(width=0.25, aes(colour=Effect_on_host)) +
  labs(y="Pangenome\nfluidity", x="Effect on host") +
  stat_summary(aes(y = pangenome_fluidity,group=1), fun=mean, geom="crossbar",width=0.4,colour="black") +  
  geom_signif(comparisons=list(c("Mutualist", "Pathogen")), annotation="***",y_position = 0.43, tip_length = 0, vjust=0.4) +
  geom_signif(comparisons=list(c("Both", "Pathogen")), annotation="***",y_position = 0.40, tip_length = 0, vjust=0.4) +
  scale_y_continuous(limits=c(0,0.48),breaks=c(0,0.1,0.2,0.3,0.4), labels=prettyZero, expand=c(0,0)) +
  theme_classic() +
  scale_colour_manual(values=c("#CC79A7", "#56B4E9", "#009E73")) +
  coord_cartesian(xlim=c(1,3)) + #Remove unknown species & non-host
  theme(plot.title = element_text(size=15, hjust=0.5), axis.title = element_text(size=13),
        axis.text = element_text(size=12), axis.title.y=element_text(angle=360,vjust=0.5),
        legend.position = "none")


## ----motility vs pangenome fluidity lm---------------------------------------------------------------------------------
summary_lm_7 <- lm_summary_table(lm(genome_fluidity~Motility, data=pangenome_lifestyles))

knitr::kable(summary_lm_7, 
             caption= "Results from a linear model with pangenome fluidity as the response variable and a species' motility as the explanatory variable.", digits=4)%>% kable_styling(latex_options = "HOLD_position")


## ----motility vs pangenome fluidity mcmcglmm, cache=TRUE---------------------------------------------------------------
mcmc_model_7 <- MCMCglmm(pangenome_fluidity ~ Motility, random=~Species, data=pangenome_lifestyles, 
                         nitt=50000,
                         prior=prior, ginverse = list(Species = Ainv),verbose = FALSE)

# R squared
mFixed <- mean(mcmc_model_7$Sol[,2]) * mcmc_model_7$X[, 2] +
  mean(mcmc_model_7$Sol[,3]) * mcmc_model_7$X[, 3]

mVarF<- var(mFixed)

# Fixed effect variance / (fixed effect variance + random effect variance)
fixed_r2_model_7 <- mVarF/(mVarF+sum(apply(mcmc_model_7$VCV,2,mean)))

#random effect variance/(random effect variance + residual variance)
random_r2_model_7 <- (sum(apply(mcmc_model_7$VCV,2,mean)[-2]))/((sum(apply(mcmc_model_7$VCV,2,mean)))+mVarF)

#Conditional R2 (total variance explained by the model)
conditional_r2_model_7 <- ((sum(apply(mcmc_model_7$VCV,2,mean)[-2]))+mVarF)/((sum(apply(mcmc_model_7$VCV,2,mean)))+mVarF)

R2_table_7 <- matrix(c(fixed_r2_model_7,random_r2_model_7,conditional_r2_model_7), ncol=1)
R2_table_7 <- format(R2_table_7, digits=4)
rownames(R2_table_7) <- c("Fixed effect", "Random effect", "Total model")
colnames(R2_table_7) <- c("R-squared value")

knitr::kable(list(summary_mcmc_glmm(mcmc_model_7),R2_table_7), caption="Results from a MCMCglmm with pangenome fluidity as the response variable, a species' motility as the fixed effect, and phylogeny as a random effect.", digits=4)%>% kable_styling(latex_options = "HOLD_position")


## ----motility vs pangenome fluidity graph, fig.cap="Non-motile species have a lower pangenome fluidity than motile species", fig.width=4.5, fig.height=3----
ggplot(pangenome_lifestyles, aes(x=Motility,y=pangenome_fluidity)) +
  geom_jitter(width=0.25, aes(colour=Motility)) +
  labs(y="Pangenome\nfluidity", x="Motility") +
  stat_summary(aes(y = pangenome_fluidity,group=1), fun=mean, geom="crossbar",width=0.4,colour="black") +  
  geom_signif(comparisons=list(c("Non-motile", "Motile")), annotation="***",y_position = 0.42, tip_length = 0, vjust=0.4) +
  scale_y_continuous(limits=c(0,0.48),breaks=c(0,0.1,0.2,0.3,0.4), labels=prettyZero, expand=c(0,0)) +
  theme_classic() +
  scale_colour_manual(values=c("#CC79A7", "#56B4E9", "#009E73")) +
  theme(plot.title = element_text(size=15, hjust=0.5), axis.title = element_text(size=13),
        axis.text = element_text(size=12), axis.title.y=element_text(angle=360,vjust=0.5),
        legend.position = "none")


## ----four traits mcmcglmm, cache=TRUE----------------------------------------------------------------------------------
pangenome_lifestyles_no_unknown_nofree <- pangenome_lifestyles %>%
  filter(Host_or_free != "Unknown") %>% filter(Host_or_free != "Free") %>%
  filter(Obligate_facultative != "Unknown") %>% filter(Effect_on_host != "Unknown")

mcmc_model_8 <- MCMCglmm(pangenome_fluidity~Obligate_facultative+
                           Intra_or_extracellular+Effect_on_host+Motility,
                         random=~Species,
                         data=pangenome_lifestyles_no_unknown_nofree,
                         nitt=50000,
                         prior=prior, ginverse = list(Species = Ainv),verbose = FALSE)

mFixed <- mean(mcmc_model_8$Sol[,2]) * mcmc_model_8$X[, 2] + 
  mean(mcmc_model_8$Sol[,3]) * mcmc_model_8$X[, 3] + 
  mean(mcmc_model_8$Sol[,4]) * mcmc_model_8$X[, 4] + 
  mean(mcmc_model_8$Sol[,5]) * mcmc_model_8$X[, 5] +
  mean(mcmc_model_8$Sol[,6]) * mcmc_model_8$X[, 6] + 
  mean(mcmc_model_8$Sol[,7]) * mcmc_model_8$X[, 7] + 
  mean(mcmc_model_8$Sol[,8]) * mcmc_model_8$X[, 8] 

mVarF<- var(mFixed)

# Fixed effect variance / (fixed effect variance + random effect variance)
fixed_r2_model_8 <- mVarF/(mVarF+sum(apply(mcmc_model_8$VCV,2,mean)))

#random effect variance/(random effect variance + residual variance)
random_r2_model_8 <- (sum(apply(mcmc_model_8$VCV,2,mean)[-2]))/((sum(apply(mcmc_model_8$VCV,2,mean)))+mVarF)

#Conditionall R2 (total variance explained by the model)
conditional_r2_model_8 <-((sum(apply(mcmc_model_8$VCV,2,mean)[-2]))+mVarF)/((sum(apply(mcmc_model_8$VCV,2,mean)))+mVarF)

R2_table_8 <- matrix(c(fixed_r2_model_8,random_r2_model_8,conditional_r2_model_8), ncol=1)
R2_table_8 <- format(R2_table_8, digits=4)
rownames(R2_table_8) <- c("Fixed effect", "Random effect", "Total model")
colnames(R2_table_8) <- c("R-squared value")

knitr::kable(list(summary_mcmc_glmm(mcmc_model_8),R2_table_8), caption="Results from a MCMCglmm with pangenome fluidity as the response variable, four lifestyle traits as fixed effects, and phylogeny as a random effect.", digits=4)%>% kable_styling(latex_options = "HOLD_position")


## ----four traits plus host mcmcglmm, cache=TRUE------------------------------------------------------------------------
pangenome_lifestyles_no_unknown_nofree_host_free_4 <- pangenome_lifestyles_host_free_4_no_unknown %>%
  filter(Host_or_free != "Unknown") %>% #filter(Host_or_free != "Free") %>%
  filter(Obligate_facultative != "Unknown") %>% filter(Effect_on_host != "Unknown")

mcmc_model_8b <- MCMCglmm(pangenome_fluidity~Host_free_4+Obligate_facultative+
                            Intra_or_extracellular+Effect_on_host+Motility,
                         random=~Species,
                         data=pangenome_lifestyles_no_unknown_nofree_host_free_4,
                         nitt=50000,
                         prior=prior, ginverse = list(Species = Ainv),verbose = FALSE)

mFixed <- mean(mcmc_model_8b$Sol[,2]) * mcmc_model_8b$X[, 2] + 
  mean(mcmc_model_8b$Sol[,3]) * mcmc_model_8b$X[, 3] + 
  mean(mcmc_model_8b$Sol[,4]) * mcmc_model_8b$X[, 4] + 
  mean(mcmc_model_8b$Sol[,5]) * mcmc_model_8b$X[, 5] +
  mean(mcmc_model_8b$Sol[,6]) * mcmc_model_8b$X[, 6] + 
  mean(mcmc_model_8b$Sol[,7]) * mcmc_model_8b$X[, 7] + 
  mean(mcmc_model_8b$Sol[,8]) * mcmc_model_8b$X[, 8] +
  mean(mcmc_model_8b$Sol[,9]) * mcmc_model_8b$X[, 9] + 
  mean(mcmc_model_8b$Sol[,10]) * mcmc_model_8b$X[, 10] +
  mean(mcmc_model_8b$Sol[,11]) * mcmc_model_8b$X[, 11]

mVarF<- var(mFixed)

# Fixed effect variance / (fixed effect variance + random effect variance)
fixed_r2_model_8b <- mVarF/(mVarF+sum(apply(mcmc_model_8b$VCV,2,mean)))

#random effect variance/(random effect variance + residual variance)
random_r2_model_8b <- (sum(apply(mcmc_model_8b$VCV,2,mean)[-2]))/((sum(apply(mcmc_model_8b$VCV,2,mean)))+mVarF)

#Conditionall R2 (total variance explained by the model)
conditional_r2_model_8b <-((sum(apply(mcmc_model_8b$VCV,2,mean)[-2]))+mVarF)/((sum(apply(mcmc_model_8b$VCV,2,mean)))+mVarF)

R2_table_8b <- matrix(c(fixed_r2_model_8b,random_r2_model_8b,conditional_r2_model_8b), ncol=1)
R2_table_8b <- format(R2_table_8b, digits=4)
rownames(R2_table_8b) <- c("Fixed effect", "Random effect", "Total model")
colnames(R2_table_8b) <- c("R-squared value")

knitr::kable(list(summary_mcmc_glmm(mcmc_model_8b),R2_table_8b), caption="Results from a MCMCglmm with pangenome fluidity as the response variable, four lifestyle traits as fixed effects, and phylogeny as a random effect.", digits=4)%>% kable_styling(latex_options = "HOLD_position")


## ----phylo lms/glms with both------------------------------------------------------------------------------------------

### Now can run some phylolms()  & phyloglms() with the original data including both. 
pangenome_lifestyles_no_unknown_nofree_with_both <- pangenome_lifestyles_no_unknown_nofree %>%
  mutate(OF=case_when(Obligate_facultative=="Obligate"~0,
                                Obligate_facultative=="Facultative"~1), 
         IE=case_when(Intra_or_extracellular=="Intracellular" ~0, 
                      Intra_or_extracellular=="Both"~0.5, 
                      Intra_or_extracellular=="Extracellular"~1),
         EH=case_when(Effect_on_host=="Pathogen" ~0, 
                      Effect_on_host=="Both"~0.5, 
                      Effect_on_host=="Mutualist"~1),
         M=case_when(Motility=="Non-motile" ~0, 
                     Motility=="Both"~0.5, 
                     Motility=="Motile"~1))


# And drop tips not in the data frame
dataTree_phylo <-drop.tip(dataTree,dataTree$tip.label[-match(pangenome_lifestyles_no_unknown_nofree_with_both$Species, dataTree$tip.label)])

## Set row names
pangenome_lifestyles_no_unknown_nofree_with_both_row <- column_to_rownames(pangenome_lifestyles_no_unknown_nofree_with_both, var="Species")

# Phylo glm because OF is binary
# AIC     logLik Pen.logLik 
# 93.40     -43.70     -42.18 
# (Intercept) -0.017219  0.480559 -0.0358 0.97142  
# IE           1.421352  0.645662  2.2014 0.02771 *
phylolm_both_OFIE <- lm_summary_table(phyloglm(OF~IE, data=pangenome_lifestyles_no_unknown_nofree_with_both_row, phy=dataTree_phylo))

knitr::kable(phylolm_both_OFIE, caption="Phylogenetic logistic regression correlated evolution model between host reliance and host location", digits=4)%>% kable_styling(latex_options = "HOLD_position")

# AIC     logLik Pen.logLik 
# 92.73     -43.37     -41.88
# (Intercept)  0.33577 0.36507  0.9197 0.35771  
# EH           1.62662 0.67932  2.3945 0.01664 * 
phylolm_both_OFEH <- lm_summary_table(phyloglm(OF~EH, data=pangenome_lifestyles_no_unknown_nofree_with_both_row, phy=dataTree_phylo))

knitr::kable(phylolm_both_OFEH, caption="Results from a phylogenetic logistic regression correlated evolution model between host reliance and effect on host", digits=4)%>% kable_styling(latex_options = "HOLD_position")

# AIC     logLik Pen.logLik 
# 91.00     -42.50     -40.78  
# (Intercept)  0.33227 0.37764  0.8799 0.378935   
# M            1.43585 0.51318  2.7979 0.005143 **
phylolm_both_OFM <- lm_summary_table(phyloglm(OF~M, data=pangenome_lifestyles_no_unknown_nofree_with_both_row, phy=dataTree_phylo))

knitr::kable(phylolm_both_OFM, caption="Results from a phylogenetic logistic regression correlated evolution model between host reliance and motility", digits=4)%>% kable_styling(latex_options = "HOLD_position")

# AIC logLik 
# 62.17 -28.08
# (Intercept) 0.551646 1.022314  0.5396    0.5905    
# EH          0.258280 0.060616  4.2609 4.235e-05 ***
# R-squared: 0.1384	Adjusted R-squared: 0.1308
phylolm_both_IEEH <- lm_summary_table(phylolm(IE~EH, data=pangenome_lifestyles_no_unknown_nofree_with_both_row, phy=dataTree_phylo))

knitr::kable(phylolm_both_IEEH, caption="Results from a phylogenetic logistic regression correlated evolution model between host location and effect on host", digits=4)%>% kable_styling(latex_options = "HOLD_position")

# AIC logLik 
# 78.84 -36.42
# (Intercept) 0.618322 1.099153  0.5625  0.5749
# M           0.027209 0.040636  0.6696  0.5045
# R-squared: 0.003952	Adjusted R-squared: -0.004863
phylolm_both_IEM <- lm_summary_table(phylolm(IE~M, data=pangenome_lifestyles_no_unknown_nofree_with_both_row, phy=dataTree_phylo))

knitr::kable(phylolm_both_IEM, caption="Results from a phylogenetic logistic regression correlated evolution model between host location and motility", digits=4)%>% kable_styling(latex_options = "HOLD_position")

# AIC logLik 
# 162.15 -78.08
# (Intercept) 0.276337 1.578934  0.1750  0.8614
# M           0.060747 0.058373  1.0407  0.3003 
phylolm_both_EHM <- lm_summary_table(phylolm(EH~M, data=pangenome_lifestyles_no_unknown_nofree_with_both_row, phy=dataTree_phylo))

knitr::kable(phylolm_both_EHM, caption="Results from a phylogenetic logistic regression correlated evolution model between effect on host and motility", digits=4)%>% kable_styling(latex_options = "HOLD_position")



## ----pagel OFIE, cache=TRUE--------------------------------------------------------------------------------------------
# We'll need to simplify the traits into binary values.
# We've merged the 'both' category with one of the two catgeories.

pangenome_lifestyles_no_unknown_nofree <- pangenome_lifestyles %>%
  filter(Host_or_free != "Unknown") %>% filter(Host_or_free != "Free") %>%
  filter(Obligate_facultative != "Unknown") %>% filter(Effect_on_host != "Unknown")

## Merge 'both' with another category 
pangenome_lifestyles_no_unknown_nofree_binary <- pangenome_lifestyles_no_unknown_nofree %>%
  mutate(OF=case_when(Obligate_facultative=="Obligate"~0,
                                Obligate_facultative=="Facultative"~1), 
         IE=case_when(Intra_or_extracellular=="Intracellular" ~0, 
                      Intra_or_extracellular=="Both"~0, 
                      Intra_or_extracellular=="Extracellular"~1),
         EH=case_when(Effect_on_host=="Pathogen" ~0, 
                      Effect_on_host=="Both"~0, 
                      Effect_on_host=="Mutualist"~1),
         M=case_when(Motility=="Non-motile" ~0, 
                     Motility=="Both"~1, 
                     Motility=="Motile"~1))

# And drop tips not in the data frame
dataTree_pagel <-drop.tip(dataTree,dataTree$tip.label[-match(pangenome_lifestyles_no_unknown_nofree_binary$Species, dataTree$tip.label)])

# First set Species as the rownames of the dataset
pangenome_lifestyles_no_unknown_nofree_binary_row <- column_to_rownames(pangenome_lifestyles_no_unknown_nofree_binary, var="Species")

# We need to make 2 vectors, one for each of our traits.
OF_vector <- setNames(pangenome_lifestyles_no_unknown_nofree_binary_row$OF,rownames(pangenome_lifestyles_no_unknown_nofree_binary_row))
IE_vector <- setNames(pangenome_lifestyles_no_unknown_nofree_binary_row$IE,rownames(pangenome_lifestyles_no_unknown_nofree_binary_row))


# Run Pagel's model of correlated evolution
OF_IE_correlated <- fitPagel(dataTree_pagel, x=OF_vector, y=IE_vector)

# likelihood-ratio:  13.045 
# p-value:  0.011058 (*)
# independent      -95.78262 199.5652
# dependent        -89.26012 194.5202
OF_IE_correlated_table_hyp <- matrix(c(13.045,0.011058), ncol=1)
rownames(OF_IE_correlated_table_hyp) <- c("Likelihood ratio","p-value")
colnames(OF_IE_correlated_table_hyp) <- c("")

OF_IE_correlated_table_fit <- matrix(c(-95.78262,199.5652,-89.26012,194.5202), ncol=2)
rownames(OF_IE_correlated_table_fit) <- c("Independent","Dependent")
colnames(OF_IE_correlated_table_fit) <- c("Log-likelihood","AIC")

knitr::kable(list(OF_IE_correlated_table_hyp,OF_IE_correlated_table_fit), caption="Results from a Pagel's correlated evolution model between host reliance and host location", digits=4)%>% kable_styling(latex_options = "HOLD_position")




## ----pagel OFEH, cache=TRUE--------------------------------------------------------------------------------------------
EH_vector <- setNames(pangenome_lifestyles_no_unknown_nofree_binary_row$EH,rownames(pangenome_lifestyles_no_unknown_nofree_binary_row))

# Run pagel model
OF_EH_correlated <- fitPagel(dataTree_pagel, x=OF_vector, y=EH_vector)

# likelihood-ratio:  8.47535 
# p-value:  0.075638 (NS)
# independent      -85.99455 179.9891
# dependent        -81.75688 179.5138

OF_EH_correlated_table_hyp <- matrix(c(8.47535,0.075638), ncol=1)
rownames(OF_EH_correlated_table_hyp) <- c("Likelihood ratio","p-value")
colnames(OF_EH_correlated_table_hyp) <- c("")

OF_EH_correlated_table_fit <- matrix(c(-85.99455,179.9891,-81.75688,179.5138), ncol=2)
rownames(OF_EH_correlated_table_fit) <- c("Independent","Dependent")
colnames(OF_EH_correlated_table_fit) <- c("Log-likelihood","AIC")

knitr::kable(list(OF_EH_correlated_table_hyp, OF_EH_correlated_table_fit), caption="Results from a Pagel's correlated evolution model between host reliance and effect on host", digits=4) %>% kable_styling(latex_options = "HOLD_position")




## ----pagel OFM, cache=TRUE---------------------------------------------------------------------------------------------
# Make M vector
M_vector <- setNames(pangenome_lifestyles_no_unknown_nofree_binary_row$M,rownames(pangenome_lifestyles_no_unknown_nofree_binary_row))

# Run pagel model
OF_M_correlated <- fitPagel(dataTree_pagel, x=OF_vector, y=M_vector)

# likelihood-ratio:  21.7542 
# p-value:  0.000224304 (***)
# independent      -116.2164 240.4328
# dependent        -105.3393 226.6785

OF_M_correlated_table_hyp <- matrix(c(21.7542,0.000224304), ncol=1)
rownames(OF_M_correlated_table_hyp) <- c("Likelihood ratio","p-value")
colnames(OF_M_correlated_table_hyp) <- c("")

OF_M_correlated_table_fit <- matrix(c(-116.2164,240.4328,-105.3393,226.6785), ncol=2)
rownames(OF_M_correlated_table_fit) <- c("Independent","Dependent")
colnames(OF_M_correlated_table_fit) <- c("Log-likelihood","AIC")

knitr::kable(list(OF_M_correlated_table_hyp,OF_M_correlated_table_fit), caption="Results from a Pagel's correlated evolution model between host reliance and motility", digits=4) %>% kable_styling(latex_options = "HOLD_position")




## ----pagel IEEH, cache=TRUE--------------------------------------------------------------------------------------------
IE_EH_correlated <- fitPagel(dataTree_pagel, x=IE_vector, y=EH_vector)

# likelihood-ratio:  1.87169 
# p-value:  0.759344 (NS)
# independent      -79.08734 166.1747
# dependent        -78.15150 172.3030

IE_EH_correlated_table_hyp <- matrix(c(1.87169,0.759344), ncol=1)
rownames(IE_EH_correlated_table_hyp) <- c("Likelihood ratio","p-value")
colnames(IE_EH_correlated_table_hyp) <- c("")

IE_EH_correlated_table_fit <- matrix(c(-79.08734,166.1747,-78.15150,172.3030), ncol=2)
rownames(IE_EH_correlated_table_fit) <- c("Independent","Dependent")
colnames(IE_EH_correlated_table_fit) <- c("Log-likelihood","AIC")

knitr::kable(list(IE_EH_correlated_table_hyp,IE_EH_correlated_table_fit), caption="Results from a Pagel's correlated evolution model between host location and effect on host", digits=4)%>%kable_styling(latex_options = "HOLD_position")




## ----pagel IEM, cache=TRUE---------------------------------------------------------------------------------------------
IE_M_correlated <- fitPagel(dataTree_pagel, x=IE_vector, y=M_vector)

# likelihood-ratio:  4.32513 
# p-value:  0.363789 
# independent      -109.3092 226.6183
# dependent        -107.1466 230.2932

IE_M_correlated_table_hyp <- matrix(c(4.32513,0.363789), ncol=1)
rownames(IE_M_correlated_table_hyp) <- c("Likelihood ratio","p-value")
colnames(IE_M_correlated_table_hyp) <- c("")

IE_M_correlated_table_fit <- matrix(c(-109.3092,226.6183,-107.1466,230.2932), ncol=2)
rownames(IE_M_correlated_table_fit) <- c("Independent","Dependent")
colnames(IE_M_correlated_table_fit) <- c("Log-likelihood","AIC")

knitr::kable(list(IE_M_correlated_table_hyp, IE_M_correlated_table_fit), caption="Results from a Pagel's correlated evolution model between host location and motility", digits=4)%>% kable_styling(latex_options = "HOLD_position")




## ----pagel EHM, cache=TRUE---------------------------------------------------------------------------------------------
EH_M_correlated <- fitPagel(dataTree_pagel, x=EH_vector, y=M_vector)

# likelihood-ratio:  6.27651  
# p-value:  0.179428 (NS)
# independent      -99.52111 207.0422
# dependent        -96.38286 208.7657

EH_M_correlated_table_hyp <- matrix(c(6.27651,0.179428), ncol=1)
rownames(EH_M_correlated_table_hyp) <- c("Likelihood ratio","p-value")
colnames(EH_M_correlated_table_hyp) <- c("")

EH_M_correlated_table_fit <- matrix(c(-99.52111,207.0422,-96.38286,208.7657), ncol=2)
rownames(EH_M_correlated_table_fit) <- c("Independent","Dependent")
colnames(EH_M_correlated_table_fit) <- c("Log-likelihood","AIC")

knitr::kable(list(EH_M_correlated_table_hyp,EH_M_correlated_table_fit), caption="Results from a Pagel's correlated evolution model between effect on host and motility", digits=4)%>% kable_styling(latex_options = "HOLD_position")




## ----four traits binomial glms, cache=TRUE-----------------------------------------------------------------------------
# We can start to explore this with some simple binomial glms

# (Intercept)  0.04082    0.28577   0.143    0.886    
# IE           2.46061    0.54595   4.507 6.57e-06 ***
summary_glm_OFIE <- lm_summary_table(glm(OF ~ IE, data=pangenome_lifestyles_no_unknown_nofree_binary, family="binomial"))

knitr::kable(summary_glm_OFIE, 
             caption= "Results from a binomial generalised linear model with species' host reliance as the response variable and species' host location as the explanatory variable.", digits=4)%>% kable_styling(latex_options = "HOLD_position")

binomial_corr_1 <- ggplot(pangenome_lifestyles_no_unknown_nofree_binary, aes(x=as.character(IE),y=as.character(OF))) +
  geom_count() +
  labs(x="Host location", y="Host Reliance", subtitle = "p<0.001 (***)") +
  scale_size_area(breaks=c(10,20,30,40,50), name="Number\nof species") +
  scale_x_discrete(labels=c("Intracellular", "Extracellular")) +
  scale_y_discrete(labels=c("Obligate","Facultative")) +
  theme_classic() +
  theme(legend.position = "none")

# (Intercept)   0.8422     0.2261   3.726 0.000195 ***
# EH            2.2023     1.0481   2.101 0.035616 *  
summary_glm_OFEH <- lm_summary_table(glm(OF ~ EH, data=pangenome_lifestyles_no_unknown_nofree_binary, family="binomial"))

knitr::kable(summary_glm_OFEH, 
             caption= "Results from a binomial generalised linear model with species' host reliance as the response variable and species' effect on host as the explanatory variable.", digits=4)%>% kable_styling(latex_options = "HOLD_position")

binomial_corr_2 <- ggplot(pangenome_lifestyles_no_unknown_nofree_binary, aes(x=as.character(EH),y=as.character(OF))) +
  geom_count() +
  labs(x="Effect on host", y="Host Reliance", subtitle = "p<0.05(*)") +
  scale_size_area(breaks=c(10,20,30,40,50), name="Number\nof species") +
  scale_x_discrete(labels=c("Pathogen", "Mutualist")) +
  scale_y_discrete(labels=c("Obligate","Facultative")) +
  theme_classic() +
  theme(legend.position = "none")

# (Intercept)   0.3075     0.2635   1.167 0.243171    
# M             2.2575     0.5819   3.879 0.000105 ***
summary_glm_OFM <- lm_summary_table(glm(OF ~ M, data=pangenome_lifestyles_no_unknown_nofree_binary, family="binomial"))

knitr::kable(summary_glm_OFM, 
             caption= "Results from a binomial generalised linear model with species' host reliance as the response variable and species' motility as the explanatory variable.", digits=4)%>% kable_styling(latex_options = "HOLD_position")

binomial_corr_3 <- ggplot(pangenome_lifestyles_no_unknown_nofree_binary, aes(x=as.character(M),y=as.character(OF))) +
  geom_count() +
  labs(x="Motility", y="Host Reliance", subtitle = "p<0.001(***)") +
  scale_size_area(breaks=c(10,20,30,40,50), name="Number\nof species") +
  scale_x_discrete(labels=c("Non-motile", "Motile")) +
  scale_y_discrete(labels=c("Obligate","Facultative")) +
  theme_classic() +
  theme(legend.position = "none")

# (Intercept)   0.1076     0.2077   0.518   0.6043  
# EH            1.1161     0.5495   2.031   0.0422 *
summary_glm_IEEH <- lm_summary_table(glm(IE~EH, data=pangenome_lifestyles_no_unknown_nofree_binary, family="binomial"))

knitr::kable(summary_glm_OFM, 
             caption= "Results from a binomial generalised linear model with species' host location as the response variable and species' effect on host as the explanatory variable.", digits=4)%>% kable_styling(latex_options = "HOLD_position")

binomial_corr_4 <- ggplot(pangenome_lifestyles_no_unknown_nofree_binary, aes(x=as.character(EH),y=as.character(IE))) +
  geom_count() +
  labs(x="Effect on host", y="Host location", subtitle = "p<0.001(***)") +
  scale_size_area(breaks=c(10,20,30,40,50), name="Number\nof species") +
  scale_x_discrete(labels=c("Pathogen", "Mutualist")) +
  scale_y_discrete(labels=c("Intracellular", "Extracellular")) +
  theme_classic()

# (Intercept)   0.0339     0.2604   0.130    0.896
# M             0.5539     0.3816   1.452    0.147
summary_glm_IEM <- lm_summary_table(glm(IE~M, data=pangenome_lifestyles_no_unknown_nofree_binary, family="binomial"))

knitr::kable(summary_glm_IEM, 
             caption= "Results from a binomial generalised linear model with species' host location as the response variable and species' motility as the explanatory variable.", digits=4)%>% kable_styling(latex_options = "HOLD_position")

binomial_corr_5 <- ggplot(pangenome_lifestyles_no_unknown_nofree_binary, aes(x=as.character(M),y=as.character(IE))) +
  geom_count() +
  labs(x="Motility", y="Host location", subtitle = "p=0.15") +
  scale_size_area(breaks=c(10,20,30,40,50), name="Number\nof species") +
  scale_x_discrete(labels=c("Non-motile", "Motile")) +
  scale_y_discrete(labels=c("Intracellular", "Extracellular")) +
  theme_classic() +
  theme(legend.position = "none")

# (Intercept)  -1.1676     0.3060  -3.815 0.000136 ***
# M            -0.6242     0.4894  -1.275 0.202160  
summary_glm_EHM <- lm_summary_table(glm(EH ~ M, data=pangenome_lifestyles_no_unknown_nofree_binary, family="binomial"))

knitr::kable(summary_glm_EHM, 
             caption= "Results from a binomial generalised linear model with species' effect on host as the response variable and species' motility as the explanatory variable.", digits=4)%>% kable_styling(latex_options = "HOLD_position")

binomial_corr_6 <- ggplot(pangenome_lifestyles_no_unknown_nofree_binary, aes(x=as.character(M),y=as.character(EH))) +
  geom_count() +
  labs(x="Motility", y="Effect on host", subtitle = "p=0.20") +
  scale_size_area(breaks=c(10,20,30,40,50), name="Number\nof species") +
  scale_x_discrete(labels=c("Non-motile", "Motile")) +  
  scale_y_discrete(labels=c("Pathogen", "Mutualist")) +
  theme_classic() +
  theme(legend.position = "none")

binomial_corr <- binomial_corr_1 + binomial_corr_2 + binomial_corr_3 + binomial_corr_4 + 
  binomial_corr_5 +binomial_corr_6 + plot_layout(ncol=2)
## Although these simple models do not take phylogeny into account, we find that:
# Host reliance (OF) is significantly correlated with the three other traits.
# And Host location (IE) is significantly correlated with effect on host (EH).


## ----four traits binomial pair graphs, fig.cap="Associations between pairs of four lifestyle traits.The size of circles correspond to the number of species in one of four combinations of lifestyles for each trait.", fig.width=7, fig.height=8----

binomial_corr



## ----phylo logistic regressions----------------------------------------------------------------------------------------

# AIC     logLik Pen.logLik 
# 91.94     -42.97     -41.25
# (Intercept)  0.26214 0.38993  0.6723 0.50140  
# IE           1.21958 0.54388  2.2424 0.02494 *
phyloglm_binary_OFIE <- lm_summary_table(phyloglm(OF~IE, data=pangenome_lifestyles_no_unknown_nofree_binary_row, phy=dataTree_pagel))

knitr::kable(phyloglm_binary_OFIE, 
             caption= "Results from a phylogenetic logistic regression model with species' host reliance as the response variable and species' host location as the explanatory variable.", digits=4)%>% kable_styling(latex_options = "HOLD_position")

# AIC     logLik Pen.logLik 
# 99.60     -46.80     -45.42 
# (Intercept)  0.68551 0.34513  1.9862 0.04701 *
# EH           1.08066 0.75212  1.4368 0.15077   
phyloglm_binary_OFEH <- lm_summary_table(phyloglm(OF~EH, data=pangenome_lifestyles_no_unknown_nofree_binary_row, phy=dataTree_pagel))

knitr::kable(phyloglm_binary_OFEH, 
             caption= "Results from a phylogenetic logistic regression model with species' host reliance as the response variable and species' effect on host as the explanatory variable.", digits=4)%>% kable_styling(latex_options = "HOLD_position")

# AIC     logLik Pen.logLik 
# 91.95     -42.98     -41.14 
# (Intercept)  0.33954 0.37784  0.8986 0.368851   
# M            1.20562 0.46422  2.5971 0.009401 **
phyloglm_binary_OFM <- lm_summary_table(phylolm(OF~M, data=pangenome_lifestyles_no_unknown_nofree_binary_row, phy=dataTree_pagel))

knitr::kable(phyloglm_binary_OFM, 
             caption= "Results from a phylogenetic logistic regression model with species' host reliance as the response variable and species' motility as the explanatory variable.", digits=4)%>% kable_styling(latex_options = "HOLD_position")

# AIC     logLik Pen.logLik 
# 89.50     -41.75     -40.00 
# (Intercept) 0.064557 0.380413  0.1697  0.8652
# EH          0.388640 0.468461  0.8296  0.4068
phyloglm_binary_IEEH <- lm_summary_table(phyloglm(IE~EH, data=pangenome_lifestyles_no_unknown_nofree_binary_row, phy=dataTree_pagel))

knitr::kable(phyloglm_binary_IEEH, 
             caption= "Results from a phylogenetic logistic regression model with species' host location as the response variable and species' effect on host as the explanatory variable.", digits=4)%>% kable_styling(latex_options = "HOLD_position")

# AIC     logLik Pen.logLik 
# 93.88     -43.94     -41.48 
# (Intercept) 0.056662 0.378757  0.1496  0.8811
# M           0.088584 0.237688  0.3727  0.7094
phyloglm_binary_IEM <- lm_summary_table(phyloglm(IE~M, data=pangenome_lifestyles_no_unknown_nofree_binary_row, phy=dataTree_pagel))

knitr::kable(phyloglm_binary_IEM, 
             caption= "Results from a phylogenetic logistic regression model with species' host location as the response variable and species' motility as the explanatory variable.", digits=4)%>% kable_styling(latex_options = "HOLD_position")

# AIC     logLik Pen.logLik 
# 72.22     -33.11     -31.42 
# (Intercept) -1.50058  0.44771 -3.3517 0.0008031 ***
# M           -0.25761  0.45060 -0.5717 0.5675141    
phyloglm_binary_EHM <- lm_summary_table(phyloglm(EH~M, data=pangenome_lifestyles_no_unknown_nofree_binary_row, phy=dataTree_pagel))

knitr::kable(phyloglm_binary_EHM, 
             caption= "Results from a phylogenetic logistic regression model with species' effect on host as the response variable and species' motility as the explanatory variable.", digits=4)%>% kable_styling(latex_options = "HOLD_position")



## ----phylopath set-up--------------------------------------------------------------------------------------------------

# Rename variables to make it easier to visualise models
pangenome_lifestyles_no_unknown_nofree<- pangenome_lifestyles_no_unknown_nofree %>%
  rename(PF = pangenome_fluidity)

#Set row names as species
rownames(pangenome_lifestyles_no_unknown_nofree) <- pangenome_lifestyles_no_unknown_nofree$Species

## Make into data frame
pangenome_lifestyles_no_unknown_nofree <- as.data.frame(pangenome_lifestyles_no_unknown_nofree)

#### Coding categorical variables ----
## Model will not work with categorical - need to code in as if discrete continuous

## Code as binary and/or discrete traits
pangenome_lifestyles_no_unknown_nofree_1 <- pangenome_lifestyles_no_unknown_nofree %>%
  mutate(OF=as.factor(case_when(Obligate_facultative=="Obligate"~0,
                                Obligate_facultative=="Facultative"~1)), 
         IE=case_when(Intra_or_extracellular=="Intracellular" ~0, 
                      Intra_or_extracellular=="Both"~0.5, 
                      Intra_or_extracellular=="Extracellular"~1),
         EH=case_when(Effect_on_host=="Pathogen" ~0, 
                      Effect_on_host=="Both"~0.5, 
                      Effect_on_host=="Mutualist"~1),
         M=case_when(Motility=="Non-motile" ~0, 
                     Motility=="Both"~0.5, 
                     Motility=="Motile"~1))



## ----phylopath simple models-------------------------------------------------------------------------------------------
# Run set of models which has each of the 4 causing PF, with all combinations of deleted

models_simple <- define_model_set(
  a = c(PF~OF+IE+EH+M),
  b = c(PF~IE+EH+M),
  c = c(PF~OF+EH+M),
  d = c(PF~OF+IE+M),
  e = c(PF~OF+IE+EH),
  f = c(PF~OF+IE),
  g = c(PF~OF+EH),
  h = c(PF~OF+M),
  i = c(PF~IE+EH),
  j = c(PF~IE+M),
  k = c(PF~EH+M),
  l = c(PF~OF),
  m = c(PF~IE),
  n = c(PF~EH),
  o = c(PF~M)
)


## ----phylopath simple models plot, fig.cap="Set of simple models, varying by which lifestyle traits cause pangenome fluidity", fig.width=8, fig.height=12----
plot_model_set(models_simple, labels = c(PF="Pangenome\nfluidity",
                OF="Host\nreliance", IE="Host\nlocation", 
                EH="Effect\non host", M="Motility"), text_size=3) # Plot all models



## ----phylopath simple models results-----------------------------------------------------------------------------------
## Run path analysis
result_1 <- phylo_path(models_simple, data=pangenome_lifestyles_no_unknown_nofree_1, tree=dataTree, 
                       model="lambda", method="logistic_MPLE")

## Path analysis results
s1 <- summary(result_1)



## ----phylopath complex models------------------------------------------------------------------------------------------
models_complex <- define_model_set(
  a = c(),
  b = c(IE~OF),
  c = c(IE~OF, M~OF),
  d = c(IE~OF, M~OF, EH~OF),
  e = c(IE~OF, M~OF, EH~OF, EH~M),
  f = c(IE~OF, M~OF+EH, EH~OF),
  g = c(IE~OF, M~OF+EH, OF~EH),
  h = c(IE~OF, M~EH, OF~EH+M),
  i = c(M~EH, OF~EH+M+IE),
  j = c(EH~M, OF~EH+M+IE),
  k = c(M~EH+OF, OF~EH+IE), 
  l = c(EH~OF+M, OF~M+IE),
  .common = c(PF~OF+IE+EH+M)
)


## ----phylopath complex models plot, fig.cap="Set of more complex models, varying by which lifestyle traits cause pangenome fluidity and how they might cause each other.", fig.width=8, fig.height=12----
plot_model_set(models_complex, labels = c(PF="Pangenome\nfluidity",
                OF="Host\nreliance", IE="Host\nlocation", 
                EH="Effect\non host", M="Motility"), text_size=2) # Plot all models


## ----phylopath complex models results----------------------------------------------------------------------------------
## Run path analysis
result_2 <- phylo_path(models_complex, data=pangenome_lifestyles_no_unknown_nofree_1, tree=dataTree, 
                       model="lambda", method="logistic_MPLE")

## Path analysis results
s2 <- summary(result_2)
s2_subset <- s2 %>% select(w,p,CICc)

knitr::kable(s2_subset, caption="Comparison of set of more complex causal models, varying by which lifestyle trait(s) directly influence pangenome fluidity and how they might cause each other.", digits = 4)%>% kable_styling(latex_options = "HOLD_position")


## ----phylopath complex models results plot, fig.cap="Comparing model support for set of more complex models of causation between four lifestyle traits and pangenome fluidity. Models are ordered by their value of `w`, a measure of model support, and numbers on the bars correspond to overall p-values of the model.", fig.width=8, fig.height=8----
plot(s2)


## ----phylopath complex models i and j, fig.cap="The two models of causation with highest support, model i (first panel) and model j (second panel)"----
i <- plot(models_complex$i) #plot single model 
j <- plot(models_complex$j) #plot single model 

i_j <- i + j

i_j


## ----phylopath complex models best model, fig.cap="Average best model of causation between four lifestyle traits and pangenome fluidiy"----
# Two methods for averaging models: "conditional" and "full".
# The methods differ in how they deal with averaging a path coefficient where the path is absent in some of the models. 
# The full method sets the coefficient (and the variance) for the missing paths to zero, 
# meaning paths that are missing in some models will shrink towards zero. 
# The conditional method only averages over models where the path appears, making it more sensitive to small effects.

# We will use conditional, since i & j only differ in the direction of the EH<->M connection.
# Any individual model cannot be circular, so only one direction can be defined at once.
# Therefore, the coefficients of each are relevant when comparing to the other models,
# And it would not make sense to set each to 0 in one of the 2 models.
average_model_2_conditional <- average(result_2, avg_method = "conditional") #default

## Plot average best models
plot(average_model_2_conditional, curvature = 0.1, algorithm="dh", type="width", 
     labels = c(IE="Host\nlocation", OF="Host\nreliance",
                PF="Pangenome\nfluidity", EH="Effect\non host", M="Motility"),
     box_x=25, box_y=15, text_size=3) 


## ----phylopath complex models coeffs, fig.cap="Path coeeficients for the average best model, with 95% confidence intervals."----
coef_plot(average_model_2_conditional, reverse_order = TRUE) + 
  ggplot2::coord_flip() + 
  ggplot2::theme_bw()


## ----phylopath model deletion------------------------------------------------------------------------------------------

## Now, to arrive at final best model, want to do deletion to the minimal model.
# This asks if removing any one path improves the model fit.
# Do deletion separately for both models i and j
models_deletion_i <- define_model_set(
  i = c(PF~OF+IE+EH+M, M~EH, OF~EH+M+IE), # original 
  a = c(PF~IE+EH+M, M~EH, OF~EH+M+IE),
  b = c(PF~OF+EH+M, M~EH, OF~EH+M+IE),
  c = c(PF~OF+IE+M, M~EH, OF~EH+M+IE),
  d = c(PF~OF+IE+EH, M~EH, OF~EH+M+IE),
  e = c(PF~OF+IE+EH+M, OF~EH+M+IE),
  f = c(PF~OF+IE+EH+M, M~EH, OF~M+IE),
  g = c(PF~OF+IE+EH+M, M~EH, OF~EH+IE),
  h = c(PF~OF+IE+EH+M, M~EH, OF~EH+M)
)

models_deletion_j <- define_model_set(
  j = c(PF~OF+IE+EH+M, EH~M, OF~EH+M+IE), # original 
  k = c(PF~IE+EH+M, EH~M, OF~EH+M+IE),
  l = c(PF~OF+EH+M, EH~M, OF~EH+M+IE),
  m = c(PF~OF+IE+M, EH~M, OF~EH+M+IE),
  n = c(PF~OF+IE+EH, EH~M, OF~EH+M+IE),
  o = c(PF~OF+IE+EH+M, OF~EH+M+IE),
  p = c(PF~OF+IE+EH+M, EH~M, OF~M+IE),
  q = c(PF~OF+IE+EH+M, EH~M, OF~EH+IE),
  r = c(PF~OF+IE+EH+M, EH~M, OF~EH+M)
)


## Run path analysis
result_3 <- phylo_path(models_deletion_i, data=pangenome_lifestyles_no_unknown_nofree_1, tree=dataTree, 
                       model="lambda", method="logistic_MPLE")

result_4 <- phylo_path(models_deletion_j, data=pangenome_lifestyles_no_unknown_nofree_1, tree=dataTree, 
                       model="lambda", method="logistic_MPLE")


## Path analysis results
s3 <- summary(result_3)

s4 <- summary(result_4)


## ----phylopath model deletion i results plot, fig.cap="Models each with a path deleted from original (model i) and comparison of support."----
i_deletion_models <- plot_model_set(models_deletion_i,labels = c(PF="Pangenome\nfluidity",
                OF="Host\nreliance", IE="Host\nlocation", 
                EH="Effect\non host", M="Motility"), text_size=3)

i_deletion_comparison <- plot(s3)

i_deletion <- i_deletion_models + i_deletion_comparison + plot_layout(ncol=1)



## ----phylopath model deletion j results plot, fig.cap="Models each with a path deleted from original (model j) and comparison of support."----
j_deletion_models <- plot_model_set(models_deletion_j, labels = c(PF="Pangenome\nfluidity",
                OF="Host\nreliance", IE="Host\nlocation", 
                EH="Effect\non host", M="Motility"), text_size=3)

j_deletion_comparison <- plot(s4)

j_deletion <- j_deletion_models + j_deletion_comparison + plot_layout(ncol=1)



## ----phylopath model deletion i and j together-------------------------------------------------------------------------
## Also combine them together - see whether i and j still the best
models_deletion_ij <- define_model_set(
  i = c(PF~OF+IE+EH+M, M~EH, OF~EH+M+IE), # original 
  a = c(PF~IE+EH+M, M~EH, OF~EH+M+IE),
  b = c(PF~OF+EH+M, M~EH, OF~EH+M+IE),
  c = c(PF~OF+IE+M, M~EH, OF~EH+M+IE),
  d = c(PF~OF+IE+EH, M~EH, OF~EH+M+IE),
  e = c(PF~OF+IE+EH+M, OF~EH+M+IE),
  f = c(PF~OF+IE+EH+M, M~EH, OF~M+IE),
  g = c(PF~OF+IE+EH+M, M~EH, OF~EH+IE),
  h = c(PF~OF+IE+EH+M, M~EH, OF~EH+M),
  j = c(PF~OF+IE+EH+M, EH~M, OF~EH+M+IE), # original 
  k = c(PF~IE+EH+M, EH~M, OF~EH+M+IE),
  l = c(PF~OF+EH+M, EH~M, OF~EH+M+IE),
  m = c(PF~OF+IE+M, EH~M, OF~EH+M+IE),
  n = c(PF~OF+IE+EH, EH~M, OF~EH+M+IE),
  o = c(PF~OF+IE+EH+M, OF~EH+M+IE),
  p = c(PF~OF+IE+EH+M, EH~M, OF~M+IE),
  q = c(PF~OF+IE+EH+M, EH~M, OF~EH+IE),
  r = c(PF~OF+IE+EH+M, EH~M, OF~EH+M)
)


# Run path analysis
result_5 <- phylo_path(models_deletion_ij, data=pangenome_lifestyles_no_unknown_nofree_1, tree=dataTree, 
                       model="lambda", method="logistic_MPLE")

## Path analysis results
s5 <- summary(result_5)


## ----phylopath model deletion i and j together plot, fig.cap="Comparison of support for models i and j, and versions each with one path deleted."----
plot(s5)


## ----phylopath model deletion i and j together average conditional, fig.cap="Average best model, combining models i and j, along with 5 models with one oath deleted."----
## Plot average best model
average_model_5_conditional <- average(result_5, avg_method = "conditional") #default

## Plot average best models
plot(average_model_5_conditional, curvature = 0.1, algorithm="dh", type="width", 
     labels = c(IE="Host\nlocation", OF="Host\nreliance",
                PF="Pangenome\nfluidity", EH="Effect\non host", M="Motility"),
     box_x=25, box_y=15, text_size=3) 



## ----additional lifestyle path binary - set up-------------------------------------------------------------------------
# Code single lifestyle variable
pangenome_lifestyles_no_unknown_nofree_2 <- pangenome_lifestyles_no_unknown_nofree %>%
  mutate(OF=case_when(Obligate_facultative=="Obligate"~0,
                                Obligate_facultative=="Facultative"~1), 
         IE=case_when(Intra_or_extracellular=="Intracellular" ~0, 
                      Intra_or_extracellular=="Both"~0, 
                      Intra_or_extracellular=="Extracellular"~1),
         EH=case_when(Effect_on_host=="Pathogen" ~0, 
                      Effect_on_host=="Both"~0, 
                      Effect_on_host=="Mutualist"~1),
         M=case_when(Motility=="Non-motile" ~0, 
                     Motility=="Both"~1, 
                     Motility=="Motile"~1))


## ----additional lifestyle path binary simple models results------------------------------------------------------------
## Run path analysis
result_1b <- phylo_path(models_simple, data=pangenome_lifestyles_no_unknown_nofree_2, tree=dataTree, 
                       model="lambda", method="logistic_MPLE")

## Path analysis results
s1b <- summary(result_1b)



## ----additional lifestyle path binary complex models-------------------------------------------------------------------
models_complex <- define_model_set(
  a = c(),
  b = c(IE~OF),
  c = c(IE~OF, M~OF),
  d = c(IE~OF, M~OF, EH~OF),
  e = c(IE~OF, M~OF, EH~OF, EH~M),
  f = c(IE~OF, M~OF+EH, EH~OF),
  g = c(IE~OF, M~OF+EH, OF~EH),
  h = c(IE~OF, M~EH, OF~EH+M),
  i = c(M~EH, OF~EH+M+IE),
  j = c(EH~M, OF~EH+M+IE),
  k = c(M~EH+OF, OF~EH+IE), 
  l = c(EH~OF+M, OF~M+IE),
  .common = c(PF~OF+IE+EH+M)
)


## ----additional lifestyle path binary complex models results-----------------------------------------------------------
## Run path analysis
result_2b <- phylo_path(models_complex, data=pangenome_lifestyles_no_unknown_nofree_2, tree=dataTree, 
                       model="lambda", method="logistic_MPLE")

## Path analysis results
s2b <- summary(result_2b)
s2b_subset <- s2b %>% select(w,p,CICc)

knitr::kable(s2b_subset, caption="Comparison of set of more complex causal models, varying by which lifestyle trait(s) directly influence pangenome fluidity and how they might cause each other.", digits=4)%>% kable_styling(latex_options = "HOLD_position")


## ----additional lifestyle path binary complex models results plot, fig.cap="Comparing model support for set of more complex models of causation between four lifestyle traits and pangenome fluidity. Models are ordered by their value of `w`, a measure of model support, and numbers on the bars correspond to overall p-values of the model.", fig.width=8, fig.height=8----
plot(s2b)


## ----additional lifestyle path binary complex models best model, fig.cap="Average best model of causation between four lifestyle traits and pangenome fluidity."----
# Two methods for averaging models: "conditional" and "full".
# The methods differ in how they deal with averaging a path coefficient where the path is absent in some of the models. 
# The full method sets the coefficient (and the variance) for the missing paths to zero, 
# meaning paths that are missing in some models will shrink towards zero. 
# The conditional method only averages over models where the path appears, making it more sensitive to small effects.

# We will use conditional, since i & j only differ in the direction of the EH<->M connection.
# Any individual model cannot be circular, so only one direction can be defined at once.
# Therefore, the coefficients of each are relevant when comparing to the other models,
# And it would not make sense to set each to 0 in one of the 2 models.
average_model_2b <- best(result_2b) #default

## Plot average best models
plot(average_model_2b, curvature = 0.1, algorithm="dh", type="width", 
     labels = c(IE="Host\nlocation", OF="Host\nreliance",
                PF="Pangenome\nfluidity", EH="Effect\non host", M="Motility"),
     box_x=25, box_y=15, text_size=3) 


## ----additional lifestyle path binary model deletion-------------------------------------------------------------------

## Now, to arrive at final best model, want to do deletion to the minimal model.
# This asks if removing any one path improves the model fit.
# Do deletion separately for both models i and j
models_deletion_d <- define_model_set(
  d = c(PF~OF+IE+EH+M, IE~OF, M~OF, EH~OF), # original 
  a = c(PF~IE+EH+M, IE~OF, M~OF, EH~OF),
  b = c(PF~OF+EH+M, IE~OF, M~OF, EH~OF),
  c = c(PF~OF+IE+M, IE~OF, M~OF, EH~OF),
  e = c(PF~OF+IE+EH, IE~OF, M~OF, EH~OF),
  f = c(PF~OF+IE+EH+M, M~OF, EH~OF),
  g = c(PF~OF+IE+EH+M, IE~OF, EH~OF),
  h = c(PF~OF+IE+EH+M, IE~OF, M~OF)
)


## Run path analysis
result_3b <- phylo_path(models_deletion_d, data=pangenome_lifestyles_no_unknown_nofree_1, tree=dataTree, 
                       model="lambda", method="logistic_MPLE")

## Path analysis results
s3b <- summary(result_3b)



## ----additional lifestyle path binary model deletion d models, fig.cap="Models each with a path deleted from original (model d).", fig.height=12----
plot_model_set(models_deletion_d, labels = c(PF="Pangenome\nfluidity",
                OF="Host\nreliance", IE="Host\nlocation", 
                EH="Effect\non host", M="Motility"), text_size=3)



## ----additional lifestyle path binary complex models results table-----------------------------------------------------

## Path analysis results
s3b_subset <- s3b %>% select(w,p,CICc)

knitr::kable(s3b_subset, caption="Comparison of support for set of models each with a path deleted from original (model d).", digits=4)%>% kable_styling(latex_options = "HOLD_position")


## ----additional lifestyle path binary model deletion d results plot, fig.cap="Comparison of support for set of models each with a path deleted from original (model d)."----
plot(s3b)



## ----additional lifestyle path binary model deletion d average best model, fig.cap="Average of three best models following model deletion, for analysis with all lifestyle traits coded as binary variables."----
## Plot average best model
average_model_3b <- average(result_3b, avg_method = "conditional") #default

## Plot average best models
plot(average_model_3b, curvature = 0.1, algorithm="dh", type="width", 
     labels = c(IE="Host\nlocation", OF="Host\nreliance",
                PF="Pangenome\nfluidity", EH="Effect\non host", M="Motility"),
     box_x=25, box_y=15, text_size=3) 



## ----Genome size vs pangenome fluidity MCMCglmm, echo=FALSE, cache=TRUE------------------------------------------------
mcmc_model_9 <- MCMCglmm(pangenome_fluidity ~ genome_size, random=~Species, data=pangenome_lifestyles,
                          nitt=50000,
                          prior=prior, ginverse = list(Species = Ainv),verbose = FALSE)

# R2 values
mFixed <- mean(mcmc_model_9$Sol[,2]) * mcmc_model_9$X[, 2]
mVarF<- var(mFixed)

# Fixed effect variance / (fixed effect variance + random effect variance)
fixed_r2_model_9 <- mVarF/(mVarF+sum(apply(mcmc_model_9$VCV,2,mean)))
#random effect variance/(random effect variance + residual variance)
random_r2_model_9 <- (sum(apply(mcmc_model_9$VCV,2,mean)[-2]))/((sum(apply(mcmc_model_9$VCV,2,mean)))+mVarF)
#Conditionall R2 (total variance explained by the model)
conditional_r2_model_9 <- ((sum(apply(mcmc_model_9$VCV,2,mean)[-2]))+mVarF)/((sum(apply(mcmc_model_9$VCV,2,mean)))+mVarF)

R2_table_9 <- matrix(c(fixed_r2_model_9,random_r2_model_9,conditional_r2_model_9), ncol=1)
R2_table_9 <- format(R2_table_9, digits=4)
rownames(R2_table_9) <- c("Fixed effect", "Random effect", "Total model")
colnames(R2_table_9) <- c("R-squared value")

knitr::kable(list(summary_mcmc_glmm(mcmc_model_9),R2_table_9), caption="Results from a MCMCglmm with pangenome fluidity as the response variable, genome size as the fixed effect, and phylogeny as a random effect.", digits=4)%>% kable_styling(latex_options = "HOLD_position")


## ----Genome size vs pangenome fluidity plot, fig.cap="Plot showing correlation between a species' average number of genes in their genomes and their pangenome fluidity. The line is the slope and intercept from the MCMCglmm analysis in the above table. N=126"----
ggplot(pangenome_lifestyles, aes(x=genome_size,y=pangenome_fluidity)) +
  geom_point(colour="black") +
  labs(y="Pangenome\nfluidity", x="Genome size \n(number of genes)") +
  geom_abline(intercept=0.1669, slope=0.00001834, colour="#0072B2", size=1) + #MCMCglmm, n=126
  scale_y_continuous(limits=c(0,0.48),breaks=c(0,0.1,0.2,0.3,0.4), labels=prettyZero, expand=c(0,0)) +
  scale_x_continuous(limits=c(0,7500),breaks=c(0,2000,4000,6000), expand=c(0,0)) +
  theme_classic() +
  theme(plot.title = element_text(size=15, hjust=0.5), axis.title = element_text(size=13),
        axis.text = element_text(size=12), axis.title.y=element_text(angle=360,vjust=0.5),
        axis.text.y = element_blank(),
        legend.position = "none")


## ----effective pop size vs pangenome fluidity MCMCglmm, echo=FALSE, cache=TRUE-----------------------------------------
## Make Ne on scale easier for stats programs - can be complain if numbers too large
# Divide Ne by 1,000,000 to get value per 1 million.
pangenome_lifestyles$Ne_small <- pangenome_lifestyles$Ne/1000000

# Data - as above, with 115 species where host or both, and know all 4 lifestyle traits
# AND have NE data
# N=75 species
pangenome_lifestyles_no_unknown_ne <- pangenome_lifestyles %>%
  filter(Ne_small != "NA")

mcmc_model_10 <- MCMCglmm(pangenome_fluidity ~ Ne_small, random=~Species,
                          data=pangenome_lifestyles_no_unknown_ne,
                          nitt=50000,
                          prior=prior, ginverse = list(Species = Ainv),verbose = FALSE)

# R2 values
mFixed <- mean(mcmc_model_10$Sol[,2]) * mcmc_model_10$X[, 2]
mVarF<- var(mFixed)

# Fixed effect variance / (fixed effect variance + random effect variance)
fixed_r2_model_10 <- mVarF/(mVarF+sum(apply(mcmc_model_10$VCV,2,mean)))
#random effect variance/(random effect variance + residual variance)
random_r2_model_10 <- (sum(apply(mcmc_model_10$VCV,2,mean)[-2]))/((sum(apply(mcmc_model_10$VCV,2,mean)))+mVarF)
#Conditionall R2 (total variance explained by the model)
conditional_r2_model_10 <- ((sum(apply(mcmc_model_10$VCV,2,mean)[-2]))+mVarF)/((sum(apply(mcmc_model_10$VCV,2,mean)))+mVarF)

R2_table_10 <- matrix(c(fixed_r2_model_10,random_r2_model_10,conditional_r2_model_10), ncol=1)
R2_table_10 <- format(R2_table_10, digits=4)
rownames(R2_table_10) <- c("Fixed effect", "Random effect", "Total model")
colnames(R2_table_10) <- c("R-squared value")

knitr::kable(list(summary_mcmc_glmm(mcmc_model_10),R2_table_10), caption="Results from a MCMCglmm with pangenome fluidity as the response variable, effective population size as the fixed effect, and phylogeny as a random effect.", digits=4)%>% kable_styling(latex_options = "HOLD_position")



## ----effective pop size vs pangenome fluidity plot, fig.cap="Plot showing correlation between a species' effective population size and their pangenome fluidity. The line is the slope and intercept of the MCMCglmm analysis in the above table. N=77."----
ggplot(pangenome_lifestyles_no_unknown_ne, aes(x=Ne_small,y=pangenome_fluidity)) +
  geom_point() +
  labs(y="Pangenome\nfluidity", x="Effective population size \n(millions)") +
  geom_abline(intercept=0.0983450, slope=0.0003523, colour="#0072B2", size=1) + #MCMCglmm, n=126
  scale_x_continuous(breaks=c(0, 100, 200, 300, 400), limits=c(0,430)) +
  scale_y_continuous(limits=c(0,0.48),breaks=c(0,0.1,0.2,0.3,0.4), labels=prettyZero, expand=c(0,0)) +
  theme_classic() +
  theme(plot.title = element_text(size=15, hjust=0.5), axis.title = element_text(size=13),
        axis.text = element_text(size=12), axis.title.y=element_text(angle=360,vjust=0.5),
        axis.text.y = element_blank(),
        legend.position = "none")


## ----NE/GS/LS path analysis set up-------------------------------------------------------------------------------------
pangenome_lifestyles_no_unknown_nofree_ne <- pangenome_lifestyles %>%
  filter(Host_or_free != "Unknown") %>% filter(Host_or_free != "Free") %>%
  filter(Obligate_facultative != "Unknown") %>% filter(Effect_on_host != "Unknown") %>%
  filter(Ne_small != "NA")

pangenome_lifestyles_no_unknown_nofree_ne <- pangenome_lifestyles_no_unknown_nofree_ne %>%
  rename(PF = pangenome_fluidity, GS = genome_size, NE = Ne_small)

#Set row names as species
rownames(pangenome_lifestyles_no_unknown_nofree_ne) <- pangenome_lifestyles_no_unknown_nofree_ne$Species

## Make into data frame
pangenome_lifestyles_no_unknown_nofree_ne <- as.data.frame(pangenome_lifestyles_no_unknown_nofree_ne)

#### Code single lifestyle variable
pangenome_lifestyles_no_unknown_nofree_ne_1 <- pangenome_lifestyles_no_unknown_nofree_ne %>%
  mutate(OF=case_when(Obligate_facultative=="Obligate"~0,
                      Obligate_facultative=="Facultative"~1), 
         IE=case_when(Intra_or_extracellular=="Intracellular" ~0, 
                      Intra_or_extracellular=="Both"~0.5, 
                      Intra_or_extracellular=="Extracellular"~1), 
         EH= case_when(Effect_on_host=="Pathogen"~0,
                       Effect_on_host=="Both"~0.5,
                       Effect_on_host=="Mutualist"~1),
         M=case_when(Motility=="Non-motile" ~0, 
                     Motility=="Both"~0.5, 
                     Motility=="Motile"~1))

pangenome_lifestyles_no_unknown_nofree_ne_1 <- pangenome_lifestyles_no_unknown_nofree_ne_1 %>%
  mutate(LS = OF+IE+EH+M)


## ----LS vs pangenome fludity MCMCglmm, echo=FALSE, cache=TRUE----------------------------------------------------------
mcmc_model_LS <- MCMCglmm(PF ~ LS, random=~Species, data=pangenome_lifestyles_no_unknown_nofree_ne_1,
                          nitt=50000,
                          prior=prior, ginverse = list(Species = Ainv),verbose = FALSE)

# (Intercept)   0.10603  0.02547  0.18939     4700 0.0157 *  
# LS            0.04594  0.02687  0.06459     4700 <2e-04 ***
# DIC: -193.8435
#summary(mcmc_model_LS)

# R2 values
mFixed <- mean(mcmc_model_LS$Sol[,2]) * mcmc_model_LS$X[, 2]
mVarF<- var(mFixed)

# Fixed effect variance / (fixed effect variance + random effect variance)
# 0.1938511
fixed_r2_model_LS <- mVarF/(mVarF+sum(apply(mcmc_model_LS$VCV,2,mean)))
#random effect variance/(random effect variance + residual variance)
# 0.4127378
random_r2_model_LS <- (sum(apply(mcmc_model_LS$VCV,2,mean)[-2]))/((sum(apply(mcmc_model_LS$VCV,2,mean)))+mVarF)
#Conditionall R2 (total variance explained by the model)
# 0.6065888
conditional_r2_model_LS <- ((sum(apply(mcmc_model_LS$VCV,2,mean)[-2]))+mVarF)/((sum(apply(mcmc_model_LS$VCV,2,mean)))+mVarF)

R2_table_LS <- matrix(c(fixed_r2_model_LS,random_r2_model_LS,conditional_r2_model_LS), ncol=1)
R2_table_LS <- format(R2_table_LS, digits=4)
rownames(R2_table_LS) <- c("Fixed effect", "Random effect", "Total model")
colnames(R2_table_LS) <- c("R-squared value")

knitr::kable(list(summary_mcmc_glmm(mcmc_model_LS),R2_table_LS), caption="Results from a MCMCglmm with pangenome fluidity as the response variable, lifestyle vairability as the fixed effect, and phylogeny as a random effect.", digits=4)%>% kable_styling(latex_options = "HOLD_position")


## ----LS vs pangenome fluidity plot, fig.cap="Plot showing correlation between lifestyle as a single variable and species' pangenome fluidity. The line is the slope and intercept from the MCMCglmm analysis in the table above.", fig.height=5----
ggplot(pangenome_lifestyles_no_unknown_nofree_ne_1, aes(x=LS,y=PF)) +
  geom_point(colour="black") +
  labs(y="Pangenome \nfluidity", x="Lifestyle") +
  geom_abline(intercept=0.10603, slope=0.04594, colour="#0072B2", size=1) + #MCMCglmm, n=115
  scale_y_continuous(limits=c(0,0.48),breaks=c(0,0.1,0.2,0.3,0.4), labels=prettyZero, expand=c(0,0)) +
  theme_classic() +
  theme(plot.title = element_text(size=15, hjust=0.5), axis.title = element_text(size=13),
        axis.text = element_text(size=12), axis.title.y=element_text(angle=360,vjust=0.5))


## ----NE/GS/LS path analysis simple models------------------------------------------------------------------------------
models_LS_GS_NE_simple <- define_model_set(
  a = c(PF~LS+GS+NE),
  b = c(PF~LS+NE),
  c = c(PF~LS+GS),
  d = c(PF~GS+NE),
  e = c(PF~LS),
  f = c(PF~GS),
  g = c(PF~NE)
)



## ----NE/GS/LS path analysis simple models plot, fig.cap="Set of simple models, varying by which of lifestyle, effective population size and genome size causes pangenome fluidity.", fig.width=8, fig.height=8----
plot_model_set(models_LS_GS_NE_simple, labels = c(LS="Lifestyle", GS="Genome\nsize",
                PF="Panenome\nfluidity", NE="Eff. pop\nsize"), text_size=4) # Plot all models


## ----NE/GS/LS path analysis simple models results----------------------------------------------------------------------
## Run path analysis
## This time Brownian Motion is the model of evolution where the models have highest support..
result_6 <- phylo_path(models_LS_GS_NE_simple, data=pangenome_lifestyles_no_unknown_nofree_ne_1,
                       tree=dataTree, 
                       model="BM", method="logistic_MPLE")

## Path analysis results
s6 <- summary(result_6)

## Plot average best model
best_model_6_conditional <- best(result_6) 



## ----NE/GS/LS path analysis simple models results table----------------------------------------------------------------
## Path analysis results
s6_subset <- s6 %>% select(w,p,CICc)

knitr::kable(s6_subset, caption="Details of of support for a simple set of models varying by which of lifestyle, effective population size and genome size causes pangenome fluidity.", digits=4)%>% kable_styling(latex_options = "HOLD_position")


## ----NE/GS/LS path analysis simple models results plot, fig.cap="Comparison of support for a simple set of models varying by which of lifestyle, effective population size and genome size causes pangenome fluidity."----
# Single model is the best - model c.
# This has GS->GF and LS->GF, but not NE->GF
plot(s6)


## ----NE/GS/LS path analysis simple models best model, fig.cap="Model with best support, out of a set of simple models. "----
## Plot average best models
plot(best_model_6_conditional, curvature = 0.1, type="width", 
     labels = c(LS="Lifestyle", GS="Genome\nsize",
                PF="Panenome\nfluidity", NE="Eff. pop\nsize"), text_size=4) 


## ----NE/GS/LS path analysis complex models-----------------------------------------------------------------------------
# The factors may themselves be connected - we can explore this.
models_LS_GS_NE_complex <- define_model_set(
  a = c(PF~LS+GS+NE),
  b = c(PF~LS+NE),
  c = c(PF~LS+GS),
  d = c(PF~GS+NE),
  e = c(PF~LS),
  f = c(PF~GS),
  g = c(PF~NE),
  h = c(PF~LS+GS+NE, NE~LS),
  i = c(PF~LS+GS, NE~LS),
  j = c(PF~LS+GS+NE, NE~LS, GS~LS),
  k = c(PF~LS+GS, NE~LS, GS~LS),
  l = c(PF~LS+GS, GS~LS),
  m = c(PF~LS, GS~LS),
  n = c(PF~LS+GS, GS~LS, GS~NE)
)



## ----NE/GS/LS path analysis complex models plot, fig.cap="Set of simple and more complex models, varying by which of lifestyle, effective population size and genome size causes pangenome fluidity, and also how those factors might influence each other.", fig.width=8, fig.height=8----
plot_model_set(models_LS_GS_NE_complex, labels = c(LS="Lifestyle", GS="Genome\nsize",
                PF="Panenome\nfluidity", NE="Eff. pop\nsize"), text_size=4) # Plot all models


## ----NE/GS/LS path analysis complex models results---------------------------------------------------------------------
## This time Brownian Motion is a better model of evolution.
result_7 <- phylo_path(models_LS_GS_NE_complex, data=pangenome_lifestyles_no_unknown_nofree_ne_1, tree=dataTree, 
                       model="BM", method="logistic_MPLE")

## Path analysis results
s7 <- summary(result_7)



## ----NE/GS/LS path analysis complex models results table---------------------------------------------------------------
## Path analysis results
s7_subset <- s7 %>% select(w,p,CICc)

knitr::kable(s7_subset, caption="Comparison of set of more complex causal models, varying by which of lifestyle, genome size and effective population size directly influence pangenome fluidity, and how they might cause each other.", digits=4)%>% kable_styling(latex_options = "HOLD_position")


## ----NE/GS/LS path analysis complex models results average models------------------------------------------------------
# Four models have high support (within 2CICs)
# Model c, from result 6, has 3rd highest support.
# Top is model n, which has LS->GS and NE->GS, in addition to LS->GF and GS->GF.
# Second is model l, which is identical to model n but with no NE->GF.
plot(s7)

## Plot average best model
average_model_7_conditional <- average(result_7, avg_method = "conditional") 
average_model_7_full <- average(result_7, avg_method = "full") 


## ----NE/GS/LS path analysis complex models best conditional, fig.cap="Best model of causal relationships between lifestyle, genome size, effective population size and pangenome fluidity; an average of four models with similar structure. This is the same model as Figure 4 in the main text."----
## Plot average best models
# Suggests lifestyle is most important - has the largest direct impact on GF an causes other two factors.
# Effective population size has minor role - and no direct causal impact on GF.
plot(average_model_7_conditional, curvature = 0.1, algorithm="dh", type="width", 
     labels = c(LS="Lifestyle", GS="Genome\nsize",
                PF="Panenome\nfluidity", NE="Eff. pop\nsize"), text_size=3) 


## ----NE/GS/LS path analysis complex models coeffs, fig.cap="Path coeficients for the average best model, with 95% confidence intervals.", fig.height=5----
coef_plot(average_model_7_conditional, reverse_order = TRUE) + 
  ggplot2::coord_flip() + 
  ggplot2::theme_bw()


## ----NE/GS/LS path analysis complex models best full,fig.cap="Best model of causal relationships between lifestyle, genome size, effective population size and pangenome fluidity; an average of four models with similar structure, with correlation coefficients weighted by setting absent paths to zero when averaging models."----
plot(average_model_7_full, curvature = 0.1, algorithm="dh", type="width", 
     labels = c(LS="Lifestyle", GS="Genome\nsize",
                PF="Pangenome\nfluidity", NE="Eff. pop\nsize"), text_size=3) 


## ----additional path - binary lifestyle set up-------------------------------------------------------------------------
#### Code single lifestyle variable
pangenome_lifestyles_no_unknown_nofree_ne_2 <- pangenome_lifestyles_no_unknown_nofree_ne %>%
  mutate(OF=case_when(Obligate_facultative=="Obligate"~0,
                                Obligate_facultative=="Facultative"~1), 
         IE=case_when(Intra_or_extracellular=="Intracellular" ~0, 
                      Intra_or_extracellular=="Both"~0, 
                      Intra_or_extracellular=="Extracellular"~1),
         EH=case_when(Effect_on_host=="Pathogen" ~0, 
                      Effect_on_host=="Both"~0, 
                      Effect_on_host=="Mutualist"~1),
         M=case_when(Motility=="Non-motile" ~0, 
                     Motility=="Both"~1, 
                     Motility=="Motile"~1))

pangenome_lifestyles_no_unknown_nofree_ne_2 <- pangenome_lifestyles_no_unknown_nofree_ne_2 %>%
  mutate(LS = OF+IE+EH+M)


## ----binary LS vs pangenome fludity MCMCglmm, echo=FALSE, cache=TRUE---------------------------------------------------
mcmc_model_LS_binary <- MCMCglmm(PF ~ LS, random=~Species, data=pangenome_lifestyles_no_unknown_nofree_ne_2,
                          nitt=50000,
                          prior=prior, ginverse = list(Species = Ainv),verbose = FALSE)

# R2 values
mFixed <- mean(mcmc_model_LS_binary$Sol[,2]) * mcmc_model_LS_binary$X[, 2]
mVarF<- var(mFixed)

# Fixed effect variance / (fixed effect variance + random effect variance)
# 0.1938511
fixed_r2_model_LS_binary <- mVarF/(mVarF+sum(apply(mcmc_model_LS_binary$VCV,2,mean)))
#random effect variance/(random effect variance + residual variance)
# 0.4127378
random_r2_model_LS_binary <- (sum(apply(mcmc_model_LS_binary$VCV,2,mean)[-2]))/((sum(apply(mcmc_model_LS_binary$VCV,2,mean)))+mVarF)
#Conditionall R2 (total variance explained by the model)
# 0.6065888
conditional_r2_model_LS_binary <- ((sum(apply(mcmc_model_LS_binary$VCV,2,mean)[-2]))+mVarF)/((sum(apply(mcmc_model_LS_binary$VCV,2,mean)))+mVarF)

R2_table_LS_binary <- matrix(c(fixed_r2_model_LS_binary,random_r2_model_LS_binary,conditional_r2_model_LS_binary), ncol=1)
R2_table_LS_binary <- format(R2_table_LS_binary, digits=4)
rownames(R2_table_LS_binary) <- c("Fixed effect", "Random effect", "Total model")
colnames(R2_table_LS_binary) <- c("R-squared value")

knitr::kable(list(summary_mcmc_glmm(mcmc_model_LS_binary),R2_table_LS_binary), caption="Results from a MCMCglmm with pangenome fluidity as the response variable, lifestyle vairability as the fixed effect (calculated with no intermediate/both category), and phylogeny as a random effect.", digits=4)%>% kable_styling(latex_options = "HOLD_position")


## ----binary LS vs pangenome fluidity plot, fig.cap="Plot showing correlation between lifestyle as a single variable (calculated with no intermediate/both category) and species' pangenome fluidity. The line is the slope and intercept from the MCMCglmm analysis in the table above."----
ggplot(pangenome_lifestyles_no_unknown_nofree_ne_2, aes(x=LS,y=PF)) +
  geom_point(colour="black") +
  labs(y="Pangenome \nfluidity", x="Lifestyle") +
  geom_abline(intercept=0.12484, slope=0.04264, colour="#0072B2", size=1) + #MCMCglmm, n=115
  scale_y_continuous(limits=c(0,0.48),breaks=c(0,0.1,0.2,0.3,0.4), labels=prettyZero, expand=c(0,0)) +
  theme_classic() +
  theme(plot.title = element_text(size=15, hjust=0.5), axis.title = element_text(size=13),
        axis.text = element_text(size=12), axis.title.y=element_text(angle=360,vjust=0.5))


## ----additional path - binary lifestyle simple models results----------------------------------------------------------
## Run path analysis
## This time Brownian Motion is the model of evolution where the models have highest support..
result_8 <- phylo_path(models_LS_GS_NE_simple, data=pangenome_lifestyles_no_unknown_nofree_ne_2,
                       tree=dataTree, 
                       model="BM", method="logistic_MPLE")

## Path analysis results
s8 <- summary(result_8)



## ----additional path - binary lifestyle simple models results table----------------------------------------------------
## Path analysis results
s8_subset <- s8 %>% select(w,p,CICc)

knitr::kable(s8_subset, caption="Details of of support for a simple set of models varying by which of lifestyle, effective population size and genome size causes pangenome fluidity.", digits=4)%>% kable_styling(latex_options = "HOLD_position")


## ----additional path - binary lifestyle simple models results plot, fig.cap="Comparison of support for a simple set of models varying by which of lifestyle (calculated with no intermediate/both category), effective population size and genome size causes pangenome fluidity."----
plot(s8)


## ----additional path - binary lifestyle simple models best model, fig.cap="Average of two models with best support, out of a set of simple models."----
average_model_8_conditional <- average(result_8, avg_method = "conditional") 
## Plot average best models
plot(average_model_8_conditional, curvature = 0.1, type="width", 
     labels = c(LS="Lifestyle", GS="Genome\nsize",
                PF="Panenome\nfluidity", NE="Eff. pop\nsize"), text_size=4) 


## ----additional path - binary lifestyle complex models results---------------------------------------------------------
## This time Brownian Motion is a better model of evolution.
result_9 <- phylo_path(models_LS_GS_NE_complex, data=pangenome_lifestyles_no_unknown_nofree_ne_2, tree=dataTree, 
                       model="BM", method="logistic_MPLE")

## Path analysis results
s9 <- summary(result_9)



## ----additional path - binary lifestyle complex models results average models, fig.cap=""------------------------------
# Four models have high support (within 2CICs)
# Model c, from result 6, has 3rd highest support.
# Top is model n, which has LS->GS and NE->GS, in addition to LS->GF and GS->GF.
# Second is model l, which is identical to model n but with no NE->GF.
plot(s9)



## ----additional path - binary lifestyle complex models best conditional, fig.cap="Best model of causal relationships between lifestyle, genome size, effective population size and pangenome fluidity; an average of four models with similar structure."----
average_model_9_conditional <- average(result_9, avg_method = "conditional") 
## Plot average best models
# Suggests lifestyle is most important - has the largest direct impact on GF an causes other two factors.
# Effective population size has minor role - and no direct causal impact on GF.
plot(average_model_9_conditional, curvature = 0.1, algorithm="dh", type="width", 
     labels = c(LS="Lifestyle", GS="Genome\nsize",
                PF="Panenome\nfluidity", NE="Eff. pop\nsize"), text_size=3) 


## ----MCA set-up--------------------------------------------------------------------------------------------------------
## Option 1 - code as binary and/or discrete traits
## Don't have Effect on host in this analysis - not related to environmental variability
## More & less indicates level of environmental variability
pangenome_lifestyles_no_unknown_nofree_MCA <- pangenome_lifestyles_no_unknown_nofree %>%
  mutate(OF=as.factor(case_when(Obligate_facultative=="Obligate"~"Low",
                      Obligate_facultative=="Facultative"~"High")), 
         IE=as.factor(case_when(Intra_or_extracellular=="Intracellular" ~"Low", 
                      Intra_or_extracellular=="Both"~"Intermediate", 
                      Intra_or_extracellular=="Extracellular"~"High")), 
         EH=as.factor(case_when(Effect_on_host=="Pathogen" ~"Low", 
                      Effect_on_host=="Both"~"Intermediate", 
                      Effect_on_host=="Mutualist"~"High")),
         M=as.factor(case_when(Motility=="Non-motile" ~"Low", 
                     Motility=="Both"~"Intermediate", 
                     Motility=="Motile"~"High")))

### Multiple Component Analysis

MCA_dataset_1 <- pangenome_lifestyles_no_unknown_nofree_MCA %>%
  select(OF,IE,EH,M)

summary_MCA_dataset_1 <- summary(MCA_dataset_1)[, 1:4]

## MCA Analysis
MCA_results_1 <- MCA(MCA_dataset_1, graph = FALSE)



## ----MCA initial results, fig.cap="Dimensions produced from MCA analysis on four lifestyle traits across 115 species. Dimensions are ordered from left to right by the percentage of total variation that they explain, which is printed above each bar."----
#print(MCA_results_1)
eig.val_1 <- get_eigenvalue(MCA_results_1)

fviz_screeplot(MCA_results_1, addlabels = TRUE, ylim = c(0, 45))



## ----MCA results exploration, fig.cap="Contribution of variables and their categories to the two dimensions explaining the most variation in lifestyle across the species. OF = host reliance (obligate/facultative); IE = host location (intra/extracellular); EH = effect on host (pathogen/mutualist); M = motility (non-motile/motile)."----
## Results
var_1 <- get_mca_var(MCA_results_1)

fviz_contrib(MCA_results_1, choice = "var", axes = 1, top = 15) +#Very uniform
fviz_contrib(MCA_results_1, choice = "var", axes = 2, top = 15) #Some not contributing


## ----MCA dimension correlations, fig.cap="Correlation between dimensions 1 and 2, with a panel for each lifestyle variable highlighting the location of the categories within that lifestyle. Dots represent possible values for a species for each dimension, and colours are to distinguish the categories. Dimension 1 groups species in a way that better reflects what we would like from a single lifestyle variable."----

fviz_ellipses(MCA_results_1, c("OF","IE", "EH","M"),geom = "point") # Dimension 1 seems to be pretty good



## ----MCA extract dimension 1-------------------------------------------------------------------------------------------
ind_1 <- get_mca_ind(MCA_results_1)
species_dimension1_1 <- as.vector(ind_1$coord[,1]) 

pangenome_lifestyles_no_unknown_nofree_MCA_dim1 <- cbind(pangenome_lifestyles_no_unknown_nofree,species_dimension1_1)

# Code original lifestyle variable from main analysis
pangenome_lifestyles_no_unknown_nofree_MCA_dim1_LS <- pangenome_lifestyles_no_unknown_nofree_MCA_dim1 %>%
  mutate(OF=case_when(Obligate_facultative=="Obligate"~0,
                      Obligate_facultative=="Facultative"~1), 
         IE=case_when(Intra_or_extracellular=="Intracellular" ~0, 
                      Intra_or_extracellular=="Both"~0.5, 
                      Intra_or_extracellular=="Extracellular"~1), 
         EH= case_when(Effect_on_host=="Pathogen"~0,
                       Effect_on_host=="Both"~0.5,
                       Effect_on_host=="Mutualist"~1),
         M=case_when(Motility=="Non-motile" ~0, 
                     Motility=="Both"~0.5, 
                     Motility=="Motile"~1))

pangenome_lifestyles_no_unknown_nofree_MCA_dim1_LS <- pangenome_lifestyles_no_unknown_nofree_MCA_dim1_LS %>%
  mutate(LS = OF+IE+EH+M, LS_MCA = (-1*species_dimension1_1))



## ----LS vs MCA LS plot, fig.cap="Correlation between our original lifestyle variable used in our main analysis and dimension 1 from the MCA."----
ggplot(pangenome_lifestyles_no_unknown_nofree_MCA_dim1_LS, aes(x=LS,y=LS_MCA)) +
  geom_point() +
  theme_classic() +
  labs(x="Lifestyle\n(summed)", y="MCA\nDimension 1") +
  theme(plot.title = element_text(size=15, hjust=0.5), axis.title = element_text(size=13),
        axis.text = element_text(size=12), axis.title.y=element_text(angle=360,vjust=0.5), 
        legend.text = element_text(size=12))


## ----MCA NE/GS/LS path analysis----------------------------------------------------------------------------------------
#Set row names as species
#rownames(pangenome_lifestyles_no_unknown_nofree_MCA_dim1_LS) <- pangenome_lifestyles_no_unknown_nofree_MCA_dim1_LS$Species

pangenome_lifestyles_no_unknown_nofree_MCA_dim1_LS$Ne_small <- pangenome_lifestyles_no_unknown_nofree_MCA_dim1_LS$Ne/1000000

pangenome_lifestyles_no_unknown_nofree_MCA_dim1_LS_ne <- pangenome_lifestyles_no_unknown_nofree_MCA_dim1_LS %>%
  filter(Ne_small != "NA")

pangenome_lifestyles_no_unknown_nofree_MCA_dim1_LS_ne <- pangenome_lifestyles_no_unknown_nofree_MCA_dim1_LS_ne %>% rename(NE = Ne_small,GS=genome_size)



## ----MCA NE/GS/LS path analysis simple models, fig.cap="Set of simple models varying by influence of lifestyle (MCA dimension 1), effective population size and genome size, on pangenome fluidity.", fig.height=6----
models_LS_MCA_GS_NE_simple <- define_model_set(
  a = c(PF~LS_MCA+GS+NE),
  b = c(PF~LS_MCA+NE),
  c = c(PF~LS_MCA+GS),
  d = c(PF~GS+NE),
  e = c(PF~LS_MCA),
  f = c(PF~GS),
  g = c(PF~NE)
)

plot_model_set(models_LS_MCA_GS_NE_simple, text_size = 4)


## ----MCA E/GS/LS path analysis simple models results-------------------------------------------------------------------

result_8 <- phylo_path(models_LS_MCA_GS_NE_simple, data=pangenome_lifestyles_no_unknown_nofree_MCA_dim1_LS_ne,
                       tree=dataTree, 
                       model="BM", method="logistic_MPLE")

## Path analysis results
s8 <- summary(result_8)

## Plot average best model
average_model_8_conditional <- average(result_8, avg_method="conditional") 

s8_subset <- s8 %>% select(w,p,CICc)

knitr::kable(s8_subset, caption="Comparison of support for a simple set of models varying by which of lifestyle (MCA dimension 1), effective population size and genome size causes pangenome fluidity.", digits=4)%>% kable_styling(latex_options = "HOLD_position")


## ----MCA E/GS/LS path analysis simple models results plot, fig.cap="Comaparison of support for a simple set of models varying by which of lifestyle (MCA dimension 1), effective population size and genome size causes pangenome fluidity."----
plot(s8)


## ----MCA NE/GS/LS path analysis simple models results d sep best model, echo=FALSE-------------------------------------
# Rejected because lifestyle not conditionally independent from 
result_8_dsep_modelc <- result_8$d_sep$c

result_8_dsep_modelc_subset <- result_8_dsep_modelc %>% select(d_sep, p)

knitr::kable(result_8_dsep_modelc_subset, caption="Specific resuls for model c, showing support for each conditional independency specified by the model. Lifestyle and genome size are not conditionally independent, indicated by the p-value of less than 0.05, meaning they are significantly correlated.", digits=4)%>% kable_styling(latex_options = "HOLD_position")


## ----MCA NE/GS/LS path analysis complex models, fig.cap="A set of more complex causal models, varying by which of lifestyle (MCA dimension 1), genome size and effective population size directly influence pangenome fluidity, and how they might influence each other.", fig.height=10----

models_LS_MCA_GS_NE_complex <- define_model_set(
  a = c(PF~LS_MCA+GS+NE),
  b = c(PF~LS_MCA+NE),
  c = c(PF~LS_MCA+GS),
  d = c(PF~GS+NE),
  e = c(PF~LS_MCA),
  f = c(PF~GS),
  g = c(PF~NE),
  h = c(PF~LS_MCA+GS+NE, NE~LS_MCA),
  i = c(PF~LS_MCA+GS, NE~LS_MCA),
  j = c(PF~LS_MCA+GS+NE, NE~LS_MCA, GS~LS_MCA),
  k = c(PF~LS_MCA+GS, NE~LS_MCA, GS~LS_MCA),
  l = c(PF~LS_MCA+GS, GS~LS_MCA),
  m = c(PF~LS_MCA, GS~LS_MCA),
  n = c(PF~LS_MCA+GS, GS~LS_MCA, GS~NE)
)

plot_model_set(models_LS_MCA_GS_NE_complex, labels = c(LS_MCA="Lifestyle\n(MCA)", GS="Genome\nsize",
                PF="Panenome\nfluidity", NE="Eff. pop\nsize"), text_size=3) # Plot all models



## ----MCA NE/GS/LS path analysis complex models results-----------------------------------------------------------------

result_9 <- phylo_path(models_LS_MCA_GS_NE_complex, data=pangenome_lifestyles_no_unknown_nofree_MCA_dim1_LS_ne, tree=dataTree, 
                       model="BM", method="logistic_MPLE")

## Path analysis results
s9 <- summary(result_9)

## Path analysis results
s9_subset <- s9 %>% select(w,p,CICc)

knitr::kable(s9_subset, caption="Comparison of set of more complex causal models, varying by which of lifestyle (MCA dimension 1), genome size and effective population size directly influence pangenome fluidity, and how they might cause each other.", digits=4)%>% kable_styling(latex_options = "HOLD_position")


## ----MCA E/GS/LS path analysis complex models results plot, fig.cap="Comparison of support for a more complex set of models varying by which of lifestyle (MCA dimension 1), genome size and effective population size directly influence pangenome fluidity, and how they might cause each other."----
plot(s9)


## ----MCA NE/GS/LS path analysis complex models results average models conditional, fig.cap="Best model of causal relationships between lifestyle, genome size, effective population size and pangenome fluidity; an average of four models."----
## Plot average best model
average_model_9_conditional <- average(result_9, avg_method = "conditional") 

plot(average_model_9_conditional, curvature = 0.1, algorithm="dh", type="width", 
     labels = c(LS_MCA="Lifestyle\n(MCA)", GS="Genome\nsize",
                PF="Panenome\nfluidity", NE="Eff. pop\nsize"), text_size=3) 


## ----MCA NE/GS/LS path analysis complex models results average models full, fig.cap="Best model of causal relationships between lifestyle, genome size, effective population size and pangenome fluidity; an average of four models with similar structure, with correlation coefficients weighted by setting absent paths to zero when averaging models."----

average_model_9_full <- average(result_9, avg_method = "full") 

plot(average_model_9_full, curvature = 0.1, algorithm="dh", type="width", 
     labels = c(LS_MCA="Lifestyle\n(MCA)", GS="Genome\nsize",
                PF="Pangenome\nfluidity", NE="Eff. pop\nsize"), text_size=3) 


## ----host/free/most and NE path set-up---------------------------------------------------------------------------------
pangenome_lifestyles_host_free_4_ne <- pangenome_lifestyles_no_unknown_nofree_ne %>%
  mutate(Host_free_4 = case_when(Host_or_free=="Host" ~ "Host",
                                 Host_or_free=="Free" ~ "Free",
                                 Host_or_free=="Both" & Category_primary_env == "Free" ~ "Mostly_free",
                                 Host_or_free=="Both" & Category_primary_env == "Host" ~ "Mostly_host",
                                 Host_or_free=="Unknown"~"Unknown"))

# Reorder factor
pangenome_lifestyles_host_free_4_ne$Host_free_4 <- factor(pangenome_lifestyles_host_free_4_ne$Host_free_4, 
                                                       levels=c("Host","Mostly_host","Mostly_free", "Free"))

pangenome_lifestyles_no_unknown_nofree_ne_2 <- pangenome_lifestyles_host_free_4_ne %>%
  mutate(OF=case_when(Obligate_facultative=="Obligate"~0,
                      Obligate_facultative=="Facultative"~1), 
         IE=case_when(Intra_or_extracellular=="Intracellular" ~0, 
                      Intra_or_extracellular=="Both"~0.5, 
                      Intra_or_extracellular=="Extracellular"~1), 
         EH= case_when(Effect_on_host=="Pathogen"~0,
                       Effect_on_host=="Both"~0.5,
                       Effect_on_host=="Mutualist"~1),
         M=case_when(Motility=="Non-motile" ~0, 
                     Motility=="Both"~0.5, 
                     Motility=="Motile"~1),
         HF=case_when(Host_or_free=="Host"~0,
                      Host_or_free=="Both"~1), # No free in 115 species
         PEnv=case_when(Category_primary_env=="Host"~0,
                        Category_primary_env=="Free"~1),
         HF_4=case_when(Host_free_4=="Host"~0,
                        Host_free_4=="Mostly_host"~0.5,
                        Host_free_4=="Mostly_free"~1)) # No free-living in 115 species

pangenome_lifestyles_no_unknown_nofree_ne_2 <- pangenome_lifestyles_no_unknown_nofree_ne_2 %>%
  mutate(LS = OF+IE+EH+M)



## ----host/free/most and NE path simple models--------------------------------------------------------------------------
##Code models
# Host/free vs NE
models_HF4_NE <- define_model_set(
  a = c(PF~NE+HF_4),
  b = c(PF~NE),
  c = c(PF~HF_4),
  d = c(PF~HF_4, NE~HF_4),
  e = c(PF~NE, HF_4~NE)
)



## ----host/free/most and NE path simple models plot, fig.cap="Set of simple models varying by which of effective population size and whether a species is mostly host-associated or free-living influences pangenome fluidity.", fig.height=7----
plot_model_set(models_HF4_NE, labels = c(HF_4="Host/free",
                PF="Pangenome\nfluidity", NE="Eff. pop\nsize"), text_size=3)


## ----host/free/most and NE path simple models results------------------------------------------------------------------
## Run path analysis
## This time Brownian Motion is a better model of evolution.
result_10 <- phylo_path(models_HF4_NE, data=pangenome_lifestyles_no_unknown_nofree_ne_2, tree=dataTree, 
                       model="BM", method="logistic_MPLE")


## Path analysis results
s10 <- summary(result_10)

## Path analysis results
s10_subset <- s10 %>% select(w,p,CICc)

knitr::kable(s10_subset, caption="Comparison of set of more complex causal models, varying by which of host-association (host, mostly host, mostly free) and effective population size influence pangenome fluidity.", digits=4)%>% kable_styling(latex_options = "HOLD_position")



## ----host/free/most and NE path simple models results plot, fig.cap="Comparison of set of more complex causal models, varying by which of host-association (host, mostly host, mostly free) and effective population size influence pangenome fluidity."----
plot(s10)


## ----host/free/most and NE path simple models results average best model conditional, fig.cap="Average best model of causal relationships between host-association (host-associated, mostly host and mostly free-living), effective population size and pangenome fluidity."----
## Plot average best model
average_model_10_conditional <- average(result_10, avg_method = "conditional") 

## Plot average best models
## Clear that host/free and/or primary env is most important for GF -> little support for NE.
plot(average_model_10_conditional, curvature = 0.1, algorithm="dh", type="width", 
     labels = c(PF="Pangenome\nfluidity", NE="Eff. pop\nsize", HF_4="Host or\nfree"), text_size=3) 


## ----host/free/most and NE path simple models results average best model full, fig.cap="Average best model of causal relationships between host-association (host-associated, mostly host and mostly free-living), effective population size and pangenome fluidity; an average of four models with similar structure, with correlation coefficients weighted by setting absent paths to zero when averaging models."----
## Plot average best model
average_model_10_full <- average(result_10, avg_method = "full") 

## Plot average best models
## Clear that host/free and/or primary env is most important for GF -> little support for NE.
plot(average_model_10_full, curvature = 0.1, algorithm="dh", type="width", 
     labels = c(PF="Pangenome\nfluidity", NE="Eff. pop\nsize", HF_4="Host or\nfree"), text_size=3) 


## ----core 100 in pangenome vs fluidity, fig.cap="Pangenome fluidity is highly correlated with the proportion of core genes in the pangenome.\nScatterplot showing how pangenome fluidity varies with the percentage of genes in a species’ pangenome which are ‘core’, defined in this plot as present in 100% of genomes. Data includes all 126 species."----
## Percentage of pangenome that is core (core=100% genomes) vs pangenome fluidity
ggplot(pangenome_lifestyles, aes(x=percentage_core_100,y=pangenome_fluidity)) +
  geom_point() +
  labs(x="Percentage of core genes in pangenome\n(core=100% genomes)", y="Pangenome\nfluidity") +
  scale_y_continuous(limits=c(0,0.43),breaks=c(0,0.1,0.2,0.3,0.4), expand=c(0,0),labels=prettyZero) +
  scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),limits=c(0,1.05),labels=prettyZero,expand=c(0,0))+  
  stat_smooth(method="lm", se=F, colour="#0072B2") +
  theme_classic() +
  theme(plot.title = element_text(size=15, hjust=0.5), axis.title = element_text(size=13),
        axis.text = element_text(size=12), axis.title.y=element_text(angle=360,vjust=0.5))


## ----histogram of pangenome size, fig.cap="Pangenome size varies considerably across species.\nHistogram showing the variation in pangenome size, meaning the number of unique genes sequenced across all genomes of a species. Data is for all 126 species."----

# Shows variation in pangenome size
ggplot(pangenome_lifestyles, aes(x=pan_size_100)) +
  geom_histogram(colour="#009E73",fill="#009E73", alpha=0.5) +
  labs(x="Pangenome size (number of genes)", y="Number of\nspecies") +
  scale_x_continuous(breaks=c(0,10000,20000,30000))+
  theme_classic() +
  theme(plot.title = element_text(size=15, hjust=0.5), axis.title = element_text(size=13),
        axis.text = element_text(size=12), axis.title.y=element_text(angle=360,vjust=0.5))



## ----stacked bar plot genome size core vs accessory, fig.cap="A stacked bar plot, with a bar for each species; the total height is the average genome size for that species, the ‘core’ are the genes present in 100% of genomes, and the ‘accessory’ genes are any remaining genes not present in 100% of genomes."----
# Show variation in pangenome size, with core & accessory at 100% threshold
pangenome_lifestyles <- pangenome_lifestyles %>%
  mutate(genome_accessory_100 = genome_size-core_100)

pangenome_lifestyles_stacked_genome <- pangenome_lifestyles %>%
  select(Species,genome_size,core_100,genome_accessory_100)

pangenome_lifestyles_stacked_genome <- pangenome_lifestyles_stacked_genome %>%
  pivot_longer(cols=c(core_100,genome_accessory_100))

pangenome_lifestyles_stacked_genome$name <- factor(pangenome_lifestyles_stacked_genome$name, levels=c("genome_accessory_100","core_100"))

ggplot(pangenome_lifestyles_stacked_genome, aes(fill=name, y=value,x=reorder(Species,genome_size))) +
  geom_bar(position="stack", stat="identity",colour="white") +
  scale_y_continuous(limits=c(0,7000), expand=c(0,0)) +
  scale_fill_manual("Type of gene", values = c("core_100" = "#0072B2","genome_accessory_100" = "#56B4E9"),
                    labels=c("core_100"="Core", "genome_accessory_100"="Accessory"))+
  labs(x="Species", y="Genome\nSize") +
  theme_classic() +
  theme(plot.title = element_text(size=15, hjust=0.5), axis.title = element_text(size=13),
        axis.text.y = element_text(size=12), axis.title.y=element_text(angle=360,vjust=0.5),
        axis.text.x = element_blank(), axis.ticks.x=element_blank())



## ----fluidity vs number of genomes, fig.height=7, fig.cap="No correlation between pangenome fluidity and the number of genomes per species. A. All species; B. Subset with species between 10-50 genomes."----
#summary(lm(pangenome_fluidity~Number_of_Strain,data=pangenome_lifestyles))


fluidity_n_genomes_a <- ggplot(pangenome_lifestyles, aes(x=Number_of_Strain, y=pangenome_fluidity)) +
  geom_point(colour="#306E80") +
  labs(x="Number of genomes", y="Pangenome \nfluidity") +
  scale_y_continuous(breaks=c(0,0.1,0.2,0.3,0.4)) +
  stat_smooth(method="lm", se=T, colour="#0072B2") +
  theme_classic() +
  theme(plot.title = element_text(size=15, hjust=0.5), axis.title = element_text(size=13),
        axis.text = element_text(size=12), axis.title.y=element_text(angle=360,vjust=0.5))

fluidity_n_genomes_b <- ggplot(pangenome_lifestyles, aes(x=Number_of_Strain, y=pangenome_fluidity)) +
  geom_point(colour="#306E80") +
  labs(x="Number of genomes", y="Pangenome \nfluidity") +
  scale_y_continuous(breaks=c(0,0.1,0.2,0.3,0.4)) +
  stat_smooth(method="lm", se=T, colour="#0072B2") +
  scale_x_continuous(limits=c(0,50)) +
  theme_classic() +
  theme(plot.title = element_text(size=15, hjust=0.5), axis.title = element_text(size=13),
        axis.text = element_text(size=12), axis.title.y=element_text(angle=360,vjust=0.5))

fluidity_n_genomes <- fluidity_n_genomes_a + fluidity_n_genomes_b + 
  plot_annotation(tag_levels = 'A') +
  plot_layout(ncol=1)

fluidity_n_genomes



## ----percentage core pangenome vs number of genomes, fig.height=7, fig.cap="The percentage of the pangenome which is core genes is negatively correlated with the number of genomes per species. A. All species; B. Subset with species between 10-50 genomes."----
#summary(lm(percentage_core_100~Number_of_Strain,data=pangenome_lifestyles))

per_core_n_genomes_a <- ggplot(pangenome_lifestyles, aes(x=Number_of_Strain, y=percentage_core_100)) +
  geom_point(colour="#306E80") +
  labs(x="Number of genomes", y="% of core\ngenes in\npangenome") +
  scale_y_continuous(limits=c(0,1),breaks=c(0,0.2,0.4,0.6,0.8,1), expand=c(0,0),labels=prettyZero) +
  stat_smooth(method="lm", se=T, colour="#0072B2") +
  theme_classic() +
  theme(plot.title = element_text(size=15, hjust=0.5), axis.title = element_text(size=13),
        axis.text = element_text(size=12), axis.title.y=element_text(angle=360,vjust=0.5))

per_core_n_genomes_b <- ggplot(pangenome_lifestyles, aes(x=Number_of_Strain, y=percentage_core_100)) +
  geom_point(colour="#306E80") +
  labs(x="Number of genomes", y="% of core\ngenes in\npangenome") +
  scale_y_continuous(limits=c(0,1),breaks=c(0,0.2,0.4,0.6,0.8,1), expand=c(0,0),labels=prettyZero) +
  stat_smooth(method="lm", se=T, colour="#0072B2") +
  scale_x_continuous(limits=c(0,50)) +
  theme_classic() +
  theme(plot.title = element_text(size=15, hjust=0.5), axis.title = element_text(size=13),
        axis.text = element_text(size=12), axis.title.y=element_text(angle=360,vjust=0.5))

per_core_n_genomes <- per_core_n_genomes_a + per_core_n_genomes_b + 
  plot_annotation(tag_levels = 'A') +
  plot_layout(ncol=1)

per_core_n_genomes


## ----pangenome size vs number of genomes, fig.height=7, fig.cap="The size of the pangenome is positively correlated with the number of genomes per species. A. All species; B. Subset with species between 10-50 genomes."----

#summary(lm(pan_size_100~Number_of_Strain,data=pangenome_lifestyles))

pan_size_n_genomes_a <- ggplot(pangenome_lifestyles, aes(x=Number_of_Strain, y=pan_size_100)) +
  geom_point(colour="#306E80") +
  labs(x="Number of genomes", y="Pangenome size\n(number of genes)") +
  stat_smooth(method="lm", se=T, colour="#0072B2") +
  theme_classic() +
  theme(plot.title = element_text(size=15, hjust=0.5), axis.title = element_text(size=13),
        axis.text = element_text(size=12), axis.title.y=element_text(angle=360,vjust=0.5))

pan_size_n_genomes_b <- ggplot(pangenome_lifestyles, aes(x=Number_of_Strain, y=pan_size_100)) +
  geom_point(colour="#306E80") +
  labs(x="Number of genomes", y="Pangenome size\n(number of genes)") +
  stat_smooth(method="lm", se=T, colour="#0072B2") +
  scale_x_continuous(limits=c(0,50)) +
  scale_y_continuous(limits=c(0,21000)) +
  theme_classic() +
  theme(plot.title = element_text(size=15, hjust=0.5), axis.title = element_text(size=13),
        axis.text = element_text(size=12), axis.title.y=element_text(angle=360,vjust=0.5))

pan_size_n_genomes <- pan_size_n_genomes_a + pan_size_n_genomes_b + 
  plot_annotation(tag_levels = 'A') +
  plot_layout(ncol=1)

pan_size_n_genomes


## ----tree plot species labels, fig.height=12, fig.width=8, fig.cap="Phylogeny of all 126 species in our dataset, which is available to download in ultrametric nexus format from github (*link*). Coloured squares indicate whether a species is host-associated, free-living, or a mixture of both, and bars indicate pangenome fluidity."----
#tree<-read.tree("tree.tre")

basic_tree <- ggtree(dataTree) + geom_tiplab(size=2.5, fontface=3) +
  geom_rootedge(T) 

tree_2 <- basic_tree + 
  geom_fruit(data=pangenome_lifestyles_host_free_4, 
             geom=geom_tile, 
             position=position_identityx(hexpand=24), #alters size of ring
             mapping=aes(y=Species,fill=Host_free_4), 
             color = "white",
             lwd = 0.2,
             linetype = 1,
             offset = 4,
             pwidth=0.25) +
  scale_fill_manual("Host-associated\n or free-living",values=c("#CC79A7", "#569AE9", "#56C8E9", "#009E73","#999999"))

tree_3 <- tree_2 +
  geom_fruit(data=pangenome_lifestyles, 
             geom=geom_bar, 
             position=position_identityx(hexpand=25),
             mapping=aes(y=Species,x=pangenome_fluidity),
             stat="Identity", 
             orientation="y",
             fill="#0072B2",
             pwidth=0.2) +
  theme(plot.margin = unit(c(0,0,0,0), "mm"))

tree_3


## ----species lifestyle table-------------------------------------------------------------------------------------------
pangenome_lifestyles_table <- pangenome_lifestyles_host_free_4 %>%
  select(Species,Host_free_4,
         Obligate_facultative,Intra_or_extracellular,
         Effect_on_host, Motility, Number_of_Strain)

pangenome_lifestyles_table <- as.data.frame(pangenome_lifestyles_table)

pangenome_lifestyles_table[is.na(pangenome_lifestyles_table)] <- "Unknown"


knitr::kable(pangenome_lifestyles_table,longtable=TRUE,booktabs=TRUE, linesep="",
             caption="Information for key lifestyle traits for all 126 species in our dataset.", 
             col.names = c("Species","Host or free","Host reliance",
                           "Host location", "Effect on host", "Motility", "# Genomes"))%>% 
  kable_styling(latex_options = c("HOLD_position","repeat_header"), font_size=8)


