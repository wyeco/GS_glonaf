##### code for genome size + glonaf data analysis  ####
#####   27 Nov 2021 ####
rm(list = ls())


library(ggplot2)
library(phylolm)
library(rr2)
library(ape)
library(geiger)
library(ggpubr)
library(dplyr)
setwd("./data_analysis")

## mported the cleaned data
GenomeSize_df_final <- read.csv("GenomeSize_glnaf.csv", header = T, sep = ',')
str(GenomeSize_df_final)
hist(GenomeSize_df_final$no_ploidy_lev)
hist(GenomeSize_df_final$holoploid)
GenomeSize_df_final$PLOIDY <- as.factor(GenomeSize_df_final$PLOIDY)
GenomeSize_df_final$taxno <- as.factor(GenomeSize_df_final$taxno)
GenomeSize_df_final$species <- as.factor(GenomeSize_df_final$species)
summary(GenomeSize_df_final)


### prepare the species taxonomy information 
library(taxize)
### Retrieve higher taxonomic names

tax_name(query = "Arceuthobium", get = c("order", 'family'), db = 'ncbi') 

## obtain the phylogeny
library("V.PhyloMaker")
library("ggtree")
library("ggplot2")

genome_taxonomy <- read.csv("./genome_taxnomy.csv", 
                            header = T, sep = ';')
genome.tre <- phylo.maker(sp.list = genome_taxonomy, 
                          tree = GBOTB.extended, nodes = nodes.info.1, scenarios="S3")
genome.phy <- genome.tre$scenario.3 

write.tree(genome.phy, file ="./genome.phylogeny.tre")
#### phylogeny ends here

library(plyr)
#####
## it does not work anymore as there are many NAs in this data.
summ_glonaf_data <- ddply(GenomeSize_df_final, "taxno", summarise,
                          N    = length(monoploid),
                          mean = mean(monoploid),
                          sd   = sd(monoploid),
                          se   = sd / sqrt(N)
)

summ_glonaf_data
##obtain the mean of each group
mean.monoploid.data <- ddply(GenomeSize_df_final, "taxno", summarise, mean = mean(monoploid))
median.monoploid.data <- ddply(genome_glonaf_data, "taxno", summarise, median = median(monoploid))
min.monoploid.data <- ddply(genome_glonaf_data, "taxno", summarise, min = min(monoploid))
max.monoploid.data <- ddply(genome_glonaf_data, "taxno", summarise, max = max(monoploid))
mean.monoploid.data 
median.monoploid.data
min.monoploid.data
max.monoploid.data

genome_glonaf_angio_data <- subset(genome_glonaf_data, TAXONOMY =="angiosperms", drop =T)
genome_glonaf_fern_data <- subset(genome_glonaf_data, TAXONOMY =="ferns", drop =T)
genome_glonaf_gymno_data <- subset(genome_glonaf_data, TAXONOMY =="gymnosp", drop =T)

mean.naturalized.angio.monoploid.data <- ddply(genome_glonaf_angio_data, "as.factor(GloNAF_incidence)", summarise, mean = mean(monoploid), sd   = sd(monoploid))
mean.naturalized.angio.monoploid.data

mean.naturalized.fern.monoploid.data <- ddply(genome_glonaf_fern_data, "as.factor(GloNAF_incidence)", summarise, mean = mean(monoploid),sd   = sd(monoploid))
mean.naturalized.fern.monoploid.data
mean.naturalized.gymno.monoploid.data <- ddply(genome_glonaf_gymno_data, "as.factor(GloNAF_incidence)", summarise, mean = mean(monoploid),sd   = sd(monoploid))
mean.naturalized.gymno.monoploid.data


anova_angio <- aov(log10(monoploid) ~ as.factor(GloNAF_incidence), genome_glonaf_angio_data)
plot(anova_angio,2)
summary(anova_angio)

anova_fern <- aov(log10(monoploid+7) ~ as.factor(GloNAF_incidence), genome_glonaf_fern_data)
plot(anova_fern)
summary(anova_fern)

anova_gymno <- aov(monoploid ~ as.factor(GloNAF_incidence), genome_glonaf_gymno_data)
plot(anova_gymno,2)
summary(anova_gymno)


####PLOT THE frequancy distribution, use the split violin plot  ####
# the functions as below
GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, draw_group = function(self, data, ..., draw_quantiles = NULL){
  data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
  grp <- data[1,'group']
  newdata <- plyr::arrange(transform(data, x = if(grp%%2==1) xminv else xmaxv), if(grp%%2==1) y else -y)
  newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
  newdata[c(1,nrow(newdata)-1,nrow(newdata)), 'x'] <- round(newdata[1, 'x']) 
  if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
    stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <= 
                                              1))
    quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
    aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
    aesthetics$alpha <- rep(1, nrow(quantiles))
    both <- cbind(quantiles, aesthetics)
    quantile_grob <- GeomPath$draw_panel(both, ...)
    ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
  }
  else {
    ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
  }
})

geom_split_violin <- function (mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, position = position, show.legend = show.legend, inherit.aes = inherit.aes, params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}

GenomeSize_df_final$taxno <- factor(GenomeSize_df_final$taxno, 
                                   levels = c("angiosperms","gymnosp", "ferns", "lycophytes"), ordered = TRUE)


library(ggplot2)
library(dplyr)
monoploid_freq <- GenomeSize_df_final %>% dplyr::filter (monoploid != "NA") %>% droplevels()  %>%
ggplot(aes(taxno, log10(monoploid), fill=as.factor(GloNAF_incidence))) + 
  geom_split_violin(trim = TRUE,scale = "area") + 
  geom_boxplot(width = 0.2, notch = F, outlier.shape = NA, coef=0, alpha = 0.5) +
  stat_summary(fun = mean, fun.min = mean, fun.max = mean,
               geom = "crossbar",
               width = 0.2, colour = "red", 
               position = position_dodge(width = .2)
  ) +
  labs(x=NULL,y="Monoploid genome size (pg)") +
  theme_classic() +
  theme(text = element_text(size = 20),
        legend.position = "none") +
  scale_fill_manual(values=c("#56B4E9", "#E69F00"), 
                    name="Group",
                    breaks=c("1", "0"),
                    labels=c("Naturalized", "Non-naturalized"))
pdf("FigS.monoploid GS violind plot updated.pdf", useDingbats=FALSE, width=5, height=5)

dev.off()


## test the significance between naturalized and non-naturalzied A: monoploid ########
## monoploid angiosperm   #######
monoploid_angio <- GenomeSize_df_final %>% dplyr::filter (monoploid != "NA" & taxno == "angiosperms") %>% droplevels() 
library(phytools)

keep.spp_monoploid_angio <-levels(monoploid_angio$species)
# keeptips_ext <- genome.phy$tip.label[match(keep.spp_ext, genome.phy$tip.label)]
remove_taxa_monoploid_angio = setdiff(genome.monoploid.phy$tip.label, keep.spp_monoploid_angio)
genome.monoploid.angio.phy <- drop.tip(genome.monoploid.phy,remove_taxa_monoploid_angio)

# Order data by tip order

monoploid_angio = monoploid_angio[(monoploid_angio$species %in% genome.monoploid.angio.phy$tip.label), ]

row.names(monoploid_angio) = monoploid_angio$species
#sort_by_tip.label
monoploid_angio <- monoploid_angio[match(genome.monoploid.angio.phy$tip.label, monoploid_angio$species),]
###check the names
library(geiger)  #for name.check
name.check(genome.monoploid.angio.phy, monoploid_angio)

monoploid_angio_anova <- phylANOVA(genome.monoploid.angio.phy, monoploid_angio$GloNAF_incidence,
                         monoploid_angio$monoploid,  nsim=1000, posthoc=TRUE, p.adj="bonferroni")

monoploid_angio_anova
## monoploid gymnosperm
monoploid_gymno <- GenomeSize_df_final %>% dplyr::filter (monoploid != "NA" & taxno == "gymnosp") %>% droplevels() 

keep.spp_monoploid_gymno <-levels(monoploid_gymno$species)
# keeptips_ext <- genome.phy$tip.label[match(keep.spp_ext, genome.phy$tip.label)]
remove_taxa_monoploid_gymno = setdiff(genome.monoploid.phy$tip.label, keep.spp_monoploid_gymno)
genome.monoploid.gymno.phy <- drop.tip(genome.monoploid.phy,remove_taxa_monoploid_gymno)

# Order data by tip order

monoploid_gymno = monoploid_gymno[(monoploid_gymno$species %in% genome.monoploid.gymno.phy$tip.label), ]

row.names(monoploid_gymno) = monoploid_gymno$species
#sort_by_tip.label
monoploid_gymno <- monoploid_gymno[match(genome.monoploid.gymno.phy$tip.label, monoploid_gymno$species),]
###check the names
library(geiger)  #for name.check
name.check(genome.monoploid.gymno.phy, monoploid_gymno)

monoploid_gymno_anova <- phylANOVA(genome.monoploid.gymno.phy, monoploid_gymno$GloNAF_incidence,
                                   monoploid_gymno$monoploid,  nsim=1000, posthoc=TRUE, p.adj="bonferroni")

monoploid_gymno_anova


####  ends of monoploid anova test

#### test the significance between naturalized and non-naturalzied B: hoploid    #####
holoploid_angio <- GenomeSize_df_final %>% dplyr::filter (holoploid != "NA" & taxno == "angiosperms") %>% droplevels() 
library(phytools)

keep.spp_holoploid_angio <-levels(holoploid_angio$species)
# keeptips_ext <- genome.phy$tip.label[match(keep.spp_ext, genome.phy$tip.label)]
remove_taxa_holoploid_angio = setdiff(genome.holoploid.phy$tip.label, keep.spp_holoploid_angio)
genome.holoploid.angio.phy <- drop.tip(genome.holoploid.phy,remove_taxa_holoploid_angio)

# Order data by tip order

holoploid_angio = holoploid_angio[(holoploid_angio$species %in% genome.holoploid.angio.phy$tip.label), ]

row.names(holoploid_angio) = holoploid_angio$species
#sort_by_tip.label
holoploid_angio <- holoploid_angio[match(genome.holoploid.angio.phy$tip.label, holoploid_angio$species),]
###check the names
library(geiger)  #for name.check
name.check(genome.holoploid.angio.phy, holoploid_angio)

holoploid_angio_anova <- phylANOVA(genome.holoploid.angio.phy, holoploid_angio$GloNAF_incidence,
                                   holoploid_angio$holoploid,  nsim=1000, posthoc=TRUE, p.adj="bonferroni")

holoploid_angio_anova
## holoploid gymnosperm
holoploid_gymno <- GenomeSize_df_final %>% dplyr::filter (holoploid != "NA" & taxno == "gymnosp") %>% droplevels() 

keep.spp_holoploid_gymno <-levels(holoploid_gymno$species)
# keeptips_ext <- genome.phy$tip.label[match(keep.spp_ext, genome.phy$tip.label)]
remove_taxa_holoploid_gymno = setdiff(genome.holoploid.phy$tip.label, keep.spp_holoploid_gymno)
genome.holoploid.gymno.phy <- drop.tip(genome.holoploid.phy,remove_taxa_holoploid_gymno)

# Order data by tip order

holoploid_gymno = holoploid_gymno[(holoploid_gymno$species %in% genome.holoploid.gymno.phy$tip.label), ]

row.names(holoploid_gymno) = holoploid_gymno$species
#sort_by_tip.label
holoploid_gymno <- holoploid_gymno[match(genome.holoploid.gymno.phy$tip.label, holoploid_gymno$species),]
###check the names
library(geiger)  #for name.check
name.check(genome.holoploid.gymno.phy, holoploid_gymno)

holoploid_gymno_anova <- phylANOVA(genome.holoploid.gymno.phy, holoploid_gymno$GloNAF_incidence,
                                   holoploid_gymno$holoploid,  nsim=1000, posthoc=TRUE, p.adj="bonferroni")

holoploid_gymno_anova

### ferns haploid
holoploid_fern <- GenomeSize_df_final %>% dplyr::filter (holoploid != "NA" & taxno == "ferns") %>% droplevels() 

keep.spp_holoploid_fern <-levels(holoploid_fern$species)
# keeptips_ext <- genome.phy$tip.label[match(keep.spp_ext, genome.phy$tip.label)]
remove_taxa_holoploid_fern = setdiff(genome.holoploid.phy$tip.label, keep.spp_holoploid_fern)
genome.holoploid.fern.phy <- drop.tip(genome.holoploid.phy,remove_taxa_holoploid_fern)

# Order data by tip order

holoploid_fern = holoploid_fern[(holoploid_fern$species %in% genome.holoploid.fern.phy$tip.label), ]

row.names(holoploid_fern) = holoploid_fern$species
#sort_by_tip.label
holoploid_fern <- holoploid_fern[match(genome.holoploid.fern.phy$tip.label, holoploid_fern$species),]
###check the names
library(geiger)  #for name.check
name.check(genome.holoploid.fern.phy, holoploid_fern)

holoploid_fern_anova <- phylANOVA(genome.holoploid.fern.phy, holoploid_fern$GloNAF_incidence,
                                  holoploid_fern$holoploid,  nsim=1000, posthoc=TRUE, p.adj="bonferroni")

holoploid_fern_anova

### lycophytes haploid
holoploid_lycophytes <- GenomeSize_df_final %>% dplyr::filter (holoploid != "NA" & taxno == "lycophytes") %>% droplevels() 

keep.spp_holoploid_lycophytes <-levels(holoploid_lycophytes$species)
# keeptips_ext <- genome.phy$tip.label[match(keep.spp_ext, genome.phy$tip.label)]
remove_taxa_holoploid_lycophytes = setdiff(genome.holoploid.phy$tip.label, keep.spp_holoploid_lycophytes)
genome.holoploid.lycophytes.phy <- drop.tip(genome.holoploid.phy,remove_taxa_holoploid_lycophytes)

# Order data by tip order

holoploid_lycophytes = holoploid_lycophytes[(holoploid_lycophytes$species %in% genome.holoploid.lycophytes.phy$tip.label), ]

row.names(holoploid_lycophytes) = holoploid_lycophytes$species
#sort_by_tip.label
holoploid_lycophytes <- holoploid_lycophytes[match(genome.holoploid.lycophytes.phy$tip.label, holoploid_lycophytes$species),]
###check the names
library(geiger)  #for name.check
name.check(genome.holoploid.lycophytes.phy, holoploid_lycophytes)

holoploid_lycophytes_anova <- phylANOVA(genome.holoploid.lycophytes.phy, holoploid_lycophytes$GloNAF_incidence,
                                        holoploid_lycophytes$holoploid,  nsim=1000, posthoc=TRUE, p.adj="bonferroni")

holoploid_lycophytes_anova

####################### frequency and violin plot end here  #######

### plot the relationship between monoploid and holoploid
ggplot(GenomeSize_df_final, aes(x= monoploid, y= holoploid)) +
  geom_smooth(method=lm, se= T) +
  geom_point(aes(color = PLOIDY, shape = PLOIDY), position = position_jitter(w = 0.01, h = 0.09)) + 
  #scale_color_manual(values =c('#999999','#E69F00', '#56B4E9'))+
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"))+
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     strip.background = element_rect(colour="white", fill="white"),
                     text = element_text(size = 16),
                     legend.position = c(0.9,0.8)) +
  labs(title="",
       x="Monoploid genome size (pg)", y = "Holoploid genome size (pg)")

pdf("Fig.S holoploid GS vs. monoploid GS-16-062021.pdf", useDingbats=FALSE, width=6, height=6)


dev.off()
summary(lm(GenomeSize_df_monoploid_final$holoploid ~ GenomeSize_df_monoploid_final$monoploid))

###########  test the monoploid gs and naturalization success  #############
GenomeSize_df_monoploid_final <- GenomeSize_df_final %>% filter (monoploid != "NA") %>% droplevels() #7252 species

## prepare the monoploid data tree
GenomeSize_df_monoploid_name<-levels(GenomeSize_df_monoploid_final$species)
remove_monoploid_taxa = setdiff(genome.phy$tip.label, GenomeSize_df_monoploid_name)
genome.monoploid.phy <- drop.tip(genome.phy,remove_monoploid_taxa)



### regression analyses only globally ######
### As suggested by Mark and others, I will run analysis using all speies together, that is,
###  not seperate them into the above three taxonomies
## running a normal glm with a quadratic term of genome size, as suggested by Mark

# Order data by tip order

GenomeSize_df_monoploid_final = GenomeSize_df_monoploid_final[(GenomeSize_df_monoploid_final$species %in% genome.monoploid.phy$tip.label), ]
row.names(GenomeSize_df_monoploid_final) = GenomeSize_df_monoploid_final$species
#sort_by_tip.label
GenomeSize_df_monoploid_final <- GenomeSize_df_monoploid_final[match(genome.monoploid.phy$tip.label, GenomeSize_df_monoploid_final$species),]
###check the names
library(geiger)  #for name.check
name.check(genome.monoploid.phy, GenomeSize_df_monoploid_final )


library(phylolm)
range_1.1 <- function(x){2*((x-min(x))/(max(x)-min(x)))-1}  ###stardardozed to -1 and 1


####### IA monoploid vs. incidence starts here ######

genome.incidence.null.model <- phyloglm(GloNAF_incidence ~ 1, phy = genome.monoploid.phy, data = GenomeSize_df_monoploid_final, 
                                        method = c("logistic_MPLE","logistic_IG10","poisson_GEE"),
                                        btol = 100, log.alpha.bound = 10,
                                        start.beta=NULL, start.alpha=NULL,
                                        boot = 0, full.matrix = TRUE)
summary(genome.incidence.null.model)

genome.global.monoploid.model <- phyloglm(GloNAF_incidence ~ range_1.1(log10(monoploid)), 
                                          phy = genome.monoploid.phy, data = GenomeSize_df_monoploid_final, 
                                          method = c("logistic_MPLE","logistic_IG10","poisson_GEE"),
                                          btol = 10, log.alpha.bound = 4,
                                          start.beta=NULL, start.alpha=NULL,
                                          boot = 0, full.matrix = TRUE)
summary(genome.global.monoploid.model)
library(rr2)
R2(genome.global.monoploid.model, genome.incidence.null.model) #0.0885



###plot the estimated slope
plot(range_1.1(log10(GenomeSize_df_monoploid_final$monoploid)),
     jitter(GenomeSize_df_monoploid_final$GloNAF_incidence,factor=0,amount=0.02),
     xlab="monoploid",ylab="Naturalization incidence",xlim=c(-1,1))
cc.global.genome_mono <- coef(genome.global.monoploid.model)
curve(plogis(cc.global.genome_mono[1]+ cc.global.genome_mono[2]*x),col="red",
      lwd=2,lty=1,add=TRUE)


## using glm model, without consideration of the phylogeny
ggplot(data=GenomeSize_df_monoploid_final,aes(range_1.1(log10(monoploid)), GloNAF_incidence)) + geom_point() + 
  geom_smooth(method = "lm", formula = y ~ poly(x, 2))

ggplot(GenomeSize_df_monoploid_final,aes(x= range_1.1(log10(monoploid)), y=GloNAF_incidence)) +
  geom_point(aes(color = taxno), shape = 21,position = position_jitter(w = 0.01, h = 0.09)) + 
  geom_smooth(method = "lm", formula = y ~ poly(x, 2))


my.model <- glm(GloNAF_incidence ~ range_1.1(log10(monoploid)),
                data=GenomeSize_df_monoploid_final, family="binomial")

summary(my.model)


my.q.model <-glm(GloNAF_incidence ~ range_1.1(log10(monoploid))+I(range_1.1(log10(monoploid))^2), 
                 data=GenomeSize_df_monoploid_final, family="binomial")
summary(my.q.model)



####quadratic  model
genome.phyloglm.q.model <- phyloglm(GloNAF_incidence ~ range_1.1(log10(monoploid)) +
                                      I(range_1.1(log10(monoploid))^2), 
                                    phy = genome.monoploid.phy, data = GenomeSize_df_monoploid_final, 
                                    method = c("logistic_MPLE","logistic_IG10","poisson_GEE"),
                                    btol = 10, log.alpha.bound = 4,
                                    start.beta=NULL, start.alpha=NULL,
                                    boot = 0, full.matrix = TRUE)

summary(genome.phyloglm.q.model)
summary(genome.global.monoploid.model)
R2(genome.phyloglm.q.model, genome.incidence.null.model) #0.081


## prepare new data to plot the prediction
New_monoploid <- data.frame( monoploid = rep(seq(from = 0.06633,to = 152.2, length = 5000),6))


X <- model.matrix(~ range_1.1(log10(monoploid)) +
                    I(range_1.1(log10(monoploid))^2), data = New_monoploid)
X1 <- model.matrix(~ range_1.1(log10(monoploid)), data = New_monoploid)


New_monoploid$Pred <- X %*% coefficients(genome.phyloglm.q.model)#model coefficients
New_monoploid$Pred.linear <- X1 %*% coefficients(genome.global.monoploid.model)#model coefficients

New_monoploid$SE <- sqrt(  diag(X %*%vcov(genome.phyloglm.q.model) %*% t(X))  )
library(viridis)
p_nat_incidence <- ggplot(GenomeSize_df_monoploid_final,aes(x= range_1.1(log10(monoploid)), y=GloNAF_incidence)) +
  geom_point(aes(color = taxno, shape = taxno), position = position_jitter(w = 0.01, h = 0.09)) + 
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"))+
   scale_shape_manual(values=c(3, 16, 17))+
   geom_line(data = New_monoploid, aes(x= range_1.1(log10(monoploid)), y = exp(Pred)/ (1+exp(Pred))),size=1, color = "red") +
   geom_line(data = New_monoploid, aes(x= range_1.1(log10(monoploid)), y = exp(Pred.linear)/ (1+exp(Pred.linear))),size=1) +
  theme_classic() + 
  coord_cartesian(ylim = c(0,1))+ # does not remove points
   labs(y="Naturalization probablity",x="Monoploid genome size (pg)") +
  theme(legend.position = c(0.85, 0.6)) +
  theme(axis.text=element_text(size=12,face="bold"),
        axis.title=element_text(size=14,face="bold")) +
  theme(plot.title = element_text(size = 12, face = "bold"),
        legend.title=element_text(size=14, face = "bold"), 
        legend.text=element_text(size=12, face = "bold")) +
  ggtitle("(a)")
####### IA monoploid vs. incidence  ends here 

##### run the naturalization extent and invasive extent 
####### IB monoploid vs. naturalization extent  starts here #########

library(aod)

##### IA:2a naturalization extent
## prepare the data
library(dplyr)
genome.monoploid.extent.data <- GenomeSize_df_monoploid_final %>% filter(reg_nat_total != "NA" ) %>% droplevels()
summary(genome.monoploid.extent.data)
str(genome.monoploid.extent.data) ##3631 species 2364 species for the new data
## prepare the extent phylogeny
#But get the tips to keep
keep.spp_monoploid_ext <-levels(genome.monoploid.extent.data$species)
# keeptips_ext <- genome.phy$tip.label[match(keep.spp_ext, genome.phy$tip.label)]
remove_taxa_monoploid_ext = setdiff(genome.phy$tip.label, keep.spp_monoploid_ext)
genome.monoploid.extent.phy <- drop.tip(genome.phy,remove_taxa_monoploid_ext)

# Order data by tip order

genome.monoploid.extent.data = genome.monoploid.extent.data[(genome.monoploid.extent.data$species %in% genome.monoploid.extent.phy$tip.label), ]

row.names(genome.monoploid.extent.data) = genome.monoploid.extent.data$species
#sort_by_tip.label
genome.monoploid.extent.data <- genome.monoploid.extent.data[match(genome.monoploid.extent.phy$tip.label, genome.monoploid.extent.data$species),]
###check the names
library(geiger)  #for name.check
name.check(genome.monoploid.extent.phy, genome.monoploid.extent.data)

str(genome.monoploid.extent.data)


## monoploid * taxomony
genome.null.extent.model <- phylolm(range_1.1(log10(reg_nat_total)) ~ 1,
                                    phy = genome.monoploid.extent.phy, data = genome.monoploid.extent.data, 
                                    method = "lambda")
summary(genome.null.extent.model)

genome.monoploid.extent.model <- phylolm(range_1.1(log10(reg_nat_total)) ~ range_1.1(log10(monoploid)),
                                         phy = genome.monoploid.extent.phy, data = genome.monoploid.extent.data,
                                         method = "lambda")

summary(genome.monoploid.extent.model)
# R2(genome.monoploid.extent.model, genome.null.extent.model, genome.extent.phy)
# R2(genome.monoploid.extent.model,  genome.extent.phy)
plot(genome.monoploid.extent.model)
## model Residual Plot
qplot(range_1.1(log10(genome.monoploid.extent.data$reg_nat_total)), residuals(genome.monoploid.extent.model))
hist(residuals(genome.monoploid.extent.model))
plot(density(resid(genome.monoploid.extent.model))) #A density plot
# A quantile normal plot - good for checking normality
qqnorm(resid(genome.monoploid.extent.model)) 
qqline(resid(genome.monoploid.extent.model))

####quadratic  model
genome.monoploid.extent.glm.model <- glm(range_1.1(log10(reg_nat_total)) ~ range_1.1(log10(monoploid)), 
                                           data = genome.monoploid.extent.data)

summary(genome.monoploid.extent.glm.model)
plot(genome.monoploid.extent.glm.model)

genome.monoploid.extent.q.glm.model <- glm(range_1.1(log10(reg_nat_total)) ~ range_1.1(log10(monoploid)) +
                                             I(range_1.1(log10(monoploid))^2), 
                                            data = genome.monoploid.extent.data)

summary(genome.monoploid.extent.q.glm.model)
plot(genome.monoploid.extent.q.glm.model)



genome.monoploid.extent.q.model <- phylolm(range_1.1(log10(reg_nat_total)) ~ range_1.1(log10(monoploid)) +
                                           I(range_1.1(log10(monoploid))^2), 
                                           phy = genome.monoploid.extent.phy, data = genome.monoploid.extent.data,
                                         method = "lambda")

summary(genome.monoploid.extent.q.model)
plot(genome.monoploid.extent.q.model)
## model Residual Plot
qplot(range_1.1(log10(genome.monoploid.extent.data$reg_nat_total)), residuals(genome.monoploid.extent.q.model))
hist(residuals(genome.monoploid.extent.q.model))
plot(density(resid(genome.monoploid.extent.q.model))) #A density plot
qqnorm(resid(genome.monoploid.extent.q.model)) 
qqline(resid(genome.monoploid.extent.q.model))


plot(range_1.1(log10(genome.monoploid.extent.data$monoploid)),
     jitter(range_1.1(log10(genome.monoploid.extent.data$reg_nat_total)),factor=0,amount=0.02),
     xlab="monoploid",ylab="no.of extent regions",xlim=c(-1,1))
cc.genome_nat_ext <- coef(genome.monoploid.extent.model)
curve(plogis(cc.genome_nat_ext[1]+ cc.genome_nat_ext[2]*x),col="red",
      lwd=2,lty=1,add=TRUE)

cc.genome_nat_ext.q <- coef(genome.monoploid.extent.q.model)
curve(plogis(cc.genome_nat_ext.q[1]+ cc.genome_nat_ext.q[2]*x + cc.genome_nat_ext.q[3]*x^2 ),col="red",
      lwd=2,lty=1,add=TRUE)

##plot

X1 <- model.matrix(~ range_1.1(log10(monoploid)), data = New_monoploid)
New_monoploid$Pred.ext <- X1 %*% coefficients(genome.monoploid.extent.model)#model coefficients

X2 <- model.matrix(~ range_1.1(log10(monoploid)) +
                     I(range_1.1(log10(monoploid))^2), data = New_monoploid)
New_monoploid$Pred.q.ext <- X2 %*% coefficients(genome.monoploid.extent.q.model)#model coefficients


New_monoploid$SE.ext <- sqrt(diag(X1 %*%vcov(genome.monoploid.extent.model) %*% t(X1)))
library(viridis)
p_naturalization_extent <- ggplot(genome.monoploid.extent.data,aes(x= range_1.1(log10(monoploid)), y= range_1.1(log10(reg_nat_total)))) +
  geom_point(aes(color = taxno, shape = taxno),position = position_jitter(w = 0.01, h = 0.09)) + 
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"))+

  scale_shape_manual(values=c(3, 16, 17))+
   geom_line(data = New_monoploid, aes(x= range_1.1(log10(monoploid)), y = Pred.ext),size=1) +
  geom_line(data = New_monoploid, aes(x= range_1.1(log10(monoploid)), y = Pred.q.ext),size=1,  color = "red") +
      theme_classic() + 
  coord_cartesian(ylim = c(-1,1))+ # does not remove points
   labs(y="Naturalization extent (Number of regions)",x="Monoploid genome size (pg)") +
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=12,face="bold"),
        axis.title=element_text(size=14,face="bold")) +
  theme(plot.title = element_text(size = 12, face = "bold"),
        legend.title=element_text(size=14, face = "bold"), 
        legend.text=element_text(size=12, face = "bold")) +
  ggtitle("(b)")




##### IC monoploid vs. invasive extent  starts here ######

genome.global.extent.invasive.data <- GenomeSize_df_monoploid_final %>%  filter (reg_no_invasive >0 ) %>% droplevels()
summary(genome.global.extent.invasive.data)
str(genome.global.extent.invasive.data) ##1060 species 680 species in the new data
## prepare the extent phylogeny
#But get the tips to keep
keep.spp_monoploid_invext <-levels(genome.global.extent.invasive.data$species)

# keeptips_invext <- genome.phy$tip.label[match(keep.spp_invext, genome.phy$tip.label)]

remove_taxa_monoploid_invext = setdiff(genome.phy$tip.label, keep.spp_monoploid_invext)
genome.extent.invasive.monoploid.phy <- drop.tip(genome.phy,remove_taxa_monoploid_invext)

# Order data by tip order

genome.global.extent.invasive.data = genome.global.extent.invasive.data[(genome.global.extent.invasive.data$species %in% genome.extent.invasive.monoploid.phy$tip.label), ]

row.names(genome.global.extent.invasive.data) = genome.global.extent.invasive.data$species
#sort_by_tip.label
genome.global.extent.invasive.data <- genome.global.extent.invasive.data[match(genome.extent.invasive.monoploid.phy$tip.label, genome.global.extent.invasive.data$species),]
##check the names
library(geiger)  #for name.check
name.check(genome.extent.invasive.monoploid.phy, genome.global.extent.invasive.data)

str(genome.global.extent.invasive.data)
# 
## monoploid * taxomony
genome.null.invasive.extent.model <- phylolm(range_1.1(log10(reg_no_invasive)) ~ 1,
                                             phy = genome.extent.invasive.monoploid.phy, data = genome.global.extent.invasive.data, 
                                             method = "lambda")
summary(genome.null.invasive.extent.model)

genome.monoploid.invasive.extent.model <- phylolm(range_1.1(log10(reg_no_invasive)) ~ range_1.1(log10(monoploid)),
                                                  phy = genome.extent.invasive.monoploid.phy, data = genome.global.extent.invasive.data,
                                                  method = "lambda")

summary(genome.monoploid.invasive.extent.model)

# R2(genome.monoploid.extent.model, phy = genome.extent.phy)
plot(genome.monoploid.invasive.extent.model)
## model Residual Plot
qplot(range_1.1(log10(genome.global.extent.invasive.data$reg_no_invasive)), residuals(genome.monoploid.invasive.extent.model))
hist(residuals(genome.monoploid.invasive.extent.model))
plot(density(resid(genome.monoploid.invasive.extent.model))) #A density plot
# A quantile normal plot - good for checking normality
qqnorm(resid(genome.monoploid.invasive.extent.model)) 
qqline(resid(genome.monoploid.invasive.extent.model))


###plot the estimated slope
plot(range_1.1(log10(genome.global.extent.invasive.data$monoploid)),
     jitter(range_1.1(log10(genome.global.extent.invasive.data$reg_no_invasive)),factor=0,amount=0.02),
     xlab="monoploid",ylab="no.of invasive regions",xlim=c(-1,1))
cc.genome_inva_ext <- coef(genome.monoploid.invasive.extent.model)
curve(plogis(cc.genome_inva_ext[1]+ cc.genome_inva_ext[2]*x),col="red",
      lwd=2,lty=1,add=TRUE)
cc.genome_inv_ext.q <- coef(genome.monoploid.extent.invasive.q.model)
curve(plogis(cc.genome_inv_ext.q[1]+ cc.genome_inv_ext.q[2]*x + cc.genome_inv_ext.q[3]*x^2 ),col="red",
      lwd=2,lty=1,add=TRUE)



genome.monoploid.extent.invasive.glm.model <- glm(range_1.1(log10(reg_no_invasive)) ~ range_1.1(log10(monoploid)),                                                      , 
                                                    data = genome.global.extent.invasive.data)
summary(genome.monoploid.extent.invasive.glm.model)
####quadratic  model

genome.monoploid.extent.invasive.q.glm.model <- glm(range_1.1(log10(reg_no_invasive)) ~ range_1.1(log10(monoploid)) +
                                             I(range_1.1(log10(monoploid))^2), 
                                           data = genome.global.extent.invasive.data)

summary(genome.monoploid.extent.invasive.q.glm.model)


genome.monoploid.extent.invasive.q.model <- phylolm(range_1.1(log10(reg_no_invasive)) ~ range_1.1(log10(monoploid)) +
                                             I(range_1.1(log10(monoploid))^2), 
                                             phy = genome.extent.invasive.monoploid.phy, data = genome.global.extent.invasive.data,
                                           method = "lambda")

summary(genome.monoploid.extent.invasive.q.model)
plot(genome.monoploid.extent.invasive.q.model)
## model Residual Plot
qplot(range_1.1(log10(genome.global.extent.invasive.data$reg_no_invasive)), residuals(genome.monoploid.extent.invasive.q.model))
hist(residuals(genome.monoploid.extent.invasive.q.model))
plot(density(resid(genome.monoploid.extent.invasive.q.model))) #A density plot
# A quantile normal plot - good for checking normality
qqnorm(resid(genome.monoploid.extent.invasive.q.model)) 
qqline(resid(genome.monoploid.extent.invasive.q.model))


###plot the estimated slope
plot(range_1.1(log10(genome.global.extent.invasive.data$monoploid)),
     jitter(range_1.1(log10(genome.global.extent.invasive.data$reg_no_invasive)),factor=0,amount=0.02),
     xlab="monoploid",ylab="no.of invasive regions",xlim=c(-1,1))
cc.genome_inva_q.ext <- coef(genome.monoploid.extent.invasive.q.model)
curve(plogis(cc.genome_inva_q.ext[1]+ cc.genome_inva_q.ext[2]*x+ cc.genome_inva_q.ext[3]*x^2 ),col="red",
      lwd=2,lty=1,add=TRUE)

##plot
X1 <- model.matrix(~ range_1.1(log10(monoploid)), data = New_monoploid)
New_monoploid$Pred.ext.inv <- X1 %*% coefficients(genome.monoploid.invasive.extent.model)#model coefficients
X2 <- model.matrix(~ range_1.1(log10(monoploid)) +
                     I(range_1.1(log10(monoploid))^2), data = New_monoploid)
New_monoploid$Pred.q.ext.inv <- X2  %*% coefficients(genome.monoploid.extent.invasive.q.model)#model coefficients


New_monoploid$SE.ext.inv <- sqrt(diag(X1 %*%vcov(genome.monoploid.invasive.extent.model) %*% t(X1)))
##plot
p_invasive_extent <- ggplot(genome.global.extent.invasive.data,
                            aes(x= range_1.1(log10(monoploid)), y= range_1.1(log10(reg_no_invasive)))) +
  geom_point(aes(color = taxno, shape = taxno),position = position_jitter(w = 0.01, h = 0.09)) + 
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"))+
  scale_shape_manual(values=c(3, 16, 17))+
  geom_line(data = New_monoploid, aes(x= range_1.1(log10(monoploid)), y = Pred.ext.inv),size=1) +
  geom_line(data = New_monoploid, aes(x= range_1.1(log10(monoploid)), y = Pred.q.ext.inv),size=1, color = "red") +
  theme_classic() + 
  coord_cartesian(ylim = c(-1,1))+ # does not remove points
  labs(y="Invasion extent (Number of regions)",x="Monoploid genome size (pg)") +
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=12,face="bold"),
        axis.title=element_text(size=14,face="bold")) +
  theme(plot.title = element_text(size = 12, face = "bold"),
        legend.title=element_text(size=14, face = "bold"), 
        legend.text=element_text(size=12, face = "bold")) +
  ggtitle("(c)")

library(ggpubr)

figure2 <-ggarrange(p_nat_incidence,                                                 # First row with incidence plot
                    ggarrange(p_naturalization_extent, p_invasive_extent, ncol = 2, labels = c("B", "C")), # Second row with two extent plot
                    nrow = 2, heights = c(0.8, 1),  ## adjust height
                    labels = "A"  )                                      # Labels of the incidence plot


pdf("Fig 2.monoploid GS incidence & extent_noSE_linear_qua.pdf", useDingbats=FALSE, width=12, height=9)
figure2
dev.off()

#######  plot the results with ploidy level shown ################
p_nat_incidence_ploidy <- ggplot(GenomeSize_df_monoploid_final,aes(x= range_1.1(log10(monoploid)), y=GloNAF_incidence)) +
  geom_point(aes(color = PLOIDY, shape = taxno), position = position_jitter(w = 0.01, h = 0.09)) + 
  scale_color_viridis(discrete=TRUE, direction = 1) +
  scale_shape_manual(values=c(3, 16, 17))+
  geom_line(data = New_monoploid, aes(x= range_1.1(log10(monoploid)), y = exp(Pred)/ (1+exp(Pred))),size=1, color = "red") +
  geom_line(data = New_monoploid, aes(x= range_1.1(log10(monoploid)), y = exp(Pred.linear)/ (1+exp(Pred.linear))),size=1) +
  theme_classic() + 
  coord_cartesian(ylim = c(0,1))+ # does not remove points
  labs(y="Naturalization probablity",x="Monoploid genome size (pg)") +
  theme(legend.position = c(0.85, 0.6)) +
  theme(axis.text=element_text(size=12,face="bold"),
        axis.title=element_text(size=14,face="bold")) +
  theme(plot.title = element_text(size = 12, face = "bold"),
        legend.title=element_text(size=14, face = "bold"), 
        legend.text=element_text(size=12, face = "bold")) +
  ggtitle("(a)")


p_naturalization_extent_ploidy <- ggplot(genome.monoploid.extent.data,aes(x= range_1.1(log10(monoploid)), y= range_1.1(log10(reg_nat_total)))) +
  geom_point(aes(color = PLOIDY, shape = taxno),position = position_jitter(w = 0.01, h = 0.09)) + 
  scale_color_viridis(discrete=TRUE, direction = 1) +
  scale_shape_manual(values=c(3, 16, 17))+
  geom_line(data = New_monoploid, aes(x= range_1.1(log10(monoploid)), y = Pred.ext),size=1) +
  geom_line(data = New_monoploid, aes(x= range_1.1(log10(monoploid)), y = Pred.q.ext),size=1,  color = "red") +
  theme_classic() + 
  coord_cartesian(ylim = c(-1,1))+ # does not remove points
  labs(y="Naturalization extent (Number of regions)",x="Monoploid genome size (pg)") +
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=12,face="bold"),
        axis.title=element_text(size=14,face="bold")) +
  theme(plot.title = element_text(size = 12, face = "bold"),
        legend.title=element_text(size=14, face = "bold"), 
        legend.text=element_text(size=12, face = "bold")) +
  ggtitle("(b)")


p_invasive_extent_ploidy <- ggplot(genome.global.extent.invasive.data,
                            aes(x= range_1.1(log10(monoploid)), y= range_1.1(log10(reg_no_invasive)))) +
  geom_point(aes(color = PLOIDY, shape = taxno),position = position_jitter(w = 0.01, h = 0.09)) + 
   scale_color_viridis(discrete=TRUE, direction = 1) +
  scale_shape_manual(values=c(3, 16, 17))+
  geom_line(data = New_monoploid, aes(x= range_1.1(log10(monoploid)), y = Pred.ext.inv),size=1) +
  geom_line(data = New_monoploid, aes(x= range_1.1(log10(monoploid)), y = Pred.q.ext.inv),size=1, color = "red") +
  theme_classic() + 
  coord_cartesian(ylim = c(-1,1))+ # does not remove points
  labs(y="Invasion extent (Number of regions)",x="Monoploid genome size (pg)") +
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=12,face="bold"),
        axis.title=element_text(size=14,face="bold")) +
  theme(plot.title = element_text(size = 12, face = "bold"),
        legend.title=element_text(size=14, face = "bold"), 
        legend.text=element_text(size=12, face = "bold")) +
  ggtitle("(c)")

figure2_ploidy <-ggarrange(p_nat_incidence_ploidy,                                                 # First row with incidence plot
                    ggarrange(p_naturalization_extent_ploidy, p_invasive_extent_ploidy, ncol = 2, labels = c("B", "C")), # Second row with two extent plot
                    nrow = 2, heights = c(0.8, 1),  ## adjust height
                    labels = "A"  )  


pdf("Fig 2.monoploid GS incidence & extent_noSE_ploidy_linear_qua_22Aug2020.pdf", useDingbats=FALSE, width=12, height=9)
figure2_ploidy
dev.off()
### monoploid ends here

#####II holoploid  #####

#### holoploid plot #####


holoploid_freq <- GenomeSize_df_final %>% filter (holoploid != "NA") %>% droplevels()  %>%
  ggplot(aes(taxno, log10(holoploid), fill=as.factor(GloNAF_incidence))) + 
  geom_split_violin(trim = TRUE,scale = "area") + 
  geom_boxplot(width = 0.2, notch = F, outlier.shape = NA, coef=0, alpha = 0.5) +
  stat_summary(fun = mean, fun.min = mean, fun.max = mean,
               geom = "crossbar",
               width = 0.2, colour = "red", 
               position = position_dodge(width = .2)
  ) +
  labs(x=NULL,y="Holoploid genome size (pg)") +
  theme_classic() +
  theme(text = element_text(size = 20),
        legend.position = c(0.6, 0.2)) +
  scale_fill_manual(values=c("#56B4E9", "#E69F00"), 
                    name="Group",
                    breaks=c("1", "0"),
                    labels=c("Naturalized", "Non-naturalized"))


monoploid_freq

library(ggpubr)

figure1_freq <-ggarrange(monoploid_freq,  holoploid_freq, ncol = 2, widths = c(1, 2), labels = c("A", "B")) # Second row with two extent plot
                         

pdf("Fig1.monoploid & holoploid GS violind plot updated 22Aug2021.pdf", useDingbats=FALSE, width=12, height=6)

dev.off()


### plots ends here

#### holoploid GS vs. naturalization success  starts here #######
GenomeSize_df_final_holoploid <- GenomeSize_df_final %>% filter (holoploid != "NA") %>% droplevels()  

## prepare the holoploid data tree
GenomeSize_df_final_name<-levels(GenomeSize_df_final_holoploid$species)
remove_holoploid_taxa = setdiff(genome.phy$tip.label, GenomeSize_df_final_name)
genome.holoploid.phy <- drop.tip(genome.phy,remove_holoploid_taxa)

# Order data by tip order
GenomeSize_df_holoploid_final <- GenomeSize_df_final_holoploid
GenomeSize_df_holoploid_final = GenomeSize_df_holoploid_final[(GenomeSize_df_holoploid_final$species %in% genome.holoploid.phy$tip.label), ]
row.names(GenomeSize_df_holoploid_final) = GenomeSize_df_holoploid_final$species
#sort_by_tip.label
GenomeSize_df_holoploid_final <- GenomeSize_df_holoploid_final[match(genome.holoploid.phy$tip.label, GenomeSize_df_holoploid_final$species),]
###check the names
library(geiger)  #for name.check
name.check(genome.holoploid.phy, GenomeSize_df_holoploid_final )

write.tree(genome.holoploid.phy, file ="./genome.holoploid.tre")
write.csv(GenomeSize_df_holoploid_final, file ="./GenomeSize_df_holoploid_final.csv")

## holoploid vs naturalization incidence  #############
genome.incidence.holoploid.null.model <- phyloglm(GloNAF_incidence ~ 1, phy = genome.holoploid.phy,
                                                  data = GenomeSize_df_holoploid_final, 
                                        method = c("logistic_MPLE","logistic_IG10","poisson_GEE"),
                                        btol = 100, log.alpha.bound = 10,
                                        start.beta=NULL, start.alpha=NULL,
                                        boot = 0, full.matrix = TRUE)
summary(genome.incidence.holoploid.null.model)

genome.global.holoploid.model <- phyloglm(GloNAF_incidence ~ range_1.1(log10(holoploid)), 
                                          phy = genome.holoploid.phy, data = GenomeSize_df_holoploid_final, 
                                          method = c("logistic_MPLE","logistic_IG10","poisson_GEE"),
                                          btol = 10, log.alpha.bound = 4,
                                          start.beta=NULL, start.alpha=NULL,
                                          boot = 0, full.matrix = TRUE)
summary(genome.global.holoploid.model)
library(rr2)
R2(genome.global.holoploid.model, genome.incidence.holoploid.null.model) #0.067



###plot the estimated slope
plot(range_1.1(log10(GenomeSize_df_holoploid_final$holoploid)),
     jitter(GenomeSize_df_holoploid_final$GloNAF_incidence,factor=0,amount=0.02),
     xlab="holoploid",ylab="Naturalization incidence",xlim=c(-1,1))
cc.global.genome_holo <- coef(genome.global.holoploid.model)
curve(plogis(cc.global.genome_holo[1]+ cc.global.genome_holo[2]*x),col="red",
      lwd=2,lty=1,add=TRUE)


####quadratic  model
genome.phyloglm.q.holo.model <- phyloglm(GloNAF_incidence ~ range_1.1(log10(holoploid)) +
                                      I(range_1.1(log10(holoploid))^2), 
                                    phy = genome.holoploid.phy, data = GenomeSize_df_holoploid_final, 
                                    method = c("logistic_MPLE","logistic_IG10","poisson_GEE"),
                                    btol = 10, log.alpha.bound = 4,
                                    start.beta=NULL, start.alpha=NULL,
                                    boot = 0, full.matrix = TRUE)

summary(genome.phyloglm.q.holo.model)
summary(genome.global.holoploid.model)
R2(genome.phyloglm.q.holo.model, genome.incidence.holoploid.null.model) #0.067


## prepare new data to plot the prediction
New_hololoid <- data.frame( holoploid = rep(seq(from = 0.1326,to = 304.40, length = 5000),6))


X.holo <- model.matrix(~ range_1.1(log10(holoploid)) +
                    I(range_1.1(log10(holoploid))^2), data = New_hololoid)
X1.holo <- model.matrix(~ range_1.1(log10(holoploid)), data = New_hololoid)


New_hololoid$Pred <- X.holo %*% coefficients(genome.phyloglm.q.holo.model)#model coefficients
New_hololoid$Pred.linear <- X1.holo %*% coefficients(genome.global.holoploid.model)#model coefficients

New_hololoid$SE <- sqrt(diag(X.holo %*%vcov(genome.phyloglm.q.holo.model) %*% t(X.holo))  )
library(viridis)
p_nat_holo_incidence <- ggplot(GenomeSize_df_holoploid_final,aes(x= range_1.1(log10(holoploid)), y=GloNAF_incidence)) +
  geom_point(aes(color = PLOIDY, shape = taxno), position = position_jitter(w = 0.01, h = 0.09)) + 
  scale_color_viridis(discrete=TRUE, direction = 1) +
  scale_shape_manual(values=c(3, 16, 17,5))+
  geom_line(data = New_hololoid, aes(x= range_1.1(log10(holoploid)), y = exp(Pred)/ (1+exp(Pred))),size=1, color = "red") +
  geom_line(data = New_hololoid, aes(x= range_1.1(log10(holoploid)), y = exp(Pred.linear)/ (1+exp(Pred.linear))),size=1) +
  theme_classic() + 
  coord_cartesian(ylim = c(0,1))+ # does not remove points
  labs(y="Naturalization probablity",x="Holoploid genome size (pg)") +
  theme(legend.position = c(0.85, 0.6)) +
  theme(axis.text=element_text(size=12,face="bold"),
        axis.title=element_text(size=14,face="bold")) +
  theme(plot.title = element_text(size = 12, face = "bold"),
        legend.title=element_text(size=14, face = "bold"), 
        legend.text=element_text(size=12, face = "bold")) +
  ggtitle("(a)")

##### holoploid vs naturalization extent  ################

genome.holoploid.extent.data <- GenomeSize_df_holoploid_final %>% filter(reg_nat_total != "NA" ) %>% droplevels()
summary(genome.holoploid.extent.data)
str(genome.holoploid.extent.data) ##3631 species 2364 species for the new data
## prepare the extent phylogeny
#But get the tips to keep
keep.spp_holoploid_ext <-levels(genome.holoploid.extent.data$species)
remove_taxa_holoploid_ext = setdiff(genome.phy$tip.label, keep.spp_holoploid_ext)
genome.holoploid.extent.phy <- drop.tip(genome.phy,remove_taxa_holoploid_ext)

# Order data by tip order

genome.holoploid.extent.data = genome.holoploid.extent.data[(genome.holoploid.extent.data$species %in% genome.holoploid.extent.phy$tip.label), ]

row.names(genome.holoploid.extent.data) = genome.holoploid.extent.data$species
#sort_by_tip.label
genome.holoploid.extent.data <- genome.holoploid.extent.data[match(genome.holoploid.extent.phy$tip.label, genome.holoploid.extent.data$species),]
###check the names
library(geiger)  #for name.check
name.check(genome.holoploid.extent.phy, genome.holoploid.extent.data)

str(genome.holoploid.extent.data)

genome.null.extent.holo.model <- phylolm(range_1.1(log10(reg_nat_total)) ~ 1,
                                    phy = genome.holoploid.extent.phy, data = genome.holoploid.extent.data, 
                                    method = "lambda")
summary(genome.null.extent.holo.model)

genome.holoploid.extent.model <- phylolm(range_1.1(log10(reg_nat_total)) ~ range_1.1(log10(holoploid)),
                                         phy = genome.holoploid.extent.phy, data = genome.holoploid.extent.data,
                                         method = "lambda")

summary(genome.holoploid.extent.model)
plot(genome.holoploid.extent.model)
## model Residual Plot
qplot(range_1.1(log10(genome.holoploid.extent.data$reg_nat_total)), residuals(genome.holoploid.extent.model))
hist(residuals(genome.holoploid.extent.model))
plot(density(resid(genome.holoploid.extent.model))) #A density plot
# A quantile normal plot - good for checking normality
qqnorm(resid(genome.holoploid.extent.model)) 
qqline(resid(genome.holoploid.extent.model))


genome.holoploid.extent.q.model <- phylolm(range_1.1(log10(reg_nat_total)) ~ range_1.1(log10(holoploid)) +
                                             I(range_1.1(log10(holoploid))^2), 
                                           phy = genome.holoploid.extent.phy, data = genome.holoploid.extent.data,
                                           method = "lambda")

summary(genome.holoploid.extent.q.model)
plot(genome.holoploid.extent.q.model)
## model Residual Plot
qplot(range_1.1(log10(genome.holoploid.extent.data$reg_nat_total)), residuals(genome.holoploid.extent.q.model))
hist(residuals(genome.holoploid.extent.q.model))
plot(density(resid(genome.holoploid.extent.q.model))) #A density plot
# A quantile normal plot - good for checking normality
qqnorm(resid(genome.holoploid.extent.q.model)) 
qqline(resid(genome.holoploid.extent.q.model))

##plot

X1.holo.nat.ext <- model.matrix(~ range_1.1(log10(holoploid)), data = New_hololoid)
New_hololoid$Pred.ext <- X1.holo.nat.ext %*% coefficients(genome.holoploid.extent.model)#model coefficients

X2.holo.nat.ext <- model.matrix(~ range_1.1(log10(holoploid)) +
                     I(range_1.1(log10(holoploid))^2), data = New_hololoid)
New_hololoid$Pred.q.ext <- X2.holo.nat.ext %*% coefficients(genome.holoploid.extent.q.model)#model coefficients


New_hololoid$SE.ext <- sqrt(diag(X1.holo.nat.ext %*%vcov(genome.holoploid.extent.model) %*% t(X1.holo.nat.ext)))

p_naturalization_holo_extent <- ggplot(genome.holoploid.extent.data,aes(x= range_1.1(log10(holoploid)), y= range_1.1(log10(reg_nat_total)))) +
  geom_point(aes(color = PLOIDY, shape = taxno),position = position_jitter(w = 0.01, h = 0.09)) + 
 # scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"))+
  scale_color_viridis(discrete=TRUE, direction = 1) +
  scale_shape_manual(values=c(3, 16, 17,5))+
  #geom_smooth(method = "lm") +
  geom_line(data = New_hololoid, aes(x= range_1.1(log10(holoploid)), y = Pred.ext),size=1) +
  geom_line(data = New_hololoid, aes(x= range_1.1(log10(holoploid)), y = Pred.q.ext),size=1,  color = "red") +
  theme_classic() + 
  #coord_cartesian(ylim = c(-1,1))+ # does not remove points
  labs(y="Naturalization extent (Number of regions)",x="holoploid genome size (pg)") +
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=12,face="bold"),
        axis.title=element_text(size=14,face="bold")) +
  theme(plot.title = element_text(size = 12, face = "bold"),
        legend.title=element_text(size=14, face = "bold"), 
        legend.text=element_text(size=12, face = "bold")) +
  ggtitle("(b)")

######   holoploid vs. invasion range   ########


genome.global.holo.extent.invasive.data <- GenomeSize_df_holoploid_final %>%  filter (reg_no_invasive >0 ) %>% droplevels()
summary(genome.global.holo.extent.invasive.data)
str(genome.global.holo.extent.invasive.data) ##1060 species 680 species in the new data
## prepare the extent phylogeny
#But get the tips to keep
keep.spp_holoploid_invext <-levels(genome.global.holo.extent.invasive.data$species)
remove_taxa_holoploid_invext = setdiff(genome.phy$tip.label, keep.spp_holoploid_invext)
genome.extent.invasive.holoploid.phy <- drop.tip(genome.phy,remove_taxa_holoploid_invext)

# Order data by tip order

genome.global.holo.extent.invasive.data = genome.global.holo.extent.invasive.data[(genome.global.holo.extent.invasive.data$species %in% genome.extent.invasive.holoploid.phy$tip.label), ]

row.names(genome.global.holo.extent.invasive.data) = genome.global.holo.extent.invasive.data$species
#sort_by_tip.label
genome.global.holo.extent.invasive.data <- genome.global.holo.extent.invasive.data[match(genome.extent.invasive.holoploid.phy$tip.label, genome.global.holo.extent.invasive.data$species),]
##check the names
library(geiger)  #for name.check
name.check(genome.extent.invasive.holoploid.phy, genome.global.holo.extent.invasive.data)

str(genome.global.holo.extent.invasive.data)
# 
## monoploid * taxomony
genome.null.holo.invasive.extent.model <- phylolm(range_1.1(log10(reg_no_invasive)) ~ 1,
                                             phy = genome.extent.invasive.holoploid.phy, data = genome.global.holo.extent.invasive.data, 
                                             method = "lambda")
summary(genome.null.holo.invasive.extent.model)

genome.holoploid.invasive.extent.model <- phylolm(range_1.1(log10(reg_no_invasive)) ~ range_1.1(log10(holoploid)),
                                                  phy = genome.extent.invasive.holoploid.phy, data = genome.global.holo.extent.invasive.data, 
                                                  method = "lambda")

summary(genome.holoploid.invasive.extent.model)

plot(genome.holoploid.invasive.extent.model)
## model Residual Plot
qplot(range_1.1(log10(genome.global.holo.extent.invasive.data$reg_no_invasive)), residuals(genome.holoploid.invasive.extent.model))
hist(residuals(genome.holoploid.invasive.extent.model))
plot(density(resid(genome.holoploid.invasive.extent.model))) #A density plot
# A quantile normal plot - good for checking normality
qqnorm(resid(genome.holoploid.invasive.extent.model)) 
qqline(resid(genome.holoploid.invasive.extent.model))


###plot the estimated slope
plot(range_1.1(log10(genome.global.holo.extent.invasive.data$holoploid)),
     jitter(range_1.1(log10(genome.global.holo.extent.invasive.data$reg_no_invasive)),factor=0,amount=0.02),
     xlab="holoploid",ylab="no.of invasive regions",xlim=c(-1,1))
cc.genome_holo_inva_ext <- coef(genome.holoploid.invasive.extent.model)
curve(plogis(cc.genome_holo_inva_ext[1]+ cc.genome_holo_inva_ext[2]*x),col="red",
      lwd=2,lty=1,add=TRUE)
cc.genome_holo_inv_ext.q <- coef(genome.holoploid.extent.invasive.q.model)
curve(plogis(cc.genome_holo_inv_ext.q[1]+ cc.genome_holo_inv_ext.q[2]*x + cc.genome_holo_inv_ext.q[3]*x^2 ),col="red",
      lwd=2,lty=1,add=TRUE)





genome.holoploid.extent.invasive.q.model <- phylolm(range_1.1(log10(reg_no_invasive)) ~ range_1.1(log10(holoploid)) +
                                                      I(range_1.1(log10(holoploid))^2), 
                                                    phy = genome.extent.invasive.holoploid.phy, data = genome.global.holo.extent.invasive.data, 
                                                    method = "lambda")

summary(genome.holoploid.extent.invasive.q.model)
plot(genome.holoploid.extent.invasive.q.model)
## model Residual Plot
qplot(range_1.1(log10(genome.global.holo.extent.invasive.data$reg_no_invasive)), residuals(genome.holoploid.extent.invasive.q.model))
hist(residuals(genome.holoploid.extent.invasive.q.model))
plot(density(resid(genome.holoploid.extent.invasive.q.model))) #A density plot
# A quantile normal plot - good for checking normality
qqnorm(resid(genome.holoploid.extent.invasive.q.model)) 
qqline(resid(genome.holoploid.extent.invasive.q.model))



##plot
X1.holo <- model.matrix(~ range_1.1(log10(holoploid)), data = New_hololoid)
New_hololoid$Pred.ext.inv <- X1.holo %*% coefficients(genome.holoploid.invasive.extent.model)#model coefficients
X.holo <- model.matrix(~ range_1.1(log10(holoploid)) +
                     I(range_1.1(log10(holoploid))^2), data = New_hololoid)
New_hololoid$Pred.q.ext.inv <- X.holo  %*% coefficients(genome.holoploid.extent.invasive.q.model)#model coefficients

##plot
p_invasive_holo_extent <- ggplot(genome.global.holo.extent.invasive.data,
                            aes(x= range_1.1(log10(holoploid)), y= range_1.1(log10(reg_no_invasive)))) +
  geom_point(aes(color = PLOIDY, shape = taxno),position = position_jitter(w = 0.01, h = 0.09)) + 
  scale_color_viridis(discrete=TRUE, direction = 1) +
  scale_shape_manual(values=c(3, 16, 17,5))+
  geom_line(data = New_hololoid, aes(x= range_1.1(log10(holoploid)), y = Pred.ext.inv),size=1) +
  geom_line(data = New_hololoid, aes(x= range_1.1(log10(holoploid)), y = Pred.q.ext.inv),size=1, color = "red") +
  theme_classic() + 
  coord_cartesian(ylim = c(-1.4,1))+ # does not remove points
  labs(y="Invasion extent (Number of regions)",x="holoploid genome size (pg)") +
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=12,face="bold"),
        axis.title=element_text(size=14,face="bold")) +
  theme(plot.title = element_text(size = 12, face = "bold"),
        legend.title=element_text(size=14, face = "bold"), 
        legend.text=element_text(size=12, face = "bold")) +
  ggtitle("(c)")



library(ggpubr)

figure2_holo <-ggarrange(p_nat_holo_incidence,                                                 # First row with incidence plot
                    ggarrange(p_naturalization_holo_extent, p_invasive_holo_extent, ncol = 2, labels = c("B", "C")), # Second row with two extent plot
                    nrow = 2, heights = c(0.8, 1),  ## adjust height
                    labels = "A"  )                                      # Labels of the incidence plot


pdf("Fig 2.holoploid GS incidence & extent_noSE_linear_qua_22Aug2021_color.pdf", useDingbats=FALSE, width=12, height=9)
figure2_holo
dev.off()

### plot the results holoploid with ploidy level shown  not used any more #######

p_nat_holo_incidence_ploidy <- ggplot(GenomeSize_df_holoploid_final,aes(x= range_1.1(log10(holoploid)), y=GloNAF_incidence)) +
  geom_point(aes(color = PLOIDY, shape = taxno), position = position_jitter(w = 0.01, h = 0.09)) + 
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"))+
  scale_shape_manual(values=c(3, 16, 17))+
  geom_line(data = New_hololoid, aes(x= range_1.1(log10(holoploid)), y = exp(Pred)/ (1+exp(Pred))),size=1, color = "red") +
  geom_line(data = New_hololoid, aes(x= range_1.1(log10(holoploid)), y = exp(Pred.linear)/ (1+exp(Pred.linear))),size=1) +
  theme_classic() + 
  coord_cartesian(ylim = c(0,1))+ # does not remove points
  labs(y="Naturalization probablity",x="Holoploid genome size (pg)") +
  theme(legend.position = c(0.95, 0.6)) +
  theme(axis.text=element_text(size=12,face="bold"),
        axis.title=element_text(size=14,face="bold")) +
  theme(plot.title = element_text(size = 12, face = "bold"),
        legend.title=element_text(size=14, face = "bold"), 
        legend.text=element_text(size=12, face = "bold")) +
  ggtitle("(a)")


p_naturalization_holo_extent_ploidy <- ggplot(genome.holoploid.extent.data,aes(x= range_1.1(log10(holoploid)), y= range_1.1(log10(reg_nat_total)))) +
  geom_point(aes(color = PLOIDY, shape = taxno),position = position_jitter(w = 0.01, h = 0.09)) + 
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"))+
  scale_shape_manual(values=c(3, 16, 17))+
  geom_line(data = New_hololoid, aes(x= range_1.1(log10(holoploid)), y = Pred.ext),size=1) +
  geom_line(data = New_hololoid, aes(x= range_1.1(log10(holoploid)), y = Pred.q.ext),size=1,  color = "red") +
  theme_classic() + 
  #coord_cartesian(ylim = c(-1,1))+ # does not remove points
  labs(y="Naturalization extent (Number of regions)",x="holoploid genome size (pg)") +
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=12,face="bold"),
        axis.title=element_text(size=14,face="bold")) +
  theme(plot.title = element_text(size = 12, face = "bold"),
        legend.title=element_text(size=14, face = "bold"), 
        legend.text=element_text(size=12, face = "bold")) +
  ggtitle("(b)")


p_invasive_holo_extent_ploidy <- ggplot(genome.global.holo.extent.invasive.data,
                                 aes(x= range_1.1(log10(holoploid)), y= range_1.1(log10(reg_no_invasive)))) +
  geom_point(aes(color = PLOIDY, shape = taxno),position = position_jitter(w = 0.01, h = 0.09)) + 
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"))+
  scale_shape_manual(values=c(3, 16, 17))+
  geom_line(data = New_hololoid, aes(x= range_1.1(log10(holoploid)), y = Pred.ext.inv),size=1) +
  geom_line(data = New_hololoid, aes(x= range_1.1(log10(holoploid)), y = Pred.q.ext.inv),size=1, color = "red") +
  theme_classic() + 
  coord_cartesian(ylim = c(-1.4,1))+ # does not remove points
  labs(y="Invasion extent (Number of regions)",x="holoploid genome size (pg)") +
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=12,face="bold"),
        axis.title=element_text(size=14,face="bold")) +
  theme(plot.title = element_text(size = 12, face = "bold"),
        legend.title=element_text(size=14, face = "bold"), 
        legend.text=element_text(size=12, face = "bold")) +
  ggtitle("(c)")




figure2_holo_ploidy <-ggarrange(p_nat_holo_incidence_ploidy,                                                 # First row with incidence plot
                         ggarrange(p_naturalization_holo_extent_ploidy, p_invasive_holo_extent_ploidy, ncol = 2, labels = c("B", "C")), # Second row with two extent plot
                         nrow = 2, heights = c(0.8, 1),  ## adjust height
                         labels = "A"  )                                      # Labels of the incidence plot


pdf("Fig 2.holoploid GS incidence & extent_noSE_ploidy_linear_qua.pdf", useDingbats=FALSE, width=12, height=9)
figure2_holo_ploidy
dev.off()

######  II NO.of ploidy level   #######

summary(GenomeSize_df_final)

hist(GenomeSize_df_final$no_ploidy_lev)

genome_ploidy_data <-   GenomeSize_df_final %>%  filter (no_ploidy_lev >0) %>% droplevels()

summary(genome_ploidy_data)

##prepare the phy
keep.ploidspp <-levels(genome_ploidy_data$species)
remove_ploidytaxa = setdiff(genome.phy$tip.label, keep.ploidspp)
genome.ploidy.phy <- drop.tip(genome.phy,remove_ploidytaxa)

# Order data by tip order

genome_ploidy_data = genome_ploidy_data[(genome_ploidy_data$species %in% genome.ploidy.phy$tip.label), ]

row.names(genome_ploidy_data) = genome_ploidy_data$species
###check the names
library(geiger)  #for name.check
name.check(genome.ploidy.phy, genome_ploidy_data)

str(genome_ploidy_data)
summary(genome_ploidy_data)

genome.ploidy.null.model <- phyloglm(GloNAF_incidence ~ 1, 
                                     phy = genome.ploidy.phy, data = genome_ploidy_data, 
                                     method = c("logistic_MPLE","logistic_IG10","poisson_GEE"),
                                     btol = 10, log.alpha.bound = 4,
                                     start.beta=NULL, start.alpha=NULL,
                                     boot = 0, full.matrix = TRUE)

genome.ploidy.model <- phyloglm(GloNAF_incidence ~ no_ploidy_lev, 
                                phy = genome.ploidy.phy, data = genome_ploidy_data, 
                                method = c("logistic_MPLE","logistic_IG10","poisson_GEE"),
                                btol = 10, log.alpha.bound = 4,
                                start.beta=NULL, start.alpha=NULL,
                                boot = 0, full.matrix = TRUE)

genome.ploidy.null.model1 <- glm(GloNAF_incidence ~ no_ploidy_lev, genome_ploidy_data,
                                 family = "binomial")
library(rr2)
summary(genome.ploidy.model)
summary(genome.ploidy.null.model1)
R2(genome.ploidy.model ,genome.ploidy.null.model1)
###plot the estimated slope
plot(genome_ploidy_data$no_ploidy_lev,
     jitter(genome_ploidy_data$GloNAF_incidence,factor=0,amount=0.02),
     xlab="no ploidy level",ylab="Naturalization incidence",xlim=c(0,6))
cc.genome_ploidy <- coef(genome.ploidy.model)
curve(plogis(cc.genome_ploidy[1]+ cc.genome_ploidy[2]*x),col="red",
      lwd=2,lty=1,add=TRUE)
pdf("Fig.no_ploidy.level & naturalization incidence.pdf", useDingbats=FALSE, width=6, height=6)


dev.off()


new_ploidy_data <-  data.frame( no_ploidy_lev = rep(seq(from = 1,to = 6, length = 5000),6))


X_ploidy <- model.matrix(~ no_ploidy_lev, data = new_ploidy_data)


new_ploidy_data$Pred <- X_ploidy %*% coefficients(genome.ploidy.model)#model coefficients
new_ploidy_data$SE <- sqrt(  diag(X_ploidy %*%vcov(genome.ploidy.model) %*% t(X_ploidy))  )


p_ploidy_nat_incidence <-ggplot(genome_ploidy_data,aes(x= no_ploidy_lev, y= GloNAF_incidence )) +
  geom_point(aes(color = taxno, shape = taxno),position = position_jitter(w = 0.4, h = 0.2)) + 
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"))+
  scale_shape_manual(values=c(3, 16, 17))+
  geom_line(data = new_ploidy_data, aes(x= no_ploidy_lev, y = exp(Pred)/ (1+exp(Pred))),size=1) +
  geom_line(data = new_ploidy_data, aes(x= no_ploidy_lev, y = exp(Pred + 1.96*SE)/ (1+exp(Pred + 1.96*SE))),size=1,linetype="dashed") +
  geom_line(data = new_ploidy_data, aes(x= no_ploidy_lev, y = exp(Pred - 1.96*SE)/ (1+exp(Pred - 1.96*SE))),size=1,linetype="dashed")+
  theme_classic() + 
  coord_cartesian(ylim = c(0,1))+ # does not remove points
  scale_x_continuous(breaks=seq(1,6,1)) +
  labs(y="Naturalization probablity",x="Number of ploidy level") +
  theme(legend.position = c(0.85,0.6)) +
  theme(axis.text=element_text(size=12,face="bold"),
        axis.title=element_text(size=14,face="bold")) +
  theme(plot.title = element_text(size = 12, face = "bold"),
        legend.title=element_text(size=14, face = "bold"), 
        legend.text=element_text(size=12, face = "bold")) +
  ggtitle("(a)")

### naturalized extent

genome_ploidy_ext_data <-   genome_ploidy_data %>%  filter (reg_nat_total >0) %>% droplevels()
summary(genome_ploidy_ext_data)
##prepare the phy
keep.ploid.extspp <-levels(genome_ploidy_ext_data$species)
remove_ploidy.exttaxa = setdiff(genome.phy$tip.label, keep.ploid.extspp)
genome.ploidy.ext.phy <- drop.tip(genome.phy,remove_ploidy.exttaxa)

# Order data by tip order

genome_ploidy_ext_data = genome_ploidy_ext_data[(genome_ploidy_ext_data$species %in% genome.ploidy.ext.phy$tip.label), ]

row.names(genome_ploidy_ext_data) = genome_ploidy_ext_data$species
###check the names
library(geiger)  #for name.check
name.check(genome.ploidy.ext.phy, genome_ploidy_ext_data)

str(genome_ploidy_ext_data)
summary(genome_ploidy_ext_data)


genome.ploidy.extent.model <- phylolm(log10(reg_nat_total) ~ no_ploidy_lev,
                                      phy = genome.ploidy.ext.phy, data = genome_ploidy_ext_data, 
                                      method = "lambda")

summary(genome.ploidy.extent.model)

genome.ploidy.glm.extent.model <- glm(log10(reg_nat_total) ~ no_ploidy_lev,
                                      data = genome_ploidy_ext_data)

summary(genome.ploidy.glm.extent.model)
plot(genome.ploidy.extent.model)
## model Residual Plot
qplot(range_1.1(log10(genome_ploidy_ext_data$reg_nat_total)), residuals(genome.ploidy.extent.model))
hist(residuals(genome.ploidy.extent.model))
plot(density(resid(genome.ploidy.extent.model))) #A density plot
# A quantile normal plot - good for checking normality
qqnorm(resid(genome.ploidy.extent.model)) 
qqline(resid(genome.ploidy.extent.model))

plot(genome_ploidy_ext_data$no_ploidy_lev,
     jitter(log10(genome_ploidy_ext_data$reg_nat_total),factor=0,amount=0.02),
     xlab="monoploid",ylab="no.of extent regions",xlim=c(1,6))
cc.genome_ploidy.nat_ext <- coef(genome.ploidy.extent.model)
curve(plogis(cc.genome_ploidy.nat_ext[1]+ cc.genome_ploidy.nat_ext[2]*x),col="red",
      lwd=2,lty=1,add=TRUE)


new_ploidy_data <-  data.frame( no_ploidy_lev = rep(seq(from = 1,to = 6, length = 5000),6))
X_ploidy <- model.matrix(~ no_ploidy_lev, data = new_ploidy_data)
new_ploidy_data$Pred.ext <- X_ploidy %*% coefficients(genome.ploidy.extent.model)#model coefficients


p_ploidy_nat_extent <- 
ggplot(genome_ploidy_ext_data,aes(x= no_ploidy_lev, y= log10(reg_nat_total))) +
 
  geom_point(aes(color = taxno, shape = taxno),position = position_jitter(w = 0.4, h = 0.2)) + 

  scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"))+
  scale_shape_manual(values=c(3, 16, 17))+

  geom_line(data = new_ploidy_data, aes(x= no_ploidy_lev, y = Pred.ext),size=1) +

  theme_classic() + 
  coord_cartesian(ylim = c(0,2.8))+ # does not remove points
  scale_x_continuous(breaks=seq(1,6,1)) +
  labs(y="Naturalization extent (Number of regions)",x="Number of ploidy level") +
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=12,face="bold"),
        axis.title=element_text(size=14,face="bold")) +
  theme(plot.title = element_text(size = 12, face = "bold"),
        legend.title=element_text(size=14, face = "bold"), 
        legend.text=element_text(size=12, face = "bold")) +
  ggtitle("(b)")


#######################################################################
### naturalization extent mainland

genome_ploidy_mainext_data <-   genome_ploidy_data %>%  filter (reg_nat_main >0) %>% droplevels()
summary(genome_ploidy_mainext_data)
##prepare the phy
keep.ploid.mainextspp <-levels(genome_ploidy_mainext_data$species)

# keepb489tips <- genome.phy$tip.label[match(keep.b489spp, genome.phy$tip.label)]

remove_ploidy.mainexttaxa = setdiff(genome.phy$tip.label, keep.ploid.mainextspp)
genome.ploidy.mainext.phy <- drop.tip(genome.phy,remove_ploidy.mainexttaxa)

# Order data by tip order

genome_ploidy_mainext_data = genome_ploidy_mainext_data[(genome_ploidy_mainext_data$species %in% genome.ploidy.mainext.phy$tip.label), ]

row.names(genome_ploidy_mainext_data) = genome_ploidy_mainext_data$species
###check the names
library(geiger)  #for name.check
name.check(genome.ploidy.mainext.phy, genome_ploidy_mainext_data)

str(genome_ploidy_mainext_data)
summary(genome_ploidy_mainext_data)
#### mainland holoploid

genome.ploidy.mainextent.model <- phylolm(range_1.1(log10(reg_nat_main)) ~ range_1.1(log10(no_ploidy_lev)),
                                          phy = genome.ploidy.mainext.phy, data = genome_ploidy_mainext_data, 
                                          method = "lambda")

summary(genome.ploidy.mainextent.model)
#R2(genome.monoploid.main.mainextent.model, phy = genome.mainextent.main.phy)
plot(genome.ploidy.mainextent.model)
## model Residual Plot
qplot(range_1.1(log10(genome_ploidy_mainext_data$reg_nat_main)), residuals(genome.ploidy.mainextent.model))
hist(residuals(genome.ploidy.mainextent.model))
plot(density(resid(genome.ploidy.mainextent.model))) #A density plot
# A quantile normal plot - good for checking normality
qqnorm(resid(genome.ploidy.mainextent.model)) 
qqline(resid(genome.ploidy.mainextent.model))

### naturalization extent island

genome_ploidy_islext_data <-   genome_ploidy_data %>%  filter (reg_nat_island >0) %>% droplevels()
summary(genome_ploidy_islext_data)
##prepare the phy
keep.ploid.islextspp <-levels(genome_ploidy_islext_data$species)
remove_ploidy.islexttaxa = setdiff(genome.phy$tip.label, keep.ploid.islextspp)
genome.ploidy.islext.phy <- drop.tip(genome.phy,remove_ploidy.islexttaxa)

# Order data by tip order

genome_ploidy_islext_data = genome_ploidy_islext_data[(genome_ploidy_islext_data$species %in% genome.ploidy.islext.phy$tip.label), ]

row.names(genome_ploidy_islext_data) = genome_ploidy_islext_data$species
###check the names
library(geiger)  #for name.check
name.check(genome.ploidy.islext.phy, genome_ploidy_islext_data)

str(genome_ploidy_islext_data)
summary(genome_ploidy_islext_data)
#### mainland holoploid

genome.ploidy.islextent.model <- phylolm(range_1.1(log10(reg_nat_island)) ~ range_1.1(log10(no_ploidy_lev)),
                                         phy = genome.ploidy.islext.phy, data = genome_ploidy_islext_data, 
                                         method = "lambda")

summary(genome.ploidy.islextent.model)
plot(genome.ploidy.islextent.model)
## model Residual Plot
qplot(range_1.1(log10(genome_ploidy_islext_data$reg_nat_island)), residuals(genome.ploidy.islextent.model))
hist(residuals(genome.ploidy.islextent.model))
plot(density(resid(genome.ploidy.islextent.model))) #A density plot
# A quantile normal plot - good for checking normality
qqnorm(resid(genome.ploidy.islextent.model)) 
qqline(resid(genome.ploidy.islextent.model))


#################### plot the results
summary.naturalized.extent.ploidyNO <- read.csv("extent.results.ploidyNO.summary.csv", header = T, sep = ';')
str(summary.naturalized.extent.ploidyNO)


# Error bars represent standard error of the mean
library(ggplot2)

# The errorbars overlapped, so use position_dodge to move them horizontally
pd <- position_dodge(0.5) # move them .05 to the left and right
ggplot(summary.naturalized.extent.ploidyNO, aes(x= group, y= estimate, shape = type)) + 
  geom_point(position= pd, size=3.5) +
  geom_errorbar(aes(ymin= estimate - SE, ymax= estimate + SE),
                width=.1, position=pd) +
  ggtitle("Naturalization extent") +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     strip.background = element_rect(colour="white", fill="white"),
                     text = element_text(size = 16),
                     legend.position = c(0.9,0.8))

pdf("Fig.4 ploidy NO naturalization invasive extent.pdf", useDingbats=FALSE, width=8, height=6)


dev.off()

###### invasive extent

genome_ploidy_invasive_data <-   genome_ploidy_data %>%  filter (reg_no_invasive >0) %>% droplevels()
summary(genome_ploidy_invasive_data)
##prepare the phy
keep.ploid.invasivespp <-levels(genome_ploidy_invasive_data$species)
remove_ploidy.invasivetaxa = setdiff(genome.phy$tip.label, keep.ploid.invasivespp)
genome.ploidy.invasive.phy <- drop.tip(genome.phy,remove_ploidy.invasivetaxa)

# Order data by tip order

genome_ploidy_invasive_data = genome_ploidy_invasive_data[(genome_ploidy_invasive_data$species %in% genome.ploidy.invasive.phy$tip.label), ]

row.names(genome_ploidy_invasive_data) = genome_ploidy_invasive_data$species
###check the names
library(geiger)  #for name.check
name.check(genome.ploidy.invasive.phy, genome_ploidy_invasive_data)

str(genome_ploidy_invasive_data)
summary(genome_ploidy_invasive_data)
#### mainland holoploid

genome.ploidy.invasiveent.model <- phylolm(log10(reg_no_invasive) ~ no_ploidy_lev,
                                           phy = genome.ploidy.invasive.phy, data = genome_ploidy_invasive_data, 
                                           method = "lambda")

summary(genome.ploidy.invasiveent.model)
plot(genome.ploidy.invasiveent.model)
## model Residual Plot
qplot(range_1.1(log10(genome_ploidy_invasive_data$reg_no_invasive)), residuals(genome.ploidy.invasiveent.model))
hist(residuals(genome.ploidy.invasiveent.model))
plot(density(resid(genome.ploidy.invasiveent.model))) #A density plot
# A quantile normal plot - good for checking normality
qqnorm(resid(genome.ploidy.invasiveent.model)) 
qqline(resid(genome.ploidy.invasiveent.model))

plot(genome_ploidy_invasive_data$no_ploidy_lev,
     jitter(log10(genome_ploidy_invasive_data$reg_no_invasive),factor=0,amount=0.02),
     xlab="monoploid",ylab="no.of extent regions",xlim=c(1,6))
cc.genome_ploidy.inv_ext <- coef(genome.ploidy.invasiveent.model)
curve(plogis(cc.genome_ploidy.inv_ext[1]+ cc.genome_ploidy.inv_ext[2]*x),col="red",
      lwd=2,lty=1,add=TRUE)


new_ploidy_data <-  data.frame( no_ploidy_lev = rep(seq(from = 1,to = 6, length = 5000),6))
X_ploidy <- model.matrix(~ no_ploidy_lev, data = new_ploidy_data)
new_ploidy_data$Pred.ext <- X_ploidy %*% coefficients(genome.ploidy.extent.model)#model coefficients
new_ploidy_data$SE.ext <- sqrt(  diag(X_ploidy %*%vcov(genome.ploidy.extent.model) %*% t(X_ploidy))  )



ggplot(genome_ploidy_invasive_data,aes(x= no_ploidy_lev, y= log10(reg_no_invasive))) +
  geom_point(aes(color = taxno, shape = taxno),position = position_jitter(w = 0.4, h = 0.2)) + 
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"))+
  scale_shape_manual(values=c(3, 16, 17))+
  geom_smooth(method = "lm") +
  geom_abline(aes(intercept=cc.genome_ploidy.inv_ext[1],slope= cc.genome_ploidy.inv_ext[2])) +
  theme_classic() + 
  # coord_cartesian(ylim = c(0,2))+ # does not remove points
  scale_x_continuous(breaks=seq(1,6,1)) +
  labs(y="Invasion extent",x="Number of ploidy level") +
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=12,face="bold"),
        axis.title=element_text(size=14,face="bold")) +
  theme(plot.title = element_text(size = 12, face = "bold"),
        legend.title=element_text(size=14, face = "bold"), 
        legend.text=element_text(size=12, face = "bold")) +
  ggtitle("(c)")



figure3_0 <-ggarrange(p_ploidy_nat_incidence,  p_ploidy_nat_extent,                                          
                      #ggarrange(p_ploidy_nat_extent, p_invasive_extent, ncol = 2, labels = c("B", "C")), # Second row with two extent plot
                      nrow = 1, widths = c(1, 0.8),  ## adjust height
                      labels = c("A", "B"))                                      # Labels of the incidence plot


pdf("Fig 3.ploidy_ incidence & extent.pdf", useDingbats=FALSE, width=12, height=9)
figure3_0
dev.off()



### II 6  diploid and ployploids

II6 <- c("diploid", "polyploid")
genome_ploidy_II6_data <-  genome_ploidy_data  %>%  filter (PLOIDY   %in% II6) %>% droplevels()
summary(genome_ploidy_II6_data)

##prepare the phy
keep.ploidII6spp <-levels(genome_ploidy_II6_data$species)
remove_ploidyII6taxa = setdiff(genome.phy$tip.label, keep.ploidII6spp)
genome.ploidy.II6.phy <- drop.tip(genome.phy,remove_ploidyII6taxa)

# Order data by tip order

genome_ploidy_II6_data = genome_ploidy_II6_data[(genome_ploidy_II6_data$species %in% genome.ploidy.II6.phy$tip.label), ]

row.names(genome_ploidy_II6_data) = genome_ploidy_II6_data$species
###check the names
library(geiger)  #for name.check
name.check(genome.ploidy.II6.phy, genome_ploidy_II6_data)

str(genome_ploidy_II6_data)
summary(genome_ploidy_II6_data)

genome.ploidy.II6.null.model <- phyloglm(GloNAF_incidence ~ 1, 
                                         phy = genome.ploidy.II6.phy, data = genome_ploidy_II6_data, 
                                         method = c("logistic_MPLE","logistic_IG10","poisson_GEE"),
                                         btol = 10, log.alpha.bound = 4,
                                         start.beta=NULL, start.alpha=NULL,
                                         boot = 0, full.matrix = TRUE)

genome.ploidy.II6.model <- phyloglm(GloNAF_incidence ~  PLOIDY-1, 
                                    phy = genome.ploidy.II6.phy, data = genome_ploidy_II6_data, 
                                    method = c("logistic_MPLE","logistic_IG10","poisson_GEE"),
                                    btol = 10, log.alpha.bound = 4,
                                    start.beta=NULL, start.alpha=NULL,
                                    boot = 0, full.matrix = TRUE)

summary(genome.ploidy.II6.model)
R2(genome.ploidy.II6.model ,genome.ploidy.II6.null.model)


pdf("Fig.no_ploidy.level & naturalization incidence.pdf", useDingbats=FALSE, width=6, height=6)


dev.off()

## II6 invasive extent

plot(genome_ploidy_II6_data$PLOIDY, genome_ploidy_II6_data$reg_nat_total)
plot(genome_ploidy_II6_data$PLOIDY, genome_ploidy_II6_data$reg_nat_main)
plot(genome_ploidy_II6_data$PLOIDY, genome_ploidy_II6_data$reg_nat_island)
plot(genome_ploidy_II6_data$PLOIDY, genome_ploidy_II6_data$reg_no_invasive)

### global
genome_ploidy_II6_ext_data <-  genome_ploidy_II6_data %>%  filter (reg_nat_total  > 0) %>% droplevels()
summary(genome_ploidy_II6_ext_data)
##prepare the phy
keep.ploid.II6extspp <-levels(genome_ploidy_II6_ext_data$species)
remove_ploidy.II6exttaxa = setdiff(genome.phy$tip.label, keep.ploid.II6extspp)
genome.ploidy.II6ext.phy <- drop.tip(genome.phy,remove_ploidy.II6exttaxa)

# Order data by tip order

genome_ploidy_II6_ext_data = genome_ploidy_II6_ext_data[(genome_ploidy_II6_ext_data$species %in% genome.ploidy.II6ext.phy$tip.label), ]

row.names(genome_ploidy_II6_ext_data) = genome_ploidy_II6_ext_data$species
###check the names
library(geiger)  #for name.check
name.check(genome.ploidy.II6ext.phy, genome_ploidy_II6_ext_data)

str(genome_ploidy_II6_ext_data)
summary(genome_ploidy_II6_ext_data)
#### mainland holoploid

genome.ploidy.II6.extent.model <- phylolm(range_1.1(log10(reg_nat_total)) ~ PLOIDY-1,
                                          phy = genome.ploidy.II6ext.phy, data = genome_ploidy_II6_ext_data, 
                                          method = "lambda")

summary(genome.ploidy.II6.extent.model)
plot(genome.ploidy.II6.extent.model)
## model Residual Plot
qplot(range_1.1(log10(genome_ploidy_II6_ext_data$reg_nat_total)), residuals(genome.ploidy.II6.extent.model))
hist(residuals(genome.ploidy.II6.extent.model))
plot(density(resid(genome.ploidy.II6.extent.model))) #A density plot
# A quantile normal plot - good for checking normality
qqnorm(resid(genome.ploidy.II6.extent.model)) 
qqline(resid(genome.ploidy.II6.extent.model))


## mainland 
genome_ploidy_II6_mainext_data <-  genome_ploidy_II6_data %>%  filter (reg_nat_main  > 0) %>% droplevels()
summary(genome_ploidy_II6_mainext_data)
##prepare the phy
keep.ploid.II6mainextspp <-levels(genome_ploidy_II6_mainext_data$species)
remove_ploidy.II6mainexttaxa = setdiff(genome.phy$tip.label, keep.ploid.II6mainextspp)
genome.ploidy.II6mainext.phy <- drop.tip(genome.phy,remove_ploidy.II6mainexttaxa)

# Order data by tip order

genome_ploidy_II6_mainext_data = genome_ploidy_II6_mainext_data[(genome_ploidy_II6_mainext_data$species %in% genome.ploidy.II6mainext.phy$tip.label), ]

row.names(genome_ploidy_II6_mainext_data) = genome_ploidy_II6_mainext_data$species
###check the names
library(geiger)  #for name.check
name.check(genome.ploidy.II6mainext.phy, genome_ploidy_II6_mainext_data)

str(genome_ploidy_II6_mainext_data)
summary(genome_ploidy_II6_mainext_data)
#### mainland 

genome.ploidy.II6.mainextent.model <- phylolm(range_1.1(log10(reg_nat_main)) ~ PLOIDY-1,
                                              phy = genome.ploidy.II6mainext.phy, data = genome_ploidy_II6_mainext_data, 
                                              method = "lambda")

summary(genome.ploidy.II6.mainextent.model)
plot(genome.ploidy.II6.mainextent.model)
## model Residual Plot
qplot(range_1.1(log10(genome_ploidy_II6_mainext_data$reg_nat_total)), residuals(genome.ploidy.II6.mainextent.model))
hist(residuals(genome.ploidy.II6.mainextent.model))
plot(density(resid(genome.ploidy.II6.mainextent.model))) #A density plot
# A quantile normal plot - good for checking normality
qqnorm(resid(genome.ploidy.II6.mainextent.model)) 
qqline(resid(genome.ploidy.II6.mainextent.model))


## Island II6
genome_ploidy_II6_islext_data <-  genome_ploidy_II6_data %>%  filter (reg_nat_island  > 0) %>% droplevels()
summary(genome_ploidy_II6_islext_data)
##prepare the phy
keep.ploid.II6islextspp <-levels(genome_ploidy_II6_islext_data$species)
remove_ploidy.II6islexttaxa = setdiff(genome.phy$tip.label, keep.ploid.II6islextspp)
genome.ploidy.II6islext.phy <- drop.tip(genome.phy,remove_ploidy.II6islexttaxa)

# Order data by tip order

genome_ploidy_II6_islext_data = genome_ploidy_II6_islext_data[(genome_ploidy_II6_islext_data$species %in% genome.ploidy.II6islext.phy$tip.label), ]

row.names(genome_ploidy_II6_islext_data) = genome_ploidy_II6_islext_data$species
###check the names
library(geiger)  #for name.check
name.check(genome.ploidy.II6islext.phy, genome_ploidy_II6_islext_data)

str(genome_ploidy_II6_islext_data)
summary(genome_ploidy_II6_islext_data)
#### islland holoploid

genome.ploidy.II6.islextent.model <- phylolm(range_1.1(log10(reg_nat_island)) ~ PLOIDY-1,
                                             phy = genome.ploidy.II6islext.phy, data = genome_ploidy_II6_islext_data, 
                                             method = "lambda")

summary(genome.ploidy.II6.islextent.model)
plot(genome.ploidy.II6.islextent.model)
## model Residual Plot
qplot(range_1.1(log10(genome_ploidy_II6_islext_data$reg_nat_island)), residuals(genome.ploidy.II6.islextent.model))
hist(residuals(genome.ploidy.II6.islextent.model))
plot(density(resid(genome.ploidy.II6.islextent.model))) #A density plot
# A quantile normal plot - good for checking normality
qqnorm(resid(genome.ploidy.II6.islextent.model)) 
qqline(resid(genome.ploidy.II6.islextent.model))


#### invasive ranges

genome_ploidy_II6_invasive_data <-  genome_ploidy_II6_data %>%  filter (reg_no_invasive  > 0) %>% droplevels()
summary(genome_ploidy_II6_invasive_data)
##prepare the phy
keep.ploid.II6invasivespp <-levels(genome_ploidy_II6_invasive_data$species)
remove_ploidy.II6invasivetaxa = setdiff(genome.phy$tip.label, keep.ploid.II6invasivespp)
genome.ploidy.II6invasive.phy <- drop.tip(genome.phy,remove_ploidy.II6invasivetaxa)

# Order data by tip order

genome_ploidy_II6_invasive_data = genome_ploidy_II6_invasive_data[(genome_ploidy_II6_invasive_data$species %in% genome.ploidy.II6invasive.phy$tip.label), ]

row.names(genome_ploidy_II6_invasive_data) = genome_ploidy_II6_invasive_data$species
###check the names
library(geiger)  #for name.check
name.check(genome.ploidy.II6invasive.phy, genome_ploidy_II6_invasive_data)

str(genome_ploidy_II6_invasive_data)
summary(genome_ploidy_II6_invasive_data)


genome.ploidy.II6.invasiveent.model <- phylolm(range_1.1(log10(reg_no_invasive)) ~ PLOIDY-1,
                                               phy = genome.ploidy.II6invasive.phy, data = genome_ploidy_II6_invasive_data, 
                                               method = "lambda")

summary(genome.ploidy.II6.invasiveent.model)
plot(genome.ploidy.II6.invasiveent.model)
## model Residual Plot
qplot(range_1.1(log10(genome_ploidy_II6_invasive_data$reg_no_invasive)), residuals(genome.ploidy.II6.invasiveent.model))
hist(residuals(genome.ploidy.II6.invasiveent.model))
plot(density(resid(genome.ploidy.II6.invasiveent.model))) #A density plot
# A quantile normal plot - good for checking normality
qqnorm(resid(genome.ploidy.II6.invasiveent.model)) 
qqline(resid(genome.ploidy.II6.invasiveent.model))

#### all the naturalized and invasive extent are not significant.

### prepare ployploids subset
genome_ploidy_II6_ploy_data <-  genome_ploidy_data  %>%  filter (PLOIDY == "polyploid") %>% droplevels()
summary(genome_ploidy_II6_ploy_data)  ## 2326 sp

##prepare the phy

keep.ploidII6_ployspp <-levels(genome_ploidy_II6_ploy_data$species)
remove_ploidyII6_ploytaxa = setdiff(genome.phy$tip.label, keep.ploidII6_ployspp)
genome.ploidy.II6_ploy.phy <- drop.tip(genome.phy,remove_ploidyII6_ploytaxa)

# Order data by tip order

genome_ploidy_II6_ploy_data = genome_ploidy_II6_ploy_data[(genome_ploidy_II6_ploy_data$species %in% genome.ploidy.II6_ploy.phy$tip.label), ]

row.names(genome_ploidy_II6_ploy_data) = genome_ploidy_II6_ploy_data$species
###check the names
library(geiger)  #for name.check
name.check(genome.ploidy.II6_ploy.phy, genome_ploidy_II6_ploy_data)

str(genome_ploidy_II6_ploy_data)
summary(genome_ploidy_II6_ploy_data)

genome.ploidy.II6_ploy.null.model <- phyloglm(GloNAF_incidence ~ 1, 
                                              phy = genome.ploidy.II6_ploy.phy, data = genome_ploidy_II6_ploy_data, 
                                              method = c("logistic_MPLE","logistic_IG10","poisson_GEE"),
                                              btol = 10, log.alpha.bound = 4,
                                              start.beta=NULL, start.alpha=NULL,
                                              boot = 0, full.matrix = TRUE)

genome.ploidy.II6_ploy.model <- phyloglm(GloNAF_incidence ~  no_ploidy_lev, 
                                         phy = genome.ploidy.II6_ploy.phy, data = genome_ploidy_II6_ploy_data, 
                                         method = c("logistic_MPLE","logistic_IG10","poisson_GEE"),
                                         btol = 10, log.alpha.bound = 4,
                                         start.beta=NULL, start.alpha=NULL,
                                         boot = 0, full.matrix = TRUE)

summary(genome.ploidy.II6_ploy.model)
R2(genome.ploidy.II6_ploy.model ,genome.ploidy.II6_ploy.null.model)


##### the interaction of genome size and ploidy level

## prepare the data
summary(genome_ploidy_data)
summary(as.factor(genome_ploidy_data$no_ploidy_lev))

genome.ploidy.gs.model <- phyloglm(GloNAF_incidence ~ range_1.1(log10(no_ploidy_lev))*range_1.1(log10(monoploid)), 
                                   phy = genome.ploidy.phy, data = genome_ploidy_data, 
                                   method = c("logistic_MPLE","logistic_IG10","poisson_GEE"),
                                   btol = 100, log.alpha.bound = 4,
                                   start.beta=NULL, start.alpha=NULL,
                                   boot = 0, full.matrix = TRUE)

summary(genome.ploidy.gs.model)
R2(genome.ploidy.gs.model ,genome.ploidy.null.model)

plot(genome.ploidy.gs.model)
hist(residuals(genome.ploidy.gs.model))
plot(density(resid(genome.ploidy.gs.model))) #A density plot
# A quantile normal plot - good for checking normality
qqnorm(resid(genome.ploidy.gs.model)) 
qqline(resid(genome.ploidy.gs.model))

###################### addad Jan2020
### plot the first two ploidy level results with genome size


### prepare ployploids subset
genome_ploidy_first_2_level_inc_data <-  genome_ploidy_data  %>%  filter (no_ploidy_lev < 3) %>% droplevels()
summary(genome_ploidy_first_2_level_inc_data)  

##prepare the phy

keep.ploid_first2_inc_spp <-levels(genome_ploidy_first_2_level_inc_data$species)

remove_ploid_first2_inc_taxa = setdiff(genome.phy$tip.label, keep.ploid_first2_inc_spp)
genomeploid_first2_inc_phy <- drop.tip(genome.phy,remove_ploid_first2_inc_taxa)

# Order data by tip order

genome_ploidy_first_2_level_inc_data = genome_ploidy_first_2_level_inc_data[(genome_ploidy_first_2_level_inc_data$species %in% genomeploid_first2_inc_phy$tip.label), ]

row.names(genome_ploidy_first_2_level_inc_data) = genome_ploidy_first_2_level_inc_data$species
###check the names
library(geiger)  #for name.check
name.check(genomeploid_first2_inc_phy, genome_ploidy_first_2_level_inc_data)

str(genome_ploidy_first_2_level_inc_data)
summary(genome_ploidy_first_2_level_inc_data)


genome_ploidy_first_2_level_inc_data$no_ploidy_lev <- as.factor(genome_ploidy_first_2_level_inc_data$no_ploidy_lev)

genome.ploidy.gs.first2_inc.model  <- phyloglm(GloNAF_incidence ~ no_ploidy_lev*range_1.1(log10(monoploid))-1, 
                                               phy = genomeploid_first2_inc_phy, data = genome_ploidy_first_2_level_inc_data, 
                                               method = c("logistic_MPLE","logistic_IG10","poisson_GEE"),
                                               btol = 100, log.alpha.bound = 4,
                                               start.beta=NULL, start.alpha=NULL,
                                               boot = 0, full.matrix = TRUE)

summary(genome.ploidy.gs.first2_inc.model)


genome.ploidy.gs.first2_inc.q.model  <- phyloglm(GloNAF_incidence ~ no_ploidy_lev + I(range_1.1(log10(monoploid))^2) +
                                                   no_ploidy_lev:range_1.1(log10(monoploid)) + no_ploidy_lev:I(range_1.1(log10(monoploid))^2)-1, 
                                                 phy = genomeploid_first2_inc_phy, data = genome_ploidy_first_2_level_inc_data, 
                                                 method = c("logistic_MPLE","logistic_IG10","poisson_GEE"),
                                                 btol = 100, log.alpha.bound = 4,
                                                 start.beta=NULL, start.alpha=NULL,
                                                 boot = 0, full.matrix = TRUE)

summary(genome.ploidy.gs.first2_inc.q.model)



### prepare ployploids subset
genome_ploidy_first_1_level_data <-  genome_ploidy_data  %>%  filter (no_ploidy_lev < 2) %>% droplevels()
summary(genome_ploidy_first_1_level_data)  

##prepare the phy

keep.ploid_first1_inc_spp <-levels(genome_ploidy_first_1_level_data$species)

remove_ploid_first1_inc_taxa = setdiff(genome.phy$tip.label, keep.ploid_first1_inc_spp)
genomeploid_first1_inc_phy <- drop.tip(genome.phy,remove_ploid_first1_inc_taxa)

# Order data by tip order

genome_ploidy_first_1_level_data = genome_ploidy_first_1_level_data[(genome_ploidy_first_1_level_data$species %in% genomeploid_first1_inc_phy$tip.label), ]

row.names(genome_ploidy_first_1_level_data) = genome_ploidy_first_1_level_data$species
###check the names
library(geiger)  #for name.check
name.check(genomeploid_first1_inc_phy, genome_ploidy_first_1_level_data)

str(genome_ploidy_first_1_level_data)
summary(genome_ploidy_first_1_level_data)

genome_ploidy_first_1_level_data$no_ploidy_lev <- as.factor(genome_ploidy_first_1_level_data$no_ploidy_lev)

genome.ploidy.gs.first1.model <- phyloglm(GloNAF_incidence ~ range_1.1(log10(monoploid)), 
                                          phy = genomeploid_first1_inc_phy, data = genome_ploidy_first_1_level_data, 
                                          method = c("logistic_MPLE","logistic_IG10","poisson_GEE"),
                                          btol = 100, log.alpha.bound = 4,
                                          start.beta=NULL, start.alpha=NULL,
                                          boot = 0, full.matrix = TRUE)

summary(genome.ploidy.gs.first1.model)

### prepare ployploids subset
genome_ploidy_first2_level_data <-  genome_ploidy_first_2_level_data  %>%  filter (no_ploidy_lev == 2) %>% droplevels()
summary(genome_ploidy_first2_level_data)  

##prepare the phy

keep.ploid_first2_inc_spp <-levels(genome_ploidy_first2_level_data$species)

remove_ploid_first2_inc_taxa = setdiff(genome.phy$tip.label, keep.ploid_first2_inc_spp)
genomeploid_first2_inc_phy <- drop.tip(genome.phy,remove_ploid_first2_inc_taxa)

# Order data by tip order

genome_ploidy_first2_level_data = genome_ploidy_first2_level_data[(genome_ploidy_first2_level_data$species %in% genomeploid_first2_inc_phy$tip.label), ]

row.names(genome_ploidy_first2_level_data) = genome_ploidy_first2_level_data$species
###check the names
library(geiger)  #for name.check
name.check(genomeploid_first2_inc_phy, genome_ploidy_first2_level_data)

str(genome_ploidy_first2_level_data)
summary(genome_ploidy_first2_level_data)

genome.ploidy.gs.first_2.model <- phyloglm(GloNAF_incidence ~ range_1.1(log10(monoploid)) , 
                                           phy = genomeploid_first2_inc_phy, data = genome_ploidy_first2_level_data, 
                                           method = c("logistic_MPLE","logistic_IG10","poisson_GEE"),
                                           btol = 100, log.alpha.bound = 4,
                                           start.beta=NULL, start.alpha=NULL,
                                           boot = 0, full.matrix = TRUE)

summary(genome.ploidy.gs.first_2.model)


#### naturalization extent global

summary(genome_ploidy_ext_data)
summary(as.factor(genome_ploidy_ext_data$no_ploidy_lev))


genome.ploidy.gs.extent.model <- phylolm(range_1.1(log10(reg_nat_total)) ~ range_1.1(log10(no_ploidy_lev))*range_1.1(log10(monoploid)),
                                         phy = genome.ploidy.ext.phy, data = genome_ploidy_ext_data, 
                                         method = "lambda")

summary(genome.ploidy.gs.extent.model)
#R2(genome.monoploid.main.extent.model, phy = genome.extent.main.phy)
plot(genome.ploidy.gs.extent.model)
## model Residual Plot
qplot(range_1.1(log10(genome_ploidy_ext_data$reg_nat_total)), residuals(genome.ploidy.gs.extent.model))
hist(residuals(genome.ploidy.gs.extent.model))
plot(density(resid(genome.ploidy.gs.extent.model))) #A density plot
# A quantile normal plot - good for checking normality
qqnorm(resid(genome.ploidy.gs.extent.model)) 
qqline(resid(genome.ploidy.gs.extent.model))

#### prepare the results for the interaction


### prepare ployploids subset
genome_ploidy_first_2_level_ext_data <-  genome_ploidy_ext_data  %>%  filter (no_ploidy_lev < 3) %>% droplevels()
summary(genome_ploidy_first_2_level_ext_data)  

##prepare the phy

keep.ploid_first2_ext_spp <-levels(genome_ploidy_first_2_level_ext_data$species)

remove_ploid_first2_ext_taxa = setdiff(genome.ploidy.ext.phy$tip.label, keep.ploid_first2_ext_spp)
genomeploid_first2_ext_phy <- drop.tip(genome.ploidy.ext.phy,remove_ploid_first2_ext_taxa)

# Order data by tip order

genome_ploidy_first_2_level_ext_data = genome_ploidy_first_2_level_ext_data[(genome_ploidy_first_2_level_ext_data$species %in% genomeploid_first2_ext_phy$tip.label), ]

row.names(genome_ploidy_first_2_level_ext_data) = genome_ploidy_first_2_level_ext_data$species
###check the names
library(geiger)  #for name.check
name.check(genomeploid_first2_ext_phy, genome_ploidy_first_2_level_ext_data)

str(genome_ploidy_first_2_level_ext_data)
summary(genome_ploidy_first_2_level_ext_data)


genome_ploidy_first_2_level_ext_data$no_ploidy_lev <- as.factor(genome_ploidy_first_2_level_ext_data$no_ploidy_lev)

genome.ploidy.gs.first2_ext.model  <- phylolm(range_1.1(log10(reg_nat_total)) ~ no_ploidy_lev*range_1.1(log10(monoploid))-1,
                                              phy = genomeploid_first2_ext_phy, data = genome_ploidy_first_2_level_ext_data, 
                                              method = "lambda")

summary(genome.ploidy.gs.first2_ext.model)


### prepare each level 
### prepare ployploids subset
genome_ploidy_first_1_level_ext_data <-  genome_ploidy_first_2_level_ext_data  %>%  filter (no_ploidy_lev == 1) %>% droplevels()
summary(genome_ploidy_first_1_level_ext_data)  

##prepare the phy

keep.ploid_first1_inc_spp <-levels(genome_ploidy_first_1_level_ext_data$species)

remove_ploid_first1_inc_taxa = setdiff(genome.ploidy.ext.phy$tip.label, keep.ploid_first1_inc_spp)
genomeploid_first1_ext_phy <- drop.tip(genome.ploidy.ext.phy,remove_ploid_first1_inc_taxa)

# Order data by tip order

genome_ploidy_first_1_level_ext_data = genome_ploidy_first_1_level_ext_data[(genome_ploidy_first_1_level_ext_data$species %in% genomeploid_first1_ext_phy$tip.label), ]

row.names(genome_ploidy_first_1_level_ext_data) = genome_ploidy_first_1_level_ext_data$species
###check the names
library(geiger)  #for name.check
name.check(genomeploid_first1_ext_phy, genome_ploidy_first_1_level_ext_data)

str(genome_ploidy_first_1_level_ext_data)
summary(genome_ploidy_first_1_level_ext_data)

genome.ploidy.gs.first_1_ext.model <- phylolm(range_1.1(log10(reg_nat_total)) ~ range_1.1(log10(monoploid)),
                                              phy = genomeploid_first1_ext_phy, data = genome_ploidy_first_1_level_ext_data, 
                                              method = "lambda")

summary(genome.ploidy.gs.first_1_ext.model)


### prepare ployploids subset:level 2
genome_ploidy_first_1_2_level_ext_data <-  genome_ploidy_first_2_level_ext_data  %>%  filter (no_ploidy_lev == 2) %>% droplevels()
summary(genome_ploidy_first_1_2_level_ext_data)  

##prepare the phy

keep.ploid_first1_2_inc_spp <-levels(genome_ploidy_first_1_2_level_ext_data$species)

remove_ploid_first1_2_inc_taxa = setdiff(genome.ploidy.ext.phy$tip.label, keep.ploid_first1_2_inc_spp)
genomeploid_first1_2_ext_phy <- drop.tip(genome.ploidy.ext.phy,remove_ploid_first1_2_inc_taxa)

# Order data by tip order

genome_ploidy_first_1_2_level_ext_data = genome_ploidy_first_1_2_level_ext_data[(genome_ploidy_first_1_2_level_ext_data$species %in% genomeploid_first1_2_ext_phy$tip.label), ]

row.names(genome_ploidy_first_1_2_level_ext_data) = genome_ploidy_first_1_2_level_ext_data$species
###check the names
library(geiger)  #for name.check
name.check(genomeploid_first1_2_ext_phy, genome_ploidy_first_1_2_level_ext_data)

str(genome_ploidy_first_1_2_level_ext_data)
summary(genome_ploidy_first_1_2_level_ext_data)

genome.ploidy.gs.first_1_2_ext.model <- phylolm(range_1.1(log10(reg_nat_total)) ~ range_1.1(log10(monoploid)),
                                                phy = genomeploid_first1_2_ext_phy, data = genome_ploidy_first_1_2_level_ext_data, 
                                                method = "lambda")

summary(genome.ploidy.gs.first_1_2_ext.model)



### naturalization extent mainland

summary(genome_ploidy_mainext_data)
summary(as.factor(genome_ploidy_mainext_data$no_ploidy_lev))
#### mainland holoploid

genome.ploidy.gs.mainextent.model <- phylolm(range_1.1(log10(reg_nat_main)) ~ no_ploidy_lev*range_1.1(log10(monoploid))-1,
                                             phy = genome.ploidy.mainext.phy, data = genome_ploidy_mainext_data, 
                                             method = "lambda")

summary(genome.ploidy.gs.mainextent.model)
plot(genome.ploidy.gs.mainextent.model)
## model Residual Plot
qplot(range_1.1(log10(genome_ploidy_mainext_data$reg_nat_main)), residuals(genome.ploidy.gs.mainextent.model))
hist(residuals(genome.ploidy.gs.mainextent.model))
plot(density(resid(genome.ploidy.gs.mainextent.model))) #A density plot
# A quantile normal plot - good for checking normality
qqnorm(resid(genome.ploidy.gs.mainextent.model)) 
qqline(resid(genome.ploidy.gs.mainextent.model))

##### check the interaction  (Mainland naturalization extent)


### prepare ployploids subset
genome_ploidy_first_2_level_main_ext_data <-  genome_ploidy_mainext_data  %>%  filter (no_ploidy_lev < 3) %>% droplevels()
summary(genome_ploidy_first_2_level_main_ext_data)  

##prepare the phy

keep.ploid_first2_main_ext_spp <-levels(genome_ploidy_first_2_level_main_ext_data$species)

remove_ploid_first2_main_ext_taxa = setdiff(genome.ploidy.ext.phy$tip.label, keep.ploid_first2_main_ext_spp)
genomeploid_first2_main_ext_phy <- drop.tip(genome.ploidy.ext.phy,remove_ploid_first2_main_ext_taxa)

# Order data by tip order

genome_ploidy_first_2_level_main_ext_data = genome_ploidy_first_2_level_main_ext_data[(genome_ploidy_first_2_level_main_ext_data$species %in% genomeploid_first2_main_ext_phy$tip.label), ]

row.names(genome_ploidy_first_2_level_main_ext_data) = genome_ploidy_first_2_level_main_ext_data$species
###check the names
library(geiger)  #for name.check
name.check(genomeploid_first2_main_ext_phy, genome_ploidy_first_2_level_main_ext_data)

str(genome_ploidy_first_2_level_main_ext_data)
summary(genome_ploidy_first_2_level_main_ext_data)


genome_ploidy_first_2_level_main_ext_data$no_ploidy_lev <- as.factor(genome_ploidy_first_2_level_main_ext_data$no_ploidy_lev)

genome.ploidy.gs.first2_main_ext.model  <- phylolm(range_1.1(log10(reg_nat_main)) ~ no_ploidy_lev*range_1.1(log10(monoploid))-1,
                                                   phy = genomeploid_first2_main_ext_phy, data = genome_ploidy_first_2_level_main_ext_data, 
                                                   method = "lambda")

summary(genome.ploidy.gs.first2_main_ext.model)


### naturalization extent island
summary(genome_ploidy_islext_data)
summary(as.factor(genome_ploidy_islext_data$no_ploidy_lev))


genome.ploidy.gs.islextent.model <- phylolm(range_1.1(log10(reg_nat_island)) ~ no_ploidy_lev*range_1.1(log10(monoploid)),
                                            phy = genome.ploidy.islext.phy, data = genome_ploidy_islext_data, 
                                            method = "lambda")

summary(genome.ploidy.gs.islextent.model)
plot(genome.ploidy.gs.islextent.model)
## model Residual Plot
qplot(range_1.1(log10(genome_ploidy_islext_data$reg_nat_island)), residuals(genome.ploidy.gs.islextent.model))
hist(residuals(genome.ploidy.gs.islextent.model))
plot(density(resid(genome.ploidy.gs.islextent.model))) #A density plot
# A quantile normal plot - good for checking normality
qqnorm(resid(genome.ploidy.gs.islextent.model)) 
qqline(resid(genome.ploidy.gs.islextent.model))

###############################
genome_ploidy_first_2_level_island_ext_data <-  genome_ploidy_islext_data  %>%  filter (no_ploidy_lev < 3) %>% droplevels()
summary(genome_ploidy_first_2_level_island_ext_data)  

##prepare the phy

keep.ploid_first2_island_ext_spp <-levels(genome_ploidy_first_2_level_island_ext_data$species)

remove_ploid_first2_island_ext_taxa = setdiff(genome.ploidy.ext.phy$tip.label, keep.ploid_first2_island_ext_spp)
genomeploid_first2_island_ext_phy <- drop.tip(genome.ploidy.ext.phy,remove_ploid_first2_island_ext_taxa)

# Order data by tip order

genome_ploidy_first_2_level_island_ext_data = genome_ploidy_first_2_level_island_ext_data[(genome_ploidy_first_2_level_island_ext_data$species %in% genomeploid_first2_island_ext_phy$tip.label), ]

row.names(genome_ploidy_first_2_level_island_ext_data) = genome_ploidy_first_2_level_island_ext_data$species
###check the names
library(geiger)  #for name.check
name.check(genomeploid_first2_island_ext_phy, genome_ploidy_first_2_level_island_ext_data)

str(genome_ploidy_first_2_level_island_ext_data)
summary(genome_ploidy_first_2_level_island_ext_data)


genome_ploidy_first_2_level_island_ext_data$no_ploidy_lev <- as.factor(genome_ploidy_first_2_level_island_ext_data$no_ploidy_lev)

genome.ploidy.gs.first2_island_ext.model  <- phylolm(range_1.1(log10(reg_nat_island)) ~ no_ploidy_lev*range_1.1(log10(monoploid))-1,
                                                     phy = genomeploid_first2_island_ext_phy, data = genome_ploidy_first_2_level_island_ext_data, 
                                                     method = "lambda")

summary(genome.ploidy.gs.first2_island_ext.model)


### invasive extent
summary(genome_ploidy_invasive_data)
summary(as.factor(genome_ploidy_invasive_data$no_ploidy_lev))


genome.ploidy.gs.invasiveent.model <- phylolm(range_1.1(log10(reg_no_invasive)) ~ no_ploidy_lev*range_1.1(log10(monoploid)),
                                              phy = genome.ploidy.invasive.phy, data = genome_ploidy_invasive_data, 
                                              method = "lambda")

summary(genome.ploidy.gs.invasiveent.model)
plot(genome.ploidy.gs.invasiveent.model)
## model Residual Plot
qplot(range_1.1(log10(genome_ploidy_invasive_data$reg_no_invasive)), residuals(genome.ploidy.gs.invasiveent.model))
hist(residuals(genome.ploidy.gs.invasiveent.model))
plot(density(resid(genome.ploidy.gs.invasiveent.model))) #A density plot
# A quantile normal plot - good for checking normality
qqnorm(resid(genome.ploidy.gs.invasiveent.model)) 
qqline(resid(genome.ploidy.gs.invasiveent.model))

##################################################
## check the interaction between monoploid and no of ploidy level
### prepare ployploids subset
genome_ploidy_first_2_level_inv_data <-  genome_ploidy_invasive_data  %>%  filter (no_ploidy_lev < 3) %>% droplevels()
summary(genome_ploidy_first_2_level_inv_data)  

##prepare the phy

keep.ploid_first2_inv_spp <-levels(genome_ploidy_first_2_level_inv_data$species)

remove_ploid_first2_inv_taxa = setdiff(genome.ploidy.invasive.phy$tip.label, keep.ploid_first2_inv_spp)
genomeploid_first2_inv_phy <- drop.tip(genome.ploidy.invasive.phy,remove_ploid_first2_inv_taxa)

# Order data by tip order

genome_ploidy_first_2_level_inv_data = genome_ploidy_first_2_level_inv_data[(genome_ploidy_first_2_level_inv_data$species %in% genomeploid_first2_inv_phy$tip.label), ]

row.names(genome_ploidy_first_2_level_inv_data) = genome_ploidy_first_2_level_inv_data$species
###check the names
library(geiger)  #for name.check
name.check(genomeploid_first2_inv_phy, genome_ploidy_first_2_level_inv_data)

str(genome_ploidy_first_2_level_inv_data)
summary(genome_ploidy_first_2_level_inv_data)


genome_ploidy_first_2_level_inv_data$no_ploidy_lev <- as.factor(genome_ploidy_first_2_level_inv_data$no_ploidy_lev)

genome.ploidy.gs.first2_inv.model  <- phylolm(range_1.1(log10(reg_no_invasive)) ~ no_ploidy_lev*range_1.1(log10(monoploid))-1,
                                              phy = genomeploid_first2_inv_phy, data = genome_ploidy_first_2_level_inv_data, 
                                              method = "lambda")

summary(genome.ploidy.gs.first2_inv.model)


### prepare each level to check the interactions
### prepare ployploids subset
genome_ploidy_first_1_level_inv_data <-  genome_ploidy_first_2_level_inv_data  %>%  filter (no_ploidy_lev == 1) %>% droplevels()
summary(genome_ploidy_first_1_level_inv_data)  

##prepare the phy

keep.ploid_first_1_inc_spp <-levels(genome_ploidy_first_1_level_inv_data$species)

remove_ploid_first_1_inc_taxa = setdiff(genome.ploidy.invasive.phy$tip.label, keep.ploid_first_1_inc_spp)
genomeploid_first_1_ext_phy <- drop.tip(genome.ploidy.invasive.phy,remove_ploid_first_1_inc_taxa)

# Order data by tip order

genome_ploidy_first_1_level_inv_data = genome_ploidy_first_1_level_inv_data[(genome_ploidy_first_1_level_inv_data$species %in% genomeploid_first_1_ext_phy$tip.label), ]

row.names(genome_ploidy_first_1_level_inv_data) = genome_ploidy_first_1_level_inv_data$species
###check the names
library(geiger)  #for name.check
name.check(genomeploid_first_1_ext_phy, genome_ploidy_first_1_level_inv_data)

str(genome_ploidy_first_1_level_inv_data)
summary(genome_ploidy_first_1_level_inv_data)

genome.ploidy.gs.first_1_inv.model <- phylolm(range_1.1(log10(reg_no_invasive)) ~ range_1.1(log10(monoploid)),
                                              phy = genomeploid_first_1_ext_phy, data = genome_ploidy_first_1_level_inv_data, 
                                              method = "lambda")

summary(genome.ploidy.gs.first_1_inv.model)


### prepare ployploids subset:level 2
genome_ploidy_first_1_2_level_inv_data <-  genome_ploidy_first_2_level_inv_data  %>%  filter (no_ploidy_lev == 2) %>% droplevels()
summary(genome_ploidy_first_1_2_level_inv_data)  

##prepare the phy

keep.ploid_first1_2_inv_spp <-levels(genome_ploidy_first_1_2_level_inv_data$species)

remove_ploid_first1_2_inv_taxa = setdiff(genome.ploidy.invasive.phy$tip.label, keep.ploid_first1_2_inv_spp)
genomeploid_first1_2_inv_phy <- drop.tip(genome.ploidy.invasive.phy,remove_ploid_first1_2_inv_taxa)

# Order data by tip order

genome_ploidy_first_1_2_level_inv_data = genome_ploidy_first_1_2_level_inv_data[(genome_ploidy_first_1_2_level_inv_data$species %in% genomeploid_first1_2_inv_phy$tip.label), ]

row.names(genome_ploidy_first_1_2_level_inv_data) = genome_ploidy_first_1_2_level_inv_data$species
###check the names
library(geiger)  #for name.check
name.check(genomeploid_first1_2_inv_phy, genome_ploidy_first_1_2_level_inv_data)

str(genome_ploidy_first_1_2_level_inv_data)
summary(genome_ploidy_first_1_2_level_inv_data)

genome.ploidy.gs.first_1_2_inv.model <- phylolm(range_1.1(log10(reg_no_invasive)) ~ range_1.1(log10(monoploid)),
                                                phy = genomeploid_first1_2_inv_phy, data = genome_ploidy_first_1_2_level_inv_data, 
                                                method = "lambda")

summary(genome.ploidy.gs.first_1_2_inv.model)





##################################################################################################
#############  check the effects of ploidy level
data_genome_polyploid_diploid <- read.csv("D:/file_EcoInfo_group_04June2021/Dropbox/postdoc2019/03Pysek_genome/data_analysis/data_genome_polyploid_diploid.csv", header = T, sep = ';')
str(data_genome_polyploid_diploid)
summary(data_genome_polyploid_diploid)

## naturalization incidence

##prepare the phy

keep.spp_ploidy_incid <-levels(data_genome_polyploid_diploid$species)

# keepb489tips <- genome.phy$tip.label[match(keep.b489spp, genome.phy$tip.label)]

remove_ploidy_incid_taxa = setdiff(genome.phy$tip.label, keep.spp_ploidy_incid)
genome_ploidy_incid_phy <- drop.tip(genome.phy,remove_ploidy_incid_taxa)

# Order data by tip order

data_genome_polyploid_diploid = data_genome_polyploid_diploid[(data_genome_polyploid_diploid$species %in% genome_ploidy_incid_phy$tip.label), ]

row.names(data_genome_polyploid_diploid) = data_genome_polyploid_diploid$species
###check the names
library(geiger)  #for name.check
name.check(genome_ploidy_incid_phy, data_genome_polyploid_diploid)

str(data_genome_polyploid_diploid)
summary(data_genome_polyploid_diploid)

library(phylolm)
genome.ploidy.di_poly_incid.model <- phyloglm(GloNAF_incidence ~  PLOIDY, 
                                              phy = genome_ploidy_incid_phy, data = data_genome_polyploid_diploid, 
                                              method = c("logistic_MPLE","logistic_IG10","poisson_GEE"),
                                              btol = 10, log.alpha.bound = 4,
                                              start.beta=NULL, start.alpha=NULL,
                                              boot = 0, full.matrix = TRUE)

summary(genome.ploidy.di_poly_incid.model)


genome.ploidy.di_poly_incid.model2 <- phyloglm(GloNAF_incidence ~  PLOIDY1-1, 
                                               phy = genome_ploidy_incid_phy, data = data_genome_polyploid_diploid, 
                                               method = c("logistic_MPLE","logistic_IG10","poisson_GEE"),
                                               btol = 10, log.alpha.bound = 4,
                                               start.beta=NULL, start.alpha=NULL,
                                               boot = 0, full.matrix = TRUE)

summary(genome.ploidy.di_poly_incid.model2)

### plot the estimates of the models
#################### plot the results
results.summary.di.polyploidy.incidence <- read.csv("results.summary.di.polyploidy.incidence.csv", header = T, sep = ';')
str(results.summary.di.polyploidy.incidence)


# Error bars represent standard error of the mean
library(ggplot2)

# The errorbars overlapped, so use position_dodge to move them horizontally
pd <- position_dodge(0.5) # move them .05 to the left and right
p_diploid_c <- ggplot(results.summary.di.polyploidy.incidence, aes(x= group, y= estimate)) + 
  geom_point(position= pd, size=3.5) +
  geom_errorbar(aes(ymin= estimate - SE, ymax= estimate + SE),
                width=.1, position=pd) +
  
  geom_hline(yintercept=0, linetype =2) +
  ggtitle("Naturalization incidence") +
  # scale_color_manual(values=c("#fdb863", "#5e3c99")) +    ###add the dots
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     strip.background = element_rect(colour="white", fill="white"),
                     text = element_text(size = 16),
                     legend.position = c(0.9,0.8))

pdf("Fig.6 di_polyploid naturalization incidence.pdf", useDingbats=FALSE, width=8, height=6)


dev.off()



### naturalization region

library(dplyr)
data_genome_polyploid_diploid_nat_reg <- data_genome_polyploid_diploid  %>%  filter (reg_nat_total > 0) %>% droplevels()
summary(data_genome_polyploid_diploid_nat_reg)  ##2585

##prepare the phy

keep.spp_ploidy_extent <-levels(data_genome_polyploid_diploid_nat_reg$species)

# keepb489tips <- genome.phy$tip.label[match(keep.b489spp, genome.phy$tip.label)]

remove_ploidy_extent_taxa = setdiff(genome.phy$tip.label, keep.spp_ploidy_extent)
genome_ploidy_extent_phy <- drop.tip(genome.phy,remove_ploidy_extent_taxa)

# Order data by tip order

data_genome_polyploid_diploid_nat_reg = data_genome_polyploid_diploid_nat_reg[(data_genome_polyploid_diploid_nat_reg$species %in% genome_ploidy_extent_phy$tip.label), ]

row.names(data_genome_polyploid_diploid_nat_reg) = data_genome_polyploid_diploid_nat_reg$species
###check the names
library(geiger)  #for name.check
name.check(genome_ploidy_extent_phy, data_genome_polyploid_diploid_nat_reg)

str(data_genome_polyploid_diploid_nat_reg)
summary(data_genome_polyploid_diploid_nat_reg)

genome.di_poly_ploidy.extent.model <- phylolm(range_1.1(log10(reg_nat_total)) ~ PLOIDY,
                                              phy = genome_ploidy_extent_phy, data = data_genome_polyploid_diploid_nat_reg, 
                                              method = "lambda")

summary(genome.di_poly_ploidy.extent.model)

genome.di_poly_ploidy.extent.model2 <- phylolm(range_1.1(log10(reg_nat_total)) ~ PLOIDY1-1,
                                               phy = genome_ploidy_extent_phy, data = data_genome_polyploid_diploid_nat_reg, 
                                               method = "lambda")

summary(genome.di_poly_ploidy.extent.model2)
anova(genome.di_poly_ploidy.extent.model_nb)

#R2(genome.monoploid.main.islextent.model, phy = genome.islextent.main.phy)
plot(genome.di_poly_ploidy.extent.model)
## model Residual Plot
qplot(range_1.1(log10(data_genome_polyploid_diploid_nat_reg$reg_nat_total)), residuals(genome.di_poly_ploidy.extent.model))
hist(residuals(genome.di_poly_ploidy.extent.model))
plot(density(resid(genome.di_poly_ploidy.extent.model))) #A density plot
# A quantile normal plot - good for checking normality
qqnorm(resid(genome.di_poly_ploidy.extent.model)) 
qqline(resid(genome.di_poly_ploidy.extent.model))

lm_di_poly_ploid <- lm(log10(data_genome_polyploid_diploid_nat_reg$reg_nat_total) ~ data_genome_polyploid_diploid_nat_reg$PLOIDY)

library(emmeans)

emmeans(lm_di_poly_ploid, pairwise ~ "PLOIDY", adjust = "bonferroni")

library(lmerTest)
####   negative binomial glm, due to the infraspecies and dual-type. 
genome.di_poly_ploidy.extent.model_nb <- glmer(range_1.1(log10(reg_nat_total)) ~ PLOIDY1-1 + (1|Genus), data_genome_polyploid_diploid_nat_reg, verbose = T)
summary(genome.di_poly_ploidy.extent.model_nb)
anova(genome.di_poly_ploidy.extent.model_nb)

##################################################
### mainland naturalization extent

library(dplyr)
data_genome_polyploid_diploid_main_nat_reg <- data_genome_polyploid_diploid  %>%  filter (reg_nat_main > 0) %>% droplevels()
summary(data_genome_polyploid_diploid_main_nat_reg)  ##2585

##prepare the phy

keep.spp_ploidy_main_extent <-levels(data_genome_polyploid_diploid_main_nat_reg$species)

# keepb489tips <- genome.phy$tip.label[match(keep.b489spp, genome.phy$tip.label)]

remove_ploidy_main_extent_taxa = setdiff(genome.phy$tip.label, keep.spp_ploidy_main_extent)
genome_ploidy_main_extent_phy <- drop.tip(genome.phy,remove_ploidy_main_extent_taxa)

# Order data by tip order

data_genome_polyploid_diploid_main_nat_reg = data_genome_polyploid_diploid_main_nat_reg[(data_genome_polyploid_diploid_main_nat_reg$species %in% genome_ploidy_main_extent_phy$tip.label), ]

row.names(data_genome_polyploid_diploid_main_nat_reg) = data_genome_polyploid_diploid_main_nat_reg$species
###check the names
library(geiger)  #for name.check
name.check(genome_ploidy_main_extent_phy, data_genome_polyploid_diploid_main_nat_reg)

str(data_genome_polyploid_diploid_main_nat_reg)
summary(data_genome_polyploid_diploid_main_nat_reg)

genome.di_poly_ploidy.main.extent.model <- phylolm(range_1.1(log10(reg_nat_main)) ~ PLOIDY,
                                                   phy = genome_ploidy_main_extent_phy, data = data_genome_polyploid_diploid_main_nat_reg, 
                                                   method = "lambda")

summary(genome.di_poly_ploidy.main.extent.model)

genome.di_poly_ploidy.main.extent.model2 <- phylolm(range_1.1(log10(reg_nat_main)) ~ PLOIDY1-1,
                                                    phy = genome_ploidy_main_extent_phy, data = data_genome_polyploid_diploid_main_nat_reg, 
                                                    method = "lambda")

summary(genome.di_poly_ploidy.main.extent.model2)
anova(genome.di_poly_ploidy.main.extent.model2)

#R2(genome.monoploid.main.islextent.model, phy = genome.islextent.main.phy)
plot(genome.di_poly_ploidy.main.extent.model)
## model Residual Plot
qplot(range_1.1(log10(data_genome_polyploid_diploid_main_nat_reg$reg_nat_main)), residuals(genome.di_poly_ploidy.main.extent.model))
hist(residuals(genome.di_poly_ploidy.main.extent.model))
plot(density(resid(genome.di_poly_ploidy.main.extent.model))) #A density plot
# A quantile normal plot - good for checking normality
qqnorm(resid(genome.di_poly_ploidy.main.extent.model)) 
qqline(resid(genome.di_poly_ploidy.main.extent.model))

### Island naturalization extent

data_genome_polyploid_diploid_island_nat_reg <- data_genome_polyploid_diploid  %>%  filter (reg_nat_island > 0) %>% droplevels()
summary(data_genome_polyploid_diploid_island_nat_reg)  ##2585

##prepare the phy

keep.spp_ploidy_island_extent <-levels(data_genome_polyploid_diploid_island_nat_reg$species)

# keepb489tips <- genome.phy$tip.label[match(keep.b489spp, genome.phy$tip.label)]

remove_ploidy_island_extent_taxa = setdiff(genome.phy$tip.label, keep.spp_ploidy_island_extent)
genome_ploidy_island_extent_phy <- drop.tip(genome.phy,remove_ploidy_island_extent_taxa)

# Order data by tip order

data_genome_polyploid_diploid_island_nat_reg = data_genome_polyploid_diploid_island_nat_reg[(data_genome_polyploid_diploid_island_nat_reg$species %in% genome_ploidy_island_extent_phy$tip.label), ]

row.names(data_genome_polyploid_diploid_island_nat_reg) = data_genome_polyploid_diploid_island_nat_reg$species
###check the names
library(geiger)  #for name.check
name.check(genome_ploidy_island_extent_phy, data_genome_polyploid_diploid_island_nat_reg)

str(data_genome_polyploid_diploid_island_nat_reg)
summary(data_genome_polyploid_diploid_island_nat_reg)

genome.di_poly_ploidy.island.extent.model <- phylolm(range_1.1(log10(reg_nat_island)) ~ PLOIDY,
                                                     phy = genome_ploidy_island_extent_phy, data = data_genome_polyploid_diploid_island_nat_reg, 
                                                     method = "lambda")

summary(genome.di_poly_ploidy.island.extent.model)

genome.di_poly_ploidy.island.extent.model2 <- phylolm(range_1.1(log10(reg_nat_island)) ~ PLOIDY1-1,
                                                      phy = genome_ploidy_island_extent_phy, data = data_genome_polyploid_diploid_island_nat_reg, 
                                                      method = "lambda")

summary(genome.di_poly_ploidy.island.extent.model2)
anova(genome.di_poly_ploidy.island.extent.model2)

#R2(genome.monoploid.island.islextent.model, phy = genome.islextent.island.phy)
plot(genome.di_poly_ploidy.island.extent.model)
## model Residual Plot
qplot(range_1.1(log10(data_genome_polyploid_diploid_island_nat_reg$reg_nat_island)), residuals(genome.di_poly_ploidy.island.extent.model))
hist(residuals(genome.di_poly_ploidy.island.extent.model))
plot(density(resid(genome.di_poly_ploidy.island.extent.model))) #A density plot
# A quantile normal plot - good for checking normality
qqnorm(resid(genome.di_poly_ploidy.island.extent.model)) 
qqline(resid(genome.di_poly_ploidy.island.extent.model))

### invasion extent
data_genome_polyploid_diploid_invasion_nat_reg <- data_genome_polyploid_diploid  %>%  filter (reg_no_invasive > 0) %>% droplevels()
summary(data_genome_polyploid_diploid_invasion_nat_reg)  ##2585

##prepare the phy

keep.spp_ploidy_invasion_extent <-levels(data_genome_polyploid_diploid_invasion_nat_reg$species)

# keepb489tips <- genome.phy$tip.label[match(keep.b489spp, genome.phy$tip.label)]

remove_ploidy_invasion_extent_taxa = setdiff(genome.phy$tip.label, keep.spp_ploidy_invasion_extent)
genome_ploidy_invasion_extent_phy <- drop.tip(genome.phy,remove_ploidy_invasion_extent_taxa)

# Order data by tip order

data_genome_polyploid_diploid_invasion_nat_reg = data_genome_polyploid_diploid_invasion_nat_reg[(data_genome_polyploid_diploid_invasion_nat_reg$species %in% genome_ploidy_invasion_extent_phy$tip.label), ]

row.names(data_genome_polyploid_diploid_invasion_nat_reg) = data_genome_polyploid_diploid_invasion_nat_reg$species
###check the names
library(geiger)  #for name.check
name.check(genome_ploidy_invasion_extent_phy, data_genome_polyploid_diploid_invasion_nat_reg)

str(data_genome_polyploid_diploid_invasion_nat_reg)
summary(data_genome_polyploid_diploid_invasion_nat_reg)

genome.di_poly_ploidy.invasion.extent.model <- phylolm(range_1.1(log10(reg_no_invasive)) ~ PLOIDY,
                                                       phy = genome_ploidy_invasion_extent_phy, data = data_genome_polyploid_diploid_invasion_nat_reg, 
                                                       method = "lambda")

summary(genome.di_poly_ploidy.invasion.extent.model)

genome.di_poly_ploidy.invasion.extent.model2 <- phylolm(range_1.1(log10(reg_no_invasive)) ~ PLOIDY1-1,
                                                        phy = genome_ploidy_invasion_extent_phy, data = data_genome_polyploid_diploid_invasion_nat_reg, 
                                                        method = "lambda")

summary(genome.di_poly_ploidy.invasion.extent.model2)
anova(genome.di_poly_ploidy.invasion.extent.model2)

#R2(genome.monoploid.invasion.islextent.model, phy = genome.islextent.invasion.phy)
plot(genome.di_poly_ploidy.invasion.extent.model)
## model Residual Plot
qplot(range_1.1(log10(data_genome_polyploid_diploid_invasion_nat_reg$reg_no_invasive)), residuals(genome.di_poly_ploidy.invasion.extent.model))
hist(residuals(genome.di_poly_ploidy.invasion.extent.model))
plot(density(resid(genome.di_poly_ploidy.invasion.extent.model))) #A density plot
# A quantile normal plot - good for checking normality
qqnorm(resid(genome.di_poly_ploidy.invasion.extent.model)) 
qqline(resid(genome.di_poly_ploidy.invasion.extent.model))


### received Petr's reply and he would like to focus on ploidy levels with enough species size, thus
### remove 5, 7, 9, and combine all bigger than 10 as one group
## I combined all the left ploidy levels in one df.
## as there are many species having more than one ploidy levels, I could not use the phylogeny in analysis if I
## compare the differences of ploidy levels.

ploidy_level_data <- read.csv("data_ploidy_level_group_31Jan2020.csv", header = T , sep = ';')

summary(ploidy_level_data)
ploidy_level_data$ploidy_level_single <- as.factor(ploidy_level_data$ploidy_level_single)

library(lmerTest)
gs.ploidy.reg <- glmer(GloNAF_incidence ~ ploidy_level_single + (1|NewAPGFamily:Genus), family = binomial(link="logit"), 
                       ploidy_level_data, glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))

summary(gs.ploidy.reg)
# confint(gs.ploidy.reg)
# exp(coef(gs.ploidy.reg)) # exponentiated coefficients
# exp(confint(gs.ploidy.reg)) # 95% CI for exponentiated coefficients
library(car)
Anova(gs.ploidy.reg, type = "II")
plot(gs.ploidy.reg)
library(effects)
plot(allEffects(gs.ploidy.reg, typical=mean))
plot(predictorEffects(gs.ploidy.reg, ~ ploidy_level_single), lines=list(multiline=TRUE),
     axes=list(grid=TRUE),
     confint=list(style="bars"))

library(multcomp)
library(emmeans)
marginal<- emmeans(gs.ploidy.reg, pairwise ~ ploidy_level_single, adjust = "bonferroni") ## bonferroni mothod for posthoc test

#### test the residual distribution
library(DHARMa)
simulationOutput.1 <- simulateResiduals(fittedModel = gs.ploidy.reg)
testDispersion(simulationOutput.1)
plotQQunif(simulationOutput.1, testUniformity = T)
plot(simulationOutput.1)
plotNormalHistogram(residuals(gs.ploidy.reg))


CLD = CLD(marginal,
          alpha   = 0.05,
          Letters = letters,         ###  Use lowercase letters for .group
          adjust  = "bonferroni") 
### Order the levels for printing

CLD$ploidy_level_single = factor(CLD$ploidy_level_single,
                                 levels=c("2", "3", "4", "6", "8", "10"))
CLD$.group=gsub(" ", "", CLD$.group)

## plot the results from emmeans
# gs.ploidy.reg.estimate <- read.csv("gs.ploidy.reg.estimate.csv", header = T, sep = ';')
# The errorbars overlapped, so use position_dodge to move them horizontally
# pd <- position_dodge(0.5) # move them .05 to the left and right
p_a <- ggplot(CLD, aes(x= as.factor(ploidy_level_single), y= emmean,
                       label = .group)) + 
  geom_point(size=3.5) +
  geom_errorbar(aes(ymin= emmean - SE, ymax= emmean + SE),
                width=.1) +
  # geom_hline(yintercept=0) +
  ggtitle("(a) Naturalization incidence") +
  # scale_color_manual(values=c("#fdb863", "#5e3c99")) +    ###add the dots
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     strip.background = element_rect(colour="white", fill="white"),
                     text = element_text(size = 16),
                     legend.position = c(0.9,0.8)) + 
  ylab("Estimated marginal means") +
  xlab("Ploidy level") +
  geom_text(nudge_x = c(0, 0, 0, 0, 0, 0),
            nudge_y = c(0.31, 0.31, 0.31, 0.31, 0.31, 0.31),
            color   = "black")

pdf("Fig.7 ploidy_level naturalization incidence.pdf", useDingbats=FALSE, width=8, height=6)

dev.off()

#### then prepare the data for extent
summary(ploidy_level_data)

## first total naturalized extent
ploidy_level_tot_naturalized_ext_data <- ploidy_level_data %>%  filter (reg_nat_total > 0) %>% droplevels()

summary(ploidy_level_tot_naturalized_ext_data)  ##2585

library(lmerTest)
ploidy_level_tot_naturalized_ext_data.reg <- glmer.nb(reg_nat_total ~ ploidy_level_single + (1|NewAPGFamily:Genus), 
                                                      ploidy_level_tot_naturalized_ext_data,verbose = T,
                                                      glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))

summary(ploidy_level_tot_naturalized_ext_data.reg)

library(car)
Anova(ploidy_level_tot_naturalized_ext_data.reg)
plot(ploidy_level_tot_naturalized_ext_data.reg)
library(effects)
plot(allEffects(ploidy_level_tot_naturalized_ext_data.reg, typical=median))
library(multcomp)
library(emmeans)
marginal2 <- emmeans(ploidy_level_tot_naturalized_ext_data.reg, pairwise ~ ploidy_level_single, adjust = "bonferroni")

#### test the residual distribution
library(DHARMa)
simulationOutput.2 <- simulateResiduals(fittedModel = ploidy_level_tot_naturalized_ext_data.reg)
testDispersion(simulationOutput.2)
plotQQunif(simulationOutput.2, testUniformity = T)
plot(simulationOutput.2)


CLD_nat_tot = CLD(marginal2,
                  alpha   = 0.05,
                  Letters = letters,         ###  Use lowercase letters for .group
                  adjust  = "bonferroni") 
### Order the levels for printing

CLD_nat_tot$ploidy_level_single = factor(CLD_nat_tot$ploidy_level_single,
                                         levels=c("2", "3", "4", "6", "8", "10"))
CLD_nat_tot$.group=gsub(" ", "", CLD_nat_tot$.group)

## plot the results from emmeans
# gs.ploidy.reg.estimate <- read.csv("gs.ploidy.reg.estimate.csv", header = T, sep = ';')
# The errorbars overlapped, so use position_dodge to move them horizontally
# pd <- position_dodge(0.5) # move them .05 to the left and right
p_b <- ggplot(CLD_nat_tot, aes(x= as.factor(ploidy_level_single), y= emmean,
                               label = .group)) + 
  geom_point(size=3.5) +
  geom_errorbar(aes(ymin= emmean - SE, ymax= emmean + SE),
                width=.1) +
  # geom_hline(yintercept=0) +
  ggtitle("(b) Naturalization extent: Global") +
  # scale_color_manual(values=c("#fdb863", "#5e3c99")) +    ###add the dots
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     strip.background = element_rect(colour="white", fill="white"),
                     text = element_text(size = 16),
                     legend.position = c(0.9,0.8)) + 
  ylab("Estimated marginal means") +
  xlab("Ploidy level") +
  geom_text(nudge_x = c(0, 0, 0, 0, 0, 0),
            nudge_y = c(0.31, 0.31, 0.31, 0.31, 0.31, 0.31),
            color   = "black")


#### mainland extent
ploidy_level_main_naturalized_ext_data <- ploidy_level_data %>%  filter (reg_nat_main > 0) %>% droplevels()

summary(ploidy_level_main_naturalized_ext_data)

library(lmerTest)
ploidy_level_main_naturalized_ext_data.reg <- glmer.nb(reg_nat_main ~ ploidy_level_single + (1|NewAPGFamily:Genus), 
                                                       ploidy_level_main_naturalized_ext_data,verbose = T,
                                                       glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))

summary(ploidy_level_main_naturalized_ext_data.reg)

library(car)
Anova(ploidy_level_main_naturalized_ext_data.reg)
plot(ploidy_level_main_naturalized_ext_data.reg)
library(effects)
plot(allEffects(ploidy_level_main_naturalized_ext_data.reg, typical=mean))
library(multcomp)
library(emmeans)
marginal3 <-emmeans(ploidy_level_main_naturalized_ext_data.reg, pairwise ~ ploidy_level_single, adjust = "bonferroni")

#### test the residual distribution
library(DHARMa)
simulationOutput.3 <- simulateResiduals(fittedModel = ploidy_level_main_naturalized_ext_data.reg)
testDispersion(simulationOutput.3)
plotQQunif(simulationOutput.3, testUniformity = T)
plot(simulationOutput.3)

################################   https://rcompanion.org/handbook/I_08.html

CLD_nat_main = CLD(marginal3,
                   alpha   = 0.05,
                   Letters = letters,         ###  Use lowercase letters for .group
                   adjust  = "bonferroni") 
### Order the levels for printing

CLD_nat_main$ploidy_level_single = factor(CLD_nat_main$ploidy_level_single,
                                          levels=c("2", "3", "4", "6", "8", "10"))
CLD_nat_main$.group=gsub(" ", "", CLD_nat_main$.group)

## plot the results from emmeans
# gs.ploidy.reg.estimate <- read.csv("gs.ploidy.reg.estimate.csv", header = T, sep = ';')
# The errorbars overlapped, so use position_dodge to move them horizontally
# pd <- position_dodge(0.5) # move them .05 to the left and right
p_c <- ggplot(CLD_nat_main, aes(x= as.factor(ploidy_level_single), y= emmean,
                                label = .group)) + 
  geom_point(size=3.5) +
  geom_errorbar(aes(ymin= emmean - SE, ymax= emmean + SE),
                width=.1) +
  # geom_hline(yintercept=0) +
  ggtitle("(c) Naturalization extent: Mainland") +
  # scale_color_manual(values=c("#fdb863", "#5e3c99")) +    ###add the dots
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     strip.background = element_rect(colour="white", fill="white"),
                     text = element_text(size = 16),
                     legend.position = c(0.9,0.8)) + 
  ylab("Estimated marginal means") +
  xlab("Ploidy level") +
  geom_text(nudge_x = c(0, 0, 0, 0, 0, 0),
            nudge_y = c(0.31, 0.31, 0.31, 0.31, 0.31, 0.31),
            color   = "black")



#### island extent
ploidy_level_island_naturalized_ext_data <- ploidy_level_data %>%  filter (reg_nat_island > 0) %>% droplevels()

summary(ploidy_level_island_naturalized_ext_data)

library(lmerTest)
ploidy_level_island_naturalized_ext_data.reg <- glmer.nb(reg_nat_island ~ ploidy_level_single + (1|NewAPGFamily:Genus), 
                                                         ploidy_level_island_naturalized_ext_data,verbose = T,
                                                         glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))

summary(ploidy_level_island_naturalized_ext_data.reg)

library(car)
Anova(ploidy_level_island_naturalized_ext_data.reg)
plot(ploidy_level_island_naturalized_ext_data.reg)
library(effects)
plot(allEffects(ploidy_level_island_naturalized_ext_data.reg, typical=mean))
library(multcomp)
library(emmeans)
marginal4 <- emmeans(ploidy_level_island_naturalized_ext_data.reg, pairwise ~ ploidy_level_single, adjust = "bonferroni")

#### test the residual distribution
library(DHARMa)
simulationOutput.4 <- simulateResiduals(fittedModel = ploidy_level_island_naturalized_ext_data.reg)
testDispersion(simulationOutput.4)
plotQQunif(simulationOutput.4, testUniformity = T)
plot(simulationOutput.4)


CLD_nat_island = CLD(marginal4,
                     alpha   = 0.05,
                     Letters = letters,         ###  Use lowercase letters for .group
                     adjust  = "bonferroni") 
### Order the levels for printing

CLD_nat_island$ploidy_level_single = factor(CLD_nat_island$ploidy_level_single,
                                            levels=c("2", "3", "4", "6", "8", "10"))
CLD_nat_island$.group=gsub(" ", "", CLD_nat_island$.group)

## plot the results from emmeans
# gs.ploidy.reg.estimate <- read.csv("gs.ploidy.reg.estimate.csv", header = T, sep = ';')
# The errorbars overlapped, so use position_dodge to move them horizontally
# pd <- position_dodge(0.5) # move them .05 to the left and right
p_d <- ggplot(CLD_nat_island, aes(x= as.factor(ploidy_level_single), y= emmean,
                                  label = .group)) + 
  geom_point(size=3.5) +
  geom_errorbar(aes(ymin= emmean - SE, ymax= emmean + SE),
                width=.1) +
  # geom_hline(yintercept=0) +
  ggtitle("(c) Naturalization extent: Island") +
  # scale_color_manual(values=c("#fdb863", "#5e3c99")) +    ###add the dots
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     strip.background = element_rect(colour="white", fill="white"),
                     text = element_text(size = 16),
                     legend.position = c(0.9,0.8)) + 
  ylab("Estimated marginal means") +
  xlab("Ploidy level") +
  geom_text(nudge_x = c(0, 0, 0, 0, 0, 0),
            nudge_y = c(0.31, 0.31, 0.31, 0.31, 0.31, 0.31),
            color   = "black")

## put them as one figure and save it
library(cowplot)
pdf("Fig.SXX ploidy_level_4in1.pdf", useDingbats=FALSE, width=12, height=8)
p_fig_3d <- plot_grid(p_a, p_b, p_c, p_d, nrow =2, align = "h")
p_fig_3d_2 <- plot_grid(p_a, p_b, nrow =2, align = "h")


dev.off()


#### invasion extent
ploidy_level_invasive_ext_data <- ploidy_level_data %>%  filter (reg_no_invasive > 0) %>% droplevels()

summary(ploidy_level_invasive_ext_data)

library(lmerTest)
ploidy_level_invasive_ext_data.reg <- glmer.nb(reg_no_invasive ~ ploidy_level_single + (1|NewAPGFamily:Genus), 
                                               ploidy_level_invasive_ext_data,verbose = T,
                                               glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))

summary(ploidy_level_invasive_ext_data.reg)

library(car)
Anova(ploidy_level_invasive_ext_data.reg)  
#Response: reg_no_invasive
# Chisq Df Pr(>Chisq)  
# ploidy_level_single 9.4908  5    0.09102 .
plot(ploidy_level_invasive_ext_data.reg)
library(effects)
plot(allEffects(ploidy_level_invasive_ext_data.reg, typical=mean))
library(multcomp)
library(emmeans)
emmeans(ploidy_level_invasive_ext_data.reg, pairwise ~ ploidy_level_single, adjust = "bonferroni")

#### test the residual distribution
library(DHARMa)
simulationOutput.5 <- simulateResiduals(fittedModel = ploidy_level_invasive_ext_data.reg)
testDispersion(simulationOutput.5)
plotQQunif(simulationOutput.5, testUniformity = T)
plot(simulationOutput.5)

#### plot fig. 3




figure3 <-ggarrange(ggarrange(p_ploidy_nat_incidence,  p_ploidy_nat_extent,  ncol = 2, labels = c("A", "B")), 
                    ggarrange(p_diploid_c, p_fig_3d_2, ncol = 2, labels = c("C", "D")), # Second row with two extent plot
                    nrow = 2, widths = c(1, 1)  ## adjust height
)                                      # Labels of the incidence plot


pdf("Fig 3.ploidy_ incidence & extent_14June2021.pdf", useDingbats=FALSE, width=10, height=9)
figure3
dev.off()





#####  III interactions of genome size and no. of ploidy level

##### then move to the last step III using the interaction of genome size and ploidy level

## prepare the data
summary(genome_ploidy_data) ## it includes 3517 species without monoploid data, thus needs
## to remove them before running model
genome_ploidy_monoploid_final <- genome_ploidy_data %>% filter (monoploid != "NA") %>% droplevels() #7252 species

summary(as.factor(genome_ploidy_monoploid_final$no_ploidy_lev))
summary(genome_ploidy_monoploid_final)
genome.ploidy.gs.model <- phyloglm(GloNAF_incidence ~ range_1.1(log10(no_ploidy_lev))*range_1.1(log10(monoploid)), 
                                   phy = genome.ploidy.phy, data = genome_ploidy_monoploid_final, 
                                   method = c("logistic_MPLE","logistic_IG10","poisson_GEE"),
                                   btol = 100, log.alpha.bound = 4,
                                   start.beta=NULL, start.alpha=NULL,
                                   boot = 0, full.matrix = TRUE)

summary(genome.ploidy.gs.model)
R2(genome.ploidy.gs.model ,genome.ploidy.null.model)

plot(genome.ploidy.gs.model)
## model Residual Plot
# qplot(range_1.1(log10(genome_ploidy_data$GloNAF_incidence)), residuals(genome.ploidy.gs.model))
hist(residuals(genome.ploidy.gs.model))
plot(density(resid(genome.ploidy.gs.model))) #A density plot
# A quantile normal plot - good for checking normality
qqnorm(resid(genome.ploidy.gs.model)) 
qqline(resid(genome.ploidy.gs.model))

###################### addad Jan2020
### plot the first two ploidy level results with genome size


### prepare ployploids subset
genome_ploidy_first_2_level_inc_data <-  genome_ploidy_monoploid_final  %>%  filter (no_ploidy_lev < 3) %>% droplevels()
summary(genome_ploidy_first_2_level_inc_data)  

##prepare the phy

keep.ploid_first2_inc_spp <-levels(genome_ploidy_first_2_level_inc_data$species)

remove_ploid_first2_inc_taxa = setdiff(genome.phy$tip.label, keep.ploid_first2_inc_spp)
genomeploid_first2_inc_phy <- drop.tip(genome.phy,remove_ploid_first2_inc_taxa)

# Order data by tip order

genome_ploidy_first_2_level_inc_data = genome_ploidy_first_2_level_inc_data[(genome_ploidy_first_2_level_inc_data$species %in% genomeploid_first2_inc_phy$tip.label), ]

row.names(genome_ploidy_first_2_level_inc_data) = genome_ploidy_first_2_level_inc_data$species
###check the names
library(geiger)  #for name.check
name.check(genomeploid_first2_inc_phy, genome_ploidy_first_2_level_inc_data)

str(genome_ploidy_first_2_level_inc_data)
summary(genome_ploidy_first_2_level_inc_data)


genome_ploidy_first_2_level_inc_data$no_ploidy_lev <- as.factor(genome_ploidy_first_2_level_inc_data$no_ploidy_lev)

genome.ploidy.gs.first2_inc.model  <- phyloglm(GloNAF_incidence ~ no_ploidy_lev*range_1.1(log10(monoploid))-1, 
                                               phy = genomeploid_first2_inc_phy, data = genome_ploidy_first_2_level_inc_data, 
                                               method = c("logistic_MPLE","logistic_IG10","poisson_GEE"),
                                               btol = 100, log.alpha.bound = 4,
                                               start.beta=NULL, start.alpha=NULL,
                                               boot = 0, full.matrix = TRUE)

summary(genome.ploidy.gs.first2_inc.model)


genome.ploidy.gs.first2_inc.q.model  <- phyloglm(GloNAF_incidence ~ no_ploidy_lev + I(range_1.1(log10(monoploid))^2) +
                                                   no_ploidy_lev:range_1.1(log10(monoploid)) + no_ploidy_lev:I(range_1.1(log10(monoploid))^2)-1, 
                                                 phy = genomeploid_first2_inc_phy, data = genome_ploidy_first_2_level_inc_data, 
                                                 method = c("logistic_MPLE","logistic_IG10","poisson_GEE"),
                                                 btol = 100, log.alpha.bound = 4,
                                                 start.beta=NULL, start.alpha=NULL,
                                                 boot = 0, full.matrix = TRUE)

summary(genome.ploidy.gs.first2_inc.q.model)



### prepare ployploids subset
genome_ploidy_first_1_level_data <-  genome_ploidy_monoploid_final  %>%  filter (no_ploidy_lev < 2) %>% droplevels()
summary(genome_ploidy_first_1_level_data)  

##prepare the phy

keep.ploid_first1_inc_spp <-levels(genome_ploidy_first_1_level_data$species)

remove_ploid_first1_inc_taxa = setdiff(genome.phy$tip.label, keep.ploid_first1_inc_spp)
genomeploid_first1_inc_phy <- drop.tip(genome.phy,remove_ploid_first1_inc_taxa)

# Order data by tip order

genome_ploidy_first_1_level_data = genome_ploidy_first_1_level_data[(genome_ploidy_first_1_level_data$species %in% genomeploid_first1_inc_phy$tip.label), ]

row.names(genome_ploidy_first_1_level_data) = genome_ploidy_first_1_level_data$species
###check the names
library(geiger)  #for name.check
name.check(genomeploid_first1_inc_phy, genome_ploidy_first_1_level_data)

str(genome_ploidy_first_1_level_data)
summary(genome_ploidy_first_1_level_data)

genome_ploidy_first_1_level_data$no_ploidy_lev <- as.factor(genome_ploidy_first_1_level_data$no_ploidy_lev)

genome.ploidy.gs.first1.model <- phyloglm(GloNAF_incidence ~ range_1.1(log10(monoploid)), 
                                          phy = genomeploid_first1_inc_phy, data = genome_ploidy_first_1_level_data, 
                                          method = c("logistic_MPLE","logistic_IG10","poisson_GEE"),
                                          btol = 100, log.alpha.bound = 4,
                                          start.beta=NULL, start.alpha=NULL,
                                          boot = 0, full.matrix = TRUE)

summary(genome.ploidy.gs.first1.model)

### prepare ployploids subset
genome_ploidy_first2_level_data <-  genome_ploidy_first_2_level_data  %>%  filter (no_ploidy_lev == 2) %>% droplevels()
summary(genome_ploidy_first2_level_data)  

##prepare the phy

keep.ploid_first2_inc_spp <-levels(genome_ploidy_first2_level_data$species)

remove_ploid_first2_inc_taxa = setdiff(genome.phy$tip.label, keep.ploid_first2_inc_spp)
genomeploid_first2_inc_phy <- drop.tip(genome.phy,remove_ploid_first2_inc_taxa)

# Order data by tip order

genome_ploidy_first2_level_data = genome_ploidy_first2_level_data[(genome_ploidy_first2_level_data$species %in% genomeploid_first2_inc_phy$tip.label), ]

row.names(genome_ploidy_first2_level_data) = genome_ploidy_first2_level_data$species
###check the names
library(geiger)  #for name.check
name.check(genomeploid_first2_inc_phy, genome_ploidy_first2_level_data)

str(genome_ploidy_first2_level_data)
summary(genome_ploidy_first2_level_data)

genome.ploidy.gs.first_2.model <- phyloglm(GloNAF_incidence ~ range_1.1(log10(monoploid)) , 
                                           phy = genomeploid_first2_inc_phy, data = genome_ploidy_first2_level_data, 
                                           method = c("logistic_MPLE","logistic_IG10","poisson_GEE"),
                                           btol = 100, log.alpha.bound = 4,
                                           start.beta=NULL, start.alpha=NULL,
                                           boot = 0, full.matrix = TRUE)

summary(genome.ploidy.gs.first_2.model)








#### naturalization extent global

summary(genome_ploidy_ext_data)
genome_ploidy_ext_monoploid_final <- genome_ploidy_ext_data %>% filter (monoploid != "NA") %>% droplevels() #2364 species

summary(as.factor(genome_ploidy_ext_monoploid_final$no_ploidy_lev))

summary(genome_ploidy_ext_monoploid_final)

genome.ploidy.gs.extent.model <- phylolm(range_1.1(log10(reg_nat_total)) ~ range_1.1(log10(no_ploidy_lev))*range_1.1(log10(monoploid)),
                                         phy = genome.ploidy.ext.phy, data = genome_ploidy_ext_monoploid_final, 
                                         method = "lambda")

summary(genome.ploidy.gs.extent.model)
#R2(genome.monoploid.main.extent.model, phy = genome.extent.main.phy)
plot(genome.ploidy.gs.extent.model)
## model Residual Plot
qplot(range_1.1(log10(genome_ploidy_ext_data$reg_nat_total)), residuals(genome.ploidy.gs.extent.model))
hist(residuals(genome.ploidy.gs.extent.model))
plot(density(resid(genome.ploidy.gs.extent.model))) #A density plot
# A quantile normal plot - good for checking normality
qqnorm(resid(genome.ploidy.gs.extent.model)) 
qqline(resid(genome.ploidy.gs.extent.model))

#### prepare the results for the interaction


### prepare ployploids subset
genome_ploidy_first_2_level_ext_data <-  genome_ploidy_ext_monoploid_final  %>%  filter (no_ploidy_lev < 3) %>% droplevels()
summary(genome_ploidy_first_2_level_ext_data)  

##prepare the phy

keep.ploid_first2_ext_spp <-levels(genome_ploidy_first_2_level_ext_data$species)

remove_ploid_first2_ext_taxa = setdiff(genome.ploidy.ext.phy$tip.label, keep.ploid_first2_ext_spp)
genomeploid_first2_ext_phy <- drop.tip(genome.ploidy.ext.phy,remove_ploid_first2_ext_taxa)

# Order data by tip order

genome_ploidy_first_2_level_ext_data = genome_ploidy_first_2_level_ext_data[(genome_ploidy_first_2_level_ext_data$species %in% genomeploid_first2_ext_phy$tip.label), ]

row.names(genome_ploidy_first_2_level_ext_data) = genome_ploidy_first_2_level_ext_data$species
###check the names
library(geiger)  #for name.check
name.check(genomeploid_first2_ext_phy, genome_ploidy_first_2_level_ext_data)

str(genome_ploidy_first_2_level_ext_data)
summary(genome_ploidy_first_2_level_ext_data)


genome_ploidy_first_2_level_ext_data$no_ploidy_lev <- as.factor(genome_ploidy_first_2_level_ext_data$no_ploidy_lev)

genome.ploidy.gs.first2_ext.model  <- phylolm(range_1.1(log10(reg_nat_total)) ~ no_ploidy_lev*range_1.1(log10(monoploid))-1,
                                              phy = genomeploid_first2_ext_phy, data = genome_ploidy_first_2_level_ext_data, 
                                              method = "lambda")

summary(genome.ploidy.gs.first2_ext.model)


### prepare each level 
### prepare ployploids subset
genome_ploidy_first_1_level_ext_data <-  genome_ploidy_first_2_level_ext_data  %>%  filter (no_ploidy_lev == 1) %>% droplevels()
summary(genome_ploidy_first_1_level_ext_data)  

##prepare the phy

keep.ploid_first1_inc_spp <-levels(genome_ploidy_first_1_level_ext_data$species)

remove_ploid_first1_inc_taxa = setdiff(genome.ploidy.ext.phy$tip.label, keep.ploid_first1_inc_spp)
genomeploid_first1_ext_phy <- drop.tip(genome.ploidy.ext.phy,remove_ploid_first1_inc_taxa)

# Order data by tip order

genome_ploidy_first_1_level_ext_data = genome_ploidy_first_1_level_ext_data[(genome_ploidy_first_1_level_ext_data$species %in% genomeploid_first1_ext_phy$tip.label), ]

row.names(genome_ploidy_first_1_level_ext_data) = genome_ploidy_first_1_level_ext_data$species
###check the names
library(geiger)  #for name.check
name.check(genomeploid_first1_ext_phy, genome_ploidy_first_1_level_ext_data)

str(genome_ploidy_first_1_level_ext_data)
summary(genome_ploidy_first_1_level_ext_data)

genome.ploidy.gs.first_1_ext.model <- phylolm(range_1.1(log10(reg_nat_total)) ~ range_1.1(log10(monoploid)),
                                              phy = genomeploid_first1_ext_phy, data = genome_ploidy_first_1_level_ext_data, 
                                              method = "lambda")

summary(genome.ploidy.gs.first_1_ext.model)


### prepare ployploids subset:level 2
genome_ploidy_first_1_2_level_ext_data <-  genome_ploidy_first_2_level_ext_data  %>%  filter (no_ploidy_lev == 2) %>% droplevels()
summary(genome_ploidy_first_1_2_level_ext_data)  

##prepare the phy

keep.ploid_first1_2_inc_spp <-levels(genome_ploidy_first_1_2_level_ext_data$species)

remove_ploid_first1_2_inc_taxa = setdiff(genome.ploidy.ext.phy$tip.label, keep.ploid_first1_2_inc_spp)
genomeploid_first1_2_ext_phy <- drop.tip(genome.ploidy.ext.phy,remove_ploid_first1_2_inc_taxa)

# Order data by tip order

genome_ploidy_first_1_2_level_ext_data = genome_ploidy_first_1_2_level_ext_data[(genome_ploidy_first_1_2_level_ext_data$species %in% genomeploid_first1_2_ext_phy$tip.label), ]

row.names(genome_ploidy_first_1_2_level_ext_data) = genome_ploidy_first_1_2_level_ext_data$species
###check the names
library(geiger)  #for name.check
name.check(genomeploid_first1_2_ext_phy, genome_ploidy_first_1_2_level_ext_data)

str(genome_ploidy_first_1_2_level_ext_data)
summary(genome_ploidy_first_1_2_level_ext_data)

genome.ploidy.gs.first_1_2_ext.model <- phylolm(range_1.1(log10(reg_nat_total)) ~ range_1.1(log10(monoploid)),
                                                phy = genomeploid_first1_2_ext_phy, data = genome_ploidy_first_1_2_level_ext_data, 
                                                method = "lambda")

summary(genome.ploidy.gs.first_1_2_ext.model)



### naturalization extent mainland

summary(genome_ploidy_mainext_data)
summary(as.factor(genome_ploidy_mainext_data$no_ploidy_lev))
#### mainland holoploid

genome.ploidy.gs.mainextent.model <- phylolm(range_1.1(log10(reg_nat_main)) ~ no_ploidy_lev*range_1.1(log10(monoploid))-1,
                                             phy = genome.ploidy.mainext.phy, data = genome_ploidy_mainext_data, 
                                             method = "lambda")

summary(genome.ploidy.gs.mainextent.model)
#R2(genome.monoploid.main.mainextent.model, phy = genome.mainextent.main.phy)
plot(genome.ploidy.gs.mainextent.model)
## model Residual Plot
qplot(range_1.1(log10(genome_ploidy_mainext_data$reg_nat_main)), residuals(genome.ploidy.gs.mainextent.model))
hist(residuals(genome.ploidy.gs.mainextent.model))
plot(density(resid(genome.ploidy.gs.mainextent.model))) #A density plot
# A quantile normal plot - good for checking normality
qqnorm(resid(genome.ploidy.gs.mainextent.model)) 
qqline(resid(genome.ploidy.gs.mainextent.model))

##### check the interaction  (Mainland naturalization extent)


### prepare ployploids subset
genome_ploidy_first_2_level_main_ext_data <-  genome_ploidy_mainext_data  %>%  filter (no_ploidy_lev < 3) %>% droplevels()
summary(genome_ploidy_first_2_level_main_ext_data)  

##prepare the phy

keep.ploid_first2_main_ext_spp <-levels(genome_ploidy_first_2_level_main_ext_data$species)

remove_ploid_first2_main_ext_taxa = setdiff(genome.ploidy.ext.phy$tip.label, keep.ploid_first2_main_ext_spp)
genomeploid_first2_main_ext_phy <- drop.tip(genome.ploidy.ext.phy,remove_ploid_first2_main_ext_taxa)

# Order data by tip order

genome_ploidy_first_2_level_main_ext_data = genome_ploidy_first_2_level_main_ext_data[(genome_ploidy_first_2_level_main_ext_data$species %in% genomeploid_first2_main_ext_phy$tip.label), ]

row.names(genome_ploidy_first_2_level_main_ext_data) = genome_ploidy_first_2_level_main_ext_data$species
###check the names
library(geiger)  #for name.check
name.check(genomeploid_first2_main_ext_phy, genome_ploidy_first_2_level_main_ext_data)

str(genome_ploidy_first_2_level_main_ext_data)
summary(genome_ploidy_first_2_level_main_ext_data)


genome_ploidy_first_2_level_main_ext_data$no_ploidy_lev <- as.factor(genome_ploidy_first_2_level_main_ext_data$no_ploidy_lev)

genome.ploidy.gs.first2_main_ext.model  <- phylolm(range_1.1(log10(reg_nat_main)) ~ no_ploidy_lev*range_1.1(log10(monoploid))-1,
                                                   phy = genomeploid_first2_main_ext_phy, data = genome_ploidy_first_2_level_main_ext_data, 
                                                   method = "lambda")

summary(genome.ploidy.gs.first2_main_ext.model)


### naturalization extent island
summary(genome_ploidy_islext_data)
summary(as.factor(genome_ploidy_islext_data$no_ploidy_lev))


genome.ploidy.gs.islextent.model <- phylolm(range_1.1(log10(reg_nat_island)) ~ no_ploidy_lev*range_1.1(log10(monoploid)),
                                            phy = genome.ploidy.islext.phy, data = genome_ploidy_islext_data, 
                                            method = "lambda")

summary(genome.ploidy.gs.islextent.model)
#R2(genome.monoploid.main.islextent.model, phy = genome.islextent.main.phy)
plot(genome.ploidy.gs.islextent.model)
## model Residual Plot
qplot(range_1.1(log10(genome_ploidy_islext_data$reg_nat_island)), residuals(genome.ploidy.gs.islextent.model))
hist(residuals(genome.ploidy.gs.islextent.model))
plot(density(resid(genome.ploidy.gs.islextent.model))) #A density plot
# A quantile normal plot - good for checking normality
qqnorm(resid(genome.ploidy.gs.islextent.model)) 
qqline(resid(genome.ploidy.gs.islextent.model))

###############################
genome_ploidy_first_2_level_island_ext_data <-  genome_ploidy_islext_data  %>%  filter (no_ploidy_lev < 3) %>% droplevels()
summary(genome_ploidy_first_2_level_island_ext_data)  

##prepare the phy

keep.ploid_first2_island_ext_spp <-levels(genome_ploidy_first_2_level_island_ext_data$species)

remove_ploid_first2_island_ext_taxa = setdiff(genome.ploidy.ext.phy$tip.label, keep.ploid_first2_island_ext_spp)
genomeploid_first2_island_ext_phy <- drop.tip(genome.ploidy.ext.phy,remove_ploid_first2_island_ext_taxa)

# Order data by tip order

genome_ploidy_first_2_level_island_ext_data = genome_ploidy_first_2_level_island_ext_data[(genome_ploidy_first_2_level_island_ext_data$species %in% genomeploid_first2_island_ext_phy$tip.label), ]

row.names(genome_ploidy_first_2_level_island_ext_data) = genome_ploidy_first_2_level_island_ext_data$species
###check the names
library(geiger)  #for name.check
name.check(genomeploid_first2_island_ext_phy, genome_ploidy_first_2_level_island_ext_data)

str(genome_ploidy_first_2_level_island_ext_data)
summary(genome_ploidy_first_2_level_island_ext_data)


genome_ploidy_first_2_level_island_ext_data$no_ploidy_lev <- as.factor(genome_ploidy_first_2_level_island_ext_data$no_ploidy_lev)

genome.ploidy.gs.first2_island_ext.model  <- phylolm(range_1.1(log10(reg_nat_island)) ~ no_ploidy_lev*range_1.1(log10(monoploid))-1,
                                                     phy = genomeploid_first2_island_ext_phy, data = genome_ploidy_first_2_level_island_ext_data, 
                                                     method = "lambda")

summary(genome.ploidy.gs.first2_island_ext.model)


### invasive extent
summary(genome_ploidy_invasive_data)
genome_ploidy_invasive_monoploid_final <- genome_ploidy_invasive_data %>% filter (monoploid != "NA") %>% droplevels() #2364 species


summary(as.factor(genome_ploidy_invasive_monoploid_final$no_ploidy_lev))


genome.ploidy.gs.invasiveent.model <- phylolm(range_1.1(log10(reg_no_invasive)) ~ no_ploidy_lev*range_1.1(log10(monoploid)),
                                              phy = genome.ploidy.invasive.phy, data = genome_ploidy_invasive_monoploid_final, 
                                              method = "lambda")

summary(genome.ploidy.gs.invasiveent.model)
#R2(genome.monoploid.main.invasiveent.model, phy = genome.invasiveent.main.phy)
plot(genome.ploidy.gs.invasiveent.model)
## model Residual Plot
qplot(range_1.1(log10(genome_ploidy_invasive_data$reg_no_invasive)), residuals(genome.ploidy.gs.invasiveent.model))
hist(residuals(genome.ploidy.gs.invasiveent.model))
plot(density(resid(genome.ploidy.gs.invasiveent.model))) #A density plot
# A quantile normal plot - good for checking normality
qqnorm(resid(genome.ploidy.gs.invasiveent.model)) 
qqline(resid(genome.ploidy.gs.invasiveent.model))

##################################################
## check the interaction between monoploid and no of ploidy level
### prepare ployploids subset
genome_ploidy_first_2_level_inv_data <-  genome_ploidy_invasive_monoploid_final  %>%  filter (no_ploidy_lev < 3) %>% droplevels()
summary(genome_ploidy_first_2_level_inv_data)  

##prepare the phy

keep.ploid_first2_inv_spp <-levels(genome_ploidy_first_2_level_inv_data$species)

remove_ploid_first2_inv_taxa = setdiff(genome.ploidy.invasive.phy$tip.label, keep.ploid_first2_inv_spp)
genomeploid_first2_inv_phy <- drop.tip(genome.ploidy.invasive.phy,remove_ploid_first2_inv_taxa)

# Order data by tip order

genome_ploidy_first_2_level_inv_data = genome_ploidy_first_2_level_inv_data[(genome_ploidy_first_2_level_inv_data$species %in% genomeploid_first2_inv_phy$tip.label), ]

row.names(genome_ploidy_first_2_level_inv_data) = genome_ploidy_first_2_level_inv_data$species
###check the names
library(geiger)  #for name.check
name.check(genomeploid_first2_inv_phy, genome_ploidy_first_2_level_inv_data)

str(genome_ploidy_first_2_level_inv_data)
summary(genome_ploidy_first_2_level_inv_data)


genome_ploidy_first_2_level_inv_data$no_ploidy_lev <- as.factor(genome_ploidy_first_2_level_inv_data$no_ploidy_lev)

genome.ploidy.gs.first2_inv.model  <- phylolm(range_1.1(log10(reg_no_invasive)) ~ no_ploidy_lev*range_1.1(log10(monoploid))-1,
                                              phy = genomeploid_first2_inv_phy, data = genome_ploidy_first_2_level_inv_data, 
                                              method = "lambda")

summary(genome.ploidy.gs.first2_inv.model)


### prepare each level to check the interactions
### prepare ployploids subset
genome_ploidy_first_1_level_inv_data <-  genome_ploidy_first_2_level_inv_data  %>%  filter (no_ploidy_lev == 1) %>% droplevels()
summary(genome_ploidy_first_1_level_inv_data)  

##prepare the phy

keep.ploid_first_1_inc_spp <-levels(genome_ploidy_first_1_level_inv_data$species)

remove_ploid_first_1_inc_taxa = setdiff(genome.ploidy.invasive.phy$tip.label, keep.ploid_first_1_inc_spp)
genomeploid_first_1_ext_phy <- drop.tip(genome.ploidy.invasive.phy,remove_ploid_first_1_inc_taxa)

# Order data by tip order

genome_ploidy_first_1_level_inv_data = genome_ploidy_first_1_level_inv_data[(genome_ploidy_first_1_level_inv_data$species %in% genomeploid_first_1_ext_phy$tip.label), ]

row.names(genome_ploidy_first_1_level_inv_data) = genome_ploidy_first_1_level_inv_data$species
###check the names
library(geiger)  #for name.check
name.check(genomeploid_first_1_ext_phy, genome_ploidy_first_1_level_inv_data)

str(genome_ploidy_first_1_level_inv_data)
summary(genome_ploidy_first_1_level_inv_data)

genome.ploidy.gs.first_1_inv.model <- phylolm(range_1.1(log10(reg_no_invasive)) ~ range_1.1(log10(monoploid)),
                                              phy = genomeploid_first_1_ext_phy, data = genome_ploidy_first_1_level_inv_data, 
                                              method = "lambda")

summary(genome.ploidy.gs.first_1_inv.model)


### prepare ployploids subset:level 2
genome_ploidy_first_1_2_level_inv_data <-  genome_ploidy_first_2_level_inv_data  %>%  filter (no_ploidy_lev == 2) %>% droplevels()
summary(genome_ploidy_first_1_2_level_inv_data)  

##prepare the phy

keep.ploid_first1_2_inv_spp <-levels(genome_ploidy_first_1_2_level_inv_data$species)

remove_ploid_first1_2_inv_taxa = setdiff(genome.ploidy.invasive.phy$tip.label, keep.ploid_first1_2_inv_spp)
genomeploid_first1_2_inv_phy <- drop.tip(genome.ploidy.invasive.phy,remove_ploid_first1_2_inv_taxa)

# Order data by tip order

genome_ploidy_first_1_2_level_inv_data = genome_ploidy_first_1_2_level_inv_data[(genome_ploidy_first_1_2_level_inv_data$species %in% genomeploid_first1_2_inv_phy$tip.label), ]

row.names(genome_ploidy_first_1_2_level_inv_data) = genome_ploidy_first_1_2_level_inv_data$species
###check the names
library(geiger)  #for name.check
name.check(genomeploid_first1_2_inv_phy, genome_ploidy_first_1_2_level_inv_data)

str(genome_ploidy_first_1_2_level_inv_data)
summary(genome_ploidy_first_1_2_level_inv_data)

genome.ploidy.gs.first_1_2_inv.model <- phylolm(range_1.1(log10(reg_no_invasive)) ~ range_1.1(log10(monoploid)),
                                                phy = genomeploid_first1_2_inv_phy, data = genome_ploidy_first_1_2_level_inv_data, 
                                                method = "lambda")

summary(genome.ploidy.gs.first_1_2_inv.model)
