#######################
### Pond Microbiome ###
#######################

# rblauste@umd.edu #
# last edited: 01/09/2023 #

###
### load packages
###

library(SCRuB)
library(ggplot2)
library(vegan)
library(reshape2)
library(gridExtra)


###
### load files
### SCRuB read counts from negative control DNA extractions and library preps
###

## metadata for all samples
meta = read.csv("metadata/sample-metadata-full.csv", h=T, row.names=1)
meta = meta[order(meta$sample.id),]
head(meta)

## metadata for scrub (contaminant removal)
meta_scrub = read.csv("metadata/scrub_list.csv", h=T, row.names=1)
meta_scrub = meta_scrub[order(rownames(meta_scrub)),]
rownames(meta_scrub) == meta$sample.id

## phylum (level 2)
taxa_phy = read.csv("taxa_tables/level-2_edit.csv", h=T, row.names=1)
colnames(taxa_phy) == meta$sample.id
taxa_phy
scr_out <- SCRuB(t(taxa_phy), meta_scrub,
                 control_order = c("control_DNA_extraction", "control_sed_processing_DNA_extract", "control_library_prep"))
taxa_phy_s = as.data.frame(t(scr_out$decontaminated_samples))

## class (level 3)
taxa_class = read.csv("taxa_tables/level-3_edit.csv", h=T, row.names=1)
colnames(taxa_class) == meta$sample.id
taxa_class
scr_out <- SCRuB(t(taxa_class), meta_scrub,
                 control_order = c("control_DNA_extraction", "control_sed_processing_DNA_extract", "control_library_prep"))
taxa_class_s = as.data.frame(t(scr_out$decontaminated_samples))

## order (level 4)
taxa_order = read.csv("taxa_tables/level-4_edit.csv", h=T, row.names=1)
colnames(taxa_order) == meta$sample.id
taxa_order
scr_out <- SCRuB(t(taxa_order), meta_scrub,
                 control_order = c("control_DNA_extraction", "control_sed_processing_DNA_extract", "control_library_prep"))
taxa_order_s = as.data.frame(t(scr_out$decontaminated_samples))

## family (level 5)
taxa_fam = read.csv("taxa_tables/level-5_edit.csv", h=T, row.names=1)
colnames(taxa_fam) == meta$sample.id
taxa_fam
scr_out <- SCRuB(t(taxa_fam), meta_scrub,
                 control_order = c("control_DNA_extraction", "control_sed_processing_DNA_extract", "control_library_prep"))
taxa_fam_s = as.data.frame(t(scr_out$decontaminated_samples))

## genus (level 6)
taxa_gen = read.csv("taxa_tables/level-6_edit.csv", h=T, row.names=1)
colnames(taxa_gen) == meta$sample.id
taxa_gen
scr_out <- SCRuB(t(taxa_gen), meta_scrub,
                 control_order = c("control_DNA_extraction", "control_sed_processing_DNA_extract", "control_library_prep"))
taxa_gen_s = as.data.frame(t(scr_out$decontaminated_samples))

## reset metadata file
meta = meta[which(meta_scrub$is_control == FALSE),]


###
### check read counts and clean tables
###

## count check per table
apply(taxa_phy, 2, sum) == apply(taxa_gen, 2, sum)

## counts stats (remove controls)
counts = apply(taxa_gen, 2, sum)
counts = counts[-c(grep("NC|PCR|Sed|Zymo", names(counts)))]

## summary stats

# all
mean(counts)
sd(counts)
sd(counts)/sqrt(length(counts))

# water (mean, sd, se)
mean(counts[grep("W", names(counts))])
sd(counts[grep("W", names(counts))])
sd(counts[grep("W", names(counts))])/sqrt(length(counts[grep("W", names(counts))]))

# sediment (mean, sd, se)
mean(counts[grep("S", names(counts))])
sd(counts[grep("S", names(counts))])
sd(counts[grep("S", names(counts))])/sqrt(length(counts[grep("S", names(counts))]))

## remove singletons
taxa_phy_s = taxa_phy_s[which(apply(taxa_phy_s, 1, sum) > 1),]
taxa_class_s = taxa_class_s[which(apply(taxa_class_s, 1, sum) > 1),]
taxa_order_s = taxa_order_s[which(apply(taxa_order_s, 1, sum) > 1),]
taxa_fam_s = taxa_fam_s[which(apply(taxa_fam_s, 1, sum) > 1),]
taxa_gen_s = taxa_gen_s[which(apply(taxa_gen_s, 1, sum) > 1),]


######################
### NORMALIZE DATA ###
######################

###
### Remove chloroplast/mitochondria counts, remove singletons, and convert to relative abundances
###

## call chloroplast (order level) + mitochondria (family level) counts
chloroplast_counts = apply(taxa_order_s[grep("Chloroplast", rownames(taxa_order_s)),], 2, sum)
mitochondria_counts = apply(taxa_fam_s[grep("Mitochondria", rownames(taxa_fam_s)),], 2, sum)

## phylum
taxa_clean = taxa_phy_s
taxa_clean[grep("Cyanob", rownames(taxa_clean)),] = taxa_clean[grep("Cyanob", rownames(taxa_clean)),] - chloroplast_counts
taxa_clean[grep("Proteo", rownames(taxa_clean)),] = taxa_clean[grep("Proteo", rownames(taxa_clean)),] - mitochondria_counts
taxa_clean = taxa_clean[which(apply(taxa_clean,1,sum)>1),]

taxa_clean = as.matrix(taxa_clean)
taxa_clean[which(taxa_clean<0)] = 0
taxa_clean = as.data.frame(taxa_clean)

taxa_clean_RA = taxa_clean
for (i in 1:dim(taxa_clean_RA)[2]) {
  taxa_clean_RA[,i] = taxa_clean[,i]/sum(taxa_clean[,i])
}
# check sums
apply(taxa_clean_RA, 2, sum)
# save table
taxa_phy_RA = taxa_clean_RA

## class
taxa_clean = taxa_class_s
taxa_clean[grep("Oxyphotobacteria", rownames(taxa_clean)),] = taxa_clean[grep("Oxyphotobacteria", rownames(taxa_clean)),] - chloroplast_counts
taxa_clean[grep("Alphaproteo", rownames(taxa_clean)),] = taxa_clean[grep("Alphaproteo", rownames(taxa_clean)),] - mitochondria_counts
taxa_clean = taxa_clean[which(apply(taxa_clean,1,sum)>1),]

taxa_clean = as.matrix(taxa_clean)
taxa_clean[which(taxa_clean<0)] = 0
taxa_clean = as.data.frame(taxa_clean)

taxa_clean_RA = taxa_clean
for (i in 1:dim(taxa_clean_RA)[2]) {
  taxa_clean_RA[,i] = taxa_clean[,i]/sum(taxa_clean[,i])
}
# check sums
apply(taxa_clean_RA, 2, sum)
# save table
taxa_class_RA = taxa_clean_RA

## order
taxa_clean = taxa_order_s
taxa_clean = taxa_clean[-c(grep("Chloroplast", rownames(taxa_clean))),]
taxa_clean[grep("Rickettsiales", rownames(taxa_clean)),] = taxa_clean[grep("Rickettsiales", rownames(taxa_clean)),] - mitochondria_counts
taxa_clean = taxa_clean[which(apply(taxa_clean,1,sum)>1),]
taxa_clean_RA = taxa_clean
for (i in 1:dim(taxa_clean_RA)[2]) {
  taxa_clean_RA[,i] = taxa_clean[,i]/sum(taxa_clean[,i])
}
# check sums
apply(taxa_clean_RA, 2, sum)
# save table
taxa_order_RA = taxa_clean_RA

## family
taxa_clean = taxa_fam_s
taxa_clean = taxa_clean[-c(grep("Chloroplast", rownames(taxa_clean))),] 
taxa_clean = taxa_clean[-c(grep("Mitochondria", rownames(taxa_clean))),] 
taxa_clean = taxa_clean[which(apply(taxa_clean,1,sum)>1),]
taxa_clean_RA = taxa_clean
for (i in 1:dim(taxa_clean_RA)[2]) {
  taxa_clean_RA[,i] = taxa_clean[,i]/sum(taxa_clean[,i])
}
# check sums
apply(taxa_clean_RA, 2, sum)
# save table
taxa_fam_RA = taxa_clean_RA

## genus
taxa_clean = taxa_gen_s
taxa_clean = taxa_clean[-c(grep("Chloroplast", rownames(taxa_clean))),] 
taxa_clean = taxa_clean[-c(grep("Mitochondria", rownames(taxa_clean))),] 
taxa_clean = taxa_clean[which(apply(taxa_clean,1,sum)>1),]
taxa_clean_RA = taxa_clean
for (i in 1:dim(taxa_clean_RA)[2]) {
  taxa_clean_RA[,i] = taxa_clean[,i]/sum(taxa_clean[,i])
}
# check sums
apply(taxa_clean_RA, 2, sum)
# save table
taxa_gen_RA = taxa_clean_RA

## save final processed RA tables
write.table(taxa_phy_RA, "taxa_tables/taxa_phy_RA.txt", sep = "\t", quote = FALSE)
write.table(taxa_class_RA, "taxa_tables/taxa_class_RA.txt", sep = "\t", quote = FALSE)
write.table(taxa_order_RA, "taxa_tables/taxa_order_RA.txt", sep = "\t", quote = FALSE)
write.table(taxa_fam_RA, "taxa_tables/taxa_fam_RA.txt", sep = "\t", quote = FALSE)
write.table(taxa_gen_RA, "taxa_tables/taxa_gen_RA.txt", sep = "\t", quote = FALSE)


#######################
### ALPHA DIVERSITY ###
#######################

## set input
taxa_table = taxa_gen_RA

## compute shannon and simpson diversity metrics
diversity_vec = matrix(nrow = dim(taxa_table)[2], ncol = 2)
diversity_vec = as.data.frame(diversity_vec)
for (a in 1:dim(taxa_table)[2]) {
  diversity_vec[a,1] = diversity(taxa_table[,a], index = "shannon")
  diversity_vec[a,2] = diversity(taxa_table[,a], index = "simpson")
}
colnames(diversity_vec) = c("Shannon", "Simpson")

## add alpha metric to metadata
meta$Shannon = diversity_vec$Shannon
meta$Simpson = diversity_vec$Simpson

## two-way ANOVA for effect of sample depth and time on alpha diversity
summary(aov(Shannon ~ sample.location.depth*time.of.day, meta[grep("Wye_water", meta$sample.type.pond),]))
TukeyHSD(aov(Shannon ~ sample.location.depth*time.of.day, meta[grep("Wye_water", meta$sample.type.pond),]))

## boxplot: Wye water by sample depth and time
# set labels
depth_names = list(
  "Bank_D1"="Bank: surface",
  "Interior_D1"="Interior: surface",
  "Interior_D2"="Interior: 1-m depth",
  "Interior_D3"="Interior: 2-m depth"
)
depth_labeller <- function(variable,value){
  return(depth_names[value])
}
# plot (FIGURE 3-A)
ggplot(meta[grep("Wye_water", meta$sample.type.pond),], 
       aes(x = time.of.day, y = Shannon, fill=time.of.day)) +
  geom_boxplot(outlier.shape = NA, size=0.8, alpha = 0.6, fill="lightgray") +
  geom_jitter(size = 2, width = 0.2, alpha = 0.6, pch=21, fill = "lightgray") +
  theme_bw() +
  ylab("Shannon Index") +
  ylim(2.4,4.3) +
  #scale_fill_manual(values=c("violet", "orange", "skyblue")) +
  facet_grid(~sample.location.depth, scales="free", space="free", 
             labeller = depth_labeller) +
  scale_x_discrete(labels = c("9:00", "12:00", "15:00")) +
  theme(legend.text = element_text(size = 11),
        legend.title = element_blank(),
        #legend.position = "none",
        axis.title.x = element_blank(),
        strip.text = element_text(size = 10),
        axis.text.x = element_text(size = 12, angle=45, hjust=1),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 12))


## ANOVA for effect of time on alpha diversity
summary(aov(Shannon ~ time.of.day, meta[grep("Wye_water", meta$sample.type.pond),]))
TukeyHSD(aov(Shannon ~ time.of.day, meta[grep("B", meta$sample.location.depth),]))
TukeyHSD(aov(Shannon ~ time.of.day, meta[grep("Wye_water", meta$sample.type.pond),]))
TukeyHSD(aov(Shannon ~ time.of.day, meta[grep("Interior_D1", meta$sample.location.depth),]))
TukeyHSD(aov(Shannon ~ time.of.day, meta[grep("D2", meta$sample.location.depth),]))
TukeyHSD(aov(Shannon ~ time.of.day, meta[grep("D3", meta$sample.location.depth),]))

## ANOVA for effect of sample depth on alpha diversity
summary(aov(Shannon ~ sample.location.depth, meta[grep("Wye_water", meta$sample.type.pond),]))
TukeyHSD(aov(Shannon ~ sample.location.depth, meta[grep("Wye_water", meta$sample.type.pond),]))

## Bank vs. Interior (D1)
summary(aov(Shannon ~ sample.location.depth, meta[grep("D1", meta$sample.location.depth),]))
summary(aov(Shannon ~ sample.location.depth, meta[intersect(grep("T1", meta$time.of.day),
                                                            grep("D1", meta$sample.location.depth)),]))
summary(aov(Shannon ~ sample.location.depth, meta[intersect(grep("T2", meta$time.of.day),
                                                            grep("D1", meta$sample.location.depth)),]))
summary(aov(Shannon ~ sample.location.depth, meta[intersect(grep("T3", meta$time.of.day),
                                                            grep("D1", meta$sample.location.depth)),]))

## Depth (Interior)
summary(aov(Shannon ~ sample.location.depth, meta[grep("Interior", meta$sample.location.depth),]))
summary(aov(Shannon ~ sample.location.depth, meta[intersect(grep("T1", meta$time.of.day),
                                                            grep("Interior", meta$sample.location.depth)),]))
summary(aov(Shannon ~ sample.location.depth, meta[intersect(grep("T2", meta$time.of.day),
                                                            grep("Interior", meta$sample.location.depth)),]))
summary(aov(Shannon ~ sample.location.depth, meta[intersect(grep("T3", meta$time.of.day),
                                                            grep("Interior", meta$sample.location.depth)),]))

## boxplot: Wye sediment by sample location (FIGURE 3-B)
ggplot(meta[grep("Wye_sed", meta$sample.type.pond),], 
       aes(x = sample.water.location, y = Shannon)) +
  geom_boxplot(outlier.shape = NA, size=0.8, alpha = 0.6, fill="lightgray") +
  geom_jitter(size = 2, width = 0.1, alpha = 0.6, pch=21, fill = "lightgray") +
  theme_bw() +
  ylab("Shannon Index") +
  ylim(2.4,4.3) +
  #scale_x_discrete(labels = c("9:00", "12:00", "15:00")) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 14, angle=45, hjust=1),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 14))

## ANOVA and t-test for effect of sediment location on alpha diversity
summary(aov(Shannon ~ sample.water.location, meta[grep("Wye_sed", meta$sample.type.pond),]))
test = meta[grep("Wye_sed", meta$sample.type.pond),]
t.test(test[grep("Bank", test$sample.water.location),]$Shannon,
       test[grep("Interior", test$sample.water.location),]$Shannon)


######################
### BETA DIVERSITY ###
######################

###
### PCoA of Wye water samples
###

# beta-diversity measure
beta <- vegdist(t(taxa_table[,grep("Wye_water", meta$sample.type.pond)]), 'bray', binary = T)

# projection
pcoa <- cmdscale(beta, k = 4, eig = TRUE)

# cleanup
ord <- as.data.frame(pcoa$points)
names(ord) <- c('pcoa1', 'pcoa2', 'pcoa3', 'pcoa4') ## rename coordinates

# Percent explained variation
eig <- eigenvals(pcoa)
100*head(eig/sum(eig))

# add metadata
ord$sample.depth = as.character(meta$sample.depth[grep("Wye_water", meta$sample.type.pond)])
ord$time.of.day = as.character(meta$time.of.day[grep("Wye_water", meta$sample.type.pond)])
ord$sample.location = as.character(meta$sample.location[grep("Wye_water", meta$sample.type.pond)])
ord$sample.location.depth = as.character(meta$sample.location.depth[grep("Wye_water", meta$sample.type.pond)])

# plot (Figure 3-C)
ggplot(data = ord, aes(x = pcoa1, y = pcoa2, 
                       fill=sample.location.depth,
                       shape=time.of.day)) +
  geom_hline(yintercept = 0, lty = 2, alpha=0.7) +
  geom_vline(xintercept = 0, lty = 2, alpha=0.7) +
  geom_point(alpha = 0.7, #pch=21, 
             size=4.5) +
  theme_bw() +
  xlab('PCoA1 (11.7%)') +
  ylab('PCoA2 (8.6%)') +
  scale_shape_manual(values=c(22, 24, 23),
                     name = "Time of Day",
                     labels = c("9:00", "12:00", "15:00")) +
  scale_fill_manual(values=c(c("magenta", "skyblue", "blue", "darkblue")),
                    labels = c("Bank: surface",
                               "Interior: surface",
                               "Interior: 1-m depth",
                               "Interior: 2-m depth"),
                    name="Sample Location") +
  guides(fill = guide_legend(override.aes=list(shape = c(21)))) +
  theme(axis.text = element_text(size=14, color = "black"),
        axis.title = element_text(size=14, color = "black"),
        legend.text = element_text(size=14, color = "black"),
        legend.title = element_text(size=14, color = "black")) 

# stats
adonis(beta ~ sample.location.depth, data = ord, permutations = 999)$aov.tab
adonis(beta ~ time.of.day, data = ord, permutations = 999)$aov.tab
adonis(beta ~ sample.location, data = ord, permutations = 999)$aov.tab

# Bank vs. Interior (D1)
beta <- vegdist(t(taxa_table[,grep("D1", meta$sample.location.depth)]), 'bray', binary = T)
adonis(beta ~ sample.location.depth, data = meta[grep("D1", meta$sample.location.depth),], permutations = 999)$aov.tab
beta <- vegdist(t(taxa_table[,intersect(grep("T1", meta$time.of.day),
                                        grep("D1", meta$sample.location.depth))]), 'bray', binary = T)
adonis(beta ~ sample.location.depth, data = meta[intersect(grep("T1", meta$time.of.day),
                                                           grep("D1", meta$sample.location.depth)),], permutations = 999)$aov.tab
beta <- vegdist(t(taxa_table[,intersect(grep("T2", meta$time.of.day),
                                        grep("D1", meta$sample.location.depth))]), 'bray', binary = T)
adonis(beta ~ sample.location.depth, data = meta[intersect(grep("T2", meta$time.of.day),
                                                           grep("D1", meta$sample.location.depth)),], permutations = 999)$aov.tab
beta <- vegdist(t(taxa_table[,intersect(grep("T3", meta$time.of.day),
                                        grep("D1", meta$sample.location.depth))]), 'bray', binary = T)
adonis(beta ~ sample.location.depth, data = meta[intersect(grep("T3", meta$time.of.day),
                                                           grep("D1", meta$sample.location.depth)),], permutations = 999)$aov.tab


###
### PCoA of Wye sediment samples
###

## PCoA of Wye water samples

# beta-diversity measure
beta <- vegdist(t(taxa_table[,grep("Wye_sed", meta$sample.type.pond)]), 'bray', binary = T)

# projection
pcoa <- cmdscale(beta, k = 4, eig = TRUE)

# cleanup
ord <- as.data.frame(pcoa$points)
names(ord) <- c('pcoa1', 'pcoa2', 'pcoa3', 'pcoa4') ## rename coordinates

# Percent explained variation
eig <- eigenvals(pcoa)
100*head(eig/sum(eig))

# add metadata
ord$sample.location = as.character(meta$sample.water.location[grep("Wye_sed", meta$sample.type.pond)])

# plot (Figure 3-D)
ggplot(data = ord, aes(x = pcoa1, y = pcoa2, 
                       fill=sample.location)) +
  geom_hline(yintercept = 0, lty = 2, alpha=0.7) +
  geom_vline(xintercept = 0, lty = 2, alpha=0.7) +
  geom_point(alpha = 0.7, pch=21, 
             size=5) +
  theme_bw() +
  xlab('PCoA1 (21.4%)') +
  ylab('PCoA2 (17.7%)') +
  scale_fill_manual(values=c(c("magenta", "darkblue")),
                    labels = c("Sediment: Bank", "Sediment: Interior"),
                    name="Sample Type") +
  guides(fill = guide_legend(override.aes=list(shape = c(21)))) +
  theme(axis.text = element_text(size=14, color = "black"),
        axis.title = element_text(size=14, color = "black"),
        legend.text = element_text(size=14, color = "black"),
        legend.title = element_text(size=14, color = "black")) 

# stats
adonis(beta ~ sample.location, data = ord, permutations = 999)$aov.tab


#################################
### RELATIVE ABUNDANCE: GENUS ###
#################################

###
### Genus
###

## set table
taxa_plot = taxa_gen_RA

## subset water/sediment samples
taxa_plot = taxa_plot[,grep("Wye", meta$sample.type.pond)]
meta_plot = meta[grep("Wye", meta$sample.type.pond),]

## add 'average_water' column and re-order table
taxa_plot$avg_water = apply(taxa_plot[,grep("Wye_water", meta$sample.type.pond)], 1, mean)
taxa_plot = taxa_plot[order(taxa_plot$avg_water, decreasing = T),]
taxa_plot[1:10, 70:71]

# generate averages for time and depth
taxa_plot$T1_9am = apply(taxa_plot[,grep("9am", meta_plot$time.of.day)], 1, mean)
taxa_plot$T2_12pm = apply(taxa_plot[,grep("12pm", meta_plot$time.of.day)], 1, mean)
taxa_plot$T3_3pm = apply(taxa_plot[,grep("3pm", meta_plot$time.of.day)], 1, mean)
taxa_plot$Bank_D1 = apply(taxa_plot[,grep("Bank_D1", meta_plot$sample.location.depth)], 1, mean)
taxa_plot$Interior_D1 = apply(taxa_plot[,grep("Interior_D1", meta_plot$sample.location.depth)], 1, mean)
taxa_plot$Interior_D2 = apply(taxa_plot[,grep("Interior_D2", meta_plot$sample.location.depth)], 1, mean)
taxa_plot$Interior_D3 = apply(taxa_plot[,grep("Interior_D3", meta_plot$sample.location.depth)], 1, mean)
# connect time and depth
taxa_plot$T1_Bank_D1 = apply(taxa_plot[,intersect(grep("9am", meta_plot$time.of.day), grep("Bank_D1", meta_plot$sample.location.depth))], 1, mean)
taxa_plot$T1_Interior_D1 = apply(taxa_plot[,intersect(grep("9am", meta_plot$time.of.day), grep("Interior_D1", meta_plot$sample.location.depth))], 1, mean)
taxa_plot$T1_Interior_D2 = apply(taxa_plot[,intersect(grep("9am", meta_plot$time.of.day), grep("Interior_D2", meta_plot$sample.location.depth))], 1, mean)
taxa_plot$T1_Interior_D3 = apply(taxa_plot[,intersect(grep("9am", meta_plot$time.of.day), grep("Interior_D3", meta_plot$sample.location.depth))], 1, mean)
taxa_plot$T2_Bank_D1 = apply(taxa_plot[,intersect(grep("12pm", meta_plot$time.of.day), grep("Bank_D1", meta_plot$sample.location.depth))], 1, mean)
taxa_plot$T2_Interior_D1 = apply(taxa_plot[,intersect(grep("12pm", meta_plot$time.of.day), grep("Interior_D1", meta_plot$sample.location.depth))], 1, mean)
taxa_plot$T2_Interior_D2 = apply(taxa_plot[,intersect(grep("12pm", meta_plot$time.of.day), grep("Interior_D2", meta_plot$sample.location.depth))], 1, mean)
taxa_plot$T2_Interior_D3 = apply(taxa_plot[,intersect(grep("12pm", meta_plot$time.of.day), grep("Interior_D3", meta_plot$sample.location.depth))], 1, mean)
taxa_plot$T3_Bank_D1 = apply(taxa_plot[,intersect(grep("3pm", meta_plot$time.of.day), grep("Bank_D1", meta_plot$sample.location.depth))], 1, mean)
taxa_plot$T3_Interior_D1 = apply(taxa_plot[,intersect(grep("3pm", meta_plot$time.of.day), grep("Interior_D1", meta_plot$sample.location.depth))], 1, mean)
taxa_plot$T3_Interior_D2 = apply(taxa_plot[,intersect(grep("3pm", meta_plot$time.of.day), grep("Interior_D2", meta_plot$sample.location.depth))], 1, mean)
taxa_plot$T3_Interior_D3 = apply(taxa_plot[,intersect(grep("3pm", meta_plot$time.of.day), grep("Interior_D3", meta_plot$sample.location.depth))], 1, mean)
# generate averages for sediment by sampling location
taxa_plot$Sed_bank = apply(taxa_plot[,grep("S3|S4|S6|S7|S9|S10", meta_plot$sample.id)], 1, mean)
taxa_plot$Sed_interior = apply(taxa_plot[,grep("S2|S5|S8|S11", meta_plot$sample.id)], 1, mean)

## check the 'top 10' taxa
colnames(taxa_plot)
sum(taxa_plot$avg_water[1:10]) # on average, these families cover XX% of the bacterial diversity
taxa_plot[1:10, 68:71]

## subset top taxa for plotting
taxa_plot_select = rbind(apply(taxa_plot, 2, sum),
                         apply(taxa_plot[grep("_Chloroflexi", rownames(taxa_plot)),], 2, sum),
                         apply(taxa_plot[grep("_Actinobacteria", rownames(taxa_plot)),], 2, sum),
                         apply(taxa_plot[grep("_Bacteroidetes", rownames(taxa_plot)),], 2, sum),
                         apply(taxa_plot[grep("_Cyanobacteria", rownames(taxa_plot)),], 2, sum),
                         apply(taxa_plot[grep("_Firmicutes", rownames(taxa_plot)),], 2, sum),
                         apply(taxa_plot[grep("_Proteobacteria", rownames(taxa_plot)),], 2, sum),
                         apply(taxa_plot[grep("_Verrucomicrobia", rownames(taxa_plot)),], 2, sum),
                         taxa_plot[c(1:10),])
rownames(taxa_plot_select) = c("D_0__Bacteria", 
                               "D_0__Bacteria;D_1__Chloroflexi",
                               "D_0__Bacteria;D_1__Actinobacteria",
                               "D_0__Bacteria;D_1__Bacteroidetes",
                               "D_0__Bacteria;D_1__Cyanobacteria",                               
                               "D_0__Bacteria;D_1__Firmicutes",                               
                               "D_0__Bacteria;D_1__Proteobacteria",
                               "D_0__Bacteria;D_1__Verrucomicrobia",
                               rownames(taxa_plot_select[c(9:18),]))
head(taxa_plot_select)

# clean taxa
rownames(taxa_plot_select)
taxa_plot_select[grep("^D_0__Bacteria;D_1__Proteobacteria$", rownames(taxa_plot_select)),] = taxa_plot_select[grep("^D_0__Bacteria;D_1__Proteobacteria$", rownames(taxa_plot_select)),] - apply(taxa_plot_select[grep("Gammapr|Alphapr", rownames(taxa_plot_select)),], 2, sum)
taxa_plot_select[grep("^D_0__Bacteria;D_1__Bacteroidetes$", rownames(taxa_plot_select)),] = taxa_plot_select[grep("^D_0__Bacteria;D_1__Bacteroidetes$", rownames(taxa_plot_select)),] - apply(taxa_plot_select[grep("Bacteroidia", rownames(taxa_plot_select)),], 2, sum)
taxa_plot_select[grep("^D_0__Bacteria;D_1__Cyanobacteria$", rownames(taxa_plot_select)),] = taxa_plot_select[grep("^D_0__Bacteria;D_1__Cyanobacteria$", rownames(taxa_plot_select)),] - apply(taxa_plot_select[grep("Cyanobiaceae|Microcystaceae", rownames(taxa_plot_select)),], 2, sum)
taxa_plot_select[grep("^D_0__Bacteria;D_1__Verrucomicrobia$", rownames(taxa_plot_select)),] = taxa_plot_select[grep("^D_0__Bacteria;D_1__Verrucomicrobia$", rownames(taxa_plot_select)),] - apply(taxa_plot_select[grep("Verrucomicrobiae", rownames(taxa_plot_select)),], 2, sum)
taxa_plot_select[grep("^D_0__Bacteria$", rownames(taxa_plot_select)),] = taxa_plot_select[grep("^D_0__Bacteria$", rownames(taxa_plot_select)),] - apply(taxa_plot_select[grep(";", rownames(taxa_plot_select)),], 2, sum)

# check sums
apply(taxa_plot_select,2,sum)
min(taxa_plot_select)

## prep plot
rownames(taxa_plot_select)
rownames(taxa_plot_select) = c("Other", "p__Chloroflexi", "p__Actinobacteria", "p__Bacteroidetes", "p__Cyanobacteria",
                               "p__Firmicutes",  "p__Proteobacteria", "p__Verrucomicrobia", 
                               "  Cyanobium PCC-6307", "  Chitinophagaceae sp.", "  Steroidobacter", "  Microcystis PCC-7914",
                               "  Fluviicola", "  Methylocystis", "  Methylococcaceae sp.", "  Limnohabitans",
                               "  Verrucomicrobiae sp.", "  Solitalea")

### WATER SAMPLES

# melt data
colnames(taxa_plot_select)
taxa_melt = melt(data.frame(taxa_plot_select[, 11:70],
                            rownames(taxa_plot_select)))
colnames(taxa_melt)[1] = "Taxon"

# factor
taxa_melt$Taxon = factor(taxa_melt$Taxon,
                         levels = c("Other", 
                                    "p__Chloroflexi", 
                                    "p__Firmicutes",  
                                    "p__Actinobacteria", 
                                    "p__Verrucomicrobia", "  Verrucomicrobiae sp.",
                                    "p__Cyanobacteria", "  Microcystis PCC-7914",  "  Cyanobium PCC-6307",
                                    "p__Bacteroidetes", "  Solitalea", "  Fluviicola", "  Chitinophagaceae sp.",
                                    "p__Proteobacteria",   "  Limnohabitans", "  Methylococcaceae sp.", "  Methylocystis", "  Steroidobacter"
                                    ))

# time variable fo full set
taxa_melt$time.of.day = "T1_9am"
taxa_melt$time.of.day[grep(paste(c("^", paste(meta_plot$sample.id[grep("12pm", meta_plot$time.of.day)], collapse="$|^"), "$"), collapse = ""), taxa_melt$variable)] = "T2_12pm"
taxa_melt$time.of.day[grep(paste(c("^", paste(meta_plot$sample.id[grep("3pm", meta_plot$time.of.day)], collapse="$|^"), "$"), collapse = ""), taxa_melt$variable)] = "T3_3pm"

# depth variable fo full set
taxa_melt$sample.location.depth = "Bank_D1"
taxa_melt$sample.location.depth[grep(paste(c("^", paste(meta_plot$sample.id[grep("Interior_D1", meta_plot$sample.location.depth)], collapse="$|^"), "$"), collapse = ""), taxa_melt$variable)] = "Interior_D1"
taxa_melt$sample.location.depth[grep(paste(c("^", paste(meta_plot$sample.id[grep("Interior_D2", meta_plot$sample.location.depth)], collapse="$|^"), "$"), collapse = ""), taxa_melt$variable)] = "Interior_D2"
taxa_melt$sample.location.depth[grep(paste(c("^", paste(meta_plot$sample.id[grep("Interior_D3", meta_plot$sample.location.depth)], collapse="$|^"), "$"), collapse = ""), taxa_melt$variable)] = "Interior_D3"

# sample location variable
taxa_melt$sample.location = "L1"
taxa_melt$sample.location[grep(paste(c("^", paste(meta_plot$sample.id[grep("^L2$", meta_plot$sample.location)], collapse="$|^"), "$"), collapse = ""), taxa_melt$variable)] = "L2"
taxa_melt$sample.location[grep(paste(c("^", paste(meta_plot$sample.id[grep("^L3$", meta_plot$sample.location)], collapse="$|^"), "$"), collapse = ""), taxa_melt$variable)] = "L3"
taxa_melt$sample.location[grep(paste(c("^", paste(meta_plot$sample.id[grep("^L4$", meta_plot$sample.location)], collapse="$|^"), "$"), collapse = ""), taxa_melt$variable)] = "L4"
taxa_melt$sample.location[grep(paste(c("^", paste(meta_plot$sample.id[grep("^L5$", meta_plot$sample.location)], collapse="$|^"), "$"), collapse = ""), taxa_melt$variable)] = "L5"
taxa_melt$sample.location[grep(paste(c("^", paste(meta_plot$sample.id[grep("^L6$", meta_plot$sample.location)], collapse="$|^"), "$"), collapse = ""), taxa_melt$variable)] = "L6"
taxa_melt$sample.location[grep(paste(c("^", paste(meta_plot$sample.id[grep("^L7$", meta_plot$sample.location)], collapse="$|^"), "$"), collapse = ""), taxa_melt$variable)] = "L7"
taxa_melt$sample.location[grep(paste(c("^", paste(meta_plot$sample.id[grep("^L8$", meta_plot$sample.location)], collapse="$|^"), "$"), collapse = ""), taxa_melt$variable)] = "L8"
taxa_melt$sample.location[grep(paste(c("^", paste(meta_plot$sample.id[grep("^L9$", meta_plot$sample.location)], collapse="$|^"), "$"), collapse = ""), taxa_melt$variable)] = "L9"
taxa_melt$sample.location[grep(paste(c("^", paste(meta_plot$sample.id[grep("^L10$", meta_plot$sample.location)], collapse="$|^"), "$"), collapse = ""), taxa_melt$variable)] = "L10"
taxa_melt$sample.location[grep(paste(c("^", paste(meta_plot$sample.id[grep("^L11$", meta_plot$sample.location)], collapse="$|^"), "$"), collapse = ""), taxa_melt$variable)] = "L11"
taxa_melt$sample.location[grep(paste(c("^", paste(meta_plot$sample.id[grep("^L12$", meta_plot$sample.location)], collapse="$|^"), "$"), collapse = ""), taxa_melt$variable)] = "L12"
# set factor
taxa_melt$sample.location = factor(taxa_melt$sample.location,
                                   levels = c("L1", "L2", "L3", "L4", "L5", "L6",
                                              "L7", "L8", "L9", "L10", "L11", "L12"))

## barplot phyla/genera -- split plots by time and depth
head(taxa_melt)
ggplot(taxa_melt,
  #taxa_melt[grep("L2|L5|L8|L11", taxa_melt$sample.location),],
  aes(x=sample.location, y=value, fill=Taxon)) +
  #aes(x=variable, y=value, fill=Taxon)) +
  geom_bar(stat="identity") +
  ylab("Relative Abundance") +
  xlab("Sample") +
  theme_bw() +
  facet_grid(sample.location.depth ~ time.of.day, scales="free", space="free") +
  scale_fill_manual(values=rev(c(
    "#FFCC99", "#FF9966", "#FF6600", "#993300", "#663300", # earth-tones
    "#66CCFF", "#3399FF", "#0000CC", "darkblue", # blues
    "#99FF99",  "#33CC33", "#006600", 
    #"#003300", # greens
    "yellow", "gold", 
    "purple",
    "pink",
    "#FF0033",
    "lightgray"))) +
  theme(#legend.position = "none",
    legend.text = element_text(size=13),
    strip.text = element_text(size=12),
    axis.text = element_text(size=12),
    axis.text.x = element_text(size=10, angle=45, hjust=1),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=13, color = "black"))


### Re-arrange plots

## barplot phyla/genera -- split plots by time and depth
head(taxa_melt)

# set labels
depth_names = list(
  "Bank_D1"="Bank: surface",
  "Interior_D1"="Interior: surface",
  "Interior_D2"="Interior: 1-m depth",
  "Interior_D3"="Interior: 2-m depth"
)
depth_labeller <- function(variable,value){
  return(depth_names[value])
}

# W-9am
W9 = ggplot(taxa_melt[grep("9", taxa_melt$time.of.day),],
            #aes(x=sample.location, y=value, fill=Taxon)) +
            aes(x=sample.location, y=value, fill=Taxon)) +
  geom_bar(stat="identity") +
  ylab("Relative Abundance") +
  xlab("Sample") +
  theme_bw() +
  facet_grid(~sample.location.depth, scales="free", space="free", labeller = depth_labeller) +
  scale_fill_manual(values=rev(c(
    "#FFCC99", "#FF9966", "#FF6600", "#993300", "#663300", # earth-tones
    "#66CCFF", "#3399FF", "#0000CC", "darkblue", # blues
    "#99FF99",  "#33CC33", "#006600", 
    #"#003300", # greens
    "yellow", "gold", 
    "purple",
    "pink",
    "#FF0033",
    "lightgray"))) +
  theme(legend.position = "none",
    legend.text = element_text(size=13),
    strip.text = element_text(size=11),
    axis.text = element_text(size=12),
    axis.text.x = element_text(size=10, angle=45, hjust=1),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=12, color = "black"))

# W-12pm
W12 = ggplot(taxa_melt[grep("12", taxa_melt$time.of.day),],
            #aes(x=sample.location, y=value, fill=Taxon)) +
            aes(x=sample.location, y=value, fill=Taxon)) +
  geom_bar(stat="identity") +
  ylab("Relative Abundance") +
  xlab("Sample") +
  theme_bw() +
  facet_grid(~sample.location.depth, scales="free", space="free", labeller = depth_labeller) +
  scale_fill_manual(values=rev(c(
    "#FFCC99", "#FF9966", "#FF6600", "#993300", "#663300", # earth-tones
    "#66CCFF", "#3399FF", "#0000CC", "darkblue", # blues
    "#99FF99",  "#33CC33", "#006600", 
    #"#003300", # greens
    "yellow", "gold", 
    "purple",
    "pink",
    "#FF0033",
    "lightgray"))) +
  theme(legend.position = "none",
    legend.text = element_text(size=13),
    strip.text = element_text(size=11),
    axis.text = element_text(size=12),
    axis.text.x = element_text(size=10, angle=45, hjust=1),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=12, color = "black"))

# W-3pm
W15 = ggplot(taxa_melt[grep("3", taxa_melt$time.of.day),],
            #aes(x=sample.location, y=value, fill=Taxon)) +
            aes(x=sample.location, y=value, fill=Taxon)) +
  geom_bar(stat="identity") +
  ylab("Relative Abundance") +
  xlab("Sample") +
  theme_bw() +
  facet_grid(~sample.location.depth, scales="free", space="free", labeller = depth_labeller) +
  scale_fill_manual(values=rev(c(
    "#FFCC99", "#FF9966", "#FF6600", "#993300", "#663300", # earth-tones
    "#66CCFF", "#3399FF", "#0000CC", "darkblue", # blues
    "#99FF99",  "#33CC33", "#006600", 
    #"#003300", # greens
    "yellow", "gold", 
    "purple",
    "pink",
    "#FF0033",
    "lightgray"))) +
  theme(legend.position = "none",
    legend.text = element_text(size=13),
    strip.text = element_text(size=11),
    axis.text = element_text(size=12),
    axis.text.x = element_text(size=10, angle=45, hjust=1),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=12, color = "black"))

# arrange plots (SUPPLEMENTAL FIGURE 2)
grid.arrange(W9, W12, W15, ncol=1)


### SEDIMENT SAMPLES

# melt data
colnames(taxa_plot_select)
taxa_melt = melt(data.frame(taxa_plot_select[, 1:10],
                            rownames(taxa_plot_select)))
colnames(taxa_melt)[1] = "Taxon"

# factor
taxa_melt$Taxon = factor(taxa_melt$Taxon,
                         levels = c("Other", 
                                    "p__Chloroflexi", 
                                    "p__Firmicutes",  
                                    "p__Actinobacteria", 
                                    "p__Verrucomicrobia", "  Verrucomicrobiae sp.",
                                    "p__Cyanobacteria", "  Microcystis PCC-7914",  "  Cyanobium PCC-6307",
                                    "p__Bacteroidetes", "  Solitalea", "  Fluviicola", "  Chitinophagaceae sp.",
                                    "p__Proteobacteria",   "  Limnohabitans", "  Methylococcaceae sp.", "  Methylocystis", "  Steroidobacter"
                         ))

# set factor
taxa_melt$variable = factor(taxa_melt$variable,
                                   levels = c("S2", "S3", "S4", "S5", "S6",
                                              "S7", "S8", "S9", "S10", "S11"))

# add sampling location
taxa_melt$location = "Bank"
taxa_melt$location[grep("S2|S5|S8|S11", taxa_melt$variable)] = "Interior"

## barplot phyla/genera -- split plots by time and depth (SUPPLEMENTAL FIGURE 2)
head(taxa_melt)
ggplot(taxa_melt,
       aes(x=variable, y=value, fill=Taxon)) +
  #aes(x=variable, y=value, fill=Taxon)) +
  geom_bar(stat="identity") +
  ylab("Relative Abundance") +
  xlab("Sample") +
  theme_bw() +
  facet_grid(~location, scales="free", space="free") +
  scale_fill_manual(values=rev(c(
    "#FFCC99", "#FF9966", "#FF6600", "#993300", "#663300", # earth-tones
    "#66CCFF", "#3399FF", "#0000CC", "darkblue", # blues
    "#99FF99",  "#33CC33", "#006600", 
    #"#003300", # greens
    "yellow", "gold", 
    "purple",
    "pink",
    "#FF0033",
    "lightgray"))) +
  theme(legend.position = "none",
    legend.text = element_text(size=13),
    strip.text = element_text(size=13),
    axis.text = element_text(size=15),
    axis.text.x = element_text(size=14, angle=45, hjust=1),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=15, color = "black"))


###
### AVERAGE BY TIME & DEPTH
###

colnames(taxa_plot_select)
# melt data
taxa_melt = melt(data.frame(taxa_plot_select[, 79:92],
                            rownames(taxa_plot_select)))
colnames(taxa_melt)[1] = "Taxon"

# add factor for time
taxa_melt$time.of.day = "T1_9am"
taxa_melt$time.of.day[grep("T2", taxa_melt$variable)] = "T2_12pm"
taxa_melt$time.of.day[grep("T3", taxa_melt$variable)] = "T3_3pm"
taxa_melt$time.of.day_N = 1
taxa_melt$time.of.day_N[grep("T2", taxa_melt$variable)] = 2
taxa_melt$time.of.day_N[grep("T3", taxa_melt$variable)] = 3

# add factor for depth
taxa_melt$sample.location.depth = "Bank_D1"
taxa_melt$sample.location.depth[grep("Interior_D1", taxa_melt$variable)] = "Interior_D1"
taxa_melt$sample.location.depth[grep("Interior_D2", taxa_melt$variable)] = "Interior_D2"
taxa_melt$sample.location.depth[grep("Interior_D3", taxa_melt$variable)] = "Interior_D3"

# add factor for sediment
taxa_melt$sample.location.depth[grep("Sed_bank", taxa_melt$variable)] = "Bank"
taxa_melt$sample.location.depth[grep("Sed_interior", taxa_melt$variable)] = "Interior"

# set factor for taxa
taxa_melt$Taxon = factor(taxa_melt$Taxon,
                         levels = c("Other", 
                                    "p__Chloroflexi", 
                                    "p__Firmicutes",  
                                    "p__Actinobacteria", 
                                    "p__Verrucomicrobia", "  Verrucomicrobiae sp.",
                                    "p__Cyanobacteria", "  Microcystis PCC-7914",  "  Cyanobium PCC-6307",
                                    "p__Bacteroidetes", "  Solitalea", "  Fluviicola", "  Chitinophagaceae sp.",
                                    "p__Proteobacteria",   "  Limnohabitans", "  Methylococcaceae sp.", "  Methylocystis", "  Steroidobacter"
                         ))

## WATER -- barplot phyla/genera (FIGURE 2-A)
ggplot(taxa_melt[grep("T", taxa_melt$variable),],
       aes(x=as.numeric(time.of.day_N), y=value, fill=Taxon)) +
  #geom_bar(stat="identity") +
  geom_area(alpha=0.9, size = 0.2, colour="black") +
  ylab("Avg. Relative Abundance") +
  xlab("Sample") +
  theme_bw() +
  scale_fill_manual(values=rev(c(
    "#FFCC99", "#FF9966", "#FF6600", "#993300", "#663300", # earth-tones
    "#66CCFF", "#3399FF", "#0000CC", "darkblue", # blues
    "#99FF99",  "#33CC33", "#006600", 
    #"#003300", # greens
    "yellow", "gold", 
    "purple",
    "pink",
    "#FF0033",
    "lightgray"))) +
  facet_grid(~sample.location.depth, scales="free", space="free", labeller = depth_labeller) +
  scale_x_continuous(breaks = c(1,2,3), limits = c(0.8,3.2), 
                     labels = c("9:00", "12:00", "15:00")) +
  theme(legend.text = element_text(size=13),
        legend.title = element_blank(),
        strip.text = element_text(size=10),
        axis.text = element_text(size=15),
        axis.text.x = element_text(size=15, angle=45, hjust=1),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=15, color = "black"))

## Sediment -- barplot phyla/genera (FIGURE 2-B)
ggplot(taxa_melt[grep("Sed", taxa_melt$variable),],
       aes(x=time.of.day, y=value, fill=Taxon)) +
  geom_bar(stat="identity") +
  #geom_area(alpha=0.9, size = 0.2, colour="black") +
  ylab("Avg. Relative Abundance") +
  xlab("Sample") +
  theme_bw() +
  facet_grid(~location, scales="free", space="free") +
  scale_fill_manual(values=rev(c(
    "#FFCC99", "#FF9966", "#FF6600", "#993300", "#663300", # earth-tones
    "#66CCFF", "#3399FF", "#0000CC", "darkblue", # blues
    "#99FF99",  "#33CC33", "#006600", 
    #"#003300", # greens
    "yellow", "gold", 
    "purple",
    "pink",
    "#FF0033",
    "lightgray"))) +
  facet_grid(~sample.location.depth, scales="free", space="free") +
  scale_x_discrete(labels = c("15:00")) +
  theme(legend.text = element_text(size=13),
        legend.title = element_blank(),
        strip.text = element_text(size=13),
        axis.text = element_text(size=15),
        axis.text.x = element_text(size=15, angle=45, hjust=1),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=15, color = "black")) +
  theme(legend.position = "none")


#################################
### RELATIVE ABUNDANCE: PHYLA ###
#################################

###
### Phylum
###

## set table
taxa_plot = taxa_phy_RA

## subset water/sediment samples
taxa_plot = taxa_plot[,grep("Wye", meta$sample.type.pond)]
meta_plot = meta[grep("Wye", meta$sample.type.pond),]

## add 'average_water' column and re-order table
taxa_plot$avg_water = apply(taxa_plot[,grep("Wye_water", meta$sample.type.pond)], 1, mean)
taxa_plot = taxa_plot[order(taxa_plot$avg_water, decreasing = T),]
head(taxa_plot)

## add 'average_sediment' column and re-order table
taxa_plot$avg_sed = apply(taxa_plot[,grep("Wye_sed", meta$sample.type.pond)], 1, mean)
taxa_plot = taxa_plot[order(taxa_plot$avg_sed, decreasing = T),]
head(taxa_plot)


######################
### STATS -- WATER ###
######################

###
### Phylum
###

## set table
taxa_plot = taxa_phy_RA

## subset water/sediment samples
taxa_plot = taxa_plot[,grep("Wye_water", meta$sample.type.pond)]
meta_plot = meta[grep("Wye_water", meta$sample.type.pond),]

## add 'average_water' column and re-order table
taxa_plot$avg_water = apply(taxa_plot, 1, mean)
taxa_plot = taxa_plot[order(taxa_plot$avg_water, decreasing = T),]
taxa_plot = taxa_plot[,-c(61)]
taxa_plot = taxa_plot[which(apply(taxa_plot, 1, sum) > 0),]
head(taxa_plot, 10)

## stats
p_phyla = as.data.frame(matrix(nrow = nrow(taxa_plot), ncol = 6))
for (i in 1:nrow(taxa_plot)) {
  p_phyla[i,] = unlist(summary(aov(as.numeric(taxa_plot[i,]) ~ sample.location.depth*time.of.day, meta_plot)))[17:19]
}
colnames(p_phyla) = c("p_depth", "p_time", "p_interaction",
                      "q_depth", "q_time", "q_interaction")
rownames(p_phyla) = rownames(taxa_plot)
p_phyla[,4] = p.adjust(p_phyla[,1], method = "BH")
p_phyla[,5] = p.adjust(p_phyla[,2], method = "BH")
p_phyla[,6] = p.adjust(p_phyla[,3], method = "BH")
# check significance
p_phyla[which(p_phyla[,4] < 0.05),]
p_phyla[which(p_phyla[,5] < 0.05),]
p_phyla[which(p_phyla[,6] < 0.05),]

## combine and export (SUPPLEMENTARY TABLE 3)
write.table(data.frame(taxa_plot, p_phyla), "stat_water_taxa_phy.txt", sep = "\t", quote = FALSE)


###
### Genus
###

## set table
taxa_plot = taxa_gen_RA

## subset water/sediment samples
taxa_plot = taxa_plot[,grep("Wye_water", meta$sample.type.pond)]
meta_plot = meta[grep("Wye_water", meta$sample.type.pond),]

## add 'average_water' column and re-order table
taxa_plot$avg_water = apply(taxa_plot, 1, mean)
taxa_plot = taxa_plot[order(taxa_plot$avg_water, decreasing = T),]
taxa_plot = taxa_plot[,-c(61)]
taxa_plot = taxa_plot[which(apply(taxa_plot, 1, sum) > 0),]
head(taxa_plot)

## stats
p_gen = as.data.frame(matrix(nrow = nrow(taxa_plot), ncol = 6))
for (i in 1:nrow(taxa_plot)) {
  p_gen[i,] = unlist(summary(aov(as.numeric(taxa_plot[i,]) ~ sample.location.depth*time.of.day, meta_plot)))[17:19]
}
colnames(p_gen) = c("p_depth", "p_time", "p_interaction",
                      "q_depth", "q_time", "q_interaction")
rownames(p_gen) = rownames(taxa_plot)
p_gen[,4] = p.adjust(p_gen[,1], method = "BH")
p_gen[,5] = p.adjust(p_gen[,2], method = "BH")
p_gen[,6] = p.adjust(p_gen[,3], method = "BH")
# check significance
p_gen[which(p_gen[,4] < 0.05),]
p_gen[which(p_gen[,5] < 0.05),]
p_gen[which(p_gen[,6] < 0.05),]

## combine and export (SUPPLEMENTARY TABLE 3)
write.table(data.frame(taxa_plot, p_gen), "stat_water_taxa_gen.txt", sep = "\t", quote = FALSE)


#########################
### STATS -- SEDIMENT ###
#########################

###
### Phylum
###

## set table
taxa_plot = taxa_phy_RA

## subset water/sediment samples
taxa_plot = taxa_plot[,grep("Wye_sed", meta$sample.type.pond)]
meta_plot = meta[grep("Wye_sed", meta$sample.type.pond),]

## add 'average_water' column and re-order table
taxa_plot$avg_sed = apply(taxa_plot, 1, mean)
taxa_plot = taxa_plot[order(taxa_plot$avg_sed, decreasing = T),]
taxa_plot = taxa_plot[,-c(11)]
taxa_plot = taxa_plot[which(apply(taxa_plot, 1, sum) > 0),]
head(taxa_plot)

## stats
p_phyla = as.data.frame(matrix(nrow = nrow(taxa_plot), ncol = 2))
for (i in 1:nrow(taxa_plot)) {
  p_phyla[i,] = t.test(as.numeric(taxa_plot[i, grep("Bank", meta_plot$sample.water.location)]),
                       as.numeric(taxa_plot[i, grep("Interior", meta_plot$sample.water.location)]))$p.val
}
colnames(p_phyla) = c("p_location",
                      "q_location")
rownames(p_phyla) = rownames(taxa_plot)
p_phyla[,2] = p.adjust(p_phyla[,1], method = "BH")
# check significance
p_phyla[which(p_phyla[,1] < 0.05),]
p_phyla[which(p_phyla[,2] < 0.05),]

## combine and export (SUPPLEMENTARY TABLE 4)
write.table(data.frame(taxa_plot, p_phyla), "stat_sediment_taxa_phy.txt", sep = "\t", quote = FALSE)


###
### Genus
###

## set table
taxa_plot = taxa_gen_RA

## subset water/sediment samples
taxa_plot = taxa_plot[,grep("Wye_sed", meta$sample.type.pond)]
meta_plot = meta[grep("Wye_sed", meta$sample.type.pond),]

## add 'average_water' column and re-order table
taxa_plot$avg_sed = apply(taxa_plot, 1, mean)
taxa_plot = taxa_plot[order(taxa_plot$avg_sed, decreasing = T),]
taxa_plot = taxa_plot[,-c(11)]
taxa_plot = taxa_plot[which(apply(taxa_plot, 1, sum) > 0),]
head(taxa_plot)

## stats
p_gen = as.data.frame(matrix(nrow = nrow(taxa_plot), ncol = 2))
for (i in 1:nrow(taxa_plot)) {
  p_gen[i,] = t.test(as.numeric(taxa_plot[i, grep("Bank", meta_plot$sample.water.location)]),
                     as.numeric(taxa_plot[i, grep("Interior", meta_plot$sample.water.location)]))$p.val
}
colnames(p_gen) = c("p_location",
                    "q_location")
rownames(p_gen) = rownames(taxa_plot)
p_gen[,2] = p.adjust(p_gen[,1], method = "BH")
# check significance
p_gen[which(p_gen[,1] < 0.05),]
p_gen[which(p_gen[,2] < 0.05),]

## combine and export (SUPPLEMENTARY TABLE 4)
write.table(data.frame(taxa_plot, p_gen), "stat_sediment_taxa_gen.txt", sep = "\t", quote = FALSE)


#############################
### Metadata correlations ###
#############################

###
### Alpha
###

alpha_cor = as.data.frame(matrix(nrow = length(colnames(meta)[c(19:47)]), 
                                 ncol=3))
for (i in c(19:47)) {
  alpha_cor[i,1] = colnames(meta)[i]
  alpha_cor[i,2] = cor.test(meta$Shannon, meta[,i], method = "spearman")$estimate
  alpha_cor[i,3] = cor.test(meta$Shannon, meta[,i], method = "spearman")$p.val
}
alpha_cor # (TABLE-1)


###
### Beta (TABLE-1)
###

## distinguish metadata with numeric values (avoid samples with NAs)
# inspect
colnames(meta)
colnames(meta)[19:47]
# all values present
colnames(meta)[c(35:42,44:47)]
# remove sample 'W18'
colnames(meta)[c(19:30)]

## beta-diversity (all samples)
beta <- vegdist(t(taxa_table[,grep("Wye_water", meta$sample.type.pond)]), 'bray', binary = T)

# stats
adonis(beta ~ PO4, data = meta[grep("Wye_water", meta$sample.type.pond),], permutations = 999)$aov.tab
adonis(beta ~ NH3, data = meta[grep("Wye_water", meta$sample.type.pond),], permutations = 999)$aov.tab
adonis(beta ~ NO3, data = meta[grep("Wye_water", meta$sample.type.pond),], permutations = 999)$aov.tab
adonis(beta ~ ECCFU, data = meta[grep("Wye_water", meta$sample.type.pond),], permutations = 999)$aov.tab
adonis(beta ~ ECMPN, data = meta[grep("Wye_water", meta$sample.type.pond),], permutations = 999)$aov.tab
adonis(beta ~ TC, data = meta[grep("Wye_water", meta$sample.type.pond),], permutations = 999)$aov.tab
adonis(beta ~ TIC, data = meta[grep("Wye_water", meta$sample.type.pond),], permutations = 999)$aov.tab
adonis(beta ~ TN, data = meta[grep("Wye_water", meta$sample.type.pond),], permutations = 999)$aov.tab
adonis(beta ~ TOC, data = meta[grep("Wye_water", meta$sample.type.pond),], permutations = 999)$aov.tab
adonis(beta ~ CDOM, data = meta[grep("Wye_water", meta$sample.type.pond),], permutations = 999)$aov.tab
adonis(beta ~ INV, data = meta[grep("Wye_water", meta$sample.type.pond),], permutations = 999)$aov.tab
adonis(beta ~ PHYC, data = meta[grep("Wye_water", meta$sample.type.pond),], permutations = 999)$aov.tab


## beta-diversity (without W18)
taxa_table[1,grep("Wye_water", meta$sample.type.pond)][,-c(10)]
beta <- vegdist(t(taxa_table[,grep("Wye_water", meta$sample.type.pond)][,-c(10)]), 'bray', binary = T)

# stats
adonis(beta ~ Temperature, data = meta[grep("Wye_water", meta$sample.type.pond),][-c(10),], permutations = 999)$aov.tab
adonis(beta ~ DO.mg.L, data = meta[grep("Wye_water", meta$sample.type.pond),][-c(10),], permutations = 999)$aov.tab
adonis(beta ~ SPC.uS.cm, data = meta[grep("Wye_water", meta$sample.type.pond),][-c(10),], permutations = 999)$aov.tab
adonis(beta ~ pH, data = meta[grep("Wye_water", meta$sample.type.pond),][-c(10),], permutations = 999)$aov.tab
adonis(beta ~ UVNOx.mg.L, data = meta[grep("Wye_water", meta$sample.type.pond),][-c(10),], permutations = 999)$aov.tab
adonis(beta ~ NTU, data = meta[grep("Wye_water", meta$sample.type.pond),][-c(10),], permutations = 999)$aov.tab
adonis(beta ~ BGA.PC.RFU, data = meta[grep("Wye_water", meta$sample.type.pond),][-c(10),], permutations = 999)$aov.tab
adonis(beta ~ BGA.PC.ug.L, data = meta[grep("Wye_water", meta$sample.type.pond),][-c(10),], permutations = 999)$aov.tab
adonis(beta ~ Chl.RFU, data = meta[grep("Wye_water", meta$sample.type.pond),][-c(10),], permutations = 999)$aov.tab
adonis(beta ~ Chl.ug.L, data = meta[grep("Wye_water", meta$sample.type.pond),][-c(10),], permutations = 999)$aov.tab
adonis(beta ~ fDOM.RFU, data = meta[grep("Wye_water", meta$sample.type.pond),][-c(10),], permutations = 999)$aov.tab
adonis(beta ~ fDOM.ppb, data = meta[grep("Wye_water", meta$sample.type.pond),][-c(10),], permutations = 999)$aov.tab


## beta-diversity (NN)
taxa_table[1,grep("Wye_water", meta$sample.type.pond)][,-c(which(is.na(meta[grep("Wye_water", meta$sample.type.pond),]$NN) == TRUE))]
beta <- vegdist(t(taxa_table[,grep("Wye_water", meta$sample.type.pond)][,-c(which(is.na(meta[grep("Wye_water", meta$sample.type.pond),]$NN) == TRUE))]), 'bray', binary = T)

# stats
adonis(beta ~ NN, data = meta[grep("Wye_water", meta$sample.type.pond),][-c(which(is.na(meta[grep("Wye_water", meta$sample.type.pond),]$NN) == TRUE)),], permutations = 999)$aov.tab

