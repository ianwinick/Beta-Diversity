library(tidyverse)
library(vegan)
library(FD)
library(devtools)
install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)

data <- read_csv("CommunityMatrix2023.csv") %>%
  column_to_rownames(var="plot")

severity <- data %>%
  mutate(severity = case_when(
    row.names(.) %in% 1:20 ~ "U",
    row.names(.) %in% 21:40 ~ "L",
    row.names(.) %in% 41:60 ~ "H"
   )) %>%
  select(severity)

tax_nmds <- metaMDS(data, distance="bray", k=2, trymax=100)
tax_nmds$stress
stressplot(tax_nmds)

scores <- as_tibble(scores(tax_nmds, display="sites"), rownames="sites") %>%
  mutate(severity=case_when(
    sites %in% c(1:20) ~ "U",
    sites %in% c(21:40) ~ "L",
    sites %in% c(41:60) ~ "H"
  ))
avg <- scores %>%
  mutate(severity=factor(severity, levels=c("U", "L", "H"))) %>%
  group_by(severity) %>%
  dplyr::summarize(NMDS1=mean(NMDS1), NMDS2=mean(NMDS2))
tnmds <- scores %>%
  mutate(severity=factor(severity, levels=c("U", "L", "H"))) %>%
  ggplot(aes(x=NMDS1, y=NMDS2, linetype=severity, color=severity)) +
  stat_ellipse(geom="polygon", aes(group=severity, fill=severity), alpha = 0.3, size=.75) +
  geom_point() +
  scale_linetype_manual(values=c("solid", "dotdash", "dotted")) +
  scale_color_manual(values=c("cadetblue3", "darkgoldenrod3", "black")) +
  scale_fill_manual(values=c("cadetblue3", "darkgoldenrod3", "black")) +
  theme_classic() +
  theme(axis.text.x=element_text(color="black", size=8, family="Arial"), 
        axis.text.y=element_text(color="black", size=8, family="Arial"),
        axis.title.x=element_text(color="black", size=10, family="Arial"),
        axis.title.y=element_text(color="black", size=10, family="Arial"),
        legend.position=0)

tax_permanova <- adonis2(vegdist(data, method="bray") ~ severity, permutations=9999, data=severity)
write.csv(pairwise.adonis(vegdist(data, method="bray"), severity$severity, perm=9999), "taxPERM.csv")
anova(betadisper(vegdist(data, method="bray"), severity$severity, type="centroid"))
TukeyHSD(betadisper(vegdist(data, method="bray"), severity$severity, type="centroid"))
permutest(betadisper(vegdist(data, method="bray"), severity$severity, type="centroid"))


################################################################################
# Community Weighted Means #####################################################
################################################################################
traitMatrix <- read_csv("traitMatrix.csv") %>%
  mutate(ldmc=log(ldmc)) %>%
  mutate(height=log(height)) %>%
  mutate(sla=log(sla)) %>%
  filter(spp %in% colnames(data)) %>%
  column_to_rownames(var="spp")

hist(traitMatrix$ldmc)
hist(traitMatrix$height)
hist(traitMatrix$sla)
cor.test(traitMatrix$sla, traitMatrix$ldmc, method="spearman")
ggplot(traitMatrix, aes(sla, ldmc)) +
  geom_point() +
  geom_smooth(method="lm", se=F) +
  theme_classic()
summary(lm(sla~ldmc, traitMatrix))

cwm <- dbFD(traitMatrix, data)$CWM

trait_nmds <- metaMDS(cwm, distance="euclidean", k=2, trymax=100)
trait_nmds$stress
stressplot(trait_nmds)

scores <- as_tibble(scores(trait_nmds, display="sites"), rownames="sites") %>%
  mutate(severity=case_when(
    sites %in% c(1:20) ~ "U",
    sites %in% c(21:40) ~ "L",
    sites %in% c(41:60) ~ "H"
  ))
avg <- scores %>%
  mutate(severity=factor(severity, levels=c("L", "H", "U"))) %>%
  group_by(severity) %>%
  summarise(NMDS1=mean(NMDS1), NMDS2=mean(NMDS2))
fnmds <- scores %>%
  mutate(severity=case_when(
    severity=="U" ~ "Unburned",
    severity=="L" ~ "Low",
    severity=="H" ~ "High"
  ))%>%
  mutate(severity=factor(severity, levels=c("Unburned", "Low", "High"))) %>%
  ggplot(aes(x=NMDS1, y=NMDS2, linetype=severity, color=severity)) +
  stat_ellipse(geom="polygon", aes(group=severity, fill=severity), alpha = 0.3, size=.75) +
  geom_point() +
  scale_linetype_manual(values=c("solid", "dotdash", "dotted")) +
  scale_color_manual(values=c("cadetblue3", "darkgoldenrod3", "black")) +
  scale_fill_manual(values=c("cadetblue3", "darkgoldenrod3", "black")) +
  labs(color="", fill="", linetype="") +
  theme_classic() +
  theme(axis.text.x=element_text(color="black", size=8, family="Arial"), 
        axis.text.y=element_text(color="black", size=8, family="Arial"),
        axis.title.x=element_text(color="black", size=10, family="Arial"),
        axis.title.y=element_text(color="black", size=10, family="Arial"),
        legend.text=element_text(color="black", size=8, family="Arial"),
        legend.title=element_text(color="black", size=10, family="Arial"),
        legend.position="bottom")

tiff("Fig2.tiff", width=8.5, height=14, units="cm", res=600)
(tnmds/fnmds)
dev.off()

fun_permanova <- adonis2(vegdist(cwm, method="euclidean") ~ severity, permutations=9999, data=severity)
write.csv(pairwise.adonis(vegdist(cwm, method="euclidean"), severity$severity, perm=9999), "funPERM.csv")
anova(betadisper(vegdist(cwm, method="euclidean"), severity$severity, type="centroid"))
TukeyHSD(betadisper(vegdist(cwm, method="euclidean"), severity$severity, type="centroid"))

# Removing plot outliar changes results