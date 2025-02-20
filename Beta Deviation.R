library(tidyverse)
library(lme4)
library(lmerTest)
library(DHARMa)
library(performance)
library(lsmeans)
library(mgcv)
library(rmcorr)
library(patchwork)
library(broom)

fonttable()
loadfonts(device = "win")

obs <- read_csv("Bootstrap Observed Beta Values/bootstrap_beta.csv")

null <- read_csv("Null Model/nullBeta_allyears.csv")

null <- null %>%
  pivot_longer(cols=c(tax, fun), names_to="type", values_to="null") %>%
  group_by(year, severity, type) %>%
  summarise(sd=sd(null), null=mean(null)) %>%
  select(year, type, severity, null, sd)

data <- obs %>%
  left_join(null) %>%
  mutate(dev=NA, across(dev, ~(beta-null)/sd)) %>%
  mutate(severity=factor(severity, levels=c("U", "L", "H")))

deviation <- data %>%
  group_by(year, type, severity) %>%
  summarise(sd=sd(dev), dev=mean(dev))

# Plot taxonomic beta deviation by severity
t <- deviation %>%
  filter(type=="tax") %>%
  ggplot(aes(year, dev, color=severity)) +
  geom_point(size=2) +
  geom_errorbar(aes(ymin=dev-sd, ymax=dev+sd), width=0.05) +
  geom_smooth(aes(linetype=severity), method="lm", formula=y~poly(x,2), se=FALSE) +
  scale_linetype_manual(values=c("solid", "dotdash", "dotted")) +
  scale_color_manual(values=c("cadetblue3", "darkgoldenrod3", "black")) +
  ylim(-41, -0) +
  xlab("Years since fire") +
  ylab("β-deviation") +
  theme_minimal() +
  theme(panel.background = element_rect(fill="ivory2", color="white"),
        panel.grid.major=element_line(color="seashell3"),
        axis.text.x=element_text(color="black", size=8, family="Arial"), 
        axis.text.y=element_text(color="black", size=8, family="Arial"),
        axis.title.x=element_text(color="black", size=10, family="Arial"),
        axis.title.y=element_text(color="black", size=10, family="Arial"),
        legend.position=0)

 # Plot functional beta deviation by severity
f <- deviation %>%
  filter(type=="fun") %>%
  mutate(severity=case_when(
    severity=="U" ~ "Unburned",
    severity=="L" ~ "Low",
    severity=="H" ~ "High"
  ))%>%
  mutate(severity=factor(severity, levels=c("Unburned", "Low", "High"))) %>%
  ggplot(aes(year, dev, color=severity)) +
  geom_point(size=2, show.legend=F) +
  geom_errorbar(aes(ymin=dev-sd, ymax=dev+sd), width=0.05, show.legend=F) +
  geom_smooth(aes(linetype=severity), method="lm", formula=y~poly(x,2), se=FALSE) +
  scale_linetype_manual(values=c("solid", "dotdash", "dotted")) + 
  scale_color_manual(values=c("cadetblue3", "darkgoldenrod3", "black"))+
  ylim(-2.9, 6.3) +
  xlab("Years since fire") +
  ylab("β-deviation") +
  labs(color="", linetype="") +
  theme_classic() +
  theme(panel.background = element_rect(fill="ivory2", color="white"),
        panel.grid.major=element_line(color="seashell3"),
        axis.text.x=element_text(color="black", size=8, family="Arial"), 
        axis.text.y=element_text(color="black", size=8, family="Arial"),
        axis.title.x=element_text(color="black", size=10, family="Arial"),
        axis.title.y=element_text(color="black", size=10, family="Arial"),
        legend.text=element_text(color="black", size=8, family="Arial"),
        legend.title=element_text(color="black", size=10, family="Arial"),
        legend.key.width=unit(1.2,"cm"),
        legend.position="bottom") +
  guides(color=guide_legend(nrow=3)) +
  guides(linetype=guide_legend(nrow=3))

#fig1 <- (t/f) +
 # plot_annotation(tag_levels = 'A')

tiff("Fig1.tiff", width=8.5, height=14, units="cm", res=600)
(t/f)
dev.off()

# Polynomial regressions
lm_tax <- lm(dev~year*severity*I(year^2), data[data$type=="tax",])
lm_fun <- lm(dev~year*severity*I(year^2), data[data$type=="fun",])
lm_tax <- lm(dev~year*severity, data[data$type=="tax",])
lm_fun <- lm(dev~year*severity, data[data$type=="fun",])

summary(lm_tax)
summary(lm_fun)

AIC(lm_fun)

write.csv(tidy(summary(lm_tax)), "lm_tax.csv")
write.csv(tidy(summary(lm_fun)), "lm_fun.csv")
