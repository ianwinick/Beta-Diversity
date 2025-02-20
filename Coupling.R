library(tidyverse)
library(vegan)
library(FD)
library(lme4)
library(lmerTest)
library(extrafont)
library(Hmisc)

df20 <- read_csv("CommunityMatrix2020.csv") %>%
  mutate(year=1)

df21 <- read_csv("CommunityMatrix2021.csv") %>%
  mutate(year=2)

df22 <- read_csv("CommunityMatrix2022.csv") %>%
  mutate(year=3)

df23 <- read_csv("CommunityMatrix2023.csv") %>%
  mutate(year=4)

traitMatrix <- read_csv("traitMatrix.csv") %>%
  column_to_rownames(var="spp") %>%
  arrange(row.names(.))

library(plyr)
df <- tibble(rbind.fill(df20, df21, df22, df23)) %>%
  mutate_all(~replace(., is.na(.), 0)) %>%
  select(sort(colnames(.))) %>%
  mutate(severity=factor(case_when(
    plot %in% 1:20 ~ "U",
    plot %in% 21:40 ~ "L",
    plot %in% 41:60 ~ "H"
  ), levels=c("U", "L", "H"))) %>%
  select(plot, year, severity, everything(.)) %>%
  pivot_longer(cols=c(-plot, -severity, -year), names_to="spp", values_to="pres") %>%
  mutate(pres=ifelse(is.na(pres), 0, pres)) %>%
  pivot_wider(names_from=spp, values_from=pres)
detach("package:plyr", unload = TRUE)

rich <- df %>%
  pivot_longer(cols=c(-plot, -severity, -year), names_to="spp", values_to="pres") %>%
  mutate(pres=ifelse(pres>0, 1, 0)) %>%
  group_by(plot, year) %>%
  dplyr::summarize(rich=sum(pres))

comm <- df %>%
  select(-plot, -severity, -year)

Q <- dbFD(traitMatrix, comm, scale.RaoQ=T, message=F)$RaoQ

wright <- tibble("plot"=df$plot, "year"=df$year, "Q"=Q) %>%
  left_join(rich) %>%
  mutate(rich=(rich-min(rich))/(max(rich)-min(rich))) %>%
  mutate(severity=factor(case_when(
    plot %in% 1:20 ~ "U",
    plot %in% 21:40 ~ "L",
    plot %in% 41:60 ~ "H"
  ), levels=c("U", "L", "H"))) %>%
  unite(sev_year, c(severity, year), sep="_") %>%
  filter(Q>0)

sev_year <- c("U_1", "U_2", "U_3", "U_4", "L_1", "L_2", "L_3", "L_4", "H_1", "H_2", "H_3", "H_4")
reg.e = NULL
reg.r2 = NULL
reg.ar2 = NULL
reg.p = NULL

for(i in sev_year){
  x <- wright %>%
    mutate(Q=(Q-mean(Q))/sd(Q), rich=(rich-mean(rich))/sd(rich)) %>%
    filter(sev_year==i)
  
  reg <- summary(lm(Q~rich, x))
  reg.e[i] = coef(reg)["rich", "Estimate"]
  reg.r2[i] = reg$r.squared
  reg.ar2[i] = reg$adj.r.squared
  reg.p[i] = coef(reg)["rich", "Pr(>|t|)"]
}

reg.data <- data.frame("e"=reg.e, "r2"=reg.r2, "ar2"=reg.ar2, "p"=reg.p)
reg.data <- rownames_to_column(reg.data, var="id")
reg.data <- reg.data %>%
  separate(id, c("severity", "year")) %>%
  mutate(severity=factor(severity, levels=c("U", "L", "H")))
reg.data %>%
  #filter(p<=0.05)
  filter(severity=="U")
write.csv(reg.data, "reg.data.csv", row.names=F)

tiff("Fig3.tiff", width=8.5, height=7, units="cm", res=600)

reg.data %>%
  mutate(severity=factor(case_when(
    severity=="U" ~ "Unburned",
    severity=="L" ~ "Low severity",
    severity=="H" ~ "High severity"
  ), levels=c("Unburned", "Low severity", "High severity"))) %>%
  ggplot(aes(year, e)) +
  geom_point(aes(color=severity), size=2) +
  #ylim(-0.25,1) +
  geom_hline(aes(yintercept=0), linetype="dashed", color="black") +
  xlab("Years since fire") +
  ylab("Pearson correlation coefficient (r)") +
  scale_color_manual(values=c("cadetblue3", "darkgoldenrod3", "black")) +
  facet_wrap(~severity) +
  theme_minimal() +
  theme(axis.text.x=element_text(color="black", size=8, family="Arial"), 
        axis.text.y=element_text(color="black", size=8, family="Arial"),
        axis.title.x=element_text(color="black", size=10, family="Arial"),
        axis.title.y=element_text(color="black", size=10, family="Arial"),
        strip.text.x = element_text(size = 10, family="Arial"),
        legend.position=0)
dev.off()

cor.data %>%
  filter(P<=0.05)

################################################################################
# Disentangle using FD #########################################################
################################################################################
fd <- data.frame(
  "plot"=df$plot,
  "year"=df$year,
  "severity"=df$severity,
  "FRic"=dbFD(traitMatrix, comm, scale.RaoQ=T, message=F)$FRic,
  "FEve"=dbFD(traitMatrix, comm, scale.RaoQ=T, message=F)$FEve,
  "FDis"=dbFD(traitMatrix, comm, scale.RaoQ=T, message=F)$FDis,
  "FDiv"=dbFD(traitMatrix, comm, scale.RaoQ=T, message=F)$FDiv,
  "Q"=dbFD(traitMatrix, comm, scale.RaoQ=T, message=F)$RaoQ)
fd <- fd %>%
  left_join(rich) %>%
  filter(rich>1)

# Regression plotted by year
wright %>%
  separate(sev_year, into=c("severity", "year")) %>%
  mutate(severity=factor(severity, levels=c("U", "L", "H"))) %>%
  ggplot(aes(rich, Q, color=severity)) +
  geom_point() +
  geom_smooth(method="lm", se=F) +
  facet_wrap(~year) +
  scale_color_manual(values=c("cadetblue3", "darkgoldenrod3", "black")) +
  theme_bw()
# Regression plotted by severity, colored by severity:year
wright %>%
  separate(sev_year, into=c("severity", "year")) %>%
  mutate(severity=factor(severity, levels=c("U", "L", "H"))) %>%
  mutate(year=factor(year)) %>%
  ggplot(aes(rich, Q, color=severity:year)) +
  geom_point() +
  geom_smooth(method="lm", se=F) +
  facet_wrap(~severity) +
  scale_color_manual(values=c("cadetblue1", "cadetblue2", "cadetblue3", "cadetblue4", "darkgoldenrod1", "darkgoldenrod2", "darkgoldenrod3", "darkgoldenrod4", "gray80", "gray60", "gray40", "black")) +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "black", color="white"),
        panel.grid.major=element_line(color="seashell2"),
        panel.grid.minor=element_line(color=NA))
# Regression plotted by severity, colored by year
fill <- c("cadetblue3", "darkgoldenrod3", "black")
wright %>%
  separate(sev_year, into=c("severity", "year")) %>%
  mutate(severity=factor(severity, levels=c("U", "L", "H"))) %>%
  mutate(year=factor(year)) %>%
  ggplot(aes(rich, Q, color=year)) +
  #geom_rect(aes(fill=severity), xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
  geom_point() +
  geom_smooth(method="lm", se=F) +
  facet_wrap(~severity, ncol=1) +
  #scale_color_manual(values=c("orange1", "orange2", "orange3", "orange4")) +
  #scale_color_manual(values=c("purple1", "purple2", "purple3", "purple4")) +
  scale_color_manual(values=c("paleturquoise3", "red4", "purple3", "black")) +
  theme_minimal() +
  theme(panel.background = element_rect(fill="ivory2", color="white"),
      panel.grid.major=element_line(color="seashell3"),
      panel.grid.minor=element_line(color=NA))