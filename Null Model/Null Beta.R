library(tidyverse)
library(vegan)
library(FD)

Community_longForm_2020 <- read_csv("Community_longForm_2020.csv")
Community_longForm_2021 <- read_csv("Community_longForm_2021.csv")
Community_longForm_2022 <- read_csv("Community_longForm_2022.csv")
Community_longForm_2023 <- read_csv("Community_longForm_2023.csv")

traitmatrix <- read_csv("traitMatrix.csv") %>%
  column_to_rownames(var="spp")

# 2020
Unburned2020 <- as.data.frame(nulling_me_softly(Community_longForm_2020, traits=traitmatrix, plots=1:20, Nsim=999)) %>%
  mutate(year=1, severity="U") %>%
  rownames_to_column(var="sim")
Low2020 <- as.data.frame(nulling_me_softly(Community_longForm_2020, traits=traitmatrix, plots=21:40, Nsim=999)) %>%
  mutate(year=1, severity="L") %>%
  rownames_to_column(var="sim")
High2020 <- as.data.frame(nulling_me_softly(Community_longForm_2020, traits=traitmatrix, plots=41:60, Nsim=999)) %>%
  mutate(year=1, severity="H") %>%
  rownames_to_column(var="sim")
null2020 <- rbind(Unburned2020, Low2020, High2020)
write.csv(null2020, "Null Model/nullBeta2020.csv", row.names=FALSE)

# 2021
Unburned2021 <- as.data.frame(nulling_me_softly(Community_longForm_2021, traits=traitmatrix, plots=1:20, Nsim=999)) %>%
  mutate(year=2, severity="U") %>%
  rownames_to_column(var="sim")
Low2021 <- as.data.frame(nulling_me_softly(Community_longForm_2021, traits=traitmatrix, plots=21:40, Nsim=999)) %>%
  mutate(year=2, severity="L") %>%
  rownames_to_column(var="sim")
High2021 <- as.data.frame(nulling_me_softly(Community_longForm_2021, traits=traitmatrix, plots=41:60, Nsim=999)) %>%
  mutate(year=2, severity="H") %>%
  rownames_to_column(var="sim")
null2021 <- rbind(Unburned2021, Low2021, High2021)
write.csv(null2021, "Null Model/nullBeta2021.csv", row.names=FALSE)

# 2022
Unburned2022 <- as.data.frame(nulling_me_softly(Community_longForm_2022, traits=traitmatrix, plots=1:20, Nsim=999)) %>%
  mutate(year=3, severity="U") %>%
  rownames_to_column(var="sim")
Low2022 <- as.data.frame(nulling_me_softly(Community_longForm_2022, traits=traitmatrix, plots=21:40, Nsim=999)) %>%
  mutate(year=3, severity="L") %>%
  rownames_to_column(var="sim")
High2022 <- as.data.frame(nulling_me_softly(Community_longForm_2022, traits=traitmatrix, plots=41:60, Nsim=999)) %>%
  mutate(year=3, severity="H") %>%
  rownames_to_column(var="sim")
null2022 <- rbind(Unburned2022, Low2022, High2022)
write.csv(null2022, "Null Model/nullBeta2022.csv", row.names=FALSE)

# 2023
Unburned2023 <- as.data.frame(nulling_me_softly(Community_longForm_2023, traits=traitmatrix, plots=1:20, Nsim=999)) %>%
  mutate(year=4, severity="U") %>%
  rownames_to_column(var="sim")
Low2023 <- as.data.frame(nulling_me_softly(Community_longForm_2023, traits=traitmatrix, plots=21:40, Nsim=999)) %>%
  mutate(year=4, severity="L") %>%
  rownames_to_column(var="sim")
High2023 <- as.data.frame(nulling_me_softly(Community_longForm_2023, traits=traitmatrix, plots=41:60, Nsim=999)) %>%
  mutate(year=4, severity="H") %>%
  rownames_to_column(var="sim")
null2023 <- rbind(Unburned2023, Low2023, High2023)
write.csv(null2023, "Null Model/nullBeta2023.csv", row.names=FALSE)

# Bind them

nullBeta <- rbind(null2020, null2021, null2022, null2023)
write.csv(nullBeta, "Null Model/nullBeta_allyears.csv", row.names=FALSE)
