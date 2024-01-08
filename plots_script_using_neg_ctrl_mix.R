        # ###Graphs ####

# Libraries ####
library(tidyverse)
library(survival)
library(survminer)
library(ggpubr)
library(gridExtra)
library(vegan)
library(gt)
library(dunn.test)
library(pracma)

# Proportion of NAs (survival plots) ####
surv <- read.csv("D:\\new_computer\\documents\\grad_school\\Labwork\\DI_Project\\TBSV\\qPCR_data\\CSV_sheets\\For_R\\new_analysis\\using_ct_of_neg_ctrls\\prop_na_surv.csv")
#Survival plot. In this, presence of virus = 'alive, [1]' recovery = 'dead [0]'
fit <- survfit(Surv(timepoint, status) ~ exp_group, data=surv)
ggsurvplot(fit, 
           pval = TRUE, 
           conf.int = TRUE, 
           risk.table = TRUE, 
           risk.table.col = "strata", 
           linetype = 2,
           palette = c("purple","tomato","orange")) #save as 786x690


fit2 <- survfit(Surv(timepoint, status) ~ target, data=surv)
ggsurvplot(fit2, 
           pval = TRUE, 
           conf.int = TRUE, 
           risk.table = TRUE, 
           risk.table.col = "strata", 
           linetype = 2,
           palette = c("tomato","orange","purple")) #save as 786x690


surv.diff <- survdiff(Surv(timepoint, status) ~ exp_group, data=surv)
surv.diff #Looks like there's no difference between the experimental groups, but it's borderline
pairwise_result <- pairwise_survdiff(Surv(timepoint, status) ~ exp_group, data=surv, p.adjust.method = "bonferroni")
pairwise_result #Yeah, no diff

surv.diff2 <- survdiff(Surv(timepoint, status) ~ target, data=surv)
surv.diff2 #Looks like there's differences between the targets
pairwise_result2 <- pairwise_survdiff(Surv(timepoint, status) ~ target, data=surv, p.adjust.method = "bonferroni")
pairwise_result2 #Diff bt DI72 and DI73 and DI73 and TBSV



# dCt Ratios ####
###
###
###
# Import data
ratios <- read.csv("D:\\new_computer\\documents\\grad_school\\Labwork\\DI_Project\\TBSV\\qPCR_data\\CSV_sheets\\For_R\\new_analysis\\using_ct_of_neg_ctrls\\dct_ratios.csv")
all.dct <- read.csv("D:\\new_computer\\documents\\grad_school\\Labwork\\DI_Project\\TBSV\\qPCR_data\\CSV_sheets\\For_R\\new_analysis\\using_ct_of_neg_ctrls\\all_dcts_w_na.csv")

# Set variables to correct things (e.g., factor or numeric)
# Specify the columns to be converted to factors
factor.columns <- c("timepoint", "plant", "date", "target", "exp_group")

# Specify the column to be converted to numeric
numeric.column <- "dct_ratio"
numeric.column2 <- "dct"

# Convert columns to factors and numeric
ratios.f <- ratios %>%
  mutate_at(vars(factor.columns), as.factor) %>%
  mutate_at(vars(numeric.column), as.numeric)
dct.f <- all.dct%>%
  mutate_at(vars(factor.columns), as.factor) %>%
  mutate_at(vars(numeric.column2), as.numeric)

# remove all NAs from the ratios dataframe
ratios.non.na.values <- ratios.f %>% filter(!is.na(dct_ratio))

# remove all non NAs from the dct dataframe
dct.na.values <- dct.f %>% filter(is.na(dct))

# Set target colors so they're consistent
target.colors <- c("TBSV" = "purple", "DI72" = "tomato", "DI73" = "orange")
expgrp.colors <- c("TBSV" = "purple", "TBSV_DI72" = "tomato", "TBSV_DI73" = "orange")

# Exploratory look at the means + sd of all samples
means <- ratios.non.na.values %>%
  group_by(exp_group, target, timepoint) %>%
  summarise(n=n(), # n=n() means all samples
            mean = mean(dct_ratio),
            sd = sd(dct_ratio))
print(means)

# Create ratios plot with boxplots and outliers visualized
plot2 <- ggplot(ratios.non.na.values, aes(x = timepoint, y = dct_ratio, color = target)) + 
  geom_boxplot(outlier.shape = 24) +
  geom_vline(xintercept = c(4.6), color = "darkgrey", linetype = "dashed") +
  facet_grid(. ~ exp_group) +
  labs(title = "dCt ratios by Target, Experimental Group, and Timepoint",
       subtitle = "Ratio = DI:TBSV",
       x = "Timepoint",
       y = "dCt ratios",
       color = "Target") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.subtitle = element_text(hjust = 0.5)) +
  scale_x_discrete(name = "Timepoint", labels = c("6dpi", "10dpi", "14dpi")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(strip.text.x = element_text(size = 10, face = "bold")) +
  scale_color_manual(values = target.colors)
plot2

# remove some outliers and re-run plot2
ratios.non.na.values <- ratios.non.na.values[-1,]
ratios.non.na.values <- ratios.non.na.values[-1,]
ratios.non.na.values <- ratios.non.na.values[-1,]

# Create plot with means and sd
plot3 <- ggplot(means, aes(x = timepoint, y = mean, color = target)) + 
  geom_point(position = position_dodge(width = 0.8)) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),
  position = position_dodge(width = 0.8),
  width = 0.25) +
  geom_vline(xintercept = c(4.6), color = "darkgrey", linetype = "dashed") +
  facet_grid(. ~ exp_group) +
  labs(title = "Mean dCt ratios by Target, Experimental Group, and Timepoint",
       subtitle = "Agroinfiltration (Ratio=DI:TBSV)",
       x = "Timepoint",
       y = "Mean dCt ratios",
       color = "Target") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.subtitle = element_text(hjust = 0.5)) +
  scale_x_discrete(name = "Timepoint", labels = c("6dpi", "10dpi", "14dpi")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(strip.text.x = element_text(size = 10, face = "bold")) +
  scale_color_manual(values = target.colors)
plot3

#Get plots together
grid.arrange(plot2,plot3,nrow=2) #save as width 600 height 690

#Make another plot for number of na values
prop.na <- read.csv("D:\\new_computer\\documents\\grad_school\\Labwork\\DI_Project\\TBSV\\qPCR_data\\CSV_sheets\\For_R\\new_analysis\\using_ct_of_neg_ctrls\\prop_na.csv")
factor.columns2 <- c("exp_group","target","timepoint")
prop.na.f <- prop.na %>%
  mutate_at(vars(factor.columns2), as.factor)

plot4 <- ggplot(prop.na.f, aes(x = timepoint, y = prop.na, color = target)) + 
  geom_point() +
  geom_line(aes(group = target)) +
  geom_vline(xintercept = c(4.6), color = "darkgrey", linetype = "dashed") +
  facet_grid(. ~ exp_group) +
  labs(title = "Proportion of negative plants by Target, Experimental Group, and Timepoint",
       subtitle = "",
       x = "Timepoint",
       y = "Proportion") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.subtitle = element_text(hjust = 0.5)) +
  scale_x_discrete(name = "Timepoint", labels = c("6dpi", "10dpi", "14dpi")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(strip.text.x = element_text(size = 10, face = "bold")) +
  scale_color_manual(values = target.colors)
plot4

#Get plots together
grid.arrange(plot2,plot3,nrow=2) #save as width 600 height 690


#Get plots together
grid.arrange(plot2,plot4,nrow=2) #save as width 600 height 690


# Compare mechanical inoculation with agroinfiltration ratio values ####
###
###
###
# Import data
ratios.comp <- read.csv("D:\\new_computer\\documents\\grad_school\\Labwork\\DI_Project\\TBSV\\qPCR_data\\CSV_sheets\\For_R\\new_analysis\\using_ct_of_neg_ctrls\\ratios_agro_v_leafsoup.csv")

# Filter data for this one
ratios.comp.factored <- ratios.comp %>%
  mutate_at(vars(timepoint,plant,target,exp_group),
            list(factor))
numeric.column3 <- "dct_ratio"
ratios.comp.factored <- ratios.comp.factored %>%
  mutate_at(vars(numeric.column3), as.numeric)
ratios.comp.noctrls <- ratios.comp.factored %>%
  filter(!exp_group %in% c("neg_ctrl", "pos_ctrl", "neg_plt", "rcvrd_leaf"))

# Create mechanical inoculation only dataset for graphing and remove all NAs

leaf.soup <- ratios.comp.noctrls %>%
  filter(!target %in% c("DI73"))
#remove the DI73 EXPERIMENTAL GROUP since there's no equiv in the mechanical plants
leaf.soup <- leaf.soup %>%
  filter(!exp_group %in% c("TBSV_DI73"))
# remove all NAs from the ratios dataframe
leaf.soup <- leaf.soup %>% filter(!is.na(dct_ratio))

# Create mechanical inoculation only summary stats

means2 <- leaf.soup %>%
  group_by(exp_group, method, timepoint) %>%
  summarise(n=n(), # n=n() means all samples
            mean = mean(dct_ratio),
            sd = sd(dct_ratio))
print(means2)

# Create means + sd plot for leaf soup

plot4 <- ggplot(means2, aes(x = timepoint, y = mean, color = exp_group)) + 
  geom_point(position = position_dodge(width = 0.8)) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),
                position = position_dodge(width = 0.8),
                width = 0.25) +
  geom_vline(xintercept = c(4.6), color = "darkgrey", linetype = "dashed") +
  facet_grid(. ~ method) +
  labs(title = "Mean dCt ratios by Method, Experimental Group, and Timepoint",
       subtitle = "Mechanical inoculation (ratio=DI:TBSV)",
       x = "Timepoint",
       y = "Mean dCt ratios",
       color = "Experimental Group") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.subtitle = element_text(hjust = 0.5)) +
  scale_x_discrete(name = "Timepoint", labels = c("6dpi", "10dpi", "14dpi")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(strip.text.x = element_text(size = 10, face = "bold")) +
  scale_color_manual(values = expgrp.colors)
plot4

#Combine method plots (make sure you adjust the subtitle on plot3)
grid.arrange(plot3,plot4,nrow=2) #save as width 600 height 690


# Looking at just TBSV ddCt ####
TBSV <- read.csv("D:\\new_computer\\documents\\grad_school\\Labwork\\DI_Project\\TBSV\\qPCR_data\\CSV_sheets\\For_R\\new_analysis\\using_ct_of_neg_ctrls\\tbsv_ddct.csv")
TBSV.factored <- TBSV %>%
  mutate_at(vars(timepoint,plant,date,target,exp_group),
            list(factor)) # factor appropriate columns

ggplot(TBSV.factored, aes(x = timepoint, y = log_2pwrneg_ddct, color = "TBSV")) + 
  geom_boxplot(outlier.shape = 24) +
  geom_vline(xintercept = c(4.6), color = "darkgrey", linetype = "dashed") +
  facet_grid(. ~ exp_group) +
  labs(title = "Log(2^-ddCt) Values for TBSV in all three experimental groups",
       subtitle = "Calibrator = Average TBSV dCt at 6dpi",
       x = "Timepoint",
       y = "Log(2^-ddCt) Values",
       color = "Target") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.subtitle = element_text(hjust = 0.5)) +
  scale_x_discrete(name = "Timepoint", labels = c("6dpi", "10dpi", "14dpi")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(strip.text.x = element_text(size = 10, face = "bold")) +
  scale_color_manual(values = "purple")

# AUDPC Plot ####
sev <- read.csv("D:\\new_computer\\documents\\grad_school\\Labwork\\DI_Project\\TBSV\\qPCR_data\\CSV_sheets\\For_R\\new_analysis\\using_ct_of_neg_ctrls\\disease_sev_for_R.csv")
#subsetting
severityTBSV <- sev %>%
  filter(!exp_group %in% c("TBSV_DI72", "TBSV_DI73", "neg_plt"))
severityDI72 <- sev %>%
  filter(!exp_group %in% c("TBSV", "TBSV_DI73", "neg_plt"))
severityDI73 <- sev %>%
  filter(!exp_group %in% c("TBSV_DI72", "TBSV", "neg_plt"))
severityneg <- sev %>%
  filter(!exp_group %in% c("TBSV_DI72", "TBSV_DI73", "TBSV"))
daysAfterInoculation<-c(0,6,10,14)

#Get summary stats (needed for this kind of table)
mean.tbsv0 <- 0
sd.tbsv0 <- 0
mean.tbsv1 <- mean(severityTBSV$mean.tp1)
sd.tbsv1 <- sd(severityTBSV$mean.tp1)
mean.tbsv2 <- mean(severityTBSV$mean.tp2)
sd.tbsv2 <- sd(severityTBSV$mean.tp2)
mean.tbsv3 <- mean(severityTBSV$mean.tp3)
sd.tbsv3 <- sd(severityTBSV$mean.tp3)

mean.DI720 <- 0
sd.DI720 <- 0
mean.DI721 <- mean(severityDI72$mean.tp1)
sd.DI721 <- sd(severityDI72$mean.tp1)
mean.DI722 <- mean(severityDI72$mean.tp2)
sd.DI722 <- sd(severityDI72$mean.tp2)
mean.DI723 <- mean(severityDI72$mean.tp3)
sd.DI723 <- sd(severityDI72$mean.tp0)

mean.DI730 <- 0
sd.DI730 <- 0
mean.DI731 <- mean(severityDI73$mean.tp1)
sd.DI731 <- sd(severityDI73$mean.tp1)
mean.DI732 <- mean(severityDI73$mean.tp2)
sd.DI732 <- sd(severityDI73$mean.tp2)
mean.DI733 <- mean(severityDI73$mean.tp3)
sd.DI733 <- sd(severityDI73$mean.tp3)

mean.neg0 <- 0
sd.neg0 <- 0
mean.neg1 <- 0
sd.neg1 <- 0
mean.neg2 <- 0
sd.neg2 <- 0
mean.neg3 <- 0
sd.neg3 <- 0

severityTBSV2 <- c(mean.tbsv0,mean.tbsv1,mean.tbsv2,mean.tbsv3)
severityDI72.2 <- c(mean.DI720,mean.DI721,mean.DI722,mean.DI723)
severityDI73.2 <- c(mean.DI730,mean.DI731,mean.DI732,mean.DI733)
severityneg2 <- c(mean.neg0,mean.neg1,mean.neg2,mean.neg3)
sd.TBSV.sev <- c(sd.tbsv0,sd.tbsv1,sd.tbsv2,sd.tbsv3)
sd.DI72.sev <- c(sd.DI720,sd.DI721,sd.DI722,sd.DI723)
sd.DI73.sev <- c(sd.DI730,sd.DI731,sd.DI732,sd.DI733)
sd.neg.sev <- c(sd.neg0,sd.neg1,sd.neg2,sd.neg3)

# Combine the data into a single data frame
df_combined <- rbind(
  data.frame(daysAfterInoculation, severity = severityTBSV2, treatment = "TBSV", sd_sev = sd.TBSV.sev),
  data.frame(daysAfterInoculation, severity = severityDI72.2, treatment = "TBSV+DI72", sd_sev = sd.DI72.sev),
  data.frame(daysAfterInoculation, severity = severityDI73.2, treatment = "TBSV+DI73", sd_sev = sd.DI73.sev),
  data.frame(daysAfterInoculation, severity = severityneg2, treatment = "Control plant", sd_sev = sd.DI73.sev)
)

# Create plot
ggplot(df_combined, aes(x = daysAfterInoculation, y = severity, color = treatment, group = treatment)) +
  geom_line() +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = severity - sd_sev, ymax = severity + sd_sev), width = 0.1) +
  ylim(0, 100) +
  labs(x = 'Days After Inoculation', y = '% Symptom Coverage', title = 'TBSV Disease Progress') +
  scale_color_manual(values = c("cornflowerblue", "purple", "tomato", "orange")) +
  annotate(
    "text",
    x = 5, y = max(df_combined$severity) + 5,
    label = paste("Kruskal-Wallis p =", format.pval(kw.test.AUDPC$p.value, digits = 3)),
    size = 3) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))

 





        # ###Stats ####

# ANOVA won't work for ratios or TBSV 2pwrneg_ddct, but will for prop.na ####
aov <- aov(dct_ratio ~ timepoint, data = ratios.non.na.values)
summary(aov)
hist(resid(aov))
shapiro.test(resid(aov)) #failure
bartlett.test(dct_ratio ~ timepoint, data = ratios.non.na.values) #failure

aov2 <- aov(dct_ratio ~ exp_group, data = ratios.non.na.values)
summary(aov2)
hist(resid(aov2))
shapiro.test(resid(aov2)) #failure
bartlett.test(dct_ratio ~ exp_group, data = ratios.non.na.values)#failure
#cannot log transform ratios due to negative values. Those neg values are useful information.

aov3 <- aov(log_2pwrneg_ddct ~ timepoint, data = TBSV.factored)
summary(aov3)
hist(resid(aov3))
shapiro.test(resid(aov3)) #failure
bartlett.test(log_2pwrneg_ddct ~ timepoint, data = TBSV.factored) # not a failure!

aov4 <- aov(prop.na ~ timepoint*exp_group, data = prop.na.f)
summary(aov4)
hist(resid(aov4))
shapiro.test(resid(aov4)) #failure

aov5 <- aov(prop.na ~ timepoint*target, data = prop.na.f)
summary(aov5)
hist(resid(aov5))
shapiro.test(resid(aov5)) #not a failure

aov6 <- aov(prop.na ~ target*timepoint, data = prop.na.f)
summary(aov6)
hist(resid(aov6))
shapiro.test(resid(aov6)) #failure, but barely

# Kruskal-Wallis test on ratios, TP by TP, first exp_group, then target ####

#STEP 1: experimental group
timepoints <- unique(ratios.non.na.values$timepoint) # List of timepoints

# Loop through each timepoint
for (timepoint in timepoints) {
  cat("Timepoint:", timepoint, "\n")
  
  # Subset data for the current timepoint
  subset.data <- ratios.non.na.values[ratios.non.na.values$timepoint == timepoint, ]
  
  # Kruskal-Wallis test
  kruskal.result <- kruskal.test(dct_ratio ~ exp_group, data = subset.data)
  print(kruskal.result)
  
  # Dunn's test for post hoc comparisons
  dunn.result <- dunn.test(subset.data$dct_ratio, g = subset.data$exp_group, method = "bonferroni")
  print(dunn.result)
}
# TP1 = not sig (p=0.1055)
# TP2 = not sig (p=0.1013)
# TP3 = not sig (p=0.4073)

# STEP 2: target
timepoints <- unique(ratios.non.na.values$timepoint) # List of timepoints with the new, outlier-less dataset

# Loop through each timepoint with the new dataset
for (timepoint in timepoints) {
  cat("Timepoint:", timepoint, "\n")
  
  # Subset data for the current timepoint
  subset.data <- ratios.non.na.values[ratios.non.na.values$timepoint == timepoint, ]
  
  # Kruskal-Wallis test
  kruskal.result <- kruskal.test(dct_ratio ~ target, data = subset.data)
  print(kruskal.result)
  
  # Dunn's test for post hoc comparisons
  dunn.result <- dunn.test(subset.data$dct_ratio, g = subset.data$target, method = "bonferroni")
  print(dunn.result)
}
# TP1 = not sig (p=0.1352)
# TP2 = sig (p=0.01625)
# TP3 = not sig (p=0.09348)

# STEP 3: Need to figure out which experimental group has the significance at TP2
TBSVgrp <- ratios.non.na.values %>%
  filter(!exp_group %in% c("TBSV_DI72", "TBSV_DI73"))
DI72grp <- ratios.non.na.values %>%
  filter(!exp_group %in% c("TBSV", "TBSV_DI73"))
DI73grp <- ratios.non.na.values %>%
  filter(!exp_group %in% c("TBSV_DI72", "TBSV"))
View(TBSVgrp) #worked

# Loopin' time for each. First up: TBSV
timepoints <- unique(TBSVgrp$timepoint) # List of timepoints with the new, outlier-less dataset

# Loop through each timepoint with the new dataset
for (timepoint in timepoints) {
  cat("Timepoint:", timepoint, "\n")
  
  # Subset data for the current timepoint
  subset.data <- TBSVgrp[TBSVgrp$timepoint == timepoint, ]
  
  # Kruskal-Wallis test
  kruskal.result <- kruskal.test(dct_ratio ~ target, data = subset.data)
  print(kruskal.result)
  
  # Dunn's test for post hoc comparisons
  dunn.result <- dunn.test(subset.data$dct_ratio, g = subset.data$target, method = "bonferroni")
  print(dunn.result)
}
#TP2 has p=0.02935

# Next is DI72 group
timepoints <- unique(DI72grp$timepoint) # List of timepoints with the new, outlier-less dataset

# Loop through each timepoint with the new dataset
for (timepoint in timepoints) {
  cat("Timepoint:", timepoint, "\n")
  
  # Subset data for the current timepoint
  subset.data <- DI72grp[DI72grp$timepoint == timepoint, ]
  
  # Kruskal-Wallis test
  kruskal.result <- kruskal.test(dct_ratio ~ target, data = subset.data)
  print(kruskal.result)
  
  # Dunn's test for post hoc comparisons
  dunn.result <- dunn.test(subset.data$dct_ratio, g = subset.data$target, method = "bonferroni")
  print(dunn.result)
}
# Not this group; no significance

# Next is DI73
timepoints <- unique(DI73grp$timepoint) # List of timepoints with the new, outlier-less dataset

# Loop through each timepoint with the new dataset
for (timepoint in timepoints) {
  cat("Timepoint:", timepoint, "\n")
  
  # Subset data for the current timepoint
  subset.data <- DI73grp[DI73grp$timepoint == timepoint, ]
  
  # Kruskal-Wallis test
  kruskal.result <- kruskal.test(dct_ratio ~ target, data = subset.data)
  print(kruskal.result)
  
  # Dunn's test for post hoc comparisons
  dunn.result <- dunn.test(subset.data$dct_ratio, g = subset.data$target, method = "bonferroni")
  print(dunn.result)
}
# Not this one either
# Looks like it's the TBSV group that has the significance


# 2-way ANOVA to look at differences in prop.na ####
aov6 <- aov(prop.na ~ target*exp_group, data = prop.na.f)
TukeyHSD(aov6)
#significant stuff:
  # DI73:TBSV-DI72:TBSV p=0.0335909
  # TBSV:TBSV-DI73:TBSV p=0.0041098
  # DI72:TBSV_DI72-DI73:TBSV p=0.0082113
  # TBSV:TBSV_DI72-DI73:TBSV p=0.0082113
  # DI72:TBSV_DI73-DI73:TBSV p=0.0321391
  # TBSV:TBSV_DI73-DI73:TBSV p=0.0041098
  # DI73:TBSV_DI72-TBSV:TBSV p=0.0293872
  # TBSV:TBSV_DI73-DI73:TBSV_DI72 p=0.0293872
  #### Notably, none of these things are saying that target in one experimental group is different than that SAME target in another.
     # In the TBSV group, 73 is different from 72 and TBSV (0.034,0.004), but 72 and TBSV are not different (0.978)
     # In the DI72 group, 72 is not different from 72 and 72 is not different from TBSV (0.057,1.00). Apparently, 73 and TBSV aren't different either (0.057)
     # In the DI73 group, Nobody is different from each other (TB:72 0.981, TB:73 0.117, 72:73 0.531)
     # TBSV IS NOT DIFFERENT across groups
     # DI72 IS NOT DIFFERENT across groups 
     # DI73 IS NOT DIFFERENT across groups

# Kruskal-Wallis test on just TBSV log(2^-ddCt), TP by TP ####
timepoints <- unique(TBSV.factored$timepoint) # List of timepoints

# Loop through each timepoint
for (timepoint in timepoints) {
  cat("Timepoint:", timepoint, "\n")
  
  # Subset data for the current timepoint
  subset.data <- TBSV.factored[TBSV.factored$timepoint == timepoint, ]
  
  # Kruskal-Wallis test
  kruskal.result <- kruskal.test(log_2pwrneg_ddct ~ exp_group, data = subset.data)
  print(kruskal.result)
  
  # Dunn's test for post hoc comparisons
  dunn.result <- dunn.test(subset.data$log_2pwrneg_ddct, g = subset.data$exp_group, method = "bonferroni")
  print(dunn.result)
}
# TP1 = borderline 0.04874 (difference bt 72 and 73 groups)
# TP2 = not sig 0.2028
# TP3 = sig 0.006284 (diff bt the TBSV and 72 group and the 72 and 73 group)


# KW test comparing TPs within an experimental group ####
#
#
exp_group <- unique(TBSV.factored$exp_group) # List of groups

# Loop through each timepoint
for (exp_group in exp_group) {
  cat("Experimental Group:", exp_group, "\n")
  
  # Subset data for the current timepoint
  subset.data2 <- ratios.factored[ratios.factored$exp_group == exp_group, ]
  
  # Kruskal-Wallis test
  kruskal.result <- kruskal.test(dct_ratio ~ timepoint, data = subset.data2)
  print(kruskal.result)
  
  # Dunn's test for post hoc comparisons
  dunn.result <- dunn.test(subset.data2$dct_ratio, g = subset.data2$timepoint, method = "bonferroni")
  print(dunn.result)
}

# KW test comparing the two methods ####
###
###
# Since KW can be sensitive to sample size, randomly select 29 plants from the agro dataset, 8 for each timepoint, to compare to the mechanical incoculation
# Only select from TBSV or TBSV_DI72 since the mechanical data doesn't have DI73

set.seed(938)  # Set seed for reproducibility

# Create a function to randomly select rows from each timepoint for specified experimental groups
rdm.rows <- function(ratios.non.na.values, timepoint_col, exp_group_col, target_col, groups = c("TBSV", "TBSV_DI72"), target = "DI72", num_rows = 27) {
  unique_tp <- unique(ratios.non.na.values[[timepoint_col]])
  
  selected_rows <- NULL
  
  for (tp in unique_tp) {
    tp_data <- ratios.non.na.values[ratios.non.na.values[[timepoint_col]] == tp & ratios.non.na.values[[exp_group_col]] %in% groups &
                                ratios.non.na.values[[target_col]] == target, ]
    
    if (nrow(tp_data) >= num_rows) {
      selected_rows <- rbind(selected_rows, tp_data[sample(1:nrow(tp_data), num_rows), , drop = FALSE])
    } else {
      selected_rows <- rbind(selected_rows, tp_data)
    }
  }
  
  return(selected_rows)
}

# Example usage
selected_rows <- rdm.rows(ratios.non.na.values, 'timepoint', 'exp_group', 'target', groups = c('TBSV', 'TBSV_DI72'), target = "DI72", num_rows = 8)

# Export the new random data to be pooled with the leaf.soup data and brought back in ("merge" just wasn't working)
write.csv(selected_rows, "D:\\new_computer\\documents\\grad_school\\Labwork\\DI_Project\\TBSV\\qPCR_data\\CSV_sheets\\For_R\\new_analysis\\using_ct_of_neg_ctrls\\merged6.csv")
merged <- read.csv("D:\\new_computer\\documents\\grad_school\\Labwork\\DI_Project\\TBSV\\qPCR_data\\CSV_sheets\\For_R\\new_analysis\\using_ct_of_neg_ctrls\\merged.csv")
merged2 <- read.csv("D:\\new_computer\\documents\\grad_school\\Labwork\\DI_Project\\TBSV\\qPCR_data\\CSV_sheets\\For_R\\new_analysis\\using_ct_of_neg_ctrls\\merged2.csv")
merged3 <- read.csv("D:\\new_computer\\documents\\grad_school\\Labwork\\DI_Project\\TBSV\\qPCR_data\\CSV_sheets\\For_R\\new_analysis\\using_ct_of_neg_ctrls\\merged3.csv")
merged4 <- read.csv("D:\\new_computer\\documents\\grad_school\\Labwork\\DI_Project\\TBSV\\qPCR_data\\CSV_sheets\\For_R\\new_analysis\\using_ct_of_neg_ctrls\\merged4.csv")
merged5 <- read.csv("D:\\new_computer\\documents\\grad_school\\Labwork\\DI_Project\\TBSV\\qPCR_data\\CSV_sheets\\For_R\\new_analysis\\using_ct_of_neg_ctrls\\merged5.csv")
merged6 <- read.csv("D:\\new_computer\\documents\\grad_school\\Labwork\\DI_Project\\TBSV\\qPCR_data\\CSV_sheets\\For_R\\new_analysis\\using_ct_of_neg_ctrls\\merged6.csv")

#RE-DO THE ABOVE FIVE TIMES, JUST TO SEE IF THE SIGNIFICANCE OF THE RESULT IS THE SAME
#SEEDS: merged = 72, merged2 = 545, merged3 = 888, merged4 = 326, merged5 = 200, merged6 = 938

# KW test
kruskal.test(dct_ratio ~ method, data = merged) #not significant (p=0.5211)
kruskal.test(dct_ratio ~ method, data = merged2) #not sig (p=0.7341)
kruskal.test(dct_ratio ~ method, data = merged3) #not sig (0.2575)
kruskal.test(dct_ratio ~ method, data = merged4) #not sig (0.7341)
kruskal.test(dct_ratio ~ method, data = merged5) #not sig (0.7484)
kruskal.test(dct_ratio ~ method, data = merged6) #not sig (0.88)

# AUDPC calculation and stats ####

audpc_tbsv <- trapz(daysAfterInoculation, severityTBSV2)
audpc_di72 <- trapz(daysAfterInoculation, severityDI72.2)
audpc_di73 <- trapz(daysAfterInoculation, severityDI73.2)
audpc_control <- trapz(daysAfterInoculation, severityneg2)

# Create dataframe of experimental group and AUDPC + perform ANOVA
df <- data.frame(
  Treatment = c(rep("TBSV", length(audpc_tbsv)),
                rep("DI72", length(audpc_di72)),
                rep("DI73", length(audpc_di73)),
                rep("Control", length(audpc_control))),
  AUDPC = c(audpc_tbsv, audpc_di72, audpc_di73, audpc_control)
)
# Perform Kruskal-Wallis test
kw.test.AUDPC <- kruskal.test(AUDPC ~ Treatment, data = df) #no sig diff
kw.test.AUDPC



# Timepoint-by-timepoint experimental group comparison + Dunn's test ####
timepoints <- unique(ratios.noctrls$timepoint) # List of timepoints

# Loop through each timepoint
for (timepoint in timepoints) {
  cat("Timepoint:", timepoint, "\n")
  
  # Subset data for the current timepoint
  subset.data <- ratios.noctrls[ratios.noctrls$timepoint == timepoint, ]
  
  # Kruskal-Wallis test
  kruskal.result <- kruskal.test(dct_ratio ~ exp_group, data = subset.data)
  print(kruskal.result)
  
  # Dunn's test for post hoc comparisons
  dunn.result <- dunn.test(subset.data$dct_ratio, g = subset.data$exp_group, method = "bonferroni")
  print(dunn.result)
}