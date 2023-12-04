# Libraries ####
library(tidyverse)
library(survival)
library(survminer)
library(ggpubr)
library(gridExtra)
library(vegan)


        # ###Graphs ####

# Proportion of NAs (survival plots) ####

#Import data + create the survival object looking at experimental group as a variable. In this, presence of virus = 'alive,' recovery = 'dead'
surv <- read.csv("https://raw.githubusercontent.com/skybird99-anthony-taylor/TBSV_DIs/main/prop_na_surv.csv?token=GHSAT0AAAAAACLE452YVJXMEFKUZX3FVVOCZLNXRYQ")
fit <- survfit(Surv(timepoint, status) ~ exp_group, data=surv)
summary(fit)

#Create the graph to go with the above survival object
ggsurvplot(fit, 
           pval = TRUE, 
           conf.int = TRUE, 
           risk.table = TRUE, 
           risk.table.col = "strata", 
           linetype = 2,
           palette = c("tomato","orange","purple"))

#Create the survival object looking at target as a variable. Like the above, presence of virus = 'alive,' recovery = 'dead'
fit2 <- survfit(Surv(timepoint, status) ~ target, data=surv)
summary(fit)

#Graph this second survival object
ggsurvplot(fit2, 
           pval = TRUE, 
           conf.int = TRUE, 
           risk.table = TRUE, 
           risk.table.col = "strata", 
           linetype = 2,
           palette = c("tomato","orange","purple"))

#Time for the stats! Do experimental group first
surv.diff <- survdiff(Surv(timepoint, status) ~ exp_group, data=surv)
surv.diff #Looks like there's differences between the experimental group! Let's do a post-hoc, then
pairwise.res <- pairwise_survdiff(Surv(timepoint, status) ~ exp_group, data=surv, p.adjust.method = "bonferroni")
pairwise.res

surv.diff2 <- survdiff(Surv(timepoint, status) ~ target, data=surv)
surv.diff2 #Looks like there's differences between the targets, too, so, post-hoc time
pairwise.res2 <- pairwise_survdiff(Surv(timepoint, status) ~ target, data=surv, p.adjust.method = "bonferroni")
pairwise.res2



# Ct Ratios ####
###
###
###
# Import data
ratios <- read.csv("https://raw.githubusercontent.com/skybird99-anthony-taylor/TBSV_DIs/main/ratios_both_against_self.csv?token=GHSAT0AAAAAACLE452Y32P2VK3QMXWIVE4IZLNXYFA")
# Set variables to factors and filter out controls for graphing
str(ratios)
ratios.factored <- ratios %>%
  mutate_at(vars(timepoint,plant,date,target,exp_group),
            list(factor))
ratios.noctrls <- ratios.factored %>%
  filter(!exp_group %in% c("neg_ctrl", "pos_ctrl", "neg_plt", "rcvrd_leaf"))
str(ratios.noctrls)

# Set target colors so they're consistent
target.colors <- c("TBSV" = "purple", "DI72" = "tomato", "DI73" = "orange")

# Exploratory look at the means + sd of all samples (we'll need these for a later graph, too)
means <- ratios.noctrls %>%
  group_by(exp_group, target, timepoint) %>%
  summarise(n=n(), # n=n() means all samples
            mean = mean(dct_ratio),
            sd = sd(dct_ratio))
print(means)

# Create ratios plot with boxplots and outliers visualized
plot2 <- ggplot(ratios.noctrls, aes(x = timepoint, y = dct_ratio, color = target)) + 
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

# Create plot with means and sd
plot3 <- ggplot(means, aes(x = timepoint, y = mean, color = target)) + 
  geom_point(position = position_dodge(width = 0.8)) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),
  position = position_dodge(width = 0.8),
  width = 0.25) +
  geom_vline(xintercept = c(4.6), color = "darkgrey", linetype = "dashed") +
  facet_grid(. ~ exp_group) +
  labs(title = "Mean dCt ratios by Target, Experimental Group, and Timepoint",
       subtitle = "Agroinfiltration (Ratio=DI:TBSV",
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

#remove the one super outlier of -607 and redo plot2
ratios.noctrls <- ratios.noctrls[-1,]
        
#Get plots together
grid.arrange(plot2,plot3,nrow=2) #save as width 600 height 690


# Compare mechanical inoculation with agroinfiltration ratio values
###
###
###
# Import data
ratios.comp <- read.csv("https://raw.githubusercontent.com/skybird99-anthony-taylor/TBSV_DIs/main/against_self_w_leaf_soup_ratios.csv?token=GHSAT0AAAAAACLE452Y7IWSHL7VFGGZFJQKZLNX6EQ")

# Filter data for this one
ratios.comp.factored <- ratios.comp %>%
  mutate_at(vars(timepoint,plant,date,target,exp_group),
            list(factor))
ratios.comp.noctrls <- ratios.comp.factored %>%
  filter(!exp_group %in% c("neg_ctrl", "pos_ctrl", "neg_plt", "rcvrd_leaf"))

# Create mechanical inoculation only dataset for graphing

leaf.soup <- ratios.comp.noctrls %>%
  filter(!method %in% c("agro"))

# Create mechanical inoculation only summary stats

means2 <- leaf.soup %>%
  group_by(exp_group, target, timepoint) %>%
  summarise(n=n(), # n=n() means all samples
            mean = mean(dct_ratio),
            sd = sd(dct_ratio))
print(means2)

# Create means + sd plot for leaf soup. This was not used in the publication, 
#but is useful for visualization and the data from it is used in the Kruskal-Wallis tests 
#that ARE mentioned in the publication

plot4 <- ggplot(means2, aes(x = timepoint, y = mean, color = target)) + 
  geom_point(position = position_dodge(width = 0.8)) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),
                position = position_dodge(width = 0.8),
                width = 0.25) +
  geom_vline(xintercept = c(4.6), color = "darkgrey", linetype = "dashed") +
  facet_grid(. ~ exp_group) +
  labs(title = "Mean dCt ratios by Target, Experimental Group, and Timepoint",
       subtitle = "Mechanical inoculation (ratio=DI:TBSV)",
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
plot4

#Combine method plots (make sure you adjust the subtitle on plot3)
grid.arrange(plot3,plot4,nrow=2) #save as width 600 height 690


# AUDPC Plot ####
###
###
###
#Import data
sev <- read.csv("https://raw.githubusercontent.com/skybird99-anthony-taylor/TBSV_DIs/main/disease_sev_for_R.csv?token=GHSAT0AAAAAACLE452YJ4G7IEA4JAZVBFU4ZLNYFPA")
#To make this work, we have to subset these data
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
  scale_color_manual(values = c("purple", "tomato", "orange", "cornflowerblue")) +
  annotate(
    "text",
    x = 5, y = max(df_combined$severity) + 5,
    label = paste("Kruskal-Wallis p =", format.pval(kw.test.AUDPC$p.value, digits = 3)),
    size = 3) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))

 


# Ratio/Phenotype overlay (not used in the publication) ####
###
###
###
#Import data
pheno <- read.csv("https://raw.githubusercontent.com/skybird99-anthony-taylor/TBSV_DIs/main/ratios_w_phenos.csv?token=GHSAT0AAAAAACLE452ZYLPNFYC3ZUOLSXLIZLNYHQA")

#Create summary stats
means3 <- pheno %>%
  group_by(exp_group, target, timepoint) %>%
  summarise(n=n(), # n=n() means all samples
            mean = mean(phenotype),
            sd = sd(phenotype))
print(means3)

#plot summary stats
ggplot(means3, aes(x = timepoint, y = mean, color = target)) + 
  geom_point(position = position_dodge(width = 0.8)) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),
                position = position_dodge(width = 0.8),
                width = 0.25) +
  geom_vline(xintercept = c(4.6), color = "darkgrey", linetype = "dashed") +
  facet_grid(. ~ exp_group) +
  labs(title = "Mean Disease Severity Score by Target, Experimental Group, and Timepoint",
       subtitle = "Ratio=DI:TBSV",
       x = "Timepoint",
       y = "Mean Disease Severity Score",
       color = "Target") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.subtitle = element_text(hjust = 0.5)) +
  scale_x_discrete(name = "Timepoint", labels = c("6dpi", "10dpi", "14dpi")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(strip.text.x = element_text(size = 10, face = "bold")) +
  scale_color_manual(values = target.colors)
       

         

        # ###Stats ####

# Kruskal-Wallis tests on ratios ####
###
###
# KW test one (timepoint and experimental group, alone)
kruskal.test(dct_ratio ~ timepoint, data = ratios.noctrls) #not sig (p=0.1245)
kruskal.test(dct_ratio ~ exp_group, data = ratios.noctrls) #sig (p=0.00045)
kruskal.test(dct_ratio ~ target, data = ratios.noctrls) #not sig (p = 0.3083)

#post-hoc tests
posthoc.res <- pairwise.wilcox.test(ratios.noctrls$dct_ratio, ratios.noctrls$exp_group, p.adj = "bonferroni")
print(posthoc.res)
posthoc.res2 <- pairwise.wilcox.test(ratios.noctrls$dct_ratio, ratios.noctrls$target, p.adj = "bonferroni")
print(posthoc.res2)

# KW test comparing the two methods ####
###
###
# Since KW can be sensitive to sample size, randomly select 29 plants from the agro dataset, 8 for each timepoint, to compare to the mechanical incoculation
# Only select from TBSV or TBSV_DI72 since the mechanical data doesn't have DI73

set.seed(938)  # Set seed for reproducibility (this seed will change for each of the 5 reps: merged = 72, merged2 = 545, merged3 = 888, merged4 = 326, merged5 = 200, merged6 = 938)

# Create a function to randomly select rows from each timepoint for specified experimental groups
rdm.rows <- function(ratios.noctrls, timepoint_col, exp_group_col, target_col, groups = c("TBSV", "TBSV_DI72"), target = "DI72", num_rows = 29) {
  unique_tp <- unique(ratios.noctrls[[timepoint_col]])
  
  selected_rows <- NULL
  
  for (tp in unique_tp) {
    tp_data <- ratios.noctrls[ratios.noctrls[[timepoint_col]] == tp & ratios.noctrls[[exp_group_col]] %in% groups &
                                ratios.noctrls[[target_col]] == target, ]
    
    if (nrow(tp_data) >= num_rows) {
      selected_rows <- rbind(selected_rows, tp_data[sample(1:nrow(tp_data), num_rows), , drop = FALSE])
    } else {
      selected_rows <- rbind(selected_rows, tp_data)
    }
  }
  
  return(selected_rows)
}

# Example usage
selected_rows <- rdm.rows(ratios.noctrls, 'timepoint', 'exp_group', 'target', groups = c('TBSV', 'TBSV_DI72'), target = "DI72", num_rows = 8)

# Export the new random data to be pooled with the leaf.soup data and brought back in ("merge" just wasn't working)
write.csv(selected_rows, "D:\\new_computer\\data_drive\\R_stuff\\TBSV_DIs_stuff\\random_rows.csv") #change this to the correct file path on your own computer
#Import merged data: remember to change 'merged' to each name/file path every time you repeat
merged <- read.csv("https://raw.githubusercontent.com/skybird99-anthony-taylor/TBSV_DIs/main/merged.csv?token=GHSAT0AAAAAACLE452YWGKRAOUUW2GQPFKMZLNYNIA")

# Filter out negative plants (remember to change 'merged' to 'merged2' etc. when needed)
merged6 <- merged6 %>%
  filter(!exp_group %in% c("neg_plt"))

#RE-DO THE ABOVE FIVE TIMES, JUST TO SEE IF THE SIGNIFICANCE OF THE RESULT IS THE SAME
#SEEDS: merged = 72, merged2 = 545, merged3 = 888, merged4 = 326, merged5 = 200, merged6 = 938

# KW test
kruskal.test(dct_ratio ~ method, data = merged) #not significant (p=0.07508)
kruskal.test(dct_ratio ~ method, data = merged2) #not sig (p=0.06031)
kruskal.test(dct_ratio ~ method, data = merged3) #not sig (0.6892)
kruskal.test(dct_ratio ~ method, data = merged4) #not sig (0.2713)
kruskal.test(dct_ratio ~ method, data = merged5) #not sig (0.1738)
kruskal.test(dct_ratio ~ method, data = merged6) #not sig (0.1499)

# AUDPC calculation and stats ####

audpc_tbsv <- trapz(daysAfterInoculation, severityTBSV2)
audpc_di72 <- trapz(daysAfterInoculation, severityDI72.2)
audpc_di73 <- trapz(daysAfterInoculation, severityDI73.2)
audpc_control <- trapz(daysAfterInoculation, severityneg2)

# Create dataframe of experimental group and AUDPC
df <- data.frame(
  Treatment = c(rep("TBSV", length(audpc_tbsv)),
                rep("DI72", length(audpc_di72)),
                rep("DI73", length(audpc_di73)),
                rep("Control", length(audpc_control))),
  AUDPC = c(audpc_tbsv, audpc_di72, audpc_di73, audpc_control)
)
# Perform Kruskal-Wallis test
kw.test.AUDPC <- kruskal.test(AUDPC ~ Treatment, data = df)
kw.test.AUDPC #no sig


