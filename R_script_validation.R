### Sehovic et al. 2024 - R script

### Analysis of 7-miRNA signature and
### model application from Chiorino et al. 2023; Scientific Reports

### Data used in this script can be accessed through Zenodo:
### https://doi.org/10.5281/zenodo.11234225

### Analysis was performed on R version 4.4.0 with all packages being updated 
### by 22nd May 2024

#### Loading packages ####

# In order to run the Bayesian analyses using the brms package, 
# make sure you have the "Rstan" package and its dependencies 
# installed as well as RTools 4.4

library(brms)
library(dplyr)
library(glmnet)
library(tidybayes)
library(metamisc)
library(stringr)
library(boot)
library(DescTools)
library(readxl)
library(pROC)
library(ggplot2)
library(CalibrationCurves)
library(ROCR)
library(readxl)
library(reshape2)
library(ggpubr)

#### Functions ####

### Model updating method selection - closed testing 
### (Vergouwe et al. 2016; Statistics in Medicine)

ClosedTest <- function(coefs, X, y){
  # Implement closed testing procedure (Version: 11-01-2013)
  # Arguments:
  # coefs: Vector containing the regression coefficients of the model that
  # is updated.
  # X: predictor matrix
  # y: outcome vector
  # Results:
  # coef_new: regression coefficients of chosen model
  require(rms)
  if(class(X)=="data.frame"){
    X <- data.matrix(X)
  }
  if(ncol(X)!=(length(coefs)-1)){
    stop("Number of predictors not equal to the number of coefficients")
  }
  n_coefs <- length(coefs)
  lp_old <- X %*% as.matrix(coefs[2:n_coefs])
  # Calculate updated model intercept
  intercept <- lrm.fit(y = y, offset = lp_old)$coefficients
  coefs_int <- c(intercept , coefs[2:n_coefs])
  # Calculate coefficients after recalibration
  recal <- lrm.fit(x = lp_old, y = y)$coefficients
  coefs_recal <- c(recal[1], recal[2] * coefs[2:n_coefs])
  # Calculate coefficients after model revision
  coefs_refit <- lrm.fit(x = X, y = y)$coefficients
  # Calculate the log-likelihood of the different models
  lp <- cbind(1, X) %*% coefs
  ll_original <- sum(y * lp - log(1 + exp(lp)))
  lp <- cbind(1, X) %*% coefs_int
  ll_intercept <- sum(y * lp - log(1 + exp(lp)))
  lp <- cbind(1, X) %*% coefs_recal
  ll_recalibration <- sum(y * lp - log(1 + exp(lp)))
  lp <- cbind(1, X) %*% coefs_refit
  ll_revision <- sum(y * lp - log(1 + exp(lp)))
  # Calculate difference in log-likelihood for testing of the models
  dev_original <- -2 * ll_original + 2 * ll_revision
  dev_intercept <- -2 * ll_intercept + 2 * ll_revision
  dev_recalibration <- -2 * ll_recalibration + 2 * ll_revision
  # See if difference in model fit was significant
  test1 <- (1 - pchisq(dev_original, ncol(X) + 1)) < 0.05
  test2 <- (1 - pchisq(dev_intercept, ncol(X))) < 0.05
  test3 <- (1 - pchisq(dev_recalibration, ncol(X) - 1)) < 0.05
  # See which model is chosen, index_test indicates the chosen model
  # 1. Original model
  # 2. Model with updated intercept
  # 3. Recalibrated model
  # 4. Revised model
  test_original <- 1 * (!test1)
  test_intercept <- 2 * ((test1)&(!test2))
  test_recalibration <- 3 * ((test1)&(test2)&(!test3))
  test_revision <- 4 * ((test1)&(test2)&(test3))
  index_test <- (test_original + test_intercept + test_recalibration +
                   test_revision)
  coefs_result <- rbind(coefs, coefs_int, coefs_recal, coefs_refit)
  
  # Output of the function
  new_coefs <- coefs_result[index_test, ]
  model <- c("Original Model", "Model with updated intercept",
             "Recalibrated model", "Model Revision")[index_test]
  cat("Method chosen by closed test procedure:\n", model, "\n",
      "Resulting coefficients:\n", new_coefs, "\n")
  res <- list(model = model, coefs = coefs_result)
  return(res)
}

### Function for applying the coefficients and intercept from the 
### published model - Sigmoid Logistic regression function

Prob_LR <- function(df, coefs) {
  
  preds <- data.frame()
  
  for (i in 1:nrow(df)) {
    
    sum <- 0
    for (j in 1:(nrow(coefs)-1)) {
      add <- df[coefs[j+1,1]][i,]*coefs[j+1,2]
      sum <- sum + add
    }
    sum_f <- (exp(coefs[1,2] + sum))/(1 + exp(coefs[1,2] + sum))
    preds <- rbind(preds, sum_f)
  }
  
  colnames(preds) <- "preds"
  rownames(preds) <- rownames(df)
  return(preds)
}

#### Data ####

### Set the working directory to where the "IPD_datasets.xlsx" is located

Dx_time <- read_xlsx("IPD_datasets.xlsx",
                     sheet = 2) # Data on time to diagnosis from blood sampling

db_val <- read_xlsx("IPD_datasets.xlsx",
                    sheet = 1) # Validation dataset

db_val_s <- db_val %>% dplyr::select(-c(status,class2,idandromeda)) # Predictors

# BC status of the samples in the validation dataset
Class_f <- ifelse(db_val$status == "Control",0,1)

# Final dataframe creation
db_val2 <- as.data.frame(db_val)
db_val2$outcome <- ifelse(db_val2$status == "Control",0,1)

##### Cleaned validation dataset #####
db_val3 <- db_val2 %>% dplyr::select(c(-idandromeda,-status, -class2))

# In the validation set, let-7b-5p was replaced by let-7a-5p in the ratios it was a part of, 
# as both of them are found among the 7 miRNA ratios and their Cts were found to be 
# highly correlated (ρ = 0.96). The mature sequences of the two miRNAs being very similar 
# with only two nucleotides being different and, in both miRNAs, the differing bases were purines.

##### Discovery dataset IPD #####

db_disc <- read_xlsx("IPD_datasets.xlsx",
                    sheet = 3)

#### Descriptive statistics ####

### Normality test for the continuous variables

shapiro.test(db_val2$mi_r_199a_3p_let_7a_5p)
shapiro.test(db_val2$mi_r_26b_5p_mi_r_142_5p)
shapiro.test(db_val2$let_7a_5p_mi_r_19b_3p) # let7b replaced by let-7a in qPCR plates (read explanation above)
shapiro.test(db_val2$mi_r_101_3p_mi_r_19b_3p)
shapiro.test(db_val2$mi_r_93_5p_mi_r_19b_3p)
shapiro.test(db_val2$let_7a_5p_mi_r_22_3p)
shapiro.test(db_val2$mi_r_21_5p_mi_r_23a_3p)
shapiro.test(db_val2$score_wcrf)

##### Time to diagnosis from blood sampling #####

# Stratifying the time to diagnosis of cases in validation dataset  
# according to the median
Dx_time$strat <- ifelse(Dx_time$Years >= median(Dx_time$Years),1,0)

plot(density(Dx_time$Years))

# Create a ggplot object
ggplot(Dx_time, aes(x = Years, fill = factor(strat))) +
  geom_density(alpha = 0.5) +
  labs(x = "Years",
       y = "Density") +
  scale_fill_manual(values = c("0" = "blue", "1" = "red")) +
  theme_classic()

ggplot(Dx_time, aes(x=Years, fill=factor(strat))) +
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity',
                  bins = 30) +
  scale_fill_manual(values = c("0" = "blue", "1" = "red")) +
  theme_classic()
  

##### Boxplots #####

# Boxplots on miRNA ratios

# Creating a miRNA ratio only validation dataset
mir_only <- db_val %>% dplyr::select(c(mi_r_199a_3p_let_7a_5p,
                                         mi_r_26b_5p_mi_r_142_5p,
                                         let_7a_5p_mi_r_19b_3p,
                                         mi_r_101_3p_mi_r_19b_3p,
                                         mi_r_93_5p_mi_r_19b_3p,
                                         let_7a_5p_mi_r_22_3p,
                                         mi_r_21_5p_mi_r_23a_3p,
                                         status))

colnames(mir_only) <- c("miR-199a-3p_let-7a-5p", "miR-26b-5p_miR-142-5p", "let-7a-5p_miR-19b-3p", "miR-101-3p_miR-19b-3p",   
                  "miR-93-5p_miR-19b-3p", "let-7a-5p_miR-22-3p", "miR-21-5p_miR-23a-3p", "Class")

# Melting the miRNA only dataset for plotting
db_m <- melt(mir_only)

ggplot(db_m, aes(x = variable, y = value, fill = Class)) +
  geom_boxplot() + 
  theme_classic() +
  facet_wrap(vars(variable), ncol=3, scales = "free") +
  theme(strip.background = element_blank()) +
  theme(panel.background = element_rect(fill = NA, color = "black")) +
  theme(axis.text.x=element_blank()) + 
  stat_compare_means(method = "wilcox.test",label = "p.format",
                     label.x.npc = 0.5, label.y.npc = 0.9) +
  labs(y = "Ratio value", x = "Ratio ID")

### Descriptive statistics on miRNA ratios

df0_clean_T <- mir_only %>% dplyr::filter(Class == "Case")
df0_clean_N <- mir_only %>% dplyr::filter(Class == "Control")

## BC Cases

summary_df_T <- data.frame()

for (i in 1:(ncol(df0_clean_T)-1)) {
  s <- summary(df0_clean_T[[i]])
  qrt1 <- s[[2]]
  med <- s[[3]]
  qrt2 <- s[[5]]
  comb <- cbind(qrt1, med, qrt2)
  summary_df_T <- rbind(summary_df_T, comb)
}

rownames(summary_df_T) <- colnames(df0_clean_T[-8])

summary_df_T

## Healthy controls

summary_df_N <- data.frame()

for (i in 1:(ncol(df0_clean_N)-1)) {
  s <- summary(df0_clean_N[[i]])
  qrt1 <- s[[2]]
  med <- s[[3]]
  qrt2 <- s[[5]]
  comb <- cbind(qrt1, med, qrt2)
  summary_df_N <- rbind(summary_df_N, comb)
}

rownames(summary_df_N) <- colnames(df0_clean_N[-8])

summary_df_N

##### Mann Whitney U and Pearson's chi-squared tests #####

# Performed on all predictors selected in the discovery dataset (12 predictors)

# Mann Whitney U test data frame where the results will be saved

db_val3_cont <- db_val3 %>% dplyr::select(-c(menopause, bmi_di, bmi_men, tabar))

MWU_res <- data.frame()

for (i in 1:(ncol(db_val3_cont)-1)) {
  test <- wilcox.test(db_val3_cont[[i]] ~ db_val3_cont$outcome)
  p <- test$p.value
  stat <- test$statistic
  comb <- cbind(p, stat)
  MWU_res <- rbind(MWU_res, comb)
}

rownames(MWU_res) <- colnames(db_val3_cont[-9])

MWU_res

# Pearson Chi-squared test data frame where the results will be saved

db_val3_cat <- db_val3 %>% dplyr::select(c(menopause, bmi_di, tabar, outcome))

chi_res <- data.frame()

for (i in 1:(ncol(db_val3_cat)-1)) {
  tbl <- table(db_val3_cat[[i]], db_val3_cat$outcome)
  chi_sq <- chisq.test(tbl) # Pearson's chi-squared test
  p <- chi_sq$p.value
  stat <- chi_sq$statistic
  comb <- cbind(p, stat)
  chi_res <- rbind(chi_res, comb)
}

rownames(chi_res) <- colnames(db_val3_cat[-4])

chi_res

##### Univariate logistic regression #####

# miRNA ratios and categorical/dichotomous variables

# Data frame where the univariate logistic regression results will be saved

Uni_LR <- data.frame()

for (i in 1:(ncol(db_val3)-2)) {
  glm <- glm(as.factor(db_val3$outcome) ~ db_val3[[i]],
             family = binomial)
  merge <- t(c((glm$coefficients[2]),
               exp(glm$coefficients[2]),
               exp(confint(glm)[2,]),
               summary(glm)$coefficients[2,4]))
  Uni_LR <- rbind(Uni_LR, merge)
}

rownames(Uni_LR) <- colnames(db_val3[-12:-13])
colnames(Uni_LR) <- c("b", "OR", "lb", "ub", "p")

Uni_LR

# Categorical variable - Tabar breast density scale

TABAR <- glm(as.factor(db_val3$outcome) ~ as.factor(db_val3$tabar), 
             family = binomial)

TABAR1_merge <- t(c((TABAR$coefficients[2]),
                    exp(TABAR$coefficients[2]),
                    exp(confint(TABAR)[2,]),
                    summary(TABAR)$coefficients[2,4]))

TABAR2_merge <- t(c((TABAR$coefficients[3]),
                    exp(TABAR$coefficients[3]),
                    exp(confint(TABAR)[3,]),
                    summary(TABAR)$coefficients[3,4]))

TABAR3_merge <- t(c((TABAR$coefficients[4]),
                    exp(TABAR$coefficients[4]),
                    exp(confint(TABAR)[4,]),
                    summary(TABAR)$coefficients[4,4]))

Tabar_res <- as.data.frame(rbind(TABAR1_merge, TABAR2_merge, TABAR3_merge))
colnames(Tabar_res) <- c("b", "OR", "lb", "ub", "p")
rownames(Tabar_res) <- c("Tabar_2", "Tabar_3", "Tabar_4or5")

Tabar_res

### Tabar as a continuous variable - assuming linearity

TABAR_cont <- glm(as.factor(db_val3$outcome) ~ db_val3$tabar, 
             family = binomial)

TABARcont_merge <- t(c((TABAR_cont$coefficients[2]),
                    exp(TABAR_cont$coefficients[2]),
                    exp(confint(TABAR_cont)[2,]),
                    summary(TABAR_cont)$coefficients[2,4]))

TABARcont_merge

#### Model application ####

# Applying the coefficients from the Discovery paper (Chiorino et al. 2023; Scientific Reports)

# Coefficients from the discovery paper
coefs_paper_all <- as.data.frame(c("Intercept", "mi_r_199a_3p_let_7a_5p", "mi_r_26b_5p_mi_r_142_5p",
"let_7a_5p_mi_r_19b_3p", "mi_r_101_3p_mi_r_19b_3p", "mi_r_93_5p_mi_r_19b_3p", 
"let_7a_5p_mi_r_22_3p", "mi_r_21_5p_mi_r_23a_3p",  "menopause",   
"bmi_di", "bmi_men", "score_wcrf", "tabar"))  
coefs_paper_all$coeffs <- c(1.14,0.16,-0.13,-0.27,-0.10,0.57,-0.05,0.02,-0.03,0.36,0.15,-0.16,0.27)

colnames(coefs_paper_all) <- c("var_names","coeffs")

# Obtaining probabilites of diagnosis
db_val2$pr_model <- Prob_LR(df = db_val2, coefs = coefs_paper_all)$preds

# Calculating Brier score and ROC AUC
BrierScore(db_val2$outcome, db_val2$pr_model)
pROC::auc(db_val2$outcome, db_val2$pr_model)
pROC::ci.auc(db_val2$outcome, db_val2$pr_model)

# Plotting the calibration curve
(app_cal <- CalibrationCurves::val.prob.ci.2(db_val2$pr_model, db_val2$outcome,CL.smooth = F))

app_cal_gg <- valProbggplot(db_val2$pr_model, db_val2$outcome, CL.smooth = T, smooth = "loess",
                            legendloc = FALSE, size = 5, lwd.smooth = 2, lwd.ideal = 2,
                            size.d01 = 7, length.seg = 0.9, d0lab = "Normal", d1lab = "Tumour",
                            dist.label2 = 0.09)

app_cal_gg_f <- app_cal_gg$ggPlot + 
  theme_classic() + theme(title = element_text(size = 25),
                          plot.title= element_text(hjust = 0.5),
                          axis.title = element_text(size = 20),
                          axis.text = element_text(size = 15),
                          legend.position = "none") +
  annotate("text", x = 0.5, y = 0.99, size = 12, label = "B",
           fontface = 2)

# tiff("Calibration_plot_application.tiff",
#      width = 3200, height = 1600, res = 300)

app_cal_gg_f

# dev.off()

roc_all <- roc(outcome ~ pr_model, ci = TRUE, data = db_val2)
# Sensitivity and specificity at the Youden's optimal cut-off
coords(roc_all, x="best", input="threshold", best.method="youden")

# Violin plot of the predicted probabilities stratified by BC status
# in the applied model

# tiff("Violin_probs_application.tiff",
#      width = 3200, height = 1600, res = 300)

ggplot(db_val2, aes(x=as.factor(status), y=pr_model, fill = as.factor(status))) + # fill=name allow to automatically dedicate a color for each group
  geom_violin() + 
  theme_classic() + 
  labs(x = "BC Status",
       y = "Predicted risk - category",
       fill = "Status")

# dev.off()

# ROC AUC curve

perf <- ROCR::performance(prediction(db_val2$pr_model, db_val2$outcome), 
                          'sens', 'fpr')

roc_data_all <- data.frame(
  FPR = perf@x.values[[1]],
  TPR = perf@y.values[[1]]
)


roc_plot_all <- ggplot(roc_data_all, aes(x = FPR, y = TPR)) +
  geom_line(color = "blue", size = 1.5) +  # Increase line size
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", 
              color = "gray", size = 1.2) +  # Increase line size
  geom_point(x = 1 - 0.9133858, y = 0.28125, color = "red", size = 3) +  # Add red dot
  annotate("text", x = 1 - 0.91, y = 0.3, 
           label = paste("Youden's cut-off:\n", round(0.8189862,3)), 
           color = "red", size = 5, vjust = -1) +
  annotate("text", x = 0.8, y = 0.3, 
           label = paste("AUC: ", round(auc(db_val2$outcome, db_val2$pr_model), 3)), 
           color = "blue", size = 7, vjust = -1) + 
  annotate("text", x = 0.5, y = 0.99, size = 12, label = "A",
           fontface = 2) +
  labs(x = "False Positive Rate",
       y = "True Positive Rate") +
  theme_classic() +
  theme(title = element_text(size = 25),
        plot.title= element_text(hjust = 0.5),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 15))

# tiff("ROC_application.tiff",
#      width = 2400, height = 2200, res = 300)

roc_plot_all

# dev.off()

tiff("Figure 3.tiff",
     width = 5800, height = 2400, res = 300)
ggarrange(roc_plot_all, app_cal_gg_f)

dev.off()

### Time to diagnosis and predicted probabilities

Cases_db <- merge(db_val2, Dx_time, by = "idandromeda")

# Linear regression between the predicted probabilities among cases 
# and time to diagnosis from blood sampling

linear_age_pred <- lm(pr_model ~ Years, Cases_db)
summary(linear_age_pred)

#### Healthy control subgroup analysis ####

# Checking differences between controls which underwent second level investigation
# and controls which did not

df_N <- db_val2 %>% dplyr::filter(class2 == "N") # controls which did not undergo second level investigation
df_N2 <- db_val2 %>% dplyr::filter(class2 == "N2") # controls which underwent second level investigation

# Selecting the appropriate variables in the two subgroups

x_N <- df_N %>% dplyr::select(c(starts_with("mi_r"), starts_with("let"),
                                bmi_di, menopause, bmi_men, score_wcrf, tabar)) %>% as.matrix()

x_N2 <- df_N2 %>% dplyr::select(c(starts_with("mi_r"), starts_with("let"),
                                  bmi_di, menopause, bmi_men, score_wcrf, tabar)) %>% as.matrix()

df1 <- as.data.frame(x_N)
df2 <- as.data.frame(x_N2)

# Comparing the selected predictors between the two subgroups

# Continuous variables
# Controls which did not undergo second level investigation

df1_cont <- df1 %>% dplyr::select(-c(bmi_di, menopause, tabar, bmi_men))

summary(df1_cont$mi_r_199a_3p_let_7a_5p)
sd(df1_cont$mi_r_199a_3p_let_7a_5p)

summary(df1_cont$mi_r_26b_5p_mi_r_142_5p)
sd(df1_cont$mi_r_26b_5p_mi_r_142_5p)

summary(df1_cont$mi_r_101_3p_mi_r_19b_3p)
sd(df1_cont$mi_r_101_3p_mi_r_19b_3p)

summary(df1_cont$mi_r_93_5p_mi_r_19b_3p)
sd(df1_cont$mi_r_93_5p_mi_r_19b_3p)

summary(df1_cont$mi_r_21_5p_mi_r_23a_3p)
sd(df1_cont$mi_r_21_5p_mi_r_23a_3p)

summary(df1_cont$let_7a_5p_mi_r_19b_3p)
sd(df1_cont$let_7a_5p_mi_r_19b_3p)

summary(df1_cont$let_7a_5p_mi_r_22_3p)
sd(df1_cont$let_7a_5p_mi_r_22_3p)

summary(df1_cont$score_wcrf)
sd(df1_cont$score_wcrf)

###
# Controls which underwent second level investigation

df2_cont <- df2 %>% dplyr::select(-c(bmi_di, menopause, tabar, bmi_men))

summary(df2_cont$mi_r_199a_3p_let_7a_5p)
sd(df2_cont$mi_r_199a_3p_let_7a_5p)

summary(df2_cont$mi_r_26b_5p_mi_r_142_5p)
sd(df2_cont$mi_r_26b_5p_mi_r_142_5p)

summary(df2_cont$mi_r_101_3p_mi_r_19b_3p)
sd(df2_cont$mi_r_101_3p_mi_r_19b_3p)

summary(df2_cont$mi_r_93_5p_mi_r_19b_3p)
sd(df2_cont$mi_r_93_5p_mi_r_19b_3p)

summary(df2_cont$mi_r_21_5p_mi_r_23a_3p)
sd(df2_cont$mi_r_21_5p_mi_r_23a_3p)

summary(df2_cont$let_7a_5p_mi_r_19b_3p)
sd(df2_cont$let_7a_5p_mi_r_19b_3p)

summary(df2_cont$let_7a_5p_mi_r_22_3p)
sd(df2_cont$let_7a_5p_mi_r_22_3p)

summary(df2_cont$score_wcrf)
sd(df2_cont$score_wcrf)

###

Stat_cont <- data.frame()

for (i in 1:ncol(df1_cont)) {
  wilcox <- wilcox.test(df1_cont[[i]], df2_cont[[i]]) # Wilcoxon two-sample test
  wilcox_res <- cbind(wilcox$statistic, wilcox$p.value)
  colnames(wilcox_res) <- c("wil_s", "wil_p")
  
  Stat_cont <- rbind(Stat_cont, wilcox_res)
  
}

row.names(Stat_cont) <- colnames(df1_cont)

Stat_cont

# Categorical variables

df_ctrls <- db_val2 %>% dplyr::filter(class2 != "T")

# Controls which did not undergo second level investigation
df1_cat <- df1 %>% dplyr::select(c(bmi_di, menopause, tabar))

summary(as.factor(df1_cat$tabar))
summary(as.factor(df1_cat$bmi_di))
summary(as.factor(df1_cat$menopause))

# Controls which underwent second level investigation
df2_cat <- df2 %>% dplyr::select(c(bmi_di, menopause, tabar))

summary(as.factor(df2_cat$tabar))
summary(as.factor(df2_cat$bmi_di))
summary(as.factor(df2_cat$menopause))

### Breast density

Tabar_tbl <- table(df_ctrls$tabar, df_ctrls$class2)
chisq.test(Tabar_tbl) # Pearson's chi-squared test

### BMI

BMI_tbl <- table(df_ctrls$bmi_di, df_ctrls$class2)
chisq.test(BMI_tbl)

### Menopause

Meno_tbl <- table(df_ctrls$menopause, df_ctrls$class2)
chisq.test(Meno_tbl)

# Comparing the two healthy control subtypes - model application

db_tot_f_N_all <- db_val2 %>% filter(class2 == "N" | class2 == "T") # Controls which did not undergo second level investigation
db_tot_f_N2_all <- db_val2 %>% filter(class2 == "N2" | class2 == "T") # Controls which underwent second level investigation

# Calibration curve of the model applied to the dataset with 
# Healthy controls that did not undergo second level investigation

# Calculating Brier score and ROC AUC
BrierScore(db_tot_f_N_all$outcome, db_tot_f_N_all$pr_model)
pROC::auc(db_tot_f_N_all$outcome, db_tot_f_N_all$pr_model)
pROC::ci.auc(db_tot_f_N_all$outcome, db_tot_f_N_all$pr_model)

CalibrationCurves::val.prob.ci.2(db_tot_f_N_all$pr_model, db_tot_f_N_all$outcome,CL.smooth = F)

# Calibration curve of the model applied to the dataset with 
# healthy controls that underwent second level investigation

# Calculating Brier score and ROC AUC
BrierScore(db_tot_f_N2_all$outcome, db_tot_f_N2_all$pr_model)
pROC::auc(db_tot_f_N2_all$outcome, db_tot_f_N2_all$pr_model)
pROC::ci.auc(db_tot_f_N2_all$outcome, db_tot_f_N2_all$pr_model)

CalibrationCurves::val.prob.ci.2(db_tot_f_N2_all$pr_model, db_tot_f_N2_all$outcome,CL.smooth = F)

# Removing the cases from the dataframe
db_tot_f_cont_all <- db_val2 %>% filter(class2 != "T")

# Violin plots of the predicted probabilities in the applied model
# among the compared healthy control subtypes

ggplot(db_tot_f_cont_all, aes(x=as.factor(class2), y=pr_model, fill = as.factor(class2))) +
  geom_violin() + 
  theme_classic() + 
  labs(x = "Healthy control subtype",
       y = "Predicted risk - category",
       fill = "Subtype")

#### Model Updating ####

model_update <- ClosedTest(coefs = coefs_paper_all$coeffs, 
                      X = db_val3[,-13], y = Class_f)

# Based on the closed test function the predicted probabilities need to be
# updated through a complete model revision

##### Ridge Regression #####

# Instead of using the coefficients suggested by the ClosedTest functiom
# we will perform model revision using ridge regression to 
# diminish overfitting and overoptimism in the validation dataset

x1 <- as.matrix(db_val_s)
y1 <- as.matrix(Class_f)

set.seed(1001)

cvfit = cv.glmnet(x1, y1, alpha = 0, family = "binomial")
tmp_coeffs <- coef(cvfit, s = "lambda.min")
var_set <- data.frame(name = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], 
                      coefficient = tmp_coeffs@x)

# New coefficients of the updated model
var_set

# New predicted probabilities in the validation dataset
Preds_ridge <- Prob_LR(df = db_val_s, coefs = var_set)

db_val2$preds_ridge <- Preds_ridge$preds

# ROC AUC of the updated model
BrierScore(db_val2$outcome, db_val2$preds_ridge)
pROC::auc(db_val2$outcome, db_val2$preds_ridge)
pROC::ci.auc(db_val2$outcome, db_val2$preds_ridge)

# Calibration curve of the updated model
(cal_up <- CalibrationCurves::val.prob.ci.2(db_val2$preds_ridge, 
                                            db_val2$outcome,CL.smooth = F))

cal_up_gg <- valProbggplot(db_val2$preds_ridge, 
                           db_val2$outcome, CL.smooth = T, smooth = "loess",
                            legendloc = FALSE, size = 5, lwd.smooth = 2, lwd.ideal = 2,
                            size.d01 = 7, length.seg = 0.9, d0lab = "Normal", d1lab = "Tumour",
                            dist.label2 = 0.09)

cal_up_gg_f <- cal_up_gg$ggPlot + 
  theme_classic() + theme(title = element_text(size = 25),
                          plot.title= element_text(hjust = 0.5),
                          axis.title = element_text(size = 20),
                          axis.text = element_text(size = 15),
                          legend.position = "none") +
  annotate("text", x = 0.5, y = 0.99999999999999, size = 12, label = "B",
           fontface = 2)

# tiff("Calibration_update.tiff",
#      width = 3200, height = 1600, res = 300)

cal_up_gg_f

# dev.off()

roc_all_ridge <- roc(outcome ~ preds_ridge, ci = TRUE, data = db_val2)
# Sensitivity and specificity of the updated model
# at the Youden's optimal cut-off
coords(roc_all_ridge, x="best", input="threshold", best.method="youden")

# ROC AUC curve

perf_up <- ROCR::performance(prediction(db_val2$preds_ridge, db_val2$outcome), 
                             'sens', 'fpr')

roc_data_up <- data.frame(
  FPR = perf_up@x.values[[1]],
  TPR = perf_up@y.values[[1]]
)


roc_plot_up <- ggplot(roc_data_up, aes(x = FPR, y = TPR)) +
  geom_line(color = "blue", size = 1.5) +  # Increase line size
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", 
              color = "gray", size = 1.2) +  # Increase line size
  geom_point(x = 1 - 0.7007874, y = 0.96875, color = "red", size = 3) +  # Add red dot
  annotate("text", x = 1 - 0.58, y = 0.80, 
           label = paste("Youden's cut-off:\n", round(0.1730981,3)), 
           color = "red", size = 5, vjust = -1) +
  annotate("text", x = 0.8, y = 0.3, 
           label = paste("AUC: ", round(auc(db_val2$outcome, db_val2$preds_ridge), 3)), 
           color = "blue", size = 7, vjust = -1) + 
  annotate("text", x = 0.5, y = 0.99999999999999, size = 12, label = "A",
           fontface = 2) +
  labs(x = "False Positive Rate",
       y = "True Positive Rate") +
  theme_classic() +
  theme(title = element_text(size = 25),
        plot.title= element_text(hjust = 0.5),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 15))

# tiff("ROC_update.tiff",
#      width = 2400, height = 2200, res = 300)

roc_plot_up

# dev.off()

tiff("Figure 4.tiff",
     width = 5800, height = 2400, res = 300)
ggarrange(roc_plot_up, cal_up_gg_f)

dev.off()

# Violin plot of the predicted probabilities stratified by BC status
# in the updated model

# tiff("Violin_probs_update.tiff",
#      width = 3200, height = 1600, res = 300)

ggplot(db_val2, aes(x=as.factor(status), y=preds_ridge, fill = as.factor(status))) + # fill=name allow to automatically dedicate a color for each group
  geom_violin() + 
  theme_classic() + 
  labs(x = "BC Status",
       y = "Predicted probability",
       fill = "Status")

# dev.off()

# Investigating potential association of predicted probabilities (updated model) 
# with time to diagnosis from blood sampling

Cases_db_up <- merge(db_val2, Dx_time, by = "idandromeda")
# Linear regression between the predicted probabilities (updated model) 
# among cases and time to diagnosis from blood sampling
linear_age_pred_up <- lm(preds_ridge ~ Years, Cases_db_up)
summary(linear_age_pred_up)

###### Bootstrap ######

# Bootstrap on the ridge logistic regression model updating was performed

#Creating the function
bootstrap_ridge1 <- function(data, indices){
  
  data <- data[indices,]
  y1 <- as.matrix(data$outcome)
  b <- data %>% dplyr::select(-outcome)
  x1 <- as.matrix(b)
  
  cvfit = cv.glmnet(x1, y1, alpha = 0, family = "binomial")
  tmp_coeffs <- coef(cvfit, s = "lambda.min")
  var_set <- data.frame(name = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], coefficient = tmp_coeffs@x)
  
  Preds_ridge_boot <- Prob_LR(df = b, coefs = var_set)
  
  pROC::auc(y1, Preds_ridge_boot$preds)
} 

set.seed(1001)
bootauc <- boot(data = db_val3, statistic = bootstrap_ridge1, R = 2000)
bootauc
boot.ci(boot.out = bootauc, type = c("norm", "basic", "perc", "bca"))

tiff("Supplementary Figure 1.tiff",
      width = 3200, height = 1600, res = 300)

plot(bootauc)

dev.off()

##### Bayesian model updating #####

set.seed(1001)
bay_logistic <- brm(outcome ~ ., data = as.data.frame(db_val3), family = bernoulli(logit),
                    prior = c(set_prior("normal(1.14,0.6931472)", class = "Intercept"),
                              set_prior("normal(0.16,0.6931472)", class = "b", coef = "mi_r_199a_3p_let_7a_5p"),
                              set_prior("normal(-0.13,0.6931472)", class = "b", coef = "mi_r_26b_5p_mi_r_142_5p"),
                              set_prior("normal(-0.27,0.6931472)", class = "b", coef = "let_7a_5p_mi_r_19b_3p"),
                              set_prior("normal(-0.10,0.6931472)", class = "b", coef = "mi_r_101_3p_mi_r_19b_3p"),
                              set_prior("normal(0.57,0.6931472)", class = "b", coef = "mi_r_93_5p_mi_r_19b_3p"),
                              set_prior("normal(-0.05,0.6931472)", class = "b", coef = "let_7a_5p_mi_r_22_3p"),
                              set_prior("normal(0.02,0.6931472)", class = "b", coef = "mi_r_21_5p_mi_r_23a_3p"),
                              set_prior("normal(-0.03,0.6931472)", class = "b", coef = "menopause"),
                              set_prior("normal(0.36,0.6931472)", class = "b", coef = "bmi_di"),
                              set_prior("normal(0.15,0.6931472)", class = "b", coef = "bmi_men"),
                              set_prior("normal(-0.16,0.6931472)", class = "b", coef = "score_wcrf"),
                              set_prior("normal(0.27,0.6931472)", class = "b", coef = "tabar")))

summary(bay_logistic)
bay_logistic %>% posterior_interval(prob = c(0.9))

# Bayesian model updating coefficients
estimates <- fixef(bay_logistic)
estimates_adt <- cbind(rownames(estimates), as.data.frame(estimates))

# Predicted probabilities from coefficients obtained from Bayesian modelling
Preds_bay <- Prob_LR(df = db_val3, coefs = estimates_adt)

# ROC AUC of the model updated by the Bayesian method
BrierScore(db_val3$outcome, Preds_bay$preds)
pROC::auc(db_val3$outcome, Preds_bay$preds)
pROC::ci.auc(db_val3$outcome, Preds_bay$preds)
# Calibration curve of the model updated by the Bayesian method
CalibrationCurves::val.prob.ci.2(Preds_bay$preds, db_val3$outcome,CL.smooth = F)

cal_bay_gg <- valProbggplot(Preds_bay$preds, db_val3$outcome, CL.smooth = T, smooth = "loess",
                           legendloc = FALSE, size = 5, lwd.smooth = 2, lwd.ideal = 2,
                           size.d01 = 7, length.seg = 0.9, d0lab = "Normal", d1lab = "Tumour",
                           dist.label2 = 0.09)

cal_bay_gg_f <- cal_bay_gg$ggPlot +
  annotate("text", x = 0.5, y = 0.99, size = 12, label = "B",
           fontface = 2) +
  theme_classic() + theme(title = element_text(size = 25),
                          plot.title= element_text(hjust = 0.5),
                          axis.title = element_text(size = 20),
                          axis.text = element_text(size = 15),
                          legend.position = "none")

cal_bay_gg_f

### ROC AUC

roc_Bay_res <- roc(db_val3$outcome ~ Preds_bay$preds, ci = TRUE)
# Sensitivity and specificity at the Youden's optimal cut-off
coords(roc_Bay_res, x="best", input="threshold", best.method="youden")

perf_Bay <- ROCR::performance(prediction(Preds_bay$preds, db_val3$outcome), 'sens', 'fpr')

roc_Bay_data <- data.frame(
  FPR = perf_Bay@x.values[[1]],
  TPR = perf_Bay@y.values[[1]]
)


roc_plot_Bay <- ggplot(roc_Bay_data, aes(x = FPR, y = TPR)) +
  geom_line(color = "blue", size = 1.5) +  # Increase line size
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", 
              color = "gray", size = 1.2) +
  annotate("text", x = 0.8, y = 0.3, 
           label = paste("AUC: ", round(auc(db_val3$outcome, 
                                            Preds_bay$preds), 3)), 
           color = "blue", size = 7, vjust = -1) + 
  annotate("text", x = 0.5, y = 0.99999999999, size = 12, label = "A",
           fontface = 2) +
  geom_point(x = 1 - 0.6614173, y = 0.96875, color = "red", size = 3) +  # Add red dot
  annotate("text", x = 1 - 0.57, y = 0.80, 
           label = paste("Youden's cut-off:\n", round(0.1886786,3)), 
           color = "red", size = 5, vjust = -1) +
  labs(x = "False Positive Rate",
       y = "True Positive Rate") +
  theme_classic() +
  theme(title = element_text(size = 25),
        plot.title= element_text(hjust = 0.5),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 15))

roc_plot_Bay

tiff("Supplementary Figure 2.tiff",
     width = 5800, height = 2400, res = 300)
ggarrange(roc_plot_Bay, cal_bay_gg_f)

dev.off()

##### IECV #####

# IECV stands for Internal External Cross-validation, please refer to 
# de Jong et al. 2021; Statistics in Medicine

db_disc$study <- "Disc"

db_val_IECV <- db_val3

names(db_val_IECV)[names(db_val_IECV) == 'let_7a_5p_mi_r_19b_3p'] <- 'let_7b_5p_mi_r_19b_3p'
# Naming let-7a miRNA as let-7b so that the variable can be merged with the discovery set

db_val_IECV$study <- "Val"

# merging the IPD of discovery and validation datasets
db_disc_clean <- db_disc[,-1]

comb_all_new <- as.data.frame(rbind(db_disc_clean, db_val_IECV))

###### IECV on all predictors ######

# design
f1 <- outcome ~ mi_r_199a_3p_let_7a_5p + mi_r_26b_5p_mi_r_142_5p + let_7b_5p_mi_r_19b_3p +
  mi_r_101_3p_mi_r_19b_3p + mi_r_93_5p_mi_r_19b_3p + let_7a_5p_mi_r_22_3p + mi_r_21_5p_mi_r_23a_3p + 
  bmi_di + menopause + bmi_men + score_wcrf + tabar

# IECV model
fit_all_new <- metapred(comb_all_new, strata = "study", formula = f1, scope = f1, 
                        family = binomial,  
                        meta.method = "REML", two.stage = F, 
                        estFUN = "glm")

(IECV_model_all_new <- as.data.frame(coef(fit_all_new$global.model)))

#generalisability of the full IECV model
gen(fit_all_new)

# Meta analysis of the perforamances in the full IECV model
metamisc::ma(fit_all_new, method ="REML")

###### IECV on generalisable predictors ######

fit_sel_new <- metapred(comb_all_new, strata = "study", formula = f1, family = binomial,
                        meta.method = "REML", two.stage = F,
                        estFUN = "glm")

(IECV_model_sel_new <- as.data.frame(coef(fit_sel_new$global.model)))
IECV_model_sel_f_new <- cbind(rownames(IECV_model_sel_new), IECV_model_sel_new)

#generalisability of the reduced IECV model
gen(fit_sel_new)

# Meta analysis of the perforamances in the reduced IECV model
metamisc::ma(fit_sel_new, method ="REML")

Pred_IECV_NM_new <- Prob_LR(df = comb_all_new,coefs = IECV_model_sel_f_new)

BrierScore(comb_all_new$outcome, Pred_IECV_NM_new$preds)
pROC::auc(comb_all_new$outcome, Pred_IECV_NM_new$preds)
pROC::ci.auc(comb_all_new$outcome, Pred_IECV_NM_new$preds)
CalibrationCurves::val.prob.ci.2(Pred_IECV_NM_new$preds, 
                                 comb_all_new$outcome,CL.smooth = F)

### 

(all_cal_IECV <- CalibrationCurves::val.prob.ci.2(Pred_IECV_NM_new$preds, comb_all_new$outcome,CL.smooth = T,
                                                 smooth = "loess", lwd.ideal = 2, lwd.smooth = 3,
                                                 cex = 1.5, cex.leg = 1,
                                                 cex.d01 = 1.5, cex.lab = 2, 
                                                 main = "miRNA ratios + NM predictors",
                                                 cex.main = 2, cex.axis = 1.5,
                                                 legendloc = c(0.75, 0.7)))

all_cal_IECV <- valProbggplot(Pred_IECV_NM_new$preds, comb_all_new$outcome, CL.smooth = T, smooth = "loess",
                           legendloc = FALSE, size = 5, lwd.smooth = 2, lwd.ideal = 2,
                           size.d01 = 7, length.seg = 0.9, d0lab = "Normal", d1lab = "Tumour",
                           dist.label2 = 0.09)

all_cal_IECV_gg <- all_cal_IECV$ggPlot + 
  annotate("text", x = 0.5, y = 0.99, size = 12, label = "B",
           fontface = 2) +
  theme_classic() + theme(title = element_text(size = 25),
                          plot.title= element_text(hjust = 0.5),
                          axis.title = element_text(size = 20),
                          axis.text = element_text(size = 15),
                          legend.position = "none")
all_cal_IECV_gg

### ROC AUC

roc_IECV_res <- roc(comb_all_new$outcome ~ Pred_IECV_NM_new$preds, ci = TRUE)
# Sensitivity and specificity at the Youden's optimal cut-off
coords(roc_IECV_res, x="best", input="threshold", best.method="youden")

perf_IECV <- ROCR::performance(prediction(Pred_IECV_NM_new$preds, comb_all_new$outcome), 'sens', 'fpr')

roc_IECV_data <- data.frame(
  FPR = perf_IECV@x.values[[1]],
  TPR = perf_IECV@y.values[[1]]
)


roc_plot_IECV <- ggplot(roc_IECV_data, aes(x = FPR, y = TPR)) +
  geom_line(color = "blue", size = 1.5) +  # Increase line size
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", 
              color = "gray", size = 1.2) +
  annotate("text", x = 0.8, y = 0.3, 
           label = paste("AUC: ", round(auc(comb_all_new$outcome, 
                                            Pred_IECV_NM_new$preds), 3)), 
           color = "blue", size = 7, vjust = -1) + 
  geom_point(x = 1 - 0.7708333, y = 0.7319588, color = "red", size = 3) +  # Add red dot
  annotate("text", x = 1 - 0.80, y = 0.76, 
           label = paste("Youden's cut-off:\n", round(0.3364745,3)), 
           color = "red", size = 5, vjust = -1) +
  annotate("text", x = 0.5, y = 0.99, size = 12, label = "A",
           fontface = 2) +
  labs(x = "False Positive Rate",
       y = "True Positive Rate") +
  theme_classic() +
  theme(title = element_text(size = 25),
        plot.title= element_text(hjust = 0.5),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 15))

roc_plot_IECV

###

tiff("Supplementary Figure 3.tiff",
      width = 5800, height = 2400, res = 300)

ggarrange(roc_plot_IECV,all_cal_IECV_gg)

dev.off()

#### miRNA only model ####

##### model application #####

# Intercept and coefficients from the discovery paper on miRNA ratios only
coefs_paper_mir <- as.data.frame(c("Intercept", "mi_r_199a_3p_let_7a_5p", "mi_r_26b_5p_mi_r_142_5p",
                                   "let_7a_5p_mi_r_19b_3p", "mi_r_101_3p_mi_r_19b_3p", "mi_r_93_5p_mi_r_19b_3p", 
                                   "let_7a_5p_mi_r_22_3p"))  
coefs_paper_mir$coeffs <- c(1.32,0.16,-0.16,-0.18,-0.07,0.47,-0.04)

colnames(coefs_paper_mir) <- c("var_names","coeffs")

# Obtaining probabilites of diagnosis
db_val2$pr_mirs <- Prob_LR(df = db_val2, coefs = coefs_paper_mir)$preds

# Computing the Brier Score and ROC AUC

BrierScore(db_val2$outcome, db_val2$pr_mirs)
pROC::auc(db_val2$outcome, db_val2$pr_mirs)
pROC::ci.auc(db_val2$outcome, db_val2$pr_mirs)

# Creating the calibration curve

CalibrationCurves::val.prob.ci.2(db_val2$pr_mirs, db_val2$outcome,CL.smooth = F)

# Violin plot of the predicted probabilities stratified by BC status
# in the applied miRNA-only model

ggplot(db_val2, aes(x=as.factor(outcome), y=pr_mirs, fill = as.factor(status))) + # fill=name allow to automatically dedicate a color for each group
  geom_violin() + 
  theme_classic() + 
  labs(x = "Class",
       y = "Predicted risk - category",
       fill = "Class")

(app_cal_mir <- CalibrationCurves::val.prob.ci.2(db_val2$pr_mirs, db_val2$outcome,CL.smooth = F))

app_cal_mir_gg <- valProbggplot(db_val2$pr_mirs, db_val2$outcome, CL.smooth = T, smooth = "loess",
                            legendloc = FALSE, size = 5, lwd.smooth = 2, lwd.ideal = 2,
                            size.d01 = 7, length.seg = 0.9, d0lab = "Normal", d1lab = "Tumour",
                            dist.label2 = 0.09)

app_cal_gg_mir_f <- app_cal_mir_gg$ggPlot + 
  theme_classic() + theme(title = element_text(size = 25),
                          plot.title= element_text(hjust = 0.5),
                          axis.title = element_text(size = 20),
                          axis.text = element_text(size = 15),
                          legend.position = "none") +
  annotate("text", x = 0.5, y = 0.99, size = 12, label = "B",
           fontface = 2)

app_cal_gg_mir_f

### ROC AUC

roc_mir <- roc(outcome ~ pr_mirs, ci = TRUE, data = db_val2)
coords(roc_mir, x="best", input="threshold", best.method="youden")

perf_mir <- ROCR::performance(prediction(db_val2$pr_mirs, db_val2$outcome), 
                          'sens', 'fpr')

roc_data_mir <- data.frame(
  FPR = perf_mir@x.values[[1]],
  TPR = perf_mir@y.values[[1]]
)


roc_plot_mir <- ggplot(roc_data_mir, aes(x = FPR, y = TPR)) +
  geom_line(color = "blue", size = 1.5) +  # Increase line size
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", 
              color = "gray", size = 1.2) +  # Increase line size
  geom_point(x = 1 - 0.7716535, y = 0.35, color = "red", size = 3) +  # Add red dot
  annotate("text", x = 1 - 0.82, y = 0.35, 
           label = paste("Youden's cut-off:\n", round(0.6511037,3)), 
           color = "red", size = 5, vjust = -1) +
  annotate("text", x = 0.8, y = 0.3, 
           label = paste("AUC: ", round(auc(db_val2$outcome, db_val2$pr_mirs), 3)), 
           color = "blue", size = 7, vjust = -1) + 
  annotate("text", x = 0.5, y = 0.99, size = 12, label = "A",
           fontface = 2) +
  labs(x = "False Positive Rate",
       y = "True Positive Rate") +
  theme_classic() +
  theme(title = element_text(size = 25),
        plot.title= element_text(hjust = 0.5),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 15))

roc_plot_mir

ggarrange(roc_plot_mir, app_cal_gg_mir_f)

##### Healthy control subgroup model application #####

db_tot_f_N_mir <- db_val2 %>% filter(class2 == "N" | class2 == "T") # controls which did not undergo second level investigation
db_tot_f_N2_mir <- db_val2 %>% filter(class2 == "N2" | class2 == "T") # controls which underwent second level investigation

# Calibration curve of the model applied to the dataset with 
# healthy controls that did not undergo second level investigation

# Calculating Brier score and ROC AUC
BrierScore(db_tot_f_N_mir$outcome, db_tot_f_N_mir$pr_mirs)
pROC::auc(db_tot_f_N_mir$outcome, db_tot_f_N_mir$pr_mirs)
pROC::ci.auc(db_tot_f_N_mir$outcome, db_tot_f_N_mir$pr_mirs)

CalibrationCurves::val.prob.ci.2(db_tot_f_N_mir$pr_mirs, db_tot_f_N_mir$outcome,CL.smooth = F)

# Calibration curve of the model applied to the dataset with 
# healthy controls that underwent second level investigation

# Calculating Brier score and ROC AUC
BrierScore(db_tot_f_N2_mir$outcome, db_tot_f_N2_mir$pr_mirs)
pROC::auc(db_tot_f_N2_mir$outcome, db_tot_f_N2_mir$pr_mirs)
pROC::ci.auc(db_tot_f_N2_mir$outcome, db_tot_f_N2_mir$pr_mirs)

CalibrationCurves::val.prob.ci.2(db_tot_f_N2_mir$pr_mirs, db_tot_f_N2_mir$outcome,CL.smooth = F)

### miR only model updating

a <- db_val2 %>% dplyr::select(c(mi_r_199a_3p_let_7a_5p,mi_r_26b_5p_mi_r_142_5p,
                                 let_7a_5p_mi_r_19b_3p,mi_r_101_3p_mi_r_19b_3p,
                                 mi_r_93_5p_mi_r_19b_3p,let_7a_5p_mi_r_22_3p))

##### Closed testing #####

miRs_only <- ClosedTest(coefs = coefs_paper_mir$coeffs, X = a, y = Class_f)

##### Ridge Regression #####

# Instead of classical remodeling we will perform ridge regression to 
# reduce overfitting and overoptimisim in the validation dataset and obtain
# more applicable coefficients for future studies

x1 <- as.matrix(a)
y1 <- as.matrix(Class_f)

set.seed(1001)

cvfit = cv.glmnet(x1, y1, alpha = 0, family = "binomial")
tmp_coeffs <- coef(cvfit, s = "lambda.min")
var_set <- data.frame(name = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], coefficient = tmp_coeffs@x)
var_set

db_val2$preds_ridge_mir <- Prob_LR(df = a, coefs = var_set)$preds

BrierScore(db_val2$outcome, db_val2$preds_ridge_mir)
pROC::auc(db_val2$outcome, db_val2$preds_ridge_mir)
pROC::ci.auc(db_val2$outcome, db_val2$preds_ridge_mir)

CalibrationCurves::val.prob.ci.2(db_val2$preds_ridge_mir, db_val2$outcome,CL.smooth = F)

# Violin plot of the predicted probabilities stratified by BC status
# in the updated miRNA-only model

ggplot(db_val2, aes(x=as.factor(status), y=preds_ridge_mir, fill = as.factor(status))) +
  geom_violin() + 
  theme_classic() + 
  labs(x = "BC Status",
       y = "Predicted probability",
       fill = "Status")

### Calibration curve

(app_cal_mir_up <- CalibrationCurves::val.prob.ci.2(db_val2$preds_ridge_mir, db_val2$outcome,CL.smooth = F))

app_cal_mir_up_gg <- valProbggplot(db_val2$preds_ridge_mir, db_val2$outcome, CL.smooth = T, smooth = "loess",
                                legendloc = FALSE, size = 5, lwd.smooth = 2, lwd.ideal = 2,
                                size.d01 = 7, length.seg = 0.9, d0lab = "Normal", d1lab = "Tumour",
                                dist.label2 = 0.09)

app_cal_gg_mir_up_f <- app_cal_mir_up_gg$ggPlot + 
  theme_classic() + theme(title = element_text(size = 25),
                          plot.title= element_text(hjust = 0.5),
                          axis.title = element_text(size = 20),
                          axis.text = element_text(size = 15),
                          legend.position = "none") +
  annotate("text", x = 0.5, y = 0.99, size = 12, label = "B",
           fontface = 2)

app_cal_gg_mir_up_f

### ROC AUC

roc_mir_ridge <- roc(outcome ~ preds_ridge_mir, ci = TRUE, data = db_val2)
coords(roc_mir_ridge, x="best", input="threshold", best.method="youden")

perf_mir_up <- ROCR::performance(prediction(db_val2$preds_ridge_mir, db_val2$outcome), 
                          'sens', 'fpr')

roc_data_mir_up <- data.frame(
  FPR = perf_mir_up@x.values[[1]],
  TPR = perf_mir_up@y.values[[1]]
)


roc_plot_mir_up <- ggplot(roc_data_mir_up, aes(x = FPR, y = TPR)) +
  geom_line(color = "blue", size = 1.5) +  # Increase line size
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", 
              color = "gray", size = 1.2) +  # Increase line size
  geom_point(x = 1 - 0.6062992, y = 0.84375, color = "red", size = 3) +  # Add red dot
  annotate("text", x = 1 - 0.68, y = 0.80, 
           label = paste("Youden's cut-off:\n", round(0.1853202,3)), 
           color = "red", size = 5, vjust = -1) +
  annotate("text", x = 0.8, y = 0.3, 
           label = paste("AUC: ", round(auc(db_val2$outcome, db_val2$preds_ridge_mir), 3)), 
           color = "blue", size = 7, vjust = -1) + 
  annotate("text", x = 0.5, y = 0.99, size = 12, label = "A",
           fontface = 2) +
  labs(x = "False Positive Rate",
       y = "True Positive Rate") +
  theme_classic() +
  theme(title = element_text(size = 25),
        plot.title= element_text(hjust = 0.5),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 15))

roc_plot_mir_up

ggarrange(roc_plot_mir_up, app_cal_gg_mir_up_f)

###### Bootstrap ######

mir_only_clean <- a
mir_only_clean$outcome <- Class_f

set.seed(1001)
bootauc_mir_ratio <- boot(data = mir_only_clean, statistic = bootstrap_ridge1, R = 2000)
bootauc_mir_ratio
boot.ci(boot.out = bootauc_mir_ratio, type = c("norm", "basic", "perc", "bca"))

plot(bootauc_mir_ratio)

##### Bayesian model updating ####

set.seed(1001)
bay_mir_ratio <- brm(outcome ~ ., data = as.data.frame(mir_only_clean), family = bernoulli(logit),
                    prior = c(set_prior("normal(1.32,0.6931472)", class = "Intercept"),
                              set_prior("normal(0.16,0.6931472)", class = "b", coef = "mi_r_199a_3p_let_7a_5p"),
                              set_prior("normal(-0.16,0.6931472)", class = "b", coef = "mi_r_26b_5p_mi_r_142_5p"),
                              set_prior("normal(-0.18,0.6931472)", class = "b", coef = "let_7a_5p_mi_r_19b_3p"),
                              set_prior("normal(-0.07,0.6931472)", class = "b", coef = "mi_r_101_3p_mi_r_19b_3p"),
                              set_prior("normal(0.47,0.6931472)", class = "b", coef = "mi_r_93_5p_mi_r_19b_3p"),
                              set_prior("normal(-0.04,0.6931472)", class = "b", coef = "let_7a_5p_mi_r_22_3p")))

summary(bay_mir_ratio)
bay_mir_ratio %>% posterior_interval(prob = c(0.9))

# Bayesian model updating coefficients
estimates_mir_ratio <- fixef(bay_mir_ratio)
estimates_adt_mir_ratio <- cbind(rownames(estimates_mir_ratio), as.data.frame(estimates_mir_ratio))

# Predicted probabilities from coefficients obtained from Bayesian modelling
Preds_bay_mir_ratio <- Prob_LR(df = mir_only_clean, coefs = estimates_adt_mir_ratio)

# ROC AUC of the model updated by the Bayesian method

BrierScore(mir_only_clean$outcome, Preds_bay_mir_ratio$preds)
pROC::auc(mir_only_clean$outcome, Preds_bay_mir_ratio$preds)
pROC::ci.auc(mir_only_clean$outcome, Preds_bay_mir_ratio$preds)
# Calibration curve of the model updated by the Bayesian method
CalibrationCurves::val.prob.ci.2(Preds_bay_mir_ratio$preds, mir_only_clean$outcome,CL.smooth = F)

cal_bay_gg_mir_ratio <- valProbggplot(Preds_bay_mir_ratio$preds, mir_only_clean$outcome, CL.smooth = T, smooth = "loess",
                            legendloc = FALSE, size = 5, lwd.smooth = 2, lwd.ideal = 2,
                            size.d01 = 7, length.seg = 0.9, d0lab = "Normal", d1lab = "Tumour",
                            dist.label2 = 0.09)

cal_bay_gg_mir_ratio_f <- cal_bay_gg_mir_ratio$ggPlot +
  annotate("text", x = 0.5, y = 0.99999999999, size = 12, label = "B",
           fontface = 2) +
  theme_classic() + theme(title = element_text(size = 25),
                          plot.title= element_text(hjust = 0.5),
                          axis.title = element_text(size = 20),
                          axis.text = element_text(size = 15),
                          legend.position = "none")

cal_bay_gg_mir_ratio_f

### ROC AUC

roc_Bay_mir_ratio_res <- roc(mir_only_clean$outcome ~ Preds_bay_mir_ratio$preds, ci = TRUE)
# Sensitivity and specificity at the Youden's optimal cut-off
coords(roc_Bay_mir_ratio_res, x="best", input="threshold", best.method="youden")

perf_Bay_mir_ratio <- ROCR::performance(prediction(Preds_bay_mir_ratio$preds, mir_only_clean$outcome), 'sens', 'fpr')

roc_Bay_mir_ratio_data <- data.frame(
  FPR = perf_Bay_mir_ratio@x.values[[1]],
  TPR = perf_Bay_mir_ratio@y.values[[1]]
)


roc_plot_Bay_mir_ratio <- ggplot(roc_Bay_mir_ratio_data, aes(x = FPR, y = TPR)) +
  geom_line(color = "blue", size = 1.5) +  # Increase line size
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", 
              color = "gray", size = 1.2) +
  annotate("text", x = 0.8, y = 0.3, 
           label = paste("AUC: ", round(auc(mir_only_clean$outcome, 
                                            Preds_bay_mir_ratio$preds), 3)), 
           color = "blue", size = 7, vjust = -1) + 
  annotate("text", x = 0.5, y = 0.99999999999, size = 12, label = "A",
           fontface = 2) +
  geom_point(x = 1 - 0.8661417, y = 0.53125, color = "red", size = 3) +  # Add red dot
  annotate("text", x = 1 - 0.88, y = 0.5, 
           label = paste("Youden's cut-off:\n", round(0.3158782,3)), 
           color = "red", size = 5, vjust = -1) +
  labs(x = "False Positive Rate",
       y = "True Positive Rate") +
  theme_classic() +
  theme(title = element_text(size = 25),
        plot.title= element_text(hjust = 0.5),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 15))

roc_plot_Bay_mir_ratio

ggarrange(roc_plot_Bay_mir_ratio, cal_bay_gg_mir_ratio_f)
