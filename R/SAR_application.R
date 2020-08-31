devtools::install_github("frycast/rsar")
library(rsar)
library(tidyverse)
library(colorspace)
library(plotly)
library(raster)
library(e1071)
library(caret)
source("R/SAR_application_helpers.R")
source("R/logchol_helpers.R")

###########################################################################
# EXTRACTIONS  ------------------------------------------------------------
###########################################################################

## Just looking at ssq

# kms_df <- get_kms_result(format = "df")
# 
# kms_df %>% 
#   filter(band == "cc") %>% 
#   mutate(m = factor(m), d = factor(d)) %>% 
#   ggplot(data = .) +
#   geom_line(aes(y = res, x = k, linetype = m, colour = m), size = 1) +
#   facet_grid(d~.)
# 
# kms_df %>% 
#   filter(band == "cc") %>% 
#   mutate(m = factor(m), d = factor(d)) %>% 
#   ggplot(data = .) +
#   geom_line(aes(y = res, x = k, linetype = d, colour = d), size = 1) +
#   facet_grid(m~.)

## Looking at Rand index and BIC

kms_ran <- get_kms_result(format = "df", result_fun = rand_index, labs_in_rf = T) 
kms_bic <- get_kms_result(format = "df", result_fun = kmsBIC)
# kms_bic3 <- get_kms_result(format = "df", result_fun = kmsBIC3)
# saveRDS(kms_ran, "SAR_app_data/MG/kms_ran.Rds")
# saveRDS(kms_bic, "SAR_app_data/MG/kms_bic.Rds")
# saveRDS(kms_bic3, "SAR_app_data/MG/kms_bic3.Rds")

kms_ran <- readRDS("SAR_app_data/MG/kms_ran.Rds")
kms_bic <- readRDS("SAR_app_data/MG/kms_bic.Rds")
#kms_bic3 <- readRDS("SAR_app_data/MG/kms_bic3.Rds")

## Maximising Rand index:
## These results cause us to choose d = 9, m = 1 and band = cc

kms_ran %>% 
  filter(k <= 8, d == 10) %>% 
  mutate(lag = factor(m), d = factor(d), band = toupper(band)) %>% 
  ggplot(data = .) +
  geom_line(aes(y = res, x = k, linetype = lag, colour = lag), size = 1) +
  facet_grid(band~.) +
  ylab("adjusted Rand index (R)") +
  xlab("number of clusters (k)")
kms_ran %>% 
  filter(band == "cc", k == 2) %>% 
  mutate(lag = factor(m)) %>% 
  ggplot(data = .) +
  geom_line(aes(y = res, x = d, linetype = lag, colour = lag), size = 1) +
  ylab("adjusted Rand index (R)") +
  xlab("patch size (d)")

## BIC plots

kms_bic %>% 
  filter(d == 9, k >= 2, m == 1, band == "cc") %>% 
  mutate(lag = factor(m), d = factor(d), band = toupper(band)) %>% 
  ggplot(data = .) +
  geom_line(aes(y = res, x = k), size = 1) +
  ylab("Bayes Information Criterion (BIC)") +
  xlab("number of clusters (k)")

## The highest Rand index, for all k

max_ran <- filter(kms_ran, res == max(res)); max_ran
write_kms_plot(plot_only = F, save_name = "minRand_cc.tif", 
               band = max_ran$band, k = max_ran$k, d = max_ran$d, m = max_ran$m)

## The lowest BIC in each band, over all k

min_vh <- filter(kms_bic, band == "vh", d == 9, m == 1) %>% filter(res == min(res))
min_vv <- filter(kms_bic, band == "vv", d == 9, m == 1) %>% filter(res == min(res))
min_cc <- filter(kms_bic, band == "cc", d == 9, m == 1) %>% filter(res == min(res)); min_cc
kms <- write_kms_plot(plot_only = T, save_name = "minBIC_vh.tif", 
                      band = "vh", k = min_vh$k, d = min_vh$d, m = min_vh$m)
kms <- write_kms_plot(plot_only = T, save_name = "minBIC_vv.tif", 
                      band = "vv", k = min_vv$k, d = min_vv$d, m = min_vv$m)
kms <- write_kms_plot(plot_only = F, save_name = "minBIC_cc.tif", 
                      band = "cc", k = min_cc$k, d = min_cc$d, m = min_cc$m)

## Pairs plot and 3D plotly scatter plot showing all clusters

cov_chol_plotly(kms$cluster, pairs_plot = T)

## Below is just a silly way to get k = 2, d = 9 and m = 1 because we shouldn't use BIC to choose m
# 
# min_vh2 <- filter(kms_bic, band == "vh", k == 2, d == 9) %>% filter(res == min(res))
# min_vv2 <- filter(kms_bic, band == "vv", k == 2, d == 9) %>% filter(res == min(res))
# min_cc2 <- filter(kms_bic, band == "cc", k == 2, d == 9) %>% filter(res == min(res))
# kms2 <- write_kms_plot(plot_only = T, save_name = "minBIC_vh_k2.tif", 
#                       band = "vh", k = min_vh2$k, d = min_vh2$d, m = min_vh2$m)
# kms2 <- write_kms_plot(plot_only = T, save_name = "minBIC_vv_k2.tif", 
#                       band = "vv", k = min_vv2$k, d = min_vv2$d, m = min_vv2$m,
#                       col = rev(colorspace::rainbow_hcl(min_vv2$k)))
kms2 <- write_kms_plot(plot_only = T, save_name = "minBIC_cc_k2.tif", 
                      band = "cc", k = 2, d = 9, m = 1,
                      col = colorspace::rainbow_hcl(2))

# a colour blind friendly palette gray:

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

## Plot the labels (GDVs) over the clusters

labs <- get_pix(band = "labs")
labsb <- matrix_to_brick(labs)[[1]]
dat_lab <- as.data.frame(labsb, xy=TRUE)
colnames(dat_lab) <- c("x","y","cluster")
dat_lab$cluster <- factor(dat_lab$cluster + 1)

r2 <- get_kms_classes_brick(kms2, d = 9)[[1]]
dat <- as.data.frame(r2, xy=TRUE)
colnames(dat) <- c("x","y","cluster")
dat$cluster <- factor(dat$cluster)
dat <- filter(dat, x <= 140.592, y >= -37.463)
pdf(file="min_Rand_with_GDV_map.pdf",width=5,height=4)
ggplot() +
  geom_raster(data = dat[dat$cluster == 2,], 
              aes(x, y), fill = cbPalette[1]) +
  geom_raster(data = dat_lab[dat_lab$cluster == 2,], 
              aes(x, y), fill = cbPalette[3], alpha = 0.4) +
  theme_minimal()
dev.off()

# Get the L clusters that overlap with GDVs the most

L <- 4
labsp <- get_pix_patched(band = "labs", d = 9)
GDV_table <- table(labsp > 0.8, kms$cluster) # Notice that non-GDV is clearer
# --------------------------------
#           1    2    3    4    5
# FALSE 1885 2772 3104 4203 3690
# TRUE   808  630   16    7  198
# 
#           6    7    8    9   10
# FALSE  107 3670  779 1358 1502
# TRUE     0    1    2   65    0
#
#          11   12   13   14   15
# FALSE 2844 1260 1049 3104 2978
# TRUE     0  377    0    0   71
# --------------------------------

cluster_GDV <- names(rev(sort(GDV_table[2,]))[1:L])
ep <- rev(sort(GDV_table[2,]/GDV_table[1,])) # Notice agreement with ratios
ep <- as.data.frame(ep); ep
write.csv(ep, "empirical_prob_overlap.csv")

## Get the clusters that have 0 probability of overlap

cluster_nGDV <- names(GDV_table[2,])[GDV_table[2,] < 10]

# agree <- matrix_to_brick(SAR_matrix(
#   as.matrix(kms$cluster %in% cluster_GDV), attr_src = labsp))
# raster::plot(agree[[1]])
# # raster::plot(matrix_to_brick(labsp))

## Plot the L clusters that overlap most with GDVs

r <- get_kms_classes_brick(kms, d = 9)[[1]]
dat <- as.data.frame(r, xy=TRUE)
colnames(dat) <- c("x","y","cluster")
dat <- filter(dat, x <= 140.592, y >= -37.463, cluster %in% cluster_GDV)
dat$cluster <- factor(dat$cluster)
pdf(file="test.pdf",width=5,height=4)
ggplot(data = dat) +
  geom_raster(aes(x, y, fill = cluster)) +
  theme_minimal()
dev.off()

## Plot the J clusters that overlap least with GDVs

r <- get_kms_classes_brick(kms, d = 9)[[1]]
dat <- as.data.frame(r, xy=TRUE)
colnames(dat) <- c("x","y","cluster")
dat <- filter(dat, x <= 140.592, y >= -37.463, cluster %in% cluster_nGDV)
dat$cluster <- factor(dat$cluster)
pdf(file="test.pdf",width=5,height=4)
ggplot(data = dat) +
  geom_raster(aes(x, y, fill = cluster)) +
  theme_minimal()
dev.off()

## Pairs plot and 3D plotly scatter plot showing just top L clusters

samp_GDV <- which(kms$cluster %in% cluster_GDV)
cov_chol_plotly(kms$cluster, subset = samp_GDV, pairs_plot = T, sample_size = 3000)

## Plot the above along with the other clusters as well

kmst <- kms
samp_GDV <- which(kmst$cluster %in% cluster_GDV)
kmst$cluster[samp_GDV] <- 1 #kmst$cluster[samp_GDV] + 100 # create separation
kmst$cluster[-samp_GDV] <- 2
cov_chol_plotly(kmst$cluster, pairs_plot = T, sample_size = 3000)


## Now do the opposite: highlight the clusters that are clearly non-GDV

kmst <- kms
samp_GDV <- which(kmst$cluster %in% cluster_nGDV)
kmst$cluster[samp_GDV] <- 1
kmst$cluster[-samp_GDV] <- 2
cov_chol_plotly(kmst$cluster, pairs_plot = T)



## Now clip out only pixels that correspond to GDVs exactly, 

labsp <- get_pix_patched(band = "labs", d = 9)
kmst <- kms
kmst$cluster[as.logical(labsp)] <- 1
kmst$cluster[!as.logical(labsp)] <- 2
cov_chol_plotly(kmst$cluster)

###########################################################################
# COMPARING WITH SARGDE ---------------------------------------------------
###########################################################################

sargde <- get_pix(band = "sargde")
sargdep <- get_pix_patched(band = "sargde", d = 9)

## Plotly scatter plot for high, med and low SARGDE

qnt <- quantile(sargdep)
sargdeq <- sargdep
sargdeq[T] <- 3 # Upper 25%
sargdeq[sargdep <= qnt[4]] <- 2 # Middle 50% 
sargdeq[sargdep < qnt[2]] <- 1 # Lower 25%
cov_chol_plotly(sargdeq, sample_size = length(sargdeq), pairs_plot = T)

## After seeing this above, refer back to this plot below, which shows all
## the clearly non-GDV classes. From this, it's easy to argue that
## the SARGDE index above shouldn't move into the upper or lower
## parts of the space too much. So, we can argue that SARGDE should be treated
## more as a range than a threshold. However, not that high values of
## SARGDE do correspond to the middle of the range.

kmst <- kms
samp_GDV <- which(kmst$cluster %in% cluster_nGDV)
kmst$cluster[samp_GDV] <- 1
kmst$cluster[-samp_GDV] <- 2
cov_chol_plotly(kmst$cluster, sample_size = length(kmst$cluster), pairs_plot = T)

## So based on the conclusion above, I can produce an upper triangle
## matrix of accuracies based on a 2D grid of upper and lower bound
## values for the SARGDE index

range <- seq(5, max(sargdep), by = 0.1)
obj <- vector(mode = "numeric", length = length(range)^2); i <- 1
for (lower in range) {
  for (upper in range) {
    if (lower >= upper) {acc[i] <- 0} else {
      class_GDV <- (sargdep >= lower) & (sargdep <= upper)
      obj[i] <- sum((labsp[class_GDV] > 0.8) + (labsp[!class_GDV] < 0.2))
    }
    i <- i + 1
  }
}
obj <- obj/length(labsp)
obj_mat <- matrix(obj, ncol = length(range))

fig <- plot_ly(z = ~obj_mat)
fig <- fig %>% add_surface()
fig

maxim <- which(obj_mat == max(obj_mat), arr.ind = TRUE)

# Based on the above plot:
lower_thresh <- range[maxim[2]]; lower_thresh
upper_thresh <- range[maxim[1]]; upper_thresh

sargde_pos <- as.matrix(as.integer(sargdep >= lower_thresh & sargdep <= upper_thresh))
sargde_posb <- matrix_to_brick(SAR_matrix(sargde_pos, attr_src = labsp))
raster::plot(sargde_posb)

sum(as.integer(sargde_pos) == (labsp > 0.8))/length(labsp) #93.3% accuracy

## Find SARGDE threshold with highest Rand index with k = 2 kms
## and also compute all the confusion matrices for this comparison

thr_kms <- SARGDE_threshold_grid(kms2$cluster)
thr_kms$threshold # 10.3

## That was also the threshold that maximised the accuracy

acc_kms <- sapply(thr_kms$results, function(x){x$conf$overall["Accuracy"]})
names(thr_kms$results)[which.max(acc_kms)]

## Find SARGDE threshold with highest Rand index on GDE labels
## and also compute all the confusion matrices for this comparison

labsp <- get_pix_patched(band = "labs", d = 9)
thr_gde <- SARGDE_threshold_grid(labsp)
thr_gde$threshold # 10.4

## But the accuracy is actually threshold is higher than that

acc_gde <- sapply(thr_gde$results, function(x){x$conf$overall["Accuracy"]})
acc_gde[which.max(acc_gde)] # 79% accuracy
thr_acc_gde <- names(thr_gde$results)[which.max(acc_gde)]

## Comparing the above two Rand maximising thresholds -- very similar visually

thr_kms$threshold # 10.3
thr_gde$threshold # 10.4

sargde_high <- as.matrix(as.integer(sargdep > thr_kms$threshold))
sargde_highb <- matrix_to_brick(SAR_matrix(sargde_high, attr_src = sargdep))
raster::plot(sargde_highb[[1]], main = thr_kms$threshold)

sargde_high <- as.matrix(as.integer(sargdep > thr_gde$threshold))
sargde_highb <- matrix_to_brick(SAR_matrix(sargde_high, attr_src = sargdep))
raster::plot(sargde_highb[[1]], main = thr_gde$threshold)

## Looking at the accuracy threshold -- much more aggressive

sargde_high <- as.matrix(as.integer(sargdep > thr_acc_gde))
sargde_highb <- matrix_to_brick(SAR_matrix(sargde_high, attr_src = sargdep))
raster::plot(sargde_highb[[1]], main = thr_acc_gde)

## Confusion matrices for SARGDE and kms, for k = 2 and higher k cases

mclust::adjustedRandIndex(sargde_high, kms2$cluster) # High
mclust::adjustedRandIndex(sargde_high, kms$cluster)  # Low
caret::confusionMatrix(factor(sargde_high), factor(kms2$cluster-1)) # Good

## ANOVA justifies using higher k to capture more variability
## and it also shows that SARGDE agrees more with kms2 than labs

summary(lm(sargdep ~ factor(kms$cluster))) # 53% (!)
summary(lm(sargdep ~ factor(kms2$cluster))) # 20%
summary(lm(sargdep ~ factor(labsp))) # 10%

## The following table indicates a lot of entropy in class 5,
## but it also shows that classes 1, 2, 9 and 12 are very SARGDE-high 

table(sargde_high, kms$cluster)

## I wonder where classes 5 and 9 sit in the log-chol space



