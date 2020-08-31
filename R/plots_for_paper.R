library(rsar)
library(tidyverse)
library(colorspace)
library(plotly)
library(raster)
library(e1071)
library(caret)
library(ggpubr)
source("R/SAR_application_helpers.R")
source("R/logchol_helpers.R")

# remotes::install_github("coolbutuseless/ggpattern")
# library(ggpattern)

## Suggestion to consider for updating / better figures -----------------

#library(GGally)
#ggpairs(iris, aes(colour = Species, alpha = 0.4))

## Figure "gdv_overlap.pdf" {fig:overlap-gdv} ---------------------------

kms_bic <- get_kms_result(format = "df", D = 9, result_fun = kmsBIC)
min_cc <- filter(kms_bic, band == "cc", d == 9, m == 1) %>% 
  filter(res == min(res))
kms <- get_kms(band = "cc", d = 9, m = 1)[[min_cc$k]]
L <- 4
labsp <- get_pix_patched(band = "labs", d = 9)
GDV_table <- table(labsp > 0.8, kms$cluster)
cluster_GDV <- names(rev(sort(GDV_table[2,]))[1:L])
classes <- kms$cluster
samp_GDV <- which(classes %in% cluster_GDV)
classes[samp_GDV] <- 1
classes[-samp_GDV] <- 2
cov_chol <- get_cov_chol(band = "cc", d = 9, m = 1)
plotlydf <- data.frame(cbind(cov_chol, classes))
samp <- sample(1:nrow(plotlydf), 5000)
plotlydf <- plotlydf[samp,,drop = F]
colnames(plotlydf) <- c("c_1","c_2","c_3","class")
cls <- as.vector(unique(plotlydf[,"class"]))
pdf(file="results/gdv_overlap.pdf",width=5,height=4)
pairs(plotlydf[,c("c_1","c_2","c_3")], 
      col = plotlydf[,"class"],
      oma=c(3,3,3,10))
par(xpd=TRUE)
legend("right", legend = c("< 5%", "> 5%"), fill= rev(cls))
dev.off()

## Figure "sargde_qnt.pdf" {fig:overlap-gdv} --------------------------------------------

sargdep <- get_pix_patched(band = "sargde", d = 9)
qnt <- quantile(sargdep)
sargdeq <- sargdep
sargdeq[T] <- 3 # Upper 25%
sargdeq[sargdep <= qnt[4]] <- 2 # Middle 50% 
sargdeq[sargdep < qnt[2]] <- 1 # Lower 25%

plotlydf2 <- data.frame(cbind(cov_chol, sargdeq))
# samp <- sample(1:nrow(plotlydf2), 3000)
# plotlydf2 <- plotlydf2[samp,,drop = F]
colnames(plotlydf2) <- c("c_1","c_2","c_3","class")
plotlydf2 <- plotlydf2[order(plotlydf2$class),]
samp <- sample(1:nrow(plotlydf2), 7000)
plotlydf2 <- plotlydf2[samp,,drop = F]
cls <- as.vector(unique(plotlydf2[,"class"]))
pdf(file="results/sargde_qnt.pdf",width=5,height=4)
pairs(plotlydf2[,c("c_1","c_2","c_3")], 
      col = plotlydf2[,"class"],
      oma=c(3,3,3,10))
par(xpd=TRUE)
legend("right", legend = c("low","med","high"), fill= c(1,2,3))
dev.off()

## Figure "ran.pdf" -----------------------------------------------------

kms_ran <- readRDS("SAR_app_data/MG/kms_ran.Rds")

p1 <- kms_ran %>% 
  filter(k <= 8, d == 10) %>% 
  mutate(lag = factor(m), d = factor(d), band = toupper(band)) %>% 
  ggplot(data = .) +
  geom_line(aes(y = res, x = k, linetype = lag, colour = lag), size = 1) +
  facet_grid(band~.) +
  ylab("adjusted Rand index (R)") +
  xlab("number of clusters (k)") +
  theme(legend.position = "none")
p2 <- kms_ran %>% 
  filter(band == "cc", k == 2) %>% 
  mutate(lag = factor(m)) %>% 
  ggplot(data = .) +
  geom_line(aes(y = res, x = d, linetype = lag, colour = lag), size = 1) +
  ylab("adjusted Rand index (R)") +
  xlab("patch size (d)")

#pdf(file="results/ran.pdf",width=5,height=4)
ggarrange(p1,p2, ncol = 2, nrow = 1)
#dev.off()

## Figure "BIC.pdf" ------------------------------------------------------
kms_bic <- readRDS("SAR_app_data/MG/kms_bic.Rds")

#pdf(file="results/BIC.pdf",width=5,height=4)
kms_bic %>% 
  filter(d == 9, k >= 2, m == 1, band == "cc") %>% 
  mutate(lag = factor(m), d = factor(d), band = toupper(band)) %>% 
  ggplot(data = .) +
  geom_line(aes(y = res, x = k), size = 1) +
  ylab("BIC inspired measure (D)") +
  xlab("number of clusters (k)")
#dev.off()

# kms 15 class map, SARGDE colour map, GDV map -------------------------------

# kms 15 class map 

r <- get_kms_classes_brick(kms, d = 9)[[1]]
dat <- as.data.frame(r, xy=TRUE)
colnames(dat) <- c("x","y","cluster")
dat <- filter(dat, x <= 140.592, y >= -37.463)
dat$cluster <- factor(dat$cluster)

labsp <- get_pix_patched(band = "labs", d = 9)
# labsb <- matrix_to_brick(labsp)[[1]]
# dat_lab <- as.data.frame(labsb, xy=TRUE)
# colnames(dat_lab) <- c("x","y","lab")
# dat_lab$lab <- factor(dat_lab$lab + 1)
# cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
GDV_table <- table(labsp > 0.5, kms$cluster)
ep <- rev(sort(GDV_table[2,]/GDV_table[1,]))
dat$P <- ep[as.character(dat$cluster)]
g1 <- ggplot(data = dat) +
  geom_raster(aes(x, y, fill = P)) +
  theme_minimal() +
  scale_fill_gradientn(colours = c(1,3))
  #geom_polygon_pattern(...) # would need to import vector polygons for this
  #geom_raster(data = dat_lab[dat_lab$lab == 2,], 
  #            aes(x, y), fill = cbPalette[3], alpha = 0.4)


# SARGDE colour map

sargdep <- matrix_to_brick(get_pix_patched(band = "sargde", d = 9))
dat2 <- as.data.frame(sargdep, xy=TRUE)
colnames(dat2) <- c("x","y","SARGDE")
dat2 <- filter(dat2, x <= 140.592, y >= -37.463)
#pdf(file="test.pdf",width=5,height=4)
g2 <- ggplot(data = dat2) +
  geom_raster(aes(x, y, fill = SARGDE)) +
  theme_minimal() +
  scale_fill_gradientn(colours = c(1,3))
#dev.off()

# GDV map

labs <- get_pix(band = "labs")
labsb <- matrix_to_brick(labs)
dat_lab <- as.data.frame(labsb, xy=TRUE)
dat_lab <- filter(dat_lab, x <= 140.592, y >= -37.463)
colnames(dat_lab) <- c("x","y","GDV")
dat_lab$GDV <- factor(dat_lab$GDV)

g3 <- ggplot(data = dat_lab) +
  geom_raster(aes(x, y, fill = GDV)) +
  theme_minimal() +
  scale_fill_manual(values = c("1","3"))

pdf(file="results/three_maps.pdf", width=5, height=10)
ggarrange(g1, g3, g2, ncol = 1, nrow = 3, labels= "AUTO")
dev.off()


