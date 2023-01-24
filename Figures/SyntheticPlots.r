# mod1 = readRDS("/Volumes/michais/Documents/PaperGithub/Codes/Synthetic/model1.rds")
# mod2 = readRDS("/Volumes/michais/Documents/PaperGithub/Codes/Synthetic/model2.rds")
# mod3 = readRDS("/Volumes/michais/Documents/PaperGithub/Codes/Synthetic/model3.rds")
# mod4 = readRDS("/Volumes/michais/Documents/PaperGithub/Codes/Synthetic/model4.rds")

library(coda); library(bayestestR); library(ggplot2); library(ggpubr)
theme_set(theme_classic())
setwd("/Volumes/michais/Documents/PaperGithub/Codes/Figures/")
results = readRDS("results.rds")
dsynth = results$synthOpt
mod1 = results$synthetic$setup1
mod2 = results$synthetic$setup2
mod3 = results$synthetic$setup3
mod4 = results$synthetic$setup4

poster1 = data.frame(mod1$stage2@mcmcsample)
poster2 = data.frame(mod2$stage2@mcmcsample)
poster3 = data.frame(mod3$stage2@mcmcsample)
poster4 = data.frame(mod4$stage2@mcmcsample)

ns = nrow(poster1)

# Remove burn in samples
post1 = poster1[3000:ns,]
post2 = poster2[3000:ns,]
post3 = poster3[3000:ns,]
post4 = poster4[3000:ns,]

# find effective sample sizes
ef1 = effectiveSize(mcmc(post1))
ef2 = effectiveSize(mcmc(post2))
ef3 = effectiveSize(mcmc(post3))
ef4 = effectiveSize(mcmc(post4))
min(c(ef1, ef2, ef3, ef4))

# Define new efss for the marginal posterior densities plot
ef1 = ef2 = ef3 = ef4 = c(1e3,1e3)*7
nr = nrow(post1)
names(mod4)[4] = "Rtrue"
names(mod4)[5] = "Ctrue"
names(mod4)[6] = "Ztrue"

# resultsAll = results
# resultsAll$synthetic$setup1 = mod1
# resultsAll$synthetic$setup2 = mod2
# resultsAll$synthetic$setup3 = mod3
#resultsAll$synthetic$setup4 = mod4
#saveRDS(resultsAll, "/Volumes/michais/Documents/PaperGithub/Codes/Figures/results.rds")
# R = R1 (= Ztrue) + R2 (= R)
Rr1 = mod1$Rtrue + mod1$Ztrue
Rr2 = mod2$Rtrue + mod2$Ztrue
Rr3 = mod3$Rtrue + mod3$Ztrue
Rr4 = mod4$Rtrue + mod4$Ztrue
Rrall = c(Rr1, Rr2, Rr3,Rr4)

indR1 = round(seq(1,nr, length.out = ef1[1]))
indR2 = round(seq(1,nr, length.out = ef2[1]))
indR3 = round(seq(1,nr, length.out = ef3[1]))
indR4 = round(seq(1,nr, length.out = ef4[1]))

dfR1 = data.frame(R = post1[indR1, 1], id = paste("Rtrue1=",Rr1, sep = ""))
dfR2 = data.frame(R = post2[indR2, 1], id = paste("Rtrue2=",Rr2, sep = ""))
dfR3 = data.frame(R = post3[indR3, 1], id = paste("Rtrue3=",Rr3, sep = ""))
dfR4 = data.frame(R = post4[indR4, 1], id = paste("Rtrue4=",Rr4, sep = ""))

mapR1 = map_estimate(dfR1$R)
mapR2 = map_estimate(dfR2$R)
mapR3 = map_estimate(dfR3$R)
mapR4 = map_estimate(dfR4$R)
mapR = c(mapR1, mapR2, mapR3, mapR4)

meanR1 = mean(dfR1$R)
meanR2 = mean(dfR2$R)
meanR3 = mean(dfR3$R)
meanR4 = mean(dfR4$R)
meanR = c(meanR1, meanR2, meanR3, meanR4)

dfRall = rbind(dfR1, dfR2, dfR3, dfR4)
colnames(dfRall)[2] = "posteriors"

cols = c("#1B9E77","#E41A1C", "#E69F00", "#377EB8")
Cr1 = mod1$Ctrue
Cr2 = mod2$Ctrue
Cr3 = mod3$Ctrue
Cr4 = mod4$Ctrue
Crall = c(Cr1, Cr2, Cr3, Cr4)

indC1 = round(seq(1,nr, length.out = ef1[2]))
indC2 = round(seq(1,nr, length.out = ef2[2]))
indC3 = round(seq(1,nr, length.out = ef3[2]))
indC4 = round(seq(1,nr, length.out = ef4[2]))

dfC1 = data.frame(C = post1[indC1, 2], id = paste("Ctrue1=",Cr1, sep = ""))
dfC2 = data.frame(C = post2[indC2, 2], id = paste("Ctrue2=",Cr2, sep = ""))
dfC3 = data.frame(C = post3[indC3, 2], id = paste("Ctrue3=",Cr3, sep = ""))
dfC4 = data.frame(C = post4[indC4, 2], id = paste("Ctrue4=",Cr4, sep = ""))

dfCall = rbind(dfC1, dfC2, dfC3, dfC4)
colnames(dfCall)[2] = "posteriors"

mapC1 = map_estimate(dfC1$C)
mapC2 = map_estimate(dfC2$C)
mapC3 = map_estimate(dfC3$C)
mapC4 = map_estimate(dfC4$C)

meanC1 = mean(dfC1$C)
meanC2 = mean(dfC2$C)
meanC3 = mean(dfC3$C)
meanC4 = mean(dfC4$C)

mapC = c(mapC1, mapC2, mapC3, mapC4)
meanC = c(meanC1, meanC2, meanC3, meanC4)

# Marginal Densities

Rest1 = data.frame(R = c(mapR1, meanR1,
                         dsynth$mod1optWK2$Rest[1],
                         dsynth$mod1optWK3$Rest[1] + dsynth$mod1optWK3$Zest[1]),
                   y = rep(0,4),
                   estimate = c("map", "mean", "WK2", "WK3"))
Rest2 = data.frame(R = c(mapR2, meanR2,
                         dsynth$mod2optWK2$Rest[1],
                         dsynth$mod2optWK3$Rest[1] + dsynth$mod2optWK3$Zest[1]),
                   y = rep(0,4),
                   estimate = c("map", "mean", "WK2", "WK3"))
Rest3 = data.frame(R = c(mapR3, meanR3,
                         dsynth$mod3optWK2$Rest[1],
                         dsynth$mod3optWK3$Rest[1] + dsynth$mod3optWK3$Zest[1]),
                   y = rep(0,4),
                   estimate = c("map", "mean", "WK2", "WK3"))
Rest4 = data.frame(R = c(mapR4, meanR4,
                         dsynth$mod4optWK2$Rest[1],
                         dsynth$mod4optWK3$Rest[1] + dsynth$mod4optWK3$Zest[1]),
                   y = rep(0,4),
                   estimate = c("map", "mean", "WK2", "WK3"))

shape = c("map" = 10, "mean" = 12, "WK2" = 3, "WK3" = 2)
sh = c(10,12,3,2)

plpl = ggplot(dfR1, aes(x = R, y = ..density..))+
  geom_density(col = "black", size = 1, alpha = 0.6, fill = cols[1]) + 
  geom_point(data = Rest1, aes(x=R, y=y, shape = estimate), size = 2)+
  geom_vline(aes(xintercept = Rr1, linetype  = "true"),size=0.5, colour = "black") +
  scale_shape_manual(values=shape)+
  scale_linetype_manual(name = "true value", values = c(true = 1))+
  theme(legend.title = element_blank())+ 
  xlim(0.5,3)+ ylim(0,16) 

plR1 = ggplot(dfR1, aes(x = R, y = ..density..))+
  geom_density(col = "black", size = 1, alpha = 0.6, fill = cols[1]) + 
  geom_point(data = Rest1, aes(x=R, y=y, shape = estimate), size = 2)+
  geom_vline(aes(xintercept = Rr1, linetype  = "true"),size=0.5, colour = "black") +
  scale_shape_manual(values=shape)+
  scale_linetype_manual(name = "true value", values = c(true = 1))+
  xlim(0.5,3)+ ylim(0,16) +
  theme(legend.position = "none",
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank()
        ,axis.ticks.y = element_blank()
  ) +  ggtitle("R posterior")

plR2 = ggplot(dfR2, aes(x = R, y = ..density..))+
  geom_density(col = "black", size = 1, alpha = 0.6, fill = cols[2]) + 
  geom_point(data = Rest2, aes(x=R, y=y), size = 2, shape = sh)+
  geom_vline(aes(xintercept = Rr2),size=0.5, colour = "black") +
  #scale_shape_manual(values=shape)+
  #scale_linetype_manual(name = "real value", values = c(real = 1))+
  xlim(0.5,3)+ ylim(0,16) + 
  theme(legend.position = "none",
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank()) 
#+ggtitle("Setup 2")

plR3 = ggplot(dfR3, aes(x = R, y = ..density..))+
  geom_density(col = "black", size = 1, alpha = 0.6, fill = cols[3]) + 
  geom_point(data = Rest3, aes(x=R, y=y), size = 2, shape = sh)+
  geom_vline(aes(xintercept = Rr3),size=0.5, colour = "black") +
  #scale_shape_manual(values=shape)+
  #scale_linetype_manual(name = "real value", values = c(real = 1))+
  xlim(0.5,3)+ ylim(0,16) +
  theme(legend.position = "none",
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank())
#+ggtitle("Setup 3")

plR4 = ggplot(dfR4, aes(x = R, y = ..density..))+
  geom_density(col = "black", size = 1, alpha = 0.6, fill = cols[4]) + 
  geom_point(data = Rest4, aes(x=R, y=y), size = 2, shape = sh)+
  geom_vline(aes(xintercept = Rr4),size=0.5, colour = "black") +
  #scale_shape_manual(values=shape)+
  #scale_linetype_manual(name = "real value", values = c(real = 1))+
  xlim(0.5,3) + ylim(0,16) +
  theme(#legend.position = "none",
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank())
#+ggtitle("Setup 4")

Rmarg = ggarrange(plR1, plR2, plR3, plR4, ncol = 1)


Cest1 = data.frame(C = c(mapC1, meanC1,
                         dsynth$mod1optWK2$Cest[1],
                         dsynth$mod1optWK3$Cest[1]),
                   y = rep(0,4),
                   estimate = c("map", "mean", "WK2", "WK3"))
Cest2 = data.frame(C = c(mapC2, meanC2,
                         dsynth$mod2optWK2$Cest[1],
                         dsynth$mod2optWK3$Cest[1]),
                   y = rep(0,4),
                   estimate = c("map", "mean", "WK2", "WK3"))
Cest3 = data.frame(C = c(mapC3, meanC3,
                         dsynth$mod3optWK2$Cest[1],
                         dsynth$mod3optWK3$Cest[1]),
                   y = rep(0,4),
                   estimate = c("map", "mean", "WK2", "WK3"))
Cest4 = data.frame(C = c(mapC4, meanC4,
                         dsynth$mod4optWK2$Cest[1],
                         dsynth$mod4optWK3$Cest[1]),
                   y = rep(0,4),
                   estimate = c("map", "mean", "WK2", "WK3"))

shape = c("map" = 10, "mean" = 12, "WK2" = 3, "WK3" = 2)

plC1 = ggplot(dfC1, aes(x = C, y = ..density..))+
  geom_density(col = "black", size = 1, alpha = 0.6, fill = cols[1]) + 
  geom_point(data = Cest1, aes(x=C, y=y), size = 2, shape = sh)+
  geom_vline(aes(xintercept = Cr1),size=0.5, colour = "black") +
  #scale_shape_manual(values=shape)+
  #scale_linetype_manual(name = "real value", values = c(real = 1))+
  xlim(0.5,3)+ ylim(0,16) +
  theme(legend.position = "none",
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank()
        ,axis.ticks.y = element_blank()
  ) + 
  ggtitle("C posterior") +
  annotate("text", x = 2.8, y = 14, label = "Setup 1")

plC2 = ggplot(dfC2, aes(x = C, y = ..density..))+
  geom_density(col = "black", size = 1, alpha = 0.6, fill = cols[2]) + 
  geom_point(data = Cest2, aes(x=C, y=y), size = 2, shape = sh)+
  geom_vline(aes(xintercept = Cr2),size=0.5, colour = "black") +
  #scale_shape_manual(values=shape)+
  #scale_linetype_manual(name = "real value", values = c(real = 1))+
  xlim(0.5,3)+ ylim(0,16) + 
  theme(legend.position = "none",
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank())+
  annotate("text", x = 2.8, y = 14, label = "Setup 2")

plC3 = ggplot(dfC3, aes(x = C, y = ..density..))+
  geom_density(col = "black", size = 1, alpha = 0.6, fill = cols[3]) + 
  geom_point(data = Cest3, aes(x=C, y=y), size = 2, shape = sh)+
  geom_vline(aes(xintercept = Cr3),size=0.5, colour = "black") +
  #scale_shape_manual(values=shape)+
  #scale_linetype_manual(name = "real value", values = c(real = 1))+
  xlim(0.5,3)+ ylim(0,16) +
  theme(legend.position = "none",
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank()) +
  annotate("text", x = 2.8, y = 14, label = "Setup 3")

plC4 = ggplot(dfC4, aes(x = C, y = ..density..))+
  geom_density(col = "black", size = 1, alpha = 0.6, fill = cols[4]) + 
  geom_point(data = Cest4, aes(x=C, y=y), size = 2, shape = sh)+
  geom_vline(aes(xintercept = Cr4),size=0.5, colour = "black") +
  #scale_shape_manual(values=shape)+
  #scale_linetype_manual(name = "real value", values = c(real = 1))+
  xlim(0.5,3) + ylim(0,16) +
  theme(#legend.position = "none",
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank()) +
  annotate("text", x = 2.8, y = 14, label = "Setup 4")
Cmarg = ggarrange(plC1, plC2, plC3, plC4, ncol = 1)
RCmarg = ggarrange(Rmarg, Cmarg, ncol = 2)

#ggarrange(Rmarg, Cmarg, nrow=2)
# Joint densities
get_density <- function(x, y, ...) {# https://slowkow.com/notes/ggplot2-color-by-density/
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}
#library(viridis)
setup1Joint = data.frame(R = dfR1$R, C = dfC1$C)
setup1Joint$density = get_density(setup1Joint$R, setup1Joint$C, n = 2*1e3)

setup2Joint = data.frame(R = dfR2$R, C = dfC2$C)
setup2Joint$density = get_density(setup2Joint$R, setup2Joint$C, n = 2*1e3)

setup3Joint = data.frame(R = dfR3$R, C = dfC3$C)
setup3Joint$density = get_density(setup3Joint$R, setup3Joint$C, n = 2*1e3)


setup4Joint = data.frame(R = dfR4$R, C = dfC4$C)
setup4Joint$density = get_density(setup4Joint$R, setup4Joint$C, n = 2*1e3)

map1 = data.frame(R = mapR1, C = mapC1)
map2 = data.frame(R = mapR2, C = mapC2)
map3 = data.frame(R = mapR3, C = mapC3)
map4 = data.frame(R = mapR4, C = mapC4)

mean1 = data.frame(R = meanR1, C = meanC1)
mean2 = data.frame(R = meanR2, C = meanC2)
mean3 = data.frame(R = meanR3, C = meanC3)
mean4 = data.frame(R = meanR4, C = meanC4)
realRC1 = data.frame(R = Rr1, C = Cr1)
realRC2 = data.frame(R = Rr2, C = Cr2)
realRC3 = data.frame(R = Rr3, C = Cr3)
realRC4 = data.frame(R = Rr4, C = Cr4)

#joint1
joint1 = ggplot(setup1Joint) + 
  geom_point(aes(x = R, y = C, color = density), size = 0.05) + 
  scale_color_distiller(palette= "Spectral")+
  geom_point(data = realRC1, aes(x = R, y=C), col = "red", size = 3,shape=23)+
  #geom_point(data = map1,  aes(x = R, y = C), shape = 1)+
  #geom_point(data = mean1,  aes(x = R, y = C), shape = 2)+
  theme(legend.position = c(0.8,0.8)
        ,legend.key.height = unit(0.35, 'cm')
        ,legend.key.width  = unit(0.35, 'cm')
        ,legend.title = element_text(size = 9)
        ,legend.text = element_text(size = 5)
  )+
  xlim(0.5,3)+
  ylim(0.5,3) +annotate("text", x = 2.5, y = 1, label = "Setup 1")

#joint2
joint2 = ggplot(setup2Joint) + 
  geom_point(aes(x = R, y = C, color = density), size = 0.05) + 
  scale_color_distiller(palette= "Spectral")+
  geom_point(data = realRC2, aes(x = R, y=C), col = "red", size = 3,shape=23)+
  theme(legend.position = c(0.8,0.8)
        ,legend.key.height = unit(0.35, 'cm')
        ,legend.key.width  = unit(0.35, 'cm')
        ,legend.title = element_text(size = 9)
        ,legend.text = element_text(size = 5)
  )+
  xlim(0.5,3)+
  ylim(0.5,3) +annotate("text", x = 2.5, y = 1, label = "Setup 2")

#col2 = c("real" = "red")
#joint3
joint3 = ggplot(setup3Joint) + 
  geom_point(aes(x = R, y = C, color = density), size = 0.05) + 
  scale_color_distiller(palette= "Spectral")+
  geom_point(data = realRC3, aes(x = R, y=C), col = "red", size = 3,shape=23)+
  theme(legend.position = c(0.8,0.8)
        ,legend.key.height = unit(0.35, 'cm')
        ,legend.key.width  = unit(0.35, 'cm')
        ,legend.title = element_text(size = 9)
        ,legend.text = element_text(size = 5)
  )+
  xlim(0.5,3)+
  ylim(0.5,3) +annotate("text", x = 2.5, y = 1, label = "Setup 3")

#joint4
joint4 = ggplot(setup4Joint) + 
  geom_point(aes(x = R, y = C, color = density), size = 0.05) + 
  scale_color_distiller(palette= "Spectral")+
  geom_point(data = realRC4, aes(x = R, y=C),col = "red", size = 3,shape=23, show.legend = TRUE)+
  theme(legend.position = c(0.8,0.8)
        ,legend.key.height = unit(0.35, 'cm')
        ,legend.key.width  = unit(0.35, 'cm')
        ,legend.title = element_text(size = 9)
        ,legend.text = element_text(size = 5)
  )+
  xlim(0.5,3)+
  ylim(0.5,3)+ annotate("text", x = 2.5, y = 1, label = "Setup 4")


Joint = ggarrange(joint1, joint2, joint3, joint4, nrow = 2, ncol = 2)
Joint
ggsave("Joint.pdf", plot = Joint, width = 12, height = 12, units = "cm")
leg = get_legend(plpl, position = "bottom")
RCmarg = ggarrange(Rmarg, Cmarg, ncol = 2, legend = "bottom", common.legend = TRUE
                   , legend.grob = get_legend(plpl, position = "bottom"))
#Joint = ggarrange(joint1, joint2, joint3, joint4, nrow = 1)
ggsave("RCmarg.pdf", plot = RCmarg, width = 18, height = 18, units = "cm")
postSynthMat = ggarrange(RCmarg,  Joint, nrow = 2, legend = "bottom", common.legend = TRUE
                         , legend.grob = get_legend(plpl, position = "bottom")
                         , heights = c(4,1.3))
postSynthMat

#ggsave("postSynthMat2.pdf", plot = postSynthMat, width = 18, height = 20, units = "cm")

#---------------------------------
# All marginal posteriors in one plot
dfR1 = data.frame(R = post1[indR1, 1], id = "Setup 1")
dfR2 = data.frame(R = post2[indR2, 1], id = "Setup 2")
dfR3 = data.frame(R = post3[indR3, 1], id = "Setup 3")
dfR4 = data.frame(R = post4[indR4, 1], id = "Setup 4")

dfRall = rbind(dfR1, dfR2, dfR3, dfR4)
colnames(dfRall)[2] = "posteriors"

dfC1 = data.frame(C = post1[indC1, 2], id = "Setup 1")
dfC2 = data.frame(C = post2[indC2, 2], id = "Setup 2")
dfC3 = data.frame(C = post3[indC3, 2], id = "Setup 3")
dfC4 = data.frame(C = post4[indC4, 2], id = "Setup 4")

dfCall = rbind(dfC1, dfC2, dfC3, dfC4)
colnames(dfCall)[2] = "posteriors"

plR = ggplot(dfRall, aes(x = R, y=..density.., fill=posteriors)) +
  geom_density(color = "black", size = 1, alpha = 0.6) + 
  geom_vline(aes(xintercept = mapR1, colour="map1"), linetype="dashed", size=1) +
  geom_vline(aes(xintercept = mapR2, colour="map2"), linetype="dashed", size=1) +
  geom_vline(aes(xintercept = mapR3, colour="map3"), linetype="dashed", size=1) +
  geom_vline(aes(xintercept = mapR4, colour="map4"), linetype="dashed", size=1) +
  scale_fill_manual(values = cols)+ 
  #scale_x_continuous("R", labels = as.character(c(Rrall, seq(0.5,3,by=0.5)[-2])), breaks = c(Rrall, seq(0.5,3,by=0.5)[-2]), limits = c(0.5,3))+
  xlim(0.5,3) +
  scale_y_continuous(breaks=NULL, limits = c(0,16))+
  scale_color_manual(name = "estimates", values = c("map1" = cols[1], "map2" = cols[2], "map3" = cols[3], "map4" = cols[4]))+
  theme_classic(base_size = 8) +
  theme(legend.position = c(0.8, 0.6),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 6),
        axis.text.x = element_text(face = c('plain','bold', 'bold', 'bold', rep('plain',4))),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.key.size = unit(0.5, "cm"),
        legend.key.width = unit(0.5,"cm")
  ) 

plC = ggplot(dfCall, aes(x = C, y=..density.., fill=posteriors)) + 
  geom_density(color = "black", size = 1, alpha = 0.6) + 
  geom_vline(aes(xintercept = mapC1, colour="map1"), linetype="dashed", size=1) +
  geom_vline(aes(xintercept = mapC2, colour="map2"), linetype="dashed", size=1) +
  geom_vline(aes(xintercept = mapC3, colour="map3"), linetype="dashed", size=1) +
  geom_vline(aes(xintercept = mapC4, colour="map4"), linetype="dashed", size=1) +
  scale_fill_manual(values = cols) +
  #scale_x_continuous("C", labels = as.character(c(Crall, seq(0.5,3,by=0.5)[-2])), breaks = c(Crall, seq(0.5,3,by=0.5)[-2]), limits = c(0.5,3))+
  xlim(0.5,3) +
  scale_y_continuous(breaks=NULL, limits = c(0,16)) +
  scale_color_manual(name = "estimates", values = c("map1" = cols[1], "map2" = cols[2], "map3" = cols[3], "map4" = cols[4]))+
  theme_classic(base_size = 8) +
  theme(legend.position = c(0.8, 0.6),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 6),
        axis.text.x = element_text(face = c('plain','bold', 'bold', 'plain', 'bold', rep('plain',4))),
        axis.title.x = element_text(size = 16),
        # axis.title.y = element_text(size = 16),
        axis.title.y = element_blank(),
        legend.key.size = unit(0.5, "cm"),
        legend.key.width = unit(0.5,"cm")
  ) 

#length(as.character(c(Crall, seq(0.5,3,by=0.5)[-2])))
postRC = ggarrange(plR,plC)
ggsave("postRC.pdf", plot = postRC, width = 16, height = 8, units = "cm")

Opt100 = results$Opt100
WK2_1 = Opt100$setUp1$WK2
names(WK2_1)[1:2] = c("R", "C")
WK2_1$setup = rep("1",100)

WK2_2 = Opt100$setUp2$WK2
names(WK2_2)[1:2] = c("R", "C")
WK2_2$setup = rep("2",100)

WK2_3 = Opt100$setUp3$WK2
names(WK2_3)[1:2] = c("R", "C")
WK2_3$setup = rep("3",100)

WK2_4 = Opt100$setUp4$WK2
names(WK2_4)[1:2] = c("R", "C")
WK2_4$setup = rep("4",100)

WK2all = rbind(WK2_1, WK2_2, WK2_3, WK2_4)
WK2BoxPlot = reshape2::melt(WK2all, id = "setup")
R1rall = c(mod1$Ztrue, mod2$Ztrue, mod3$Ztrue, mod4$Ztrue)
Points = data.frame(variable = c(rep("R", 4), rep("C", 4)), setup = factor(rep(1:4, 2)), value = c(Rrall + R1rall, Crall))

boxWK2R = ggplot(data = WK2BoxPlot[WK2BoxPlot$variable == "R",], aes(y = value, x = setup, fill = setup)) + 
  geom_jitter(aes(color = setup), size = 0.5)+
  geom_boxplot(fill = NA) +
  ylim(0,2) +
  scale_colour_manual(values = cols, aesthetics = "colour")+
  theme_bw()+
  ylab("R")+
  annotate("text", x = 2.5, y = 2, label = "WK2 optimization", fontface = 2)+
  geom_point(data = Points[Points$variable == "R",],
             aes(y = value, x = setup), shape=23, col ="red", size = 3, fill = NA)+
  theme(panel.grid = element_blank()
        ,legend.position = "none")

boxWK2C = ggplot(data = WK2BoxPlot[WK2BoxPlot$variable == "C",], aes(y = value, x = setup, fill = setup)) + 
  geom_jitter(aes(color = setup), size = 0.5)+
  geom_boxplot(fill = NA) +
  ylim(0,2) +
  scale_colour_manual(values = cols, aesthetics = "colour")+
  theme_bw()+
  ylab("C")+
  geom_point(data = Points[Points$variable == "C",],
             aes(y = value, x = setup), shape=23, col ="red", size = 3, fill = NA)+
  theme(panel.grid = element_blank()
        ,legend.position = "none"
  ) 
boxWK2 = ggarrange(boxWK2R, boxWK2C)

#------------------------------
# Wk3 boxplots
WK3_1 = Opt100$setUp1$WK3
names(WK3_1)[1:3] = c("R", "C","Z")
WK3_1$setup = rep("1",100)

WK3_2 = Opt100$setUp2$WK3
names(WK3_2)[1:3] = c("R", "C","Z")
WK3_2$setup = rep("2",100)

WK3_3 = Opt100$setUp3$WK3
names(WK3_3)[1:3] = c("R", "C","Z")
WK3_3$setup = rep("3",100)

WK3_4 = Opt100$setUp4$WK3
names(WK3_4)[1:3] = c("R", "C","Z")
WK3_4$setup = rep("4",100)

WK3all = rbind(WK3_1, WK3_2, WK3_3, WK3_4)
WK3BoxPlot = reshape2::melt(WK3all, id = "setup")
R1rall = c(mod1$Ztrue, mod2$Ztrue, mod3$Ztrue, mod4$Ztrue)
R2rall = c(mod1$Rtrue, mod2$Rtrue, mod3$Rtrue, mod4$Rtrue)
Points = data.frame(variable = c(rep("R", 4), rep("C", 4), rep("Z", 4)), setup = factor(rep(1:4, 3)), value = c(R2rall, Crall, R1rall))

boxWK3R2 = ggplot(data = WK3BoxPlot[WK3BoxPlot$variable == "R",], aes(y = value, x = setup, fill = setup)) + 
  geom_jitter(aes(color = setup), size = 0.5)+
  geom_boxplot(fill = NA) +
  ylim(0,2) +
  scale_colour_manual(values = cols, aesthetics = "colour")+
  theme_bw()+
  ylab(expression(R[2]))+
  geom_point(data = Points[Points$variable == "R",],
             aes(y = value, x = setup), shape=23, col ="red", size = 3, fill = NA)+
  theme(panel.grid = element_blank()
        ,legend.position = "none") 

boxWK3R2
boxWK3C = ggplot(data = WK3BoxPlot[WK3BoxPlot$variable == "C",], aes(y = value, x = setup, fill = setup)) + 
  geom_jitter(aes(color = setup), size = 0.5)+
  geom_boxplot(fill = NA) +
  ylim(0,2) +
  scale_colour_manual(values = cols, aesthetics = "colour")+
  theme_bw()+
  ylab("C")+
  geom_point(data = Points[Points$variable == "C",],
             aes(y = value, x = setup), shape=23, col ="red", size = 3, fill = NA)+
  theme(panel.grid = element_blank()
        ,legend.position = "none") 

boxWK3R1 = ggplot(data = WK3BoxPlot[WK3BoxPlot$variable == "Z",], aes(y = value, x = setup, fill = setup)) + 
  geom_jitter(aes(color = setup), size = 0.5)+
  geom_boxplot(fill = NA) +
  ylim(0,2) +
  scale_colour_manual(values = cols, aesthetics = "colour")+
  theme_bw()+
  ylab(expression(R[1]))+
  geom_point(data = Points[Points$variable == "Z",],
             aes(y = value, x = setup), shape=23, col ="red", size = 3, fill = NA)+
  theme(panel.grid = element_blank()
        ,legend.position = "none"
  )+ annotate("text", x = 2.5, y = 2, label = "WK3 optimization", fontface = 2)


boxWK3 = ggarrange(boxWK3R1, boxWK3R2, boxWK3C, nrow = 1)
ggsave("boxPlotWK3.pdf", plot = boxWK3, units = "cm", width = 18, height = 8)
ggsave("boxPlotWK2.pdf", plot = boxWK2, units = "cm", width = 12, height = 8)

### Plot in Section 2
WK2 = results$WK2and3$WK2
WK3 = results$WK2and3$WK3
colsWK = c("WK2" = "red", "WK3" = "grey")
dftext = data.frame(text = "Z", x = 0.45, y = 158)
WK23 = ggplot() +
  geom_line(data = WK3, aes(x = time, y = pressure, group = variable, color = "WK3")) +
  geom_line(data = WK2, aes(x = time, y = pressure, color = "WK2"), size = 1) +
  geom_segment(aes(x = 0.4, y = 165, xend = 0.4, yend = 150),
               arrow = arrow(length = unit(0.5, "cm")))+
  geom_text(data = dftext, aes(x=x, y=y, label = "R[1]"),parse = T) +
  scale_color_manual(name = "model", values = colsWK)+ 
  theme(legend.position = c(0.8, 0.6),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 10),
        #axis.text.x = element_text(face = c('plain','bold', 'bold', rep('plain',5))),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.key.size = unit(1, "cm"),
        legend.key.width = unit(1,"cm")
  ) 
WK23 

ggsave("WK23.pdf", plot = WK23, units = "cm", width = 16, height = 8)

### All R marginals and C marrginals in single plots
R1_ci = round(quantile(dfR1$R, probs = c(0.05,0.95)),3)
R2_ci = round(quantile(dfR2$R, probs = c(0.05,0.95)),3)
R3_ci = round(quantile(dfR3$R, probs = c(0.05,0.95)),3)

C1_ci = round(quantile(dfC1$C, probs = c(0.05,0.95)),3)
C2_ci = round(quantile(dfC2$C, probs = c(0.05,0.95)),3)
C3_ci = round(quantile(dfC3$C, probs = c(0.05,0.95)),3)

CIs = data.frame(lower = c( R1_ci[1], R2_ci[1], R3_ci[1] 
                            ,C1_ci[1], C2_ci[1], C3_ci[1]),
                 upper = c( R1_ci[2], R2_ci[2], R3_ci[2] 
                            ,C1_ci[2], C2_ci[2], C3_ci[2])
                 
)
rownames(CIs) = c(paste0("R", 1:3),paste0("C", 1:3))
#CIs
mapEst = data.frame(mapR = mapR, mapC = mapC)
meanEst = data.frame(meanR = meanR, meanC = meanC)
#--------------------------------
#--------------------------------
### Bias-corrected, Pure model & Bias predictions
mm1 = mod1
mm2 = mod2
mm3 = mod3
mm4 = mod4

data1 = subset(mm1$field_data, select = -Q)
colnames(data1) = c("pressure", "time")
data2 = subset(mm2$field_data, select = -Q)
colnames(data2) = c("pressure", "time")
data3 = subset(mm3$field_data, select = -Q)
colnames(data3) = c("pressure", "time")
data4 = subset(mm4$field_data, select = -Q)
colnames(data4) = c("pressure", "time")

tt = mm1$valOOS@newdesign$t
valOOS1 = as.data.frame(mm1$valOOS@validate)
valOOS1$time = tt
names(valOOS1) = c("bias_corrected", "tau_bc", "pure_model", "tau_pm", "bias", "bias_lower", "bias_upper" , "time")
bc1 = valOOS1[,c(1,2,8)] # OOS bias corrected predictions
bc1$id = 'Setup 1'
bc1$lower = qnorm(0.05, mean = bc1$bias_corrected, sd = bc1$tau_bc)
bc1$upper = qnorm(0.95, mean = bc1$bias_corrected, sd = bc1$tau_bc)
pm1 = valOOS1[,c(3,4,8)] # OOS pure model predictions
pm1$id = 'Setup 1'
pm1$lower = qnorm(0.05, mean = pm1$pure_model, sd = pm1$tau_pm)
pm1$upper = qnorm(0.95,  mean = pm1$pure_model, sd = pm1$tau_pm)
bias1 = valOOS1[,c(5,6,7,8)] # OOS bias predictions
bias1$id = 'Setup 1'

tt = mm2$valOOS@newdesign$t
valOOS2 = as.data.frame(mm2$valOOS@validate) 
valOOS2$time = tt
names(valOOS2) = c("bias_corrected", "tau_bc", "pure_model", "tau_pm", "bias", "bias_lower", "bias_upper" , "time")
bc2 = valOOS2[,c(1,2,8)] # OOS bias corrected predictions
bc2$id = 'Setup 2'
bc2$lower = qnorm(0.05, mean = bc2$bias_corrected, sd = bc2$tau_bc)
bc2$upper = qnorm(0.95, mean = bc2$bias_corrected, sd = bc2$tau_bc)
pm2 = valOOS2[,c(3,4,8)] # OOS pure model predictions
pm2$id = 'Setup 2'
pm2$lower = qnorm(0.05, mean = pm2$pure_model, sd = pm2$tau_pm)
pm2$upper = qnorm(0.95,  mean = pm2$pure_model, sd = pm2$tau_pm)
bias2 = valOOS2[,c(5,6,7,8)] # OOS bias predictions
bias2$id = 'Setup 2'

valOOS3 = as.data.frame(mm3$valOOS@validate) 
valOOS3$time = tt
names(valOOS3) = c("bias_corrected", "tau_bc", "pure_model", "tau_pm", "bias", "bias_lower", "bias_upper" , "time")

bc3 = valOOS3[,c(1,2,8)] # OOS bias corrected predictions
bc3$id = 'Setup 3'
bc3$lower = qnorm(0.05, mean = bc3$bias_corrected, sd = bc3$tau_bc)
bc3$upper = qnorm(0.95, mean = bc3$bias_corrected, sd = bc3$tau_bc)
pm3 = valOOS3[,c(3,4,8)] # OOS pure model predictions
pm3$id = 'Setup 3'
pm3$lower = qnorm(0.05, mean = pm3$pure_model, sd = pm3$tau_pm)
pm3$upper = qnorm(0.95,  mean = pm3$pure_model, sd = pm3$tau_pm)
bias3 = valOOS3[,c(5,6,7,8)] # OOS bias predictions
bias3$id = 'Setup 3'

tt = mm4$valOOS@newdesign$t
valOOS4 = as.data.frame(mm4$valOOS@validate)
valOOS4$time = tt
names(valOOS4) = c("bias_corrected", "tau_bc", "pure_model", "tau_pm", "bias", "bias_lower", "bias_upper" , "time")
bc4 = valOOS4[,c(1,2,8)] # OOS bias corrected predictions
bc4$id = 'Setup 4'
bc4$lower = qnorm(0.05, mean = bc4$bias_corrected, sd = bc4$tau_bc)
bc4$upper = qnorm(0.95, mean = bc4$bias_corrected, sd = bc4$tau_bc)
pm4 = valOOS4[,c(3,4,8)] # OOS pure model predictions
pm4$id = 'Setup 4'
pm4$lower = qnorm(0.05, mean = pm4$pure_model, sd = pm4$tau_pm)
pm4$upper = qnorm(0.95,  mean = pm4$pure_model, sd = pm4$tau_pm)
bias4 = valOOS4[,c(5,6,7,8)] # OOS bias predictions
bias4$id = 'Setup 4'

bcAll = rbind(bc1, bc2, bc3, bc4)
pmAll = rbind(pm1, pm2, pm3, pm4)
biasAll = rbind(bias1, bias2, bias3, bias4)

colsWK = c("WK2" = "red", "WK3" = "grey")
cols = c("Setup 1" = "#1B9E77","Setup 2" = "#E41A1C","Setup 3" = "#E69F00","Setup 4" = "#377EB8")

BiasCor_plot= ggplot()+
  geom_ribbon(data = bc1, aes(ymin = lower,  ymax = upper, x = time), fill = cols[1], alpha = 0.3)+
  geom_ribbon(data = bc2, aes(ymin = lower,  ymax = upper, x = time), fill = cols[2], alpha = 0.3)+
  geom_ribbon(data = bc3, aes(ymin = lower,  ymax = upper, x = time), fill = cols[3], alpha = 0.3)+
  geom_ribbon(data = bc4, aes(ymin = lower,  ymax = upper, x = time), fill = cols[4], alpha = 0.3)+
  geom_line(data = bc1, aes(x = time, y = bias_corrected, colour = "Setup 1"), size =1)+
  geom_line(data = bc2, aes(x = time, y = bias_corrected, colour = "Setup 2"), size =1)+
  geom_line(data = bc3, aes(x = time, y = bias_corrected, colour = "Setup 3"), size =1)+
  geom_line(data = bc4, aes(x = time, y = bias_corrected, colour = "Setup 4"), size =1)+
  geom_point(data = data1, aes(x = time, y = pressure), colour = "#1B9E77", size = 0.8)+
  geom_point(data = data2, aes(x = time, y = pressure), colour = "#E41A1C", size = 0.8)+
  geom_point(data = data3, aes(x = time, y = pressure), colour = "#E69F00", size = 0.8)+
  geom_point(data = data4, aes(x = time, y = pressure), colour = "#377EB8", size = 0.8)+
  scale_color_manual(name = "Bias-corrected\npredictions", values = cols)+
  theme_classic(base_size = 8) +
  ylab("pressure")+
  theme(legend.position = c(0.8, 0.8),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        # axis.title.y = element_blank(),
        legend.key.size = unit(0.5, "cm"),
        legend.key.width = unit(0.5,"cm")
  )+ ylim(-25,210)

PureModel_plot= ggplot()+
  geom_ribbon(data = pm1, aes(ymin = lower,  ymax = upper, x = time), fill = cols[1], alpha = 0.3)+
  geom_ribbon(data = pm2, aes(ymin = lower,  ymax = upper, x = time), fill = cols[2], alpha = 0.3)+
  geom_ribbon(data = pm3, aes(ymin = lower,  ymax = upper, x = time), fill = cols[3], alpha = 0.3)+
  geom_ribbon(data = pm4, aes(ymin = lower,  ymax = upper, x = time), fill = cols[4], alpha = 0.3)+
  geom_line(data = pm1, aes(x = time, y = pure_model, colour = "Setup 1"), size =1)+
  geom_line(data = pm2, aes(x = time, y = pure_model, colour = "Setup 2"), size =1)+
  geom_line(data = pm3, aes(x = time, y = pure_model, colour = "Setup 3"), size =1)+
  geom_line(data = pm4, aes(x = time, y = pure_model, colour = "Setup 4"), size =1)+
  geom_point(data = data1, aes(x = time, y = pressure), colour = "#1B9E77", size = 0.8)+
  geom_point(data = data2, aes(x = time, y = pressure), colour = "#E41A1C", size = 0.8)+
  geom_point(data = data3, aes(x = time, y = pressure), colour = "#E69F00", size = 0.8)+
  geom_point(data = data4, aes(x = time, y = pressure), colour = "#377EB8", size = 0.8)+
  scale_color_manual(name = "Pure-model\npredictions", values = cols)+
  theme_classic(base_size = 8) +
  ylab("pressure")+
  theme(legend.position = c(0.8, 0.8),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        # axis.title.y = element_blank(),
        legend.key.size = unit(0.3, "cm"),
        legend.key.width = unit(0.3,"cm"),
        legend.background=element_blank()
  )+ ylim(-25,210)

Bias_plot= ggplot()+
  geom_ribbon(data = bias1, aes(ymin = bias_lower,  ymax = bias_upper, x = time), fill = cols[1], alpha = 0.3)+
  geom_ribbon(data = bias2, aes(ymin = bias_lower,  ymax = bias_upper, x = time), fill = cols[2], alpha = 0.3)+
  geom_ribbon(data = bias3, aes(ymin = bias_lower,  ymax = bias_upper, x = time), fill = cols[3], alpha = 0.3)+
  geom_ribbon(data = bias4, aes(ymin = bias_lower,  ymax = bias_upper, x = time), fill = cols[4], alpha = 0.3)+
  geom_line(data = bias1, aes(x = time, y = bias, colour = "Setup 1"), size =1)+
  geom_line(data = bias2, aes(x = time, y = bias, colour = "Setup 2"), size =1)+
  geom_line(data = bias3, aes(x = time, y = bias, colour = "Setup 3"), size =1)+
  geom_line(data = bias4, aes(x = time, y = bias, colour = "Setup 4"), size =1)+
  #geom_point(data = data1, aes(x = time, y = pressure), colour = "#1B9E77")+
  #geom_point(data = data2, aes(x = time, y = pressure), colour = "#E41A1C")+
  #geom_point(data = data3, aes(x = time, y = pressure), colour = "#377EB8")+
  scale_color_manual(name = "Bias", values = cols)+
  theme_classic(base_size = 8) +
  ylab("bias")+
  # ylim(...) +
  theme(legend.position = c(0.8, 0.8),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        # axis.title.y = element_blank(),
        legend.key.size = unit(0.3, "cm"),
        legend.key.width = unit(0.3,"cm")
  ) + ylim(-25,210)

BiasCor_plot
ggsave("BiasCor.pdf", plot = BiasCor_plot, width = 14, height = 8, units = "cm")
PureModel_plot
ggsave("PureModel.pdf", plot = PureModel_plot, width = 9, height = 5, units = "cm")
Bias_plot
ggsave("Bias.pdf", plot = Bias_plot, width = 9, height = 5, units = "cm")

PureModelandBias = ggarrange(PureModel_plot, Bias_plot, nrow = 1)
ggsave("PureModelandBias.pdf", plot = PureModelandBias, width = 19, height = 8, units = "cm")


BiasCor_plot2= ggplot()+
  geom_ribbon(data = bc1, aes(ymin = lower,  ymax = upper, x = time), fill = cols[1], alpha = 0.3)+
  geom_ribbon(data = bc2, aes(ymin = lower,  ymax = upper, x = time), fill = cols[2], alpha = 0.3)+
  geom_ribbon(data = bc3, aes(ymin = lower,  ymax = upper, x = time), fill = cols[3], alpha = 0.3)+
  geom_ribbon(data = bc4, aes(ymin = lower,  ymax = upper, x = time), fill = cols[4], alpha = 0.3)+
  geom_line(data = bc1, aes(x = time, y = bias_corrected, colour = "Setup 1"), size =1)+
  geom_line(data = bc2, aes(x = time, y = bias_corrected, colour = "Setup 2"), size =1)+
  geom_line(data = bc3, aes(x = time, y = bias_corrected, colour = "Setup 3"), size =1)+
  geom_line(data = bc4, aes(x = time, y = bias_corrected, colour = "Setup 4"), size =1)+
  geom_point(data = data1, aes(x = time, y = pressure), colour = "#1B9E77", size = 0.8)+
  geom_point(data = data2, aes(x = time, y = pressure), colour = "#E41A1C", size = 0.8)+
  geom_point(data = data3, aes(x = time, y = pressure), colour = "#E69F00", size = 0.8)+
  geom_point(data = data4, aes(x = time, y = pressure), colour = "#377EB8", size = 0.8)+
  scale_color_manual(name = "Bias-corrected\npredictions", values = cols)+
  theme_classic(base_size = 8) +
  ylab("pressure")+
  theme(legend.position = c(0.8, 0.8),
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        axis.title.x = element_text(size = 12, face = "plain"),
        axis.title.y = element_text(size = 12, face = "plain"),
        # axis.title.y = element_blank(),
        legend.key.size = unit(0.5, "cm"),
        legend.key.width = unit(0.5,"cm"),
        plot.title = element_text( face = "bold", size = 14)
  )+ ylim(-25,210) + ggtitle("Bias corrected")

PureModel_plot2= ggplot()+
  geom_ribbon(data = pm1, aes(ymin = lower,  ymax = upper, x = time), fill = cols[1], alpha = 0.3)+
  geom_ribbon(data = pm2, aes(ymin = lower,  ymax = upper, x = time), fill = cols[2], alpha = 0.3)+
  geom_ribbon(data = pm3, aes(ymin = lower,  ymax = upper, x = time), fill = cols[3], alpha = 0.3)+
  geom_ribbon(data = pm4, aes(ymin = lower,  ymax = upper, x = time), fill = cols[4], alpha = 0.3)+
  geom_line(data = pm1, aes(x = time, y = pure_model, colour = "Setup 1"), size =1)+
  geom_line(data = pm2, aes(x = time, y = pure_model, colour = "Setup 2"), size =1)+
  geom_line(data = pm3, aes(x = time, y = pure_model, colour = "Setup 3"), size =1)+
  geom_line(data = pm4, aes(x = time, y = pure_model, colour = "Setup 4"), size =1)+
  geom_point(data = data1, aes(x = time, y = pressure), colour = "#1B9E77", size = 0.8)+
  geom_point(data = data2, aes(x = time, y = pressure), colour = "#E41A1C", size = 0.8)+
  geom_point(data = data3, aes(x = time, y = pressure), colour = "#E69F00", size = 0.8)+
  geom_point(data = data4, aes(x = time, y = pressure), colour = "#377EB8", size = 0.8)+
  scale_color_manual(name = "Pure-model\npredictions", values = cols)+
  theme_classic(base_size = 8) +
  ylab("pressure")+
  theme(legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        axis.title.x = element_text(size = 12, face = "plain"),
        axis.title.y = element_text(size = 12, face = "plain"),
        # axis.title.y = element_blank(),
        legend.key.size = unit(0.5, "cm"),
        legend.key.width = unit(0.5,"cm"),
        plot.title = element_text( face = "bold", size = 14)
  )+ ylim(-25,210) + ggtitle("Pure model")

Bias_plot2= ggplot()+
  geom_ribbon(data = bias1, aes(ymin = bias_lower,  ymax = bias_upper, x = time), fill = cols[1], alpha = 0.3)+
  geom_ribbon(data = bias2, aes(ymin = bias_lower,  ymax = bias_upper, x = time), fill = cols[2], alpha = 0.3)+
  geom_ribbon(data = bias3, aes(ymin = bias_lower,  ymax = bias_upper, x = time), fill = cols[3], alpha = 0.3)+
  geom_ribbon(data = bias4, aes(ymin = bias_lower,  ymax = bias_upper, x = time), fill = cols[4], alpha = 0.3)+
  geom_line(data = bias1, aes(x = time, y = bias, colour = "Setup 1"), size =1)+
  geom_line(data = bias2, aes(x = time, y = bias, colour = "Setup 2"), size =1)+
  geom_line(data = bias3, aes(x = time, y = bias, colour = "Setup 3"), size =1)+
  geom_line(data = bias4, aes(x = time, y = bias, colour = "Setup 4"), size =1)+
  #geom_point(data = data1, aes(x = time, y = pressure), colour = "#1B9E77")+
  #geom_point(data = data2, aes(x = time, y = pressure), colour = "#E41A1C")+
  #geom_point(data = data3, aes(x = time, y = pressure), colour = "#377EB8")+
  scale_color_manual(name = "Bias", values = cols)+
  theme_classic(base_size = 8) +
  # ylab("bias")+
  # ylim(...) +
  theme(legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        axis.title.x = element_text(size = 12, face = "plain"),
        axis.title.y = element_blank(),
        # axis.title.y = element_blank(),
        legend.key.size = unit(0.5, "cm"),
        legend.key.width = unit(0.5,"cm"),
        plot.title = element_text( face = "bold", size = 14)
  )+ ylim(-25,210) + ggtitle("Bias")

BiasCor_plot
ggsave("BiasCor2.pdf", plot = BiasCor_plot2, width = 14, height = 8, units = "cm")
PureModel_plot
ggsave("PureMode2l.pdf", plot = PureModel_plot2, width = 9, height = 5, units = "cm")
Bias_plot
ggsave("Bias2.pdf", plot = Bias_plot2, width = 9, height = 5, units = "cm")

PureModelandBias2 = ggarrange(PureModel_plot2, Bias_plot2, nrow = 1)
ggsave("PureModelandBias2.pdf", plot = PureModelandBias2, width = 18, height = 6.5, units = "cm")
