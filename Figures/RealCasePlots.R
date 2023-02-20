library(ggplot2);library(coda); library(bayestestR); library(ggplot2); library(ggpubr)
theme_set(theme_bw())

results = readRDS("results.rds")

Br1 = results$real$individual1
Br2 = results$real$individual2

poster1 = data.frame(Br1$stage2@mcmcsample)

post = poster1

# Remove burn in samples
ns = nrow(poster1)
post = post[3000:ns,]
post1 = post

# effective sample size
ef = effectiveSize(mcmc(post))

nr = nrow(post)
indR = round(seq(1,nr, length.out = ef[1]))
indC = round(seq(1,nr, length.out = ef[2]))


dfR = data.frame(R = post[indR, 1])
dfC = data.frame(C = post[indC, 2])

mapR = map_estimate(dfR$R)
mapC = map_estimate(dfC$C)
#mapR = c(mapR1, mapR2, mapR3)
meanR = mean(dfR$R)
meanC = mean(dfC$C)
#meanR = c(meanR1, meanR2, meanR3)
colhist = "slategray4"

# Marginal posteriors densities
histR = ggplot(dfR, aes(x = R))+
  geom_histogram(bins = 50, color=colhist, fill=colhist)+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(), 
        panel.grid = element_blank()) +
  scale_y_reverse()

histC = ggplot(dfC, aes(x = C))+
  geom_histogram(bins = 50, color=colhist, fill=colhist)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), 
        panel.grid = element_blank())+
  rotate()

# Joint density of calibration parameters
dfRC = data.frame(post[indC,1:2])

ID = c(0.8, 1.2, 1.6, 2.0)

joint_den = ggplot(dfRC, aes(x=R, y=C) ) +
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
  scale_fill_distiller(palette= "Spectral", direction=1) +
  scale_x_continuous(expand = c(0, 0)) +
  #scale_y_continuous(expand = c(0, 0),"ID", labels = c(as.character(ID[1:3]),"2.0"), breaks = ID) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(
    legend.position='none',
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    axis.text.x=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.x=element_blank(),
    axis.ticks.y=element_blank(), 
    panel.grid = element_blank()
  )


# Predictions
data1 = subset(Br1$field_data, select = -Q)
colnames(data1) = c("time","pressure")

tt = Br1$valOOS@newdesign$t
valOOS = as.data.frame(Br1$valOOS@validate) # out of sample predictions
valOOS$time = tt
names(valOOS) = c("bias_corrected", "tau_bc", "pure_model", "tau_pm", "bias", "bias_lower", "bias_upper" , "time")
bc = valOOS[,c(1,2,8)] # OOS bias corrected predictions
#bc$id = 'model1'
bc$lower = qnorm(0.05, mean = bc$bias_corrected, sd = bc$tau_bc)
bc$upper = qnorm(0.95, mean = bc$bias_corrected, sd = bc$tau_bc)
pm = valOOS[,c(3,4,8)] # OOS pure model predictions
#pm$id = 'model1'
pm$lower = qnorm(0.05, mean = pm$pure_model, sd = pm$tau_pm)
pm$upper = qnorm(0.95,  mean = pm$pure_model, sd = pm$tau_pm)
bias = valOOS[,c(5,6,7,8)] # OOS bias predictions
#bias$id = 'model1'

PrVal = c(60,90,120) 
colPred = "red"
Pred_plot = ggplot()+
  geom_point(data = data1, aes(x = time, y = pressure), colour = "black", size = 0.2, alpha = 0.7)+
  geom_ribbon(data = bc, aes(ymin = lower,  ymax = upper, x = time), fill = colPred, alpha = 0.3) +
  geom_line(data = bc, aes(x = time, y = bias_corrected), colour = colPred, size =0.4) +
  scale_y_continuous("pressure", labels = c(as.character(PrVal)), breaks = PrVal) +
  theme( panel.grid = element_blank())
pl1 = ggarrange(joint_den, histC, histR, Pred_plot ,  nrow = 2, ncol = 2)

#--------------------------------------------
#--------------------------------------------

# Posterior sample of R,C, marginal and noise precisions
poster2 = data.frame(Br2$stage2@mcmcsample)
post = poster2

# Remove burn in samples
ns = nrow(poster1)
post = post[3000:ns,]
post2 = post

# effective sample size
ef = effectiveSize(mcmc(post))

nr = nrow(post)
indR = round(seq(1,nr, length.out = ef[1]))
indC = round(seq(1,nr, length.out = ef[2]))

dfR = data.frame(R = post[indR, 1])
dfC = data.frame(C = post[indC, 2])

mapR = map_estimate(dfR$R)
mapC = map_estimate(dfC$C)
#mapR = c(mapR1, mapR2, mapR3)
meanR = mean(dfR$R)
meanC = mean(dfC$C)
#meanR = c(meanR1, meanR2, meanR3)

colhist = "slategray4"

# Marginal posteriors densities
histR = ggplot(dfR, aes(x = R))+
  geom_histogram(bins = 50, color=colhist, fill=colhist)+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(), 
        panel.grid = element_blank()) +
  scale_y_reverse()

histC = ggplot(dfC, aes(x = C))+
  geom_histogram(bins = 50, color=colhist, fill=colhist)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), 
        panel.grid = element_blank())+
  rotate()

# Joint density of calibration parameters
dfRC = data.frame(post[indC,1:2])

joint_den = ggplot(dfRC, aes(x=R, y=C) ) +
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
  scale_fill_distiller(palette= "Spectral", direction=1) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(
    legend.position='none',
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    axis.text.x=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.x=element_blank(),
    axis.ticks.y=element_blank()
  )

data1 = subset(Br2$field_data, select = -Q)
colnames(data1) = c("time","pressure")

tt = Br2$valOOS@newdesign$t
valOOS = as.data.frame(Br2$valOOS@validate)
valOOS$time = tt
names(valOOS) = c("bias_corrected", "tau_bc", "pure_model", "tau_pm", "bias", "bias_lower", "bias_upper" , "time")
bc = valOOS[,c(1,2,8)] # OOS bias corrected predictions
#bc$id = 'model1'
bc$lower = qnorm(0.05, mean = bc$bias_corrected, sd = bc$tau_bc)
bc$upper = qnorm(0.95, mean = bc$bias_corrected, sd = bc$tau_bc)
pm = valOOS[,c(3,4,8)] # OOS pure model predictions
#pm$id = 'model1'
pm$lower = qnorm(0.05, mean = pm$pure_model, sd = pm$tau_pm)
pm$upper = qnorm(0.95,  mean = pm$pure_model, sd = pm$tau_pm)
bias = valOOS[,c(5,6,7,8)] # OOS bias predictions
#bias$id = 'model1'

PrVal = c(60,90,120) 
colPred = "red"
Pred_plot = ggplot()+
  geom_point(data = data1, aes(x = time, y = pressure), colour = "black", size = 0.2, alpha = 0.7)+
  geom_ribbon(data = bc, aes(ymin = lower,  ymax = upper, x = time), fill = colPred, alpha = 0.3) +
  geom_line(data = bc, aes(x = time, y = bias_corrected), colour = colPred, size =0.4) +
  scale_y_continuous("pressure", labels = c(as.character(PrVal)), breaks = PrVal) +
  theme( panel.grid = element_blank())

pl2 =  ggarrange(joint_den, histC, histR, Pred_plot ,  nrow = 2, ncol = 2)
pl2



plAll = ggarrange(pl2, pl1, ncol = 2)
plAll
RealCase = plAll
#ggsave("RealCase.pdf", width = 15, height = 6, units = "cm")

# Create tables
histR + scale_x_discrete(position = "top", breaks = c(1.00,1.25,1.50,1.75))
# library(cowplot)
# ggdraw(switch_axis_position(histR + axis = 'R'))
# hR = histR+ scale_x_continuous(sec.axis = dup_axis()) 

#ggsave("RealCase.pdf", width = 15, height = 6, units = "cm")
RealCase1 = pl2
RealCase1
#ggsave("RealCase1.pdf", width = 12, height = 10, units = "cm")
RealCase2 = pl1
RealCase2
#ggsave("RealCase2.pdf", width = 10, height = 10, units = "cm")

