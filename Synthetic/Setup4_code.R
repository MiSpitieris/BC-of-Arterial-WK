library(SAVE)
# model and field data 
# modelField  = readRDS("modelField_data.rds")
# model_data = modelField$synthetic$setup4$model_data
# field_data = modelField$synthetic$setup4$field_data
modelField = readRDS("modelField_data.rds") # model and field data 
model_data = modelField$synthetic$setup4$model_data # model data from lhs 
field_data = modelField$synthetic$setup4$field_data # observed data with 3 replicates
pred_data  = modelField$synthetic$setup4$pred_data # prediction data 
true_par = modelField$synthetic$setup4$true # true parameter values

RCuniq = unique(model_data[,c("R", "C")])
par(mfrow = c(2,1))
plot(1, type="n", xlim=c(0.01, 0.99), ylim=c(0, 500), xlab = "time", ylab = "Pressure (mmHg)", main = "Model data")
for (i in 1:nrow(RCuniq)) {
  ind = model_data[, "R"] == RCuniq[i,1] & model_data[, "C"] == RCuniq[i,2]
  lines(model_data[ind, "t"], model_data[ind, "y"], xlab = "time", ylab = "Pressure (mmHg)", lwd=1, col=i)
}
plot(field_data$t,field_data$y, xlab = "time", ylab = "Pressure (mmHg)", main = "Field data with 3 replicates")
par(mfrow = c(1,1))

# *********************************
#     Run Stage 1 to obtain       *
#   Bias Hyperparameters MLEs     *   
# *********************************

kc = list(multistart = 10^2) # multistart optimization for the emulator model using DiceKrig 
tic1 = Sys.time()
set.seed(123)
stage1 = SAVE(response.name="y", controllable.names=c("Q","t")
              , calibration.names=c("R", "C")
              , field.data=field_data, model.data=model_data
              , bestguess=list(R=1.2,C=1)
              , kriging.controls=kc
              , verbose = TRUE
)
toc1 = Sys.time()
(run_stage1 = toc1-tic1)
stage1


# *********************************
#                MCMC             *
# *********************************

# uniform priors on [0.5, 3] for both R and C
pr = c(uniform(var.name = "R", lower = 0.5,upper = 3)
       ,uniform(var.name = "C", lower = 0.5,upper = 3))

nsamples = 10^6
tic  = Sys.time()
stage2 = bayesfit(stage1, prior = pr, n.iter =  nsamples, nMH=20
                  , prob.prop = c(0.1,0.5)
                  , mcmcMultmle = c(5,1)
                  , verbose = TRUE
)
toc  = Sys.time()
(run.time = toc - tic)

# In sample  predictions
newDes = unique(field_data[,2:3])

BurnIn = 3000

Thin = length(BurnIn:nrow(stage2@mcmcsample))/1000  #thining

ticIS = Sys.time()
valIS =  validate(object = stage2, newdesign = newDes, n.burnin = BurnIn, n.thin = Thin) 
tocIS = Sys.time()
runIS = tocIS - ticIS
predRealIS = predictreality(object = stage2, newdesign = newDes, n.burnin = BurnIn, n.thin = Thin)

# newInd = round(seq(1,101,  length.out = 40))
# timedomains = model_data$t; Q = model_data$Q
# newDesOOS =  data.frame(t = timedomains, Q =Q)  #out of sample design
# newDesOOS = newDesOOS[newInd,]
# newDesOOS$Q = newDesOOS$Q/diff(range(Q)) 
# 
# ticOOS = Sys.time()
# valOOS = validate(object = stage2, newdesign = newDesOOS, n.burnin = BurnIn, n.thin = Thin)
# predRealOOS = predictreality(object = stage2, newdesign = newDesOOS, n.burnin = BurnIn, n.thin = Thin)
# tocOOS = Sys.time()
# runOOS = tocOOS - ticOOS

# out of sample predictions
ticOOS = Sys.time()
valOOS = validate(object = stage2, newdesign = pred_data, n.burnin = BurnIn, n.thin = Thin)
predRealOOS = predictreality(object = stage2, newdesign = pred_data, n.burnin = BurnIn, n.thin = Thin)
tocOOS = Sys.time()
runOOS = tocOOS - ticOOS


post = data.frame(stage2@mcmcsample)

#true_par = modelField$synthetic$setup4$true

model4 = list(stage1 = stage1, stage2 = stage2, nsamples = nsamples, 
              Rtrue = true_par$Rtrue, Ctrue = true_par$Ctrue, Ztrue = true_par$Ztrue,  
              model_data = model_data, field_data = field_data,
              valIS = valIS, predRealIS = predRealIS,
              valOOS = valOOS, predRealOOS = predRealOOS,
              run.time = run.time, runOOS = runOOS)

saveRDS(model4, "model4.rds")
#rm(list = ls())
