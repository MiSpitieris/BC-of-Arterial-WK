library(SAVE)
# model and field data 
modelField = readRDS("modelField_data.rds")
model_data = modelField$real$individual2$model_data
field_data = modelField$real$individual2$field_data
pred_data  = modelField$real$individual2$pred_data

RCuniq = unique(model_data[,c("R", "C")])

par(mfrow = c(2,1))
plot(1, type="n", xlim=c(0.01, 0.99), ylim=c(0, 500), xlab = "time", ylab = "Pressure (mmHg)", main = "Model data")
for (i in 1:nrow(RCuniq)) {
  ind = model_data[, "R"] == RCuniq[i,1] & model_data[, "C"] == RCuniq[i,2]
  lines(model_data[ind, "t"], model_data[ind, "y"], xlab = "time", ylab = "Pressure (mmHg)", lwd=1, col=i)
}

plot(field_data$t,field_data$y, xlab = "time", ylab = "Pressure (mmHg)", main = "Field data with 3 replicates")
par(mfrow = c(1,1))


kc = list(multistart = 10^2,nugget = sqrt(.Machine$double.eps))

tic1 = Sys.time()
set.seed(123)
stage1 = SAVE(response.name="y", controllable.names=c("Q","t")
              ,calibration.names=c("R", "C")
              ,field.data=field_data, model.data=model_data
              ,bestguess=list(R=1.2,C=1.2)
              ,kriging.controls=kc
)
toc1 = Sys.time()
(run_stage1 = toc1-tic1)
stage1
pr = c(uniform(var.name = "R", lower = 0.5,upper = 2.4)
       ,uniform(var.name = "C", lower = 0.5,upper = 2.4))

nsamples = 10^6
tic  = Sys.time()
stage2 = bayesfit(stage1, prior = pr, n.iter =  nsamples, nMH=20
                  , prob.prop = c(0.1,0.5)
                  , mcmcMultmle = c(5,1)
                  ,verbose = TRUE
)
toc  = Sys.time()
(run.time = toc - tic)

# In sample  predictions
newDes = unique(field_data[,1:2])

BurnIn = 3000

Thin = length(BurnIn:nrow(stage2@mcmcsample))/1000

ticIS = Sys.time()
valIS =  validate(object = stage2, newdesign = newDes, n.burnin = BurnIn, n.thin = Thin)
tocIS = Sys.time()
runIS = tocIS - ticIS
predRealIS = predictreality(object = stage2, newdesign = newDes, n.burnin = BurnIn, n.thin = Thin)
# nrow(d)
# newInd = round(seq(1,nrow(d),  length.out = 80))
# nrow(newDes)
# newInd = round(seq(1,nrow(newDes),  length.out = 80))
# newDesOOS =  data.frame(t = timedomains, Q = inflow) #out of sample design
# newDesOOS = newDesOOS[newInd,]
# newDesOOS$Q = newDesOOS$Q/diff(range(Q)) 

ticOOS = Sys.time()
valOOS = validate(object = stage2, newdesign = pred_data, n.burnin = BurnIn, n.thin = Thin)
predRealOOS = predictreality(object = stage2, newdesign = pred_data, n.burnin = BurnIn, n.thin = Thin)
tocOOS = Sys.time()
runOOS = tocOOS - ticOOS


real_dat_results = list(stage1 = stage1, stage2 = stage2, nsamples = nsamples, 
                        model_data = model_data, field_data = field_data, run_time = run.time,
                        valIS = valIS, predRealIS = predRealIS,
                        valOOS = valOOS, predRealOOS = predRealOOS,
                        run.time = run.time, runOOS = runOOS)
saveRDS(real_dat_results, "BrachialRes2.rds")
# #model1 = readRDS("model1.rds")
# rm(list = ls())
# 
#model1  = readRDS("/Volumes/michais/Documents/Rcode/model1.rds")
