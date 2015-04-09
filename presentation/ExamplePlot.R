# setwd('~/Google Drive/Longcode/presentation/')
setwd('C:/Users/msbcr563/Google Drive/Longcode/presentation/')

# Historical production and consumption -----------------------------------
Data = read.csv('ExampleDataSubset.csv',header = TRUE)

prod = Data$Production + Data$Import - Data$Export 
storage = Data$Storage
# netwithdraw = Data$NetWithdraw
prodAstorage = prod + Data$NetStorage
con = Data$Consumption
price = Data$Price
# capacity = Data$Capacity
Date = Data$Date
base = Data$Base.Gas
Future1 = Data$Future1

EveryJan = Date[seq(from = 1,to = length(Date),by = 12)]



# Production and Consumption of Natural Gas -------------------------------
plot(con,type = 'l',col = 1,xaxt = 'n',ylim = c(min(con),max(con)),xlab = NA, ylab = NA,lwd = 5)
title(main = 'Production and Consumption of Natural Gas')
axis(1,seq(from = 1,to = length(Date),by = 12),labels = EveryJan)
par(new = TRUE)
plot(prod,type = 'l', col = 2, axes = FALSE,xlab = NA, ylab = NA,,ylim = c(min(con),max(con)),lwd = 5)
legend(x = 1, y = max(con),legend = c('Consumption','Production'),col = 1:2,lty = 1,lwd = 5,bty = 'n')



# Production, Consumption and Price of Natural Gas ------------------------
plot(con,type = 'l',col = 1,xaxt = 'n',ylim = c(min(con),max(con)),xlab = NA, ylab = NA,lwd = 5)
title(main = 'Production, Consumption and Price of Natural Gas')
axis(1,seq(from = 1,to = length(Date),by = 12),labels = EveryJan)
par(new = TRUE)
plot(prod,type = 'l', col = 2, axes = FALSE,xlab = NA, ylab = NA,,ylim = c(min(con),max(con)),lwd = 5)
par(new = TRUE)
plot(price,type = 'l', col = 3,axes = FALSE, xlab = NA, ylab = NA,lwd = 5)
axis(side = 4)
legend(x = 5, y = max(price),legend = c('Consumption','Production','Price'),col = 1:3,lty = 1)

# Production, Consumption, Price and storage of Natural Gas ------------------------
plot(con,type = 'l',col = 1,xaxt = 'n',ylim = c(min(con),max(con)),xlab = NA, ylab = NA,lwd = 5)
title(main = 'Production, Consumption and Price of Natural Gas')
axis(1,seq(from = 1,to = length(Date),by = 12),labels = EveryJan)
par(new = TRUE)
plot(prod,type = 'l', col = 2, axes = FALSE,xlab = NA, ylab = NA,,ylim = c(min(con),max(con)),lwd = 5)
par(new = TRUE)
plot(price,type = 'l', col = 3,axes = FALSE, xlab = NA, ylab = NA,lwd = 5)
axis(side = 4)
legend(x = 5, y = max(price),legend = c('Consumption','Production','Price'),col = 1:3,lty = 1)
par(new = TRUE)
plot(storage,type = 'l', col = 4,axes = FALSE, xlab = NA, ylab = NA,lwd = 5)
axis(side = 4)


# Spot Price of Natural Gas ----------------------------------------------------
plot(price,type = 'l', col = 3, xlab = NA, ylab = NA,lwd = 5,xaxt = 'n')
title(main = 'Price of Natural Gas')
axis(1,seq(from = 1,to = length(Date),by = 12),labels = EveryJan)


# Future 1 and spot Price of Natural Gas -------------------------------------------
plot(price,type = 'l', col = 3, xlab = NA, ylab = NA,lwd = 5,xaxt = 'n')
lines(Future1,type = 'l',col = 5,lwd = 5)
title(main = 'Future 1 and Spot Price of Natural Gas')
axis(1,seq(from = 1,to = length(Date),by = 12),labels = EveryJan)




# Price and Capacity of Storage of Natural Gas ----------------------------
plot(price,type = 'l', col = 3, xlab = NA, ylab = NA,lwd = 5,xaxt = 'n')
title(main = 'Price and Capacity of Storage of Natural Gas')
# par(new = TRUE)
# plot(netwithdraw,type = 'l', col = 4,axes = FALSE, xlab = NA, ylab = NA,lwd = 5,)
# axis(side = 4)
# axis(1,seq(from = 1,to = length(Date),by = 12),labels = EveryJan)


# Future Prices -----------------------------------------------------------
FutureData = read.csv('price.csv',header = FALSE)

FutureData = FutureData[,1]

DateFuture = paste0('Jan-',15:24)

plot(FutureData, type = 'l', xaxt = 'n', ylab = 'Price $', xlab = NA, lwd =5)
axis(1,seq(from = 1,to = length(FutureData),by = 12),labels = DateFuture)
title(main = 'Future Prices of Natural Gas')
abline(v = 14,col = 2,lwd = 5)

# Future volume ------------------------------------------------------------------

Vdata = read.csv('VolumeFuture.csv')
volume = Vdata$Total.Volume[-c(1,nrow(Vdata))]
plot(volume,type='l',xlab = NA, lwd = 5, xaxt = 'n')
title(main = 'The volumes of futures')
axis(1,seq(from = 1,to = length(volume),by = 12),labels = DateFuture[-c(9,10)])


# Future price first 14 
plot(FutureData[1:14], type = 'l', xaxt = 'n', ylab = 'Price $', xlab = NA, lwd =5)
axis(1,at = 1:14, labels =Vdata$Month[1:14])
title(main = 'Future Prices of Natural Gas')

# Prices and storage level ------------------------------------------------

# ratio = (storage-base)/(capacity-base)
# plot(price,type = 'l', col = 3, xlab = NA, ylab = NA,lwd = 5,xaxt = 'n')
# title(main = 'Price and Storage Level of Natural Gas')
# par(new = TRUE)
# plot(ratio,type = 'l', col = 6,axes = FALSE, xlab = NA, ylab = NA,lwd = 5)
# axis(side = 4)
# axis(1,seq(from = 1,to = length(Date),by = 12),labels = EveryJan)




